#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "comm_brick.h"
#include "comm_tiled.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "element_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "math_extra.h"
#include "error.h"
#include "memory.h"

using namespace CAC_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define BIG 1.0e20

enum{SINGLE,MULTI};               // same as in Comm
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

CommBrick::CommBrick(CAC *cac) : 
  Comm(cac),
  atom_sendnum(NULL), atom_recvnum(NULL), elem_sendnum(NULL), elem_recvnum(NULL),
  sendproc(NULL), recvproc(NULL), 
  atom_size_forward_recv(NULL), elem_size_forward_recv(NULL),
  atom_size_reverse_send(NULL), elem_size_reverse_send(NULL),  
  atom_size_reverse_recv(NULL), elem_size_reverse_recv(NULL),  
  slablo(NULL), slabhi(NULL), multilo(NULL), multihi(NULL),
  cutghostmulti(NULL), pbc_flag(NULL), pbc(NULL), 
  atom_firstrecv(NULL), elem_firstrecv(NULL), 
  atom_sendlist(NULL), atom_maxsendlist(NULL), 
  elem_sendlist(NULL), elem_maxsendlist(NULL), 
  buf_send(NULL), buf_recv(NULL)
{
  style = 0;
  layout = LAYOUT_UNIFORM;
  pbc_flag = NULL;
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommBrick::~CommBrick()
{
  free_swap();
  if (mode == MULTI) {
    free_multi();
    memory->destroy(cutghostmulti);
  }

  if (atom_sendlist) for (int i = 0; i < maxswap; i++) memory->destroy(atom_sendlist[i]);
  memory->sfree(atom_sendlist);
  if (elem_sendlist) for (int i = 0; i < maxswap; i++) memory->destroy(elem_sendlist[i]);
  memory->sfree(elem_sendlist);
  memory->destroy(atom_maxsendlist);
  memory->destroy(elem_maxsendlist);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommBrick::CommBrick(CAC *cac, Comm *oldcomm) : Comm(*oldcomm)
{
  if (oldcomm->layout == LAYOUT_TILED)
    error->all(FLERR,"Cannot change to comm_style brick from tiled layout");

  style = 0;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_buffers();
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommBrick
   ------------------------------------------------------------------------- */

void CommBrick::init_buffers()
{
  multilo = multihi = NULL;
  cutghostmulti = NULL;

  // bufextra = max size of one exchanged atom
  //          = allowed overflow of sendbuf in exchange()
  // atomvec, fix reset these 2 maxexchange values if needed
  // only necessary if their size > BUFEXTRA

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");

  nswap = 0;
  maxswap = 6;
  allocate_swap(maxswap);

  atom_sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"comm:atom_sendlist");
  elem_sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"comm:elem_sendlist");

  memory->create(atom_maxsendlist,maxswap,"comm:atom_maxsendlist");
  memory->create(elem_maxsendlist,maxswap,"comm:elem_maxsendlist");

  for (int i = 0; i < maxswap; i++) {
    atom_maxsendlist[i] = BUFMIN;
    elem_maxsendlist[i] = BUFMIN;
    memory->create(atom_sendlist[i],BUFMIN,"comm:atom_sendlist[i]");
    memory->create(elem_sendlist[i],BUFMIN,"comm:elem_sendlist[i]");
  }
}

/* ---------------------------------------------------------------------- */

void CommBrick::init()
{
  Comm::init();

  // memory for multi-style communication

  if (mode == MULTI && multilo == NULL) {
    allocate_multi(maxswap);
    memory->create(cutghostmulti,atom->ntypes+1,3,"comm:cutghostmulti");
  }
  if (mode == SINGLE && multilo) {
    free_multi();
    memory->destroy(cutghostmulti);
  }
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single mode sets slab boundaries (slablo,slabhi) based on max cutoff
   multi mode sets type-dependent slab boundaries (multilo,multihi)
   ------------------------------------------------------------------------- */

void CommBrick::setup()
{
  // cutghost[] = max distance at which ghost atoms need to be acquired
  // for orthogonal:
  //   cutghost is in box coords = neigh->cutghost in all 3 dims
  //   add ghostskin from compute in necessary
  // for triclinic:
  //   neigh->cutghost = distance between tilted planes in box coords
  //   cutghost is in lamda coords = distance between those planes
  // for multi:
  //   cutghostmulti = same as cutghost, only for each atom type

  int i,j;
  int ntypes = atom->ntypes;
  double *prd,*sublo,*subhi;

  double cut = MAX(neighbor->cutneighmax,cutghostuser);
 

  // check if additional skin from compute is larger

  for (i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->ghostskinflag) 
      cut = MAX(modify->compute[i]->ghostskin,cut);

  // if pair exists, check for threebody pair
  // double cut since neighbor of ghost is required

  if (force->pair)
    if (force->pair->threebody_flag)
      cut *= 2.0;

  // additional skin from element size



  if (triclinic == 0) {
    prd = domain->prd;
    sublo = domain->sublo;
    subhi = domain->subhi;
    cutghost[0] = cutghost[1] = cutghost[2] = cut;
    if (element->nelements) {
      cutghost[0] += element->max_size_array[0];
      cutghost[1] += element->max_size_array[1];
      cutghost[2] += element->max_size_array[2];
    }

    if (mode == MULTI) {
      error->all(FLERR,"Multi ghost mode has not been considered yet");
      //double *cuttype = neighbor->cuttype;
      //for (i = 1; i <= ntypes; i++) {
      //  cut = 0.0;
      //  if (cutusermulti) cut = cutusermulti[i];
      //  cutghostmulti[i][0] = MAX(cut,cuttype[i]);
      //  cutghostmulti[i][1] = MAX(cut,cuttype[i]);
      //  cutghostmulti[i][2] = MAX(cut,cuttype[i]);
      //}
    }
  } else {
    prd = domain->prd_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
    double *h_inv = domain->h_inv;
    double length0,length1,length2;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    length2 = h_inv[2];

    cutghost[0] = cut * length0;
    cutghost[1] = cut * length1;
    cutghost[2] = cut * length2;
    if (element->nelements) {
      double element_size_lamda[3];
      double element_size_x[3];
      element_size_x[0] = element->max_size_array[0] + domain->boxlo[0];
      element_size_x[1] = element->max_size_array[1] + domain->boxlo[1];
      element_size_x[2] = element->max_size_array[2] + domain->boxlo[2];
      domain->x2lamda(element_size_x,element_size_lamda);
      cutghost[0] += element_size_lamda[0];
      cutghost[1] += element_size_lamda[1];
      cutghost[2] += element_size_lamda[2];
    }



    if (mode == MULTI) {
      error->all(FLERR, "Multi ghost mode has not been considered yet");
      //   double *cuttype = neighbor->cuttype;
      //   for (i = 1; i <= ntypes; i++) {
      //     cut = 0.0;
      //     if (cutusermulti) cut = cutusermulti[i];
      //     cutghostmulti[i][0] = length0 * MAX(cut,cuttype[i]);
      //     cutghostmulti[i][1] = length1 * MAX(cut,cuttype[i]);
      //     cutghostmulti[i][2] = length2 * MAX(cut,cuttype[i]);
      //   }
    }
  }

  // recvneed[idim][0/1] = # of procs away I recv atoms from, within cutghost
  //   0 = from left, 1 = from right
  //   do not cross non-periodic boundaries, need[2] = 0 for 2d
  // sendneed[idim][0/1] = # of procs away I send atoms to
  //   0 = to left, 1 = to right
  //   set equal to recvneed[idim][1/0] of neighbor proc
  // maxneed[idim] = max procs away any proc recvs atoms in either direction
  // layout = UNIFORM = uniform sized sub-domains:
  //   maxneed is directly computable from sub-domain size
  //     limit to procgrid-1 for non-PBC
  //   recvneed = maxneed except for procs near non-PBC
  //   sendneed = recvneed of neighbor on each side
  // layout = NONUNIFORM = non-uniform sized sub-domains:
  //   compute recvneed via updown() which accounts for non-PBC
  //   sendneed = recvneed of neighbor on each side
  //   maxneed via Allreduce() of recvneed

  int *periodicity = domain->periodicity;
  int left,right;

  if (layout == LAYOUT_UNIFORM) {
    maxneed[0] = static_cast<int> (cutghost[0] * procgrid[0] / prd[0]) + 1;
    maxneed[1] = static_cast<int> (cutghost[1] * procgrid[1] / prd[1]) + 1;
    maxneed[2] = static_cast<int> (cutghost[2] * procgrid[2] / prd[2]) + 1;
    if (domain->dimension == 2) maxneed[2] = 0;
    if (!periodicity[0]) maxneed[0] = MIN(maxneed[0],procgrid[0]-1);
    if (!periodicity[1]) maxneed[1] = MIN(maxneed[1],procgrid[1]-1);
    if (!periodicity[2]) maxneed[2] = MIN(maxneed[2],procgrid[2]-1);

    if (!periodicity[0]) {
      recvneed[0][0] = MIN(maxneed[0],myloc[0]);
      recvneed[0][1] = MIN(maxneed[0],procgrid[0]-myloc[0]-1);
      left = myloc[0] - 1;
      if (left < 0) left = procgrid[0] - 1;
      sendneed[0][0] = MIN(maxneed[0],procgrid[0]-left-1);
      right = myloc[0] + 1;
      if (right == procgrid[0]) right = 0;
      sendneed[0][1] = MIN(maxneed[0],right);
    } else recvneed[0][0] = recvneed[0][1] =
      sendneed[0][0] = sendneed[0][1] = maxneed[0];

    if (!periodicity[1]) {
      recvneed[1][0] = MIN(maxneed[1],myloc[1]);
      recvneed[1][1] = MIN(maxneed[1],procgrid[1]-myloc[1]-1);
      left = myloc[1] - 1;
      if (left < 0) left = procgrid[1] - 1;
      sendneed[1][0] = MIN(maxneed[1],procgrid[1]-left-1);
      right = myloc[1] + 1;
      if (right == procgrid[1]) right = 0;
      sendneed[1][1] = MIN(maxneed[1],right);
    } else recvneed[1][0] = recvneed[1][1] =
      sendneed[1][0] = sendneed[1][1] = maxneed[1];

    if (!periodicity[2]) {
      recvneed[2][0] = MIN(maxneed[2],myloc[2]);
      recvneed[2][1] = MIN(maxneed[2],procgrid[2]-myloc[2]-1);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = MIN(maxneed[2],procgrid[2]-left-1);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = MIN(maxneed[2],right);
    } else recvneed[2][0] = recvneed[2][1] =
      sendneed[2][0] = sendneed[2][1] = maxneed[2];

  } else {
    recvneed[0][0] = updown(0,0,myloc[0],prd[0],periodicity[0],xsplit);
    recvneed[0][1] = updown(0,1,myloc[0],prd[0],periodicity[0],xsplit);
    left = myloc[0] - 1;
    if (left < 0) left = procgrid[0] - 1;
    sendneed[0][0] = updown(0,1,left,prd[0],periodicity[0],xsplit);
    right = myloc[0] + 1;
    if (right == procgrid[0]) right = 0;
    sendneed[0][1] = updown(0,0,right,prd[0],periodicity[0],xsplit);

    recvneed[1][0] = updown(1,0,myloc[1],prd[1],periodicity[1],ysplit);
    recvneed[1][1] = updown(1,1,myloc[1],prd[1],periodicity[1],ysplit);
    left = myloc[1] - 1;
    if (left < 0) left = procgrid[1] - 1;
    sendneed[1][0] = updown(1,1,left,prd[1],periodicity[1],ysplit);
    right = myloc[1] + 1;
    if (right == procgrid[1]) right = 0;
    sendneed[1][1] = updown(1,0,right,prd[1],periodicity[1],ysplit);

    if (domain->dimension == 3) {
      recvneed[2][0] = updown(2,0,myloc[2],prd[2],periodicity[2],zsplit);
      recvneed[2][1] = updown(2,1,myloc[2],prd[2],periodicity[2],zsplit);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = updown(2,1,left,prd[2],periodicity[2],zsplit);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = updown(2,0,right,prd[2],periodicity[2],zsplit);
    } else recvneed[2][0] = recvneed[2][1] =
      sendneed[2][0] = sendneed[2][1] = 0;

    int all[6];
    MPI_Allreduce(&recvneed[0][0],all,6,MPI_INT,MPI_MAX,world);
    maxneed[0] = MAX(all[0],all[1]);
    maxneed[1] = MAX(all[2],all[3]);
    maxneed[2] = MAX(all[4],all[5]);
  }

  // allocate comm memory

  nswap = 2 * (maxneed[0]+maxneed[1]+maxneed[2]);
  if (nswap > maxswap) grow_swap(nswap);

  // setup parameters for each exchange:
  // sendproc = proc to send to at each swap
  // recvproc = proc to recv from at each swap
  // for mode SINGLE:
  //   slablo/slabhi = boundaries for slab of atoms to send at each swap
  //   use -BIG/midpt/BIG to insure all atoms included even if round-off occurs
  //   if round-off, atoms recvd across PBC can be < or > than subbox boundary
  //   note that borders() only loops over subset of atoms during each swap
  //   treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
  // for mode MULTI:
  //   multilo/multihi is same, with slablo/slabhi for each atom type
  // pbc_flag: 0 = nothing across a boundary, 1 = something across a boundary
  // pbc = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
  // for triclinic, slablo/hi and pbc_border will be used in lamda (0-1) coords
  // 1st part of if statement is sending to the west/south/down
  // 2nd part of if statement is sending to the east/north/up

  int dim,ineed;

  int iswap = 0;
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2*maxneed[dim]; ineed++) {
      pbc_flag[iswap] = 0;
      pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
        pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;

      if (ineed % 2 == 0) {
        sendproc[iswap] = procneigh[dim][0];
        recvproc[iswap] = procneigh[dim][1];
        if (mode == SINGLE) {
          if (ineed < 2) slablo[iswap] = -BIG;
          else slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
          slabhi[iswap] = sublo[dim] + cutghost[dim];
        } else {
          for (i = 1; i <= ntypes; i++) {
            if (ineed < 2) multilo[iswap][i] = -BIG;
            else multilo[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
            multihi[iswap][i] = sublo[dim] + cutghostmulti[i][dim];
          }
        }
        if (myloc[dim] == 0) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = 1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = 1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = 1;
          }
        }

      } else {
        sendproc[iswap] = procneigh[dim][1];
        recvproc[iswap] = procneigh[dim][0];
        if (mode == SINGLE) {
          slablo[iswap] = subhi[dim] - cutghost[dim];
          if (ineed < 2) slabhi[iswap] = BIG;
          else slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
        } else {
          for (i = 1; i <= ntypes; i++) {
            multilo[iswap][i] = subhi[dim] - cutghostmulti[i][dim];
            if (ineed < 2) multihi[iswap][i] = BIG;
            else multihi[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
          }
        }
        if (myloc[dim] == procgrid[dim]-1) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = -1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = -1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = -1;
          }
        }
      }

      iswap++;
    }
  }
}

/* ----------------------------------------------------------------------
   walk up/down the extent of nearby processors in dim and dir
   loc = myloc of proc to start at
   dir = 0/1 = walk to left/right
   do not cross non-periodic boundaries
   is not called for z dim in 2d
   return how many procs away are needed to encompass cutghost away from loc
   ------------------------------------------------------------------------- */

int CommBrick::updown(int dim, int dir, int loc,
    double prd, int periodicity, double *split)
{
  int index,count;
  double frac,delta;

  if (dir == 0) {
    frac = cutghost[dim]/prd;
    index = loc - 1;
    delta = 0.0;
    count = 0;
    while (delta < frac) {
      if (index < 0) {
        if (!periodicity) break;
        index = procgrid[dim] - 1;
      }
      count++;
      delta += split[index+1] - split[index];
      index--;
    }

  } else {
    frac = cutghost[dim]/prd;
    index = loc + 1;
    delta = 0.0;
    count = 0;
    while (delta < frac) {
      if (index >= procgrid[dim]) {
        if (!periodicity) break;
        index = 0;
      }
      count++;
      delta += split[index+1] - split[index];
      index++;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   forward communication of atom/element coords every timestep
   other per-atom/element attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm(int dummy)
{
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **x = atom->x;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {

    // exchange atoms

    if (atom->natoms) {

      if (sendproc[iswap] != me) {
        if (comm_x_only) {
          if (atom_size_forward_recv[iswap]) {
            if (atom_size_forward_recv[iswap]) buf = x[atom_firstrecv[iswap]];
            else buf = NULL;
            MPI_Irecv(buf,atom_size_forward_recv[iswap],MPI_DOUBLE,
                recvproc[iswap],0,world,&request);
          }
          n = avec->pack_comm(atom_sendnum[iswap],atom_sendlist[iswap],
              buf_send,pbc_flag[iswap],pbc[iswap]);
          if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (atom_size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
        } else if (ghost_velocity) {
          error->all(FLERR,"ghost_velocity case has not been considered yet");
          //  if (size_forward_recv[iswap])
          //    MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
          //        recvproc[iswap],0,world,&request);
          //  n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
          //      buf_send,pbc_flag[iswap],pbc[iswap]);
          //  if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          //  if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          //  avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
        } else {
          if (atom_size_forward_recv[iswap])
            MPI_Irecv(buf_recv,atom_size_forward_recv[iswap],MPI_DOUBLE,
                recvproc[iswap],0,world,&request);
          n = avec->pack_comm(atom_sendnum[iswap],atom_sendlist[iswap],
              buf_send,pbc_flag[iswap],pbc[iswap]);
          if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (atom_size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          avec->unpack_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf_recv);
        }

      } else {
        if (comm_x_only) {
          if (atom_sendnum[iswap])
            avec->pack_comm(atom_sendnum[iswap],atom_sendlist[iswap],
                x[atom_firstrecv[iswap]],pbc_flag[iswap],pbc[iswap]);
        } else if (ghost_velocity) {
          error->all(FLERR,"ghost_velocity case has not been considered yet");
          //avec->pack_comm_vel(atom_sendnum[iswap],sendlist[iswap],
          //    buf_send,pbc_flag[iswap],pbc[iswap]);
          //avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
        } else {
          avec->pack_comm(atom_sendnum[iswap],atom_sendlist[iswap],
              buf_send,pbc_flag[iswap],pbc[iswap]);
          avec->unpack_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf_send);
        }
      }
    }

    // exchange elements

    if (element->nelements) {
      if (sendproc[iswap] != me) {
        if (elem_size_forward_recv[iswap]) 
          MPI_Irecv(buf_recv, elem_size_forward_recv[iswap],MPI_DOUBLE,
              recvproc[iswap],0,world,&request);
        n = evec->pack_comm(elem_sendnum[iswap],elem_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (elem_size_forward_recv[iswap]) MPI_Wait(&request, MPI_STATUS_IGNORE);
        evec->unpack_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf_recv);

      } else {
        evec->pack_comm(elem_sendnum[iswap],elem_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);
        evec->unpack_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms && element nodes every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm()
{
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **f = atom->f;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  if (atom->natoms)
    for (int iswap = nswap-1; iswap >= 0; iswap--) {
      if (sendproc[iswap] != me) {

        // exchange atoms

        if (comm_f_only) {
          if (atom_size_reverse_recv[iswap])
            MPI_Irecv(buf_recv,atom_size_reverse_recv[iswap],MPI_DOUBLE,
                sendproc[iswap],0,world,&request);
          if (atom_size_reverse_send[iswap]) {
            buf = f[atom_firstrecv[iswap]];
            MPI_Send(buf,atom_size_reverse_send[iswap],MPI_DOUBLE,
                recvproc[iswap],0,world);
          }
          if (atom_size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
        } else {
          if (atom_size_reverse_recv[iswap])
            MPI_Irecv(buf_recv,atom_size_reverse_recv[iswap],MPI_DOUBLE,
                sendproc[iswap],0,world,&request);
          n = avec->pack_reverse(atom_recvnum[iswap],atom_firstrecv[iswap],buf_send);
          if (n) MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          if (atom_size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
        avec->unpack_reverse(atom_sendnum[iswap],atom_sendlist[iswap],buf_recv);

      } else {
        if (comm_f_only) {
          if (atom_sendnum[iswap])
            avec->unpack_reverse(atom_sendnum[iswap],atom_sendlist[iswap],
                f[atom_firstrecv[iswap]]);
        } else {
          avec->pack_reverse(atom_recvnum[iswap],atom_firstrecv[iswap],buf_send);
          avec->unpack_reverse(atom_sendnum[iswap],atom_sendlist[iswap],buf_send);
        }
      }
    }

  // exchange element if element exist and threebody potential

  if (element->nelements && force->pair->threebody_flag) {
    for (int iswap = nswap-1; iswap >= 0; iswap--) {
      if (sendproc[iswap] != me) {
        if (elem_size_reverse_recv[iswap])
          MPI_Irecv(buf_recv,elem_size_reverse_recv[iswap],MPI_DOUBLE,
              sendproc[iswap],0,world,&request);
        n = evec->pack_reverse(elem_recvnum[iswap],elem_firstrecv[iswap],buf_send);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
        if (elem_size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);

        evec->unpack_reverse(elem_sendnum[iswap],elem_sendlist[iswap],buf_recv);
      } else {
        evec->pack_reverse(elem_recvnum[iswap],elem_firstrecv[iswap],buf_send);
        evec->unpack_reverse(elem_sendnum[iswap],elem_sendlist[iswap],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   exchange move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a stencil proc's box
   can happen if atom moves outside of non-periodic bounary
   or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
   ------------------------------------------------------------------------- */

void CommBrick::exchange()
{

  bigint n,sum;
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *sublo,*subhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  int dimension = domain->dimension;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (atom_map_style) atom->map_clear();
  atom->nghost = 0;
  if (elem_map_style) element->map_clear();
  element->nghost = 0;

  // insure send buf is large enough for single atom
  // bufextra = max size of one atom = allowed overflow of sendbuf
  // fixes can change per-atom size requirement on-the-fly

  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // Exchange Atoms
  // loop over dimensions

  if (atom->natoms) {
    for (int dim = 0; dim < dimension; dim++) {

      // fill buffer with atoms leaving my box, using < and >=
      // when atom is deleted, fill it in with last atom

      x = atom->x;
      nlocal = atom->nlocal;
      lo = sublo[dim];
      hi = subhi[dim];

      i = nsend = 0;

      while (i < nlocal) {
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend) grow_send(nsend,1);
          nsend += avec->pack_exchange(i,&buf_send[nsend]);
          avec->copy(nlocal-1,i,1);
          nlocal--;
        } else i++;
      }

      atom->nlocal = nlocal;

      // send/recv atoms in both directions
      // send size of message first so receiver can realloc buf_recv if needed
      // if 1 proc in dimension, no send/recv
      //   set nrecv = 0 so buf_send atoms will be lost
      // if 2 procs in dimension, single send/recv
      // if more than 2 procs in dimension, send/recv to both neighbors

      if (procgrid[dim] == 1) nrecv = 0;
      else {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
            &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,
            MPI_STATUS_IGNORE);
        nrecv = nrecv1;
        if (procgrid[dim] > 2) {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
              &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,
              MPI_STATUS_IGNORE);
          nrecv += nrecv2;
        }
        if (nrecv > maxrecv) grow_recv(nrecv);

        MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
            world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        if (procgrid[dim] > 2) {
          MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
              world,&request);
          MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
      }

      // check incoming atoms to see if they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      //   atom may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        value = buf_recv[m+dim+1];

        if (value >= lo && value < hi) m += avec->unpack_exchange(&buf_recv[m]);
        else m += static_cast<int> (buf_recv[m]);
      }
    }

    n = atom->nlocal;
    MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
    if (sum != atom->natoms) error->all(FLERR,"Atom lost");
  }

  // Exchange Elements
  // loop over dimensions

  if (element->nelements) {
    for (int dim = 0; dim < dimension; dim++) {

      // fill buffer with elements leaving my box, using < and >=
      // when element is deleted, fill it with last element

      x = element->x;
      nlocal = element->nlocal;
      lo = sublo[dim];
      hi = subhi[dim];
      i = nsend = 0;
      MPI_Barrier(world);
      while (i < nlocal) {

        // check if element i center is in my box
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend) grow_send(nsend,1);
          nsend += element->evec->pack_exchange(i,&buf_send[nsend]);
          evec->copy(nlocal-1,i,1);
          nlocal--;
        } else i++;
      }
      element->nlocal = nlocal;

      // send/recv elements in both directions
      // send size of message first so receiver can realloc buf_recv if needed
      // if 1 proc in dimension, no send/recv
      //   set nrecv = 0 so buf_send elements will be lost
      // if 2 procs in dimension, single send/recv
      // if more than 2 procs in dimension, send/recv to both neighbors

      if (procgrid[dim] == 1) nrecv = 0;
      else {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,&nrecv1,1,
            MPI_INT,procneigh[dim][1],0,world,MPI_STATUS_IGNORE);
        nrecv = nrecv1;
        if (procgrid[dim] > 2) {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,&nrecv2,1,
              MPI_INT,procneigh[dim][0],0,world,MPI_STATUS_IGNORE);
          nrecv +=nrecv2;
        }
        if (nrecv > maxrecv) grow_recv(nrecv);

        MPI_Irecv(buf_recv,nrecv1,
            MPI_DOUBLE,procneigh[dim][1],0,world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        if (procgrid[dim] > 2) {
          MPI_Irecv(&buf_recv[nrecv1],nrecv2,
              MPI_DOUBLE,procneigh[dim][0],0,world,&request);
          MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
      }

      // check incoming elements to see if they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      // element may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        value = buf_recv[m+dim+1];
        if (value >= lo && value < hi) m += evec->unpack_exchange(&buf_recv[m]);
        else m += static_cast<int> (buf_recv[m]);
      } 
    }

    n = element->nlocal;
    MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
    if (sum != element->nelements) error->all(FLERR,"Element lost");
  }

}

/* ----------------------------------------------------------------------
   borders list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a forward_comm(), so don't need to explicitly
   call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms/elements must be in lamda coords (0-1) before borders is called
   ------------------------------------------------------------------------- */

void CommBrick::borders()
{
  int i,n,itype,iswap,dim,ineed,twoneed,max,nlocal,nghost;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;

  // do swaps over all 3 dimensions

  // Swap atoms

  if (atom->natoms) {

    iswap = 0;
    smax = rmax = 0;

    for (dim = 0; dim < 3; dim++) {
      nlast = 0;
      twoneed = 2*maxneed[dim];
      for (ineed = 0; ineed < twoneed; ineed++) {

        // find atoms within slab boundaries lo/hi using <= and >=
        // check atoms between nfirst and nlast
        //   for first swaps in a dim, check owned and ghost
        //   for later swaps in a dim, only check newly arrived ghosts
        // store sent atom indices in sendlist for use in future timesteps

        x = atom->x;
        nlocal = atom->nlocal;
        nghost = atom->nghost;

        if (mode == SINGLE) {
          lo = slablo[iswap];
          hi = slabhi[iswap];
        } else {
          type = atom->type;
          mlo = multilo[iswap];
          mhi = multihi[iswap];
        }
        if (ineed % 2 == 0) {
          nfirst = nlast;
          nlast = nlocal + nghost;
        }

        nsend = 0;

        // sendflag = 0 if I do not send on this swap
        // sendneed test indicates receiver no longer requires data
        // e.g. due to non-PBC or non-uniform sub-domains

        if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
        else sendflag = 1;

        // find send atoms according to SINGLE vs MULTI
        // all atoms eligible versus only atoms in bordergroup
        // can only limit loop to bordergroup for first sends (ineed < 2)
        // on these sends, break loop in two: owned (in group) and ghost

        if (sendflag) {
          if (!bordergroup || ineed >= 2) {
            if (mode == SINGLE) {
              for (i = nfirst; i < nlast; i++)
                if (x[i][dim] >= lo && x[i][dim] <= hi) {
                  if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
                  atom_sendlist[iswap][nsend++] = i;
                }
            } else {
              for (i = nfirst; i < nlast; i++) {
                itype = type[i];
                if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                  if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
                  atom_sendlist[iswap][nsend++] = i;
                }
              }
            }

          } else {
            error->one(FLERR,"fisrt group not supported");
            //if (mode == SINGLE) {
            //  ngroup = atom->nfirst;
            //  for (i = 0; i < ngroup; i++)
            //    if (x[i][dim] >= lo && x[i][dim] <= hi) {
            //      if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
            //      atom_sendlist[iswap][nsend++] = i;
            //    }
            //  for (i = nlocal; i < nlast; i++)
            //    if (x[i][dim] >= lo && x[i][dim] <= hi) {
            //      if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
            //      atom_sendlist[iswap][nsend++] = i;
            //    }
            //} else {
            //  ngroup = atom->nfirst;
            //  for (i = 0; i < ngroup; i++) {
            //    itype = type[i];
            //    if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
            //      if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
            //      atom_sendlist[iswap][nsend++] = i;
            //    }
            //  }
            //  for (i = nlocal; i < nlast; i++) {
            //    itype = type[i];
            //    if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
            //      if (nsend == atom_maxsendlist[iswap]) grow_atom_list(iswap,nsend);
            //      atom_sendlist[iswap][nsend++] = i;
            //    }
            //  }
            //}
          }
        }

        // pack up list of border atoms

        if (nsend*atom_size_border > maxsend) grow_send(nsend*atom_size_border,0);
        if (ghost_velocity)
          error->all(FLERR,"Ghost velocity has not been considered yet");
        //n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
        //    pbc_flag[iswap],pbc[iswap]);
        else
          n = avec->pack_border(nsend,atom_sendlist[iswap],buf_send,
              pbc_flag[iswap],pbc[iswap]);

        // swap atoms with other proc
        // no MPI calls except SendRecv if nsend/nrecv = 0
        // put incoming ghosts at end of my atom arrays
        // if swapping with self, simply copy, no messages

        if (sendproc[iswap] != me) {
          MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
              &nrecv,1,MPI_INT,recvproc[iswap],0,world,
              MPI_STATUS_IGNORE);
          if (nrecv*atom_size_border > maxrecv) grow_recv(nrecv*atom_size_border);
          if (nrecv) MPI_Irecv(buf_recv,nrecv*atom_size_border,MPI_DOUBLE,
              recvproc[iswap],0,world,&request);
          if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else {
          nrecv = nsend;
          buf = buf_send;
        }

        // unpack buffer

        if (ghost_velocity)
          error->all(FLERR,"Ghost velocity has not been considered yet");
        //avec->unpack_border_vel(nrecv,nlocal+nghost,buf);
        else
          avec->unpack_border(nrecv,nlocal+nghost,buf);

        // set all pointers & counters

        smax = MAX(smax,nsend);
        rmax = MAX(rmax,nrecv);
        atom_sendnum[iswap] = nsend;
        atom_recvnum[iswap] = nrecv;
        atom_size_forward_recv[iswap] = nrecv*atom_size_forward;
        atom_size_reverse_send[iswap] = nrecv*atom_size_reverse;
        atom_size_reverse_recv[iswap] = nsend*atom_size_reverse;
        atom_firstrecv[iswap] = nlocal + nghost;
        nghost += nrecv;
        atom->nghost = nghost;
        iswap++;
      }
    }

    // insure send/recv buffers are long enough for all forward & reverse comm

    max = MAX(atom_maxforward*smax,atom_maxreverse*rmax);
    if (max > maxsend) grow_send(max,0);
    max = MAX(atom_maxforward*rmax,atom_maxreverse*smax);
    if (max > maxrecv) grow_recv(max);

    // reset global->local map

    if (atom_map_style) atom->map_set();
  }

  // Swap elements

  if (element->nelements) {

    iswap = 0;
    smax = rmax = 0;
    for (dim = 0; dim < 3; dim++) {
      nlast = 0;
      twoneed = 2*maxneed[dim];
      for (ineed = 0; ineed < twoneed; ineed++) {

        // find atoms within slab boundaries lo/hi using <= and >=
        // check atoms between nfirst and nlast
        //   for first swaps in a dim, check owned and ghost
        //   for later swaps in a dim, only check newly arrived ghosts
        // store sent atom indices in sendlist for use in future timesteps

        x = element->x;
        nlocal = element->nlocal;
        nghost = element->nghost;

        if (mode == SINGLE) {
          lo = slablo[iswap];
          hi = slabhi[iswap];
        } else {
          type = element->ctype;
          mlo = multilo[iswap];
          mhi = multihi[iswap];
        }
        if (ineed % 2 == 0) {
          nfirst = nlast;
          nlast = nlocal + nghost;
        }

        nsend = 0;

        // sendflag = 0 if I do not send on this swap
        // sendneed test indicates receiver no longer requires data
        // e.g. due to non-PBC or non-uniform sub-domains

        if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
        else sendflag = 1;

        // find send atoms according to SINGLE vs MULTI
        // all atoms eligible versus only atoms in bordergroup
        // can only limit loop to bordergroup for first sends (ineed < 2)
        // on these sends, break loop in two: owned (in group) and ghost

        if (sendflag) {
          if (!bordergroup || ineed >= 2) {
            if (mode == SINGLE) {
              for (i = nfirst; i < nlast; i++)
                if (x[i][dim] >= lo && x[i][dim] <= hi) {
                  if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
                  elem_sendlist[iswap][nsend++] = i;
                }
            } else {
              for (i = nfirst; i < nlast; i++) {
                itype = type[i];
                if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                  if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
                  elem_sendlist[iswap][nsend++] = i;
                }
              }
            }

          } else {

            error->one(FLERR,"fisrt group not supported");
            //if (mode == SINGLE) {
            //  ngroup = element->nfirst;
            //  for (i = 0; i < ngroup; i++)
            //    if (x[i][dim] >= lo && x[i][dim] <= hi) {
            //      if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
            //      elem_sendlist[iswap][nsend++] = i;
            //    }
            //  for (i = nlocal; i < nlast; i++)
            //    if (x[i][dim] >= lo && x[i][dim] <= hi) {
            //      if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
            //      elem_sendlist[iswap][nsend++] = i;
            //    }
            //} else {
            //  ngroup = element->nfirst;
            //  for (i = 0; i < ngroup; i++) {
            //    itype = type[i];
            //    if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
            //      if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
            //      elem_sendlist[iswap][nsend++] = i;
            //    }
            //  }
            //  for (i = nlocal; i < nlast; i++) {
            //    itype = type[i];
            //    if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
            //      if (nsend == elem_maxsendlist[iswap]) grow_elem_list(iswap,nsend);
            //      elem_sendlist[iswap][nsend++] = i;
            //    }
            //  }
            //}
          }
        }

        // pack up list of border atoms

        if (nsend*elem_size_border > maxsend) grow_send(nsend*elem_size_border,0);
        if (ghost_velocity)
          error->all(FLERR,"Ghost velocity has not been considered yet");
        //n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
        //    pbc_flag[iswap],pbc[iswap]);
        else
          n = evec->pack_border(nsend,elem_sendlist[iswap],buf_send,
              pbc_flag[iswap],pbc[iswap]);

        // swap atoms with other proc
        // no MPI calls except SendRecv if nsend/nrecv = 0
        // put incoming ghosts at end of my atom arrays
        // if swapping with self, simply copy, no messages

        if (sendproc[iswap] != me) {
          MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
              &nrecv,1,MPI_INT,recvproc[iswap],0,world,
              MPI_STATUS_IGNORE);
          if (nrecv*elem_size_border > maxrecv) grow_recv(nrecv*elem_size_border);
          if (nrecv) MPI_Irecv(buf_recv,nrecv*elem_size_border,MPI_DOUBLE,
              recvproc[iswap],0,world,&request);
          if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else {
          nrecv = nsend;
          buf = buf_send;
        }

        // unpack buffer

        if (ghost_velocity)
          error->all(FLERR,"Ghost velocity has not been considered yet");
        //avec->unpack_border_vel(nrecv,nlocal+nghost,buf);
        else
          evec->unpack_border(nrecv,nlocal+nghost,buf);

        // set all pointers & counters

        smax = MAX(smax,nsend);
        rmax = MAX(rmax,nrecv);
        elem_sendnum[iswap] = nsend;
        elem_recvnum[iswap] = nrecv;
        elem_size_forward_recv[iswap] = nrecv*elem_size_forward;
        elem_size_reverse_send[iswap] = nrecv*elem_size_reverse;
        elem_size_reverse_recv[iswap] = nsend*elem_size_reverse;
        elem_firstrecv[iswap] = nlocal + nghost;
        nghost += nrecv; 
        element->nghost = nghost;
        iswap++;
      }
    }

    // insure send/recv buffers are long enough for all forward & reverse comm

    max = MAX(elem_maxforward*smax,elem_maxreverse*rmax);
    if (max > maxsend) grow_send(max,0);
    max = MAX(elem_maxforward*rmax,elem_maxreverse*smax);
    if (max > maxrecv) grow_recv(max);

    // reset global->local map

    if (elem_map_style) element->map_set();
  }

}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm_pair(Pair *pair)
{
  int iswap,n,nsize;
  double *buf;
  MPI_Request request;

  // communicate atoms

  if (atom->natoms) {

    nsize = pair->comm_atom_forward;

    if (nsize) {
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = pair->pack_atom_forward_comm(atom_sendnum[iswap],atom_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (atom_recvnum[iswap])
            MPI_Irecv(buf_recv,nsize*atom_recvnum[iswap],MPI_DOUBLE,
                recvproc[iswap],0,world,&request);
          if (atom_sendnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (atom_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        pair->unpack_atom_forward_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf);
      }
    }
  }

  // communicate elements

  if (element->nelements) {

    nsize = pair->comm_elem_forward;

    if (nsize) {
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = pair->pack_elem_forward_comm(elem_sendnum[iswap],elem_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (elem_recvnum[iswap])
            MPI_Irecv(buf_recv, elem_recvnum[iswap]*nsize, MPI_DOUBLE, 
                recvproc[iswap],0,world,&request);
          if (elem_sendnum[iswap]) 
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (elem_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv; 
        } else buf = buf_send;

        // unpack buffer

        pair->unpack_elem_forward_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm_pair(Pair *pair)
{
  int iswap,n,nsize;
  double *buf;
  MPI_Request request;

  if (atom->natoms) {
    nsize = MAX(pair->comm_atom_reverse,pair->comm_atom_reverse_off);

    // only communicate if nsize > 0

    if (nsize) {
      for (iswap = nswap-1; iswap >= 0; iswap--) {

        // pack buffer

        n = pair->pack_atom_reverse_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf_send);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (atom_sendnum[iswap])
            MPI_Irecv(buf_recv,nsize*atom_sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                world,&request);
          if (atom_recvnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          if (atom_sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        pair->unpack_atom_reverse_comm(atom_sendnum[iswap],atom_sendlist[iswap],buf);
      }
    }
  }

  if (element->nelements) {
    nsize = MAX(pair->comm_elem_reverse,pair->comm_elem_reverse_off);

    // only communicate if nsize > 0

    if (nsize) {
      for (iswap = nswap-1; iswap >= 0; iswap--) {

        // pack buffer

        n = pair->pack_elem_reverse_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf_send);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (elem_sendnum[iswap])
            MPI_Irecv(buf_recv,nsize*elem_sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                world,&request);
          if (elem_recvnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          if (elem_sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        pair->unpack_elem_reverse_comm(elem_sendnum[iswap],elem_sendlist[iswap],buf);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_forward
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm_fix(Fix *fix, int size)
{
  /*
     int iswap,n,nsize;
     double *buf;
     MPI_Request request;

     if (size) nsize = size;
     else nsize = fix->comm_forward;

     for (iswap = 0; iswap < nswap; iswap++) {

  // pack buffer

  n = fix->pack_forward_comm(sendnum[iswap],sendlist[iswap],
  buf_send,pbc_flag[iswap],pbc[iswap]);

  // exchange with another proc
  // if self, set recv buffer to send buffer

  if (sendproc[iswap] != me) {
  if (recvnum[iswap])
  MPI_Irecv(buf_recv,nsize*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
  world,&request);
  if (sendnum[iswap])
  MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
  if (recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
  buf = buf_recv;
  } else buf = buf_send;

  // unpack buffer

  fix->unpack_forward_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
  */
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_forward
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm_fix(Fix *fix, int size)
{
  /*
     int iswap,n,nsize;
     double *buf;
     MPI_Request request;

     if (size) nsize = size;
     else nsize = fix->comm_reverse;

     for (iswap = nswap-1; iswap >= 0; iswap--) {

  // pack buffer

  n = fix->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

  // exchange with another proc
  // if self, set recv buffer to send buffer

  if (sendproc[iswap] != me) {
  if (sendnum[iswap])
  MPI_Irecv(buf_recv,nsize*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
  world,&request);
  if (recvnum[iswap])
  MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
  if (sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
  buf = buf_recv;
  } else buf = buf_send;

  // unpack buffer

  fix->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
  */
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for pack size to insure buf_send is big enough
   handshake sizes before each Irecv/Send to insure buf_recv is big enough
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm_fix_variable(Fix *fix)
{
  /*
     int iswap,nsend,nrecv;
     double *buf;
     MPI_Request request;

     for (iswap = nswap-1; iswap >= 0; iswap--) {

  // pack buffer

  nsend = fix->pack_reverse_comm_size(recvnum[iswap],firstrecv[iswap]);
  if (nsend > maxsend) grow_send(nsend,0);
  nsend = fix->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

  // exchange with another proc
  // if self, set recv buffer to send buffer

  if (sendproc[iswap] != me) {
  MPI_Sendrecv(&nsend,1,MPI_INT,recvproc[iswap],0,
  &nrecv,1,MPI_INT,sendproc[iswap],0,world,
  MPI_STATUS_IGNORE);
  if (sendnum[iswap]) {
  if (nrecv > maxrecv) grow_recv(nrecv);
  MPI_Irecv(buf_recv,maxrecv,MPI_DOUBLE,sendproc[iswap],0,
  world,&request);
  }
  if (recvnum[iswap])
  MPI_Send(buf_send,nsend,MPI_DOUBLE,recvproc[iswap],0,world);
  if (sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
  buf = buf_recv;
  } else buf = buf_send;

  // unpack buffer

  fix->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
  */
}
/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm_compute(Compute *compute)
{

  int iswap,n,nsize;
  double *buf;
  MPI_Request request;

  if (atom->natoms) {
    nsize = compute->comm_atom_forward;
    if (nsize)
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = compute->pack_atom_forward_comm(atom_sendnum[iswap],atom_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (atom_recvnum[iswap])
            MPI_Irecv(buf_recv,nsize*atom_recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
          if (atom_sendnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (atom_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        compute->unpack_atom_forward_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf);
      }
  }

  if (element->nelements) {
    nsize = compute->comm_elem_forward;
    if (nsize)
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = compute->pack_elem_forward_comm(elem_sendnum[iswap],elem_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (elem_recvnum[iswap])
            MPI_Irecv(buf_recv,nsize*elem_recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
          if (elem_sendnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (elem_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        compute->unpack_elem_forward_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf);
      }
  }

}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm_compute(Compute *compute)
{
  int iswap,n;
  double *buf;
  MPI_Request request;

  int nsize;

  if (atom->natoms) {
    nsize = compute->comm_atom_reverse;
    if (nsize)
      for (iswap = nswap-1; iswap >= 0; iswap--) {

        // pack buffer

        n = compute->pack_atom_reverse_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf_send);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (atom_sendnum[iswap])
            MPI_Irecv(buf_recv,nsize*atom_sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                world,&request);
          if (atom_recvnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          if (atom_sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        compute->unpack_atom_reverse_comm(atom_sendnum[iswap],atom_sendlist[iswap],buf);
      }
  }

  if (element->nelements) {
    nsize = compute->comm_elem_reverse;
    if (nsize)
      for (iswap = nswap-1; iswap >= 0; iswap--) {

        // pack buffer

        n = compute->pack_elem_reverse_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf_send);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (elem_sendnum[iswap])
            MPI_Irecv(buf_recv,nsize*elem_sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
                world,&request);
          if (elem_recvnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          if (elem_sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        compute->unpack_elem_reverse_comm(elem_sendnum[iswap],elem_sendlist[iswap],buf);
      }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm_dump(Dump *dump)
{

  int iswap,n,nsize;
  double *buf;
  MPI_Request request;

  if (atom->natoms) {

    nsize = dump->comm_atom_forward;

    if (nsize)
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = dump->pack_atom_forward_comm(atom_sendnum[iswap],atom_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (atom_recvnum[iswap])
            MPI_Irecv(buf_recv,nsize*atom_recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
          if (atom_sendnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (atom_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        dump->unpack_atom_forward_comm(atom_recvnum[iswap],atom_firstrecv[iswap],buf);
      }
  }

  if (element->nelements) {

    nsize = dump->comm_elem_forward;

    if (nsize)
      for (iswap = 0; iswap < nswap; iswap++) {

        // pack buffer

        n = dump->pack_elem_forward_comm(elem_sendnum[iswap],elem_sendlist[iswap],
            buf_send,pbc_flag[iswap],pbc[iswap]);

        // exchange with another proc
        // if self, set recv buffer to send buffer

        if (sendproc[iswap] != me) {
          if (elem_recvnum[iswap])
            MPI_Irecv(buf_recv,nsize*elem_recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
          if (elem_sendnum[iswap])
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
          if (elem_recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          buf = buf_recv;
        } else buf = buf_send;

        // unpack buffer

        dump->unpack_elem_forward_comm(elem_recvnum[iswap],elem_firstrecv[iswap],buf);
      }
  }

}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommBrick::reverse_comm_dump(Dump *dump)
{
  /*
     int iswap,n;
     double *buf;
     MPI_Request request;

     int nsize = dump->comm_reverse;

     for (iswap = nswap-1; iswap >= 0; iswap--) {

  // pack buffer

  n = dump->pack_reverse_comm(recvnum[iswap],firstrecv[iswap],buf_send);

  // exchange with another proc
  // if self, set recv buffer to send buffer

  if (sendproc[iswap] != me) {
  if (sendnum[iswap])
  MPI_Irecv(buf_recv,nsize*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,
  world,&request);
  if (recvnum[iswap])
  MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
  if (sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
  buf = buf_recv;
  } else buf = buf_send;

  // unpack buffer

  dump->unpack_reverse_comm(sendnum[iswap],sendlist[iswap],buf);
  }
  */
}

/* ----------------------------------------------------------------------
   forward communication of N values in per-atom array
   ------------------------------------------------------------------------- */

void CommBrick::forward_comm_array(int nsize, double **array)
{
  /*
     int i,j,k,m,iswap,last;
     double *buf;
     MPI_Request request;

  // insure send/recv bufs are big enough for nsize
  // based on smax/rmax from most recent borders() invocation

  if (nsize > maxforward) {
  maxforward = nsize;
  if (maxforward*smax > maxsend) grow_send(maxforward*smax,0);
  if (maxforward*rmax > maxrecv) grow_recv(maxforward*rmax);
  }

  for (iswap = 0; iswap < nswap; iswap++) {

  // pack buffer

  m = 0;
  for (i = 0; i < sendnum[iswap]; i++) {
  j = sendlist[iswap][i];
  for (k = 0; k < nsize; k++)
  buf_send[m++] = array[j][k];
  }

  // exchange with another proc
  // if self, set recv buffer to send buffer

  if (sendproc[iswap] != me) {
  if (recvnum[iswap])
  MPI_Irecv(buf_recv,nsize*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
  world,&request);
  if (sendnum[iswap])
  MPI_Send(buf_send,nsize*sendnum[iswap],MPI_DOUBLE,
  sendproc[iswap],0,world);
  if (recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
  buf = buf_recv;
  } else buf = buf_send;

  // unpack buffer

  m = 0;
  last = firstrecv[iswap] + recvnum[iswap];
  for (i = firstrecv[iswap]; i < last; i++)
  for (k = 0; k < nsize; k++)
  array[i][k] = buf[m++];
  }
  */
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
   ------------------------------------------------------------------------- */

int CommBrick::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  /*
     int nsend,nrecv,nrecv1,nrecv2;
     MPI_Request request;

     nrecv = n;
     if (nrecv > maxrecv) grow_recv(nrecv);
     memcpy(buf_recv,inbuf,nrecv*sizeof(double));

  // loop over dimensions

  for (int dim = 0; dim < 3; dim++) {

  // no exchange if only one proc in a dimension

  if (procgrid[dim] == 1) continue;

  // send/recv info in both directions using same buf_recv
  // if 2 procs in dimension, single send/recv
  // if more than 2 procs in dimension, send/recv to both neighbors

  nsend = nrecv;
  MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
  &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,MPI_STATUS_IGNORE);
  nrecv += nrecv1;
  if (procgrid[dim] > 2) {
  MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
  &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,
  MPI_STATUS_IGNORE);
  nrecv += nrecv2;
  } else nrecv2 = 0;

  if (nrecv > maxrecv) grow_recv(nrecv);

  MPI_Irecv(&buf_recv[nsend],nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
  world,&request);
  MPI_Send(buf_recv,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
  MPI_Wait(&request,MPI_STATUS_IGNORE);

  if (procgrid[dim] > 2) {
  MPI_Irecv(&buf_recv[nsend+nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
  world,&request);
  MPI_Send(buf_recv,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
  MPI_Wait(&request,MPI_STATUS_IGNORE);
  }
  }

  outbuf = buf_recv;
  return nrecv;
  */
  return 0;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
   ------------------------------------------------------------------------- */

void CommBrick::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommBrick::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap atom_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommBrick::grow_atom_list(int iswap, int n)
{
  atom_maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(atom_sendlist[iswap],atom_maxsendlist[iswap],"comm:atom_sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap elem_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommBrick::grow_elem_list(int iswap, int n)
{
  elem_maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  memory->grow(elem_sendlist[iswap],elem_maxsendlist[iswap],"comm:elem_sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
   ------------------------------------------------------------------------- */

void CommBrick::grow_swap(int n)
{
  free_swap();
  allocate_swap(n);
  if (mode == MULTI) {
    free_multi();
    allocate_multi(n);
  }

  atom_sendlist = (int **)
    memory->srealloc(atom_sendlist,n*sizeof(int *),"comm:atom_sendlist");
  elem_sendlist = (int **)
    memory->srealloc(elem_sendlist,n*sizeof(int *),"comm:elem_sendlist");
  memory->grow(atom_maxsendlist,n,"comm:atom_maxsendlist");
  memory->grow(elem_maxsendlist,n,"comm:elem_maxsendlist");
  for (int i = maxswap; i < n; i++) {
    atom_maxsendlist[i] = BUFMIN;
    elem_maxsendlist[i] = BUFMIN;
    memory->create(atom_sendlist[i],BUFMIN,"comm:atom_sendlist[i]");
    memory->create(elem_sendlist[i],BUFMIN,"comm:elem_sendlist[i]");
  }
  maxswap = n;
}

/* ----------------------------------------------------------------------
   allocation of swap info
   ------------------------------------------------------------------------- */

void CommBrick::allocate_swap(int n)
{
  memory->create(atom_sendnum,n,"comm:atom_sendnum");
  memory->create(elem_sendnum,n,"comm:elem_sendnum");
  memory->create(atom_recvnum,n,"comm:atom_recvnum");
  memory->create(elem_recvnum,n,"comm:elem_recvnum");
  memory->create(sendproc,n,"comm:sendproc");
  memory->create(recvproc,n,"comm:recvproc");
  memory->create(atom_size_forward_recv,n,"comm:atom_size");
  memory->create(elem_size_forward_recv,n,"comm:elem_size");
  memory->create(atom_size_reverse_send,n,"comm:atom_size");
  memory->create(elem_size_reverse_send,n,"comm:elem_size");
  memory->create(atom_size_reverse_recv,n,"comm:atom_size");
  memory->create(elem_size_reverse_recv,n,"comm:elem_size");
  memory->create(slablo,n,"comm:slablo");
  memory->create(slabhi,n,"comm:slabhi");
  memory->create(atom_firstrecv,n,"comm:atom_firstrecv");
  memory->create(elem_firstrecv,n,"comm:elem_firstrecv");
  memory->create(pbc_flag,n,"comm:pbc_flag");
  memory->create(pbc,n,6,"comm:pbc");
}

/* ----------------------------------------------------------------------
   allocation of multi-type swap info
   ------------------------------------------------------------------------- */

void CommBrick::allocate_multi(int n)
{
  multilo = memory->create(multilo,n,atom->ntypes+1,"comm:multilo");
  multihi = memory->create(multihi,n,atom->ntypes+1,"comm:multihi");
}

/* ----------------------------------------------------------------------
   free memory for swaps
   ------------------------------------------------------------------------- */

void CommBrick::free_swap()
{
  memory->destroy(atom_sendnum);
  memory->destroy(elem_sendnum);
  memory->destroy(atom_recvnum);
  memory->destroy(elem_recvnum);
  memory->destroy(sendproc);
  memory->destroy(recvproc);
  memory->destroy(atom_size_forward_recv);
  memory->destroy(elem_size_forward_recv);
  memory->destroy(atom_size_reverse_send);
  memory->destroy(elem_size_reverse_send);
  memory->destroy(atom_size_reverse_recv);
  memory->destroy(elem_size_reverse_recv);
  memory->destroy(slablo);
  memory->destroy(slabhi);
  memory->destroy(atom_firstrecv);
  memory->destroy(elem_firstrecv);
  memory->destroy(pbc_flag);
  memory->destroy(pbc);
}

/* ----------------------------------------------------------------------
   free memory for multi-type swaps
   ------------------------------------------------------------------------- */

void CommBrick::free_multi()
{
  memory->destroy(multilo);
  memory->destroy(multihi);
  multilo = multihi = NULL;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint CommBrick::memory_usage()
{
  bigint bytes = 0;
  bytes += nprocs * sizeof(int);    // grid2proc
  for (int i = 0; i < nswap; i++) {
    bytes += memory->usage(atom_sendlist[i],atom_maxsendlist[i]);
    bytes += memory->usage(elem_sendlist[i],elem_maxsendlist[i]);
  }
  bytes += memory->usage(buf_send,maxsend+bufextra);
  bytes += memory->usage(buf_recv,maxrecv);
  return bytes;
}
