#include <string.h>
#include "comm_tiled.h"
#include "comm_brick.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "element_vec.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"

using namespace CAC_NS;
using namespace MathExtra;

#define BUFFACTOR 1.5
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define EPSILON 1.0e-6

#define DELTA_PROCS 16

enum{SINGLE,MULTI};               // same as in Comm
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(CAC *cac) : Comm(cac)
{

  style = 1;
  layout = LAYOUT_UNIFORM;
  pbc_flag = NULL;
  init_buffers();
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommTiled::CommTiled(CAC *cac, Comm *oldcomm) : Comm(*oldcomm)
{
  style = 1;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommTiled::~CommTiled()
{
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
  memory->destroy(overlap);
  deallocate_swap(nswap);
  memory->sfree(rcbinfo);
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
   ------------------------------------------------------------------------- */

void CommTiled::init_buffers()
{
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

  maxoverlap = 0;
  overlap = NULL;

  nswap = 2 * domain->dimension;

  allocate_swap(nswap);
  rcbinfo = NULL;
}

/* ---------------------------------------------------------------------- */

void CommTiled::init()
{

  Comm::init();

  // temporary restrictions

  if (triclinic)
    error->all(FLERR,"Cannot yet use comm_style tiled with triclinic box");
  if (mode == MULTI)
    error->all(FLERR,"Cannot yet use comm_style tiled with multi-mode comm");
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size and tiling
   ------------------------------------------------------------------------- */

void CommTiled::setup()
{
  int i,j,n;

  // domain properties used in setup method and methods it calls

  dimension = domain->dimension;
  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;

  int *periodicity = domain->periodicity;

  // set function pointers

  if (layout != LAYOUT_TILED) {
    box_drop = &CommTiled::box_drop_brick;
    box_other = &CommTiled::box_other_brick;
    box_touch = &CommTiled::box_touch_brick;
    point_drop = &CommTiled::point_drop_brick;
  } else {
    box_drop = &CommTiled::box_drop_tiled;
    box_other = &CommTiled::box_other_tiled;
    box_touch = &CommTiled::box_touch_tiled;
    point_drop = &CommTiled::point_drop_tiled;
  }

  // if RCB decomp exists and just changed, gather needed global RCB info

  if (layout == LAYOUT_TILED) coord2proc_setup();

  // set cutoff for comm forward and comm reverse
  // check that cutoff < any periodic box length

  double cut = MAX(neighbor->cutneighmax,cutghostuser);

  // additional skin from computes
  
  for (i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->ghostskinflag) 
      cut = MAX(modify->compute[i]->ghostskin,cut);

  // check for threebody pair
  // double cut since neighbor of ghost is required

  if (force->pair->threebody_flag)
    cut *= 2.0;

  // additional skin from element size
  
  if (element->nelements) cut += element->max_size;



  cutghost[0] = cutghost[1] = cutghost[2] = cut;

  if ((periodicity[0] && cut > prd[0]) ||
      (periodicity[1] && cut > prd[1]) ||
      (dimension == 3 && periodicity[2] && cut > prd[2]))
    error->all(FLERR,"Communication cutoff for comm_style tiled "
        "cannot exceed periodic box length");

  // if cut = 0.0, set to epsilon to induce nearest neighbor comm
  // this is b/c sendproc is used below to infer touching exchange procs
  // exchange procs will be empty (leading to lost atoms) if sendproc = 0
  // will reset sendproc/etc to 0 after exchange is setup, down below

  int cutzero = 0;
  if (cut == 0.0) {
    cutzero = 1;
    cut = MIN(prd[0],prd[1]);
    if (dimension == 3) cut = MIN(cut,prd[2]);
    cut *= EPSILON*EPSILON;
  }

  // setup forward/reverse communication
  // loop over 6 swap directions
  // determine which procs I will send to and receive from in each swap
  // done by intersecting ghost box with all proc sub-boxes it overlaps
  // sets nsendproc, nrecvproc, sendproc, recvproc
  // sets sendother, recvother, sendself, pbc_flag, pbc, sendbox
  // resets nprocmax

  int noverlap1,indexme;
  double lo1[3],hi1[3],lo2[3],hi2[3];
  int one,two;

  int iswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    for (int idir = 0; idir < 2; idir++) {

      // one = first ghost box in same periodic image
      // two = second ghost box wrapped across periodic boundary
      // either may not exist

      one = 1;
      copy3(sublo,lo1);
      copy3(subhi,hi1);
      if (idir == 0) {
        lo1[idim] = sublo[idim] - cut;
        hi1[idim] = sublo[idim];
      } else {
        lo1[idim] = subhi[idim];
        hi1[idim] = subhi[idim] + cut;
      }

      two = 0;
      if (periodicity[idim]) 
        if ((idir == 0 && lo1[idim] < boxlo[idim]) ||
            (idir == 1 && hi1[idim] > boxhi[idim])) two = 1;

      if (two) {
        copy3(sublo,lo2);
        copy3(subhi,hi2);
        if (idir == 0) {
          lo2[idim] = lo1[idim] + prd[idim];
          hi2[idim] = boxhi[idim];
          if (sublo[idim] == boxlo[idim]) one = 0;
        } else {
          lo2[idim] = boxlo[idim];
          hi2[idim] = hi1[idim] - prd[idim];
          if (subhi[idim] == boxhi[idim]) one = 0;
        }
      }

      if (one) {
        if (idir == 0) lo1[idim] = MAX(lo1[idim],boxlo[idim]);
        else hi1[idim] = MIN(hi1[idim],boxhi[idim]);
        if (lo1[idim] == hi1[idim]) one = 0;
      }

      // noverlap = # of overlaps of box1/2 with procs via box_drop()
      // overlap = list of overlapping procs
      // if overlap with self, indexme = index of me in list

      indexme = -1;
      noverlap = 0;
      if (one) (this->*box_drop)(idim,lo1,hi1,indexme);
      noverlap1 = noverlap;
      if (two) (this->*box_drop)(idim,lo2,hi2,indexme);

      // if self is in overlap list, move it to end of list

      if (indexme >= 0) {
        int tmp = overlap[noverlap-1];
        overlap[noverlap-1] = overlap[indexme];
        overlap[indexme] = tmp;
      }

      // reallocate 2nd dimensions of all send/recv arrays, based on noverlap
      // # of sends of this swap = # of recvs of iswap +/- 1

      if (noverlap > nprocmax[iswap]) {
        int oldmax = nprocmax[iswap];
        while (nprocmax[iswap] < noverlap) nprocmax[iswap] += DELTA_PROCS;
        grow_swap_send(iswap,nprocmax[iswap],oldmax);
        if (idir == 0) grow_swap_recv(iswap+1,nprocmax[iswap]);
        else grow_swap_recv(iswap-1,nprocmax[iswap]);
      }

      // overlap how has list of noverlap procs
      // includes PBC effects

      if (noverlap && overlap[noverlap-1] == me) sendself[iswap] = 1;
      else sendself[iswap] = 0;
      if (noverlap && noverlap-sendself[iswap]) sendother[iswap] = 1;
      else sendother[iswap] = 0;

      nsendproc[iswap] = noverlap;
      for (i = 0; i < noverlap; i++) sendproc[iswap][i] = overlap[i];

      if (idir == 0) {
        recvother[iswap+1] = sendother[iswap];
        nrecvproc[iswap+1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap+1][i] = overlap[i];
      } else {
        recvother[iswap-1] = sendother[iswap];
        nrecvproc[iswap-1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap-1][i] = overlap[i];
      }

      // compute sendbox for each of my sends
      // obox = intersection of ghostbox with other proc's sub-domain
      // sbox = what I need to send to other proc
      //      = sublo to MIN(sublo+cut,subhi) in idim, for idir = 0
      //      = MIN(subhi-cut,sublo) to subhi in idim, for idir = 1
      //      = obox in other 2 dims
      // if sbox touches other proc's sub-box boundaries in lower dims,
      //   extend sbox in those lower dims to include ghost atoms

      double oboxlo[3],oboxhi[3],sbox[6];

      for (i = 0; i < noverlap; i++) {
        pbc_flag[iswap][i] = 0;
        pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] =
          pbc[iswap][i][3] = pbc[iswap][i][4] = pbc[iswap][i][5] = 0;

        (this->*box_other)(idim,idir,overlap[i],oboxlo,oboxhi);

        if (i < noverlap1) {
          sbox[0] = MAX(oboxlo[0],lo1[0]);
          sbox[1] = MAX(oboxlo[1],lo1[1]);
          sbox[2] = MAX(oboxlo[2],lo1[2]);
          sbox[3] = MIN(oboxhi[0],hi1[0]);
          sbox[4] = MIN(oboxhi[1],hi1[1]);
          sbox[5] = MIN(oboxhi[2],hi1[2]);
        } else {
          pbc_flag[iswap][i] = 1;
          if (idir == 0) pbc[iswap][i][idim] = 1;
          else pbc[iswap][i][idim] = -1;
          sbox[0] = MAX(oboxlo[0],lo2[0]);
          sbox[1] = MAX(oboxlo[1],lo2[1]);
          sbox[2] = MAX(oboxlo[2],lo2[2]);
          sbox[3] = MIN(oboxhi[0],hi2[0]);
          sbox[4] = MIN(oboxhi[1],hi2[1]);
          sbox[5] = MIN(oboxhi[2],hi2[2]);
        }

        if (idir == 0) {
          sbox[idim] = sublo[idim];
          if (i < noverlap1) sbox[3+idim] = MIN(sbox[3+idim]+cut,subhi[idim]);
          else sbox[3+idim] = MIN(sbox[3+idim]-prd[idim]+cut,subhi[idim]);
        } else {
          if (i < noverlap1) sbox[idim] = MAX(sbox[idim]-cut,sublo[idim]);
          else sbox[idim] = MAX(sbox[idim]+prd[idim]-cut,sublo[idim]);
          sbox[3+idim] = subhi[idim];
        }

        if (idim >= 1) {
          if (sbox[0] == oboxlo[0]) sbox[0] -= cut;
          if (sbox[3] == oboxhi[0]) sbox[3] += cut;
        }
        if (idim == 2) {
          if (sbox[1] == oboxlo[1]) sbox[1] -= cut;
          if (sbox[4] == oboxhi[1]) sbox[4] += cut;
        }

        memcpy(sendbox[iswap][i],sbox,6*sizeof(double));
      }

      iswap++;
    }
  }

  // setup exchange communication = subset of forward/reverse comm procs
  // loop over dimensions
  // determine which procs I will exchange with in each dimension
  // subset of procs that touch my proc in forward/reverse comm
  // sets nexchproc & exchproc, resets nexchprocmax

  int proc;

  for (int idim = 0; idim < dimension; idim++) {

    // overlap = list of procs that touch my sub-box in idim
    // proc can appear twice in list if touches in both directions
    // 2nd add-to-list checks to insure each proc appears exactly once

    noverlap = 0;
    iswap = 2*idim;
    n = nsendproc[iswap];
    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc,idim,0)) {
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap,maxoverlap,"comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }
    noverlap1 = noverlap;
    iswap = 2*idim+1;
    n = nsendproc[iswap];

    MPI_Barrier(world);

    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc,idim,1)) {
        for (j = 0; j < noverlap1; j++)
          if (overlap[j] == proc) break;
        if (j < noverlap1) continue;
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap,maxoverlap,"comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }

    MPI_Barrier(world);

    // reallocate exchproc and exchnum if needed based on noverlap

    if (noverlap > nexchprocmax[idim]) {
      while (nexchprocmax[idim] < noverlap) nexchprocmax[idim] += DELTA_PROCS;
      delete [] exchproc[idim];
      exchproc[idim] = new int[nexchprocmax[idim]];
      delete [] exchnum[idim];
      exchnum[idim] = new int[nexchprocmax[idim]];
    }

    nexchproc[idim] = noverlap;
    for (i = 0; i < noverlap; i++) exchproc[idim][i] = overlap[i];
  }

  // reset sendproc/etc to 0 if cut is really 0.0

  if (cutzero) {
    for (i = 0; i < nswap; i++) {
      nsendproc[i] = nrecvproc[i] =
        sendother[i] = recvother[i] = sendself[i] = 0;
    }
  }

  // reallocate MPI Requests and Statuses as needed

  int nmax = 0;
  for (i = 0; i < nswap; i++) nmax = MAX(nmax,nprocmax[i]);
  for (i = 0; i < dimension; i++) nmax = MAX(nmax,nexchprocmax[i]);
  if (nmax > maxreqstat) {
    maxreqstat = nmax;
    delete [] requests;
    requests = new MPI_Request[maxreqstat];
  }

}

/* ----------------------------------------------------------------------
   forward communication of atom/element coords every timestep
   other per-atom/element attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm(int dummy)
{
  int i,irecv,n,nsend,nrecv;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **x = atom->x;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    // exchange atoms

    if (atom->natoms) {

      if (comm_x_only) {
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++)
            MPI_Irecv(x[atom_firstrecv[iswap][i]],atom_size_forward_recv[iswap][i],
                MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            n = avec->pack_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
                buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          avec->pack_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
              x[atom_firstrecv[iswap][nrecv]],pbc_flag[iswap][nsend],
              pbc[iswap][nsend]);
        }
        if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
      } else if (ghost_velocity) {
        //if (recvother[iswap]) {
        //  for (i = 0; i < nrecv; i++)
        //    MPI_Irecv(&buf_recv[atom_size_forward*atom_forward_recv_offset[iswap][i]],
        //        atom_size_forward_recv[iswap][i],
        //        MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        //}
        //if (sendother[iswap]) {
        //  for (i = 0; i < nsend; i++) {
        //    n = avec->pack_comm_vel(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
        //        buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        //    MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        //  }
        //}
        //if (sendself[iswap]) {
        //  avec->pack_comm_vel(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
        //      buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        //  avec->unpack_comm_vel(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
        //      buf_send);
        //}
        //if (recvother[iswap]) {
        //  for (i = 0; i < nrecv; i++) {
        //    MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        //    avec->unpack_comm_vel(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
        //        &buf_recv[atom_size_forward*
        //        atom_forward_recv_offset[iswap][irecv]]);
        //  }
        //}

        error->all(FLERR,"ghost_velocity case has not been considered yet");
      } else {
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++)
            MPI_Irecv(&buf_recv[atom_size_forward*atom_forward_recv_offset[iswap][i]],
                atom_size_forward_recv[iswap][i],
                MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            n = avec->pack_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
                buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          avec->pack_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
          avec->unpack_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
              buf_send);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
            avec->unpack_comm(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
                &buf_recv[atom_size_forward*
                atom_forward_recv_offset[iswap][irecv]]);
          }
        }
      }
    }

    // exchange elements

    if (element->nelements) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[elem_size_forward*elem_forward_recv_offset[iswap][i]],
              elem_size_forward_recv[iswap][i],
              MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = evec->pack_comm(elem_sendnum[iswap][i],elem_sendlist[iswap][i],
              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        evec->pack_comm(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
            buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        evec->unpack_comm(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
            buf_send);
      }

      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          evec->unpack_comm(elem_recvnum[iswap][irecv],elem_firstrecv[iswap][irecv],
              &buf_recv[elem_size_forward*
              elem_forward_recv_offset[iswap][irecv]]);
        }
      }
    }
  }
}
/* ----------------------------------------------------------------------
   reverse communication of forces on atoms (only) every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm()
{

  int i,irecv,n,nsend,nrecv;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **f = atom->f;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (comm_f_only) {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Irecv(&buf_recv[atom_size_reverse*atom_reverse_recv_offset[iswap][i]],
              atom_size_reverse_recv[iswap][i],MPI_DOUBLE,
              sendproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Send(f[atom_firstrecv[iswap][i]],atom_size_reverse_send[iswap][i],
              MPI_DOUBLE,recvproc[iswap][i],0,world);
      }
      if (sendself[iswap]) {
        avec->unpack_reverse(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
            f[atom_firstrecv[iswap][nrecv]]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          avec->unpack_reverse(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
              &buf_recv[atom_size_reverse*
              atom_reverse_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++)
          MPI_Irecv(&buf_recv[atom_size_reverse*atom_reverse_recv_offset[iswap][i]],
              atom_size_reverse_recv[iswap][i],MPI_DOUBLE,
              sendproc[iswap][i],0,world,&requests[i]);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          n = avec->pack_reverse(atom_recvnum[iswap][i],atom_firstrecv[iswap][i],
              buf_send);
          MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_reverse(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
            buf_send);
        avec->unpack_reverse(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
            buf_send);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          avec->unpack_reverse(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
              &buf_recv[atom_size_reverse*
              atom_reverse_recv_offset[iswap][irecv]]);
        }
      }
    }
  }

  // exchange elements if threebody potential

  if (force->pair->threebody_flag)
    for (int iswap = nswap-1; iswap >= 0; iswap--) {
      nsend = nsendproc[iswap] - sendself[iswap];
      nrecv = nrecvproc[iswap] - sendself[iswap];

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++)
          MPI_Irecv(&buf_recv[elem_size_reverse*elem_reverse_recv_offset[iswap][i]],
              elem_size_reverse_recv[iswap][i],MPI_DOUBLE,
              sendproc[iswap][i],0,world,&requests[i]);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          n = evec->pack_reverse(elem_recvnum[iswap][i],elem_firstrecv[iswap][i],
              buf_send);
          MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        evec->pack_reverse(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
            buf_send);
        evec->unpack_reverse(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
            buf_send);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          evec->unpack_reverse(elem_sendnum[iswap][irecv],elem_sendlist[iswap][irecv],
              &buf_recv[elem_size_reverse*
              elem_reverse_recv_offset[iswap][irecv]]);
        }
      }
    }


}

/* ----------------------------------------------------------------------
   exchange move atoms to correct processors
   atoms exchanged with procs that touch sub-box in each of 3 dims
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a touching proc's box
   can happen if atom moves outside of non-periodic bounary
   or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
   ------------------------------------------------------------------------- */

void CommTiled::exchange()
{

  int i,m,nexch,nsend,nrecv,nlocal,proc,offset;
  double lo,hi,value;
  double **x;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  dimension = domain->dimension;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (atom_map_style) atom->map_clear();
  if (elem_map_style) element->map_clear();
  atom->nghost = 0;
  element->nghost = 0;

  // insure send buf is large enough for single atom
  // bufextra = max size of one atom = allowed overflow of sendbuf
  // fixes can change per-atom size requirement on-the-fly

  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");

  // domain properties used in exchange method and methods it calls
  // subbox bounds for orthogonal or triclinic

  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  // exchange atoms

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
          proc = (this->*point_drop)(dim,x[i]);
          if (proc != me) {
            buf_send[nsend++] = proc;
            nsend += avec->pack_exchange(i,&buf_send[nsend]);
          }
          avec->copy(nlocal-1,i,1);
          nlocal--;
        } else i++;
      }
      atom->nlocal = nlocal;

      // send and recv atoms from neighbor procs that touch my sub-box in dim
      // no send/recv with self
      // send size of message first
      // receiver may receive multiple messages, realloc buf_recv if needed

      nexch = nexchproc[dim];
      if (!nexch) continue;

      for (m = 0; m < nexch; m++)
        MPI_Irecv(&exchnum[dim][m],1,MPI_INT,
            exchproc[dim][m],0,world,&requests[m]);
      for (m = 0; m < nexch; m++)
        MPI_Send(&nsend,1,MPI_INT,exchproc[dim][m],0,world);
      MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);

      nrecv = 0;
      for (m = 0; m < nexch; m++) nrecv += exchnum[dim][m];
      if (nrecv > maxrecv) grow_recv(nrecv);

      offset = 0;
      for (m = 0; m < nexch; m++) {
        MPI_Irecv(&buf_recv[offset],exchnum[dim][m],
            MPI_DOUBLE,exchproc[dim][m],0,world,&requests[m]);
        offset += exchnum[dim][m];
      }
      for (m = 0; m < nexch; m++)
        MPI_Send(buf_send,nsend,MPI_DOUBLE,exchproc[dim][m],0,world);
      MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);

      // check incoming atoms to see if I own it and they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      //   atom may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        proc = static_cast<int> (buf_recv[m++]);
        if (proc == me) {
          value = buf_recv[m+dim+1];
          if (value >= lo && value < hi) {
            m += avec->unpack_exchange(&buf_recv[m]);
            continue;
          }
        }
        m += static_cast<int> (buf_recv[m]);
      }
    }
  }


  // exchange elements

  if (element->nelements) {

    for (int dim = 0; dim < dimension; dim++) {

      // fill buffer with elements leaving my box, using < and >=
      // when element is deleted, fill it in with last element 

      x = element->x;
      nlocal = element->nlocal;

      lo = sublo[dim];
      hi = subhi[dim];
      i = nsend = 0;

      while (i < nlocal) {
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend) grow_send(nsend,1);
          proc = (this->*point_drop)(dim,x[i]);
          if (proc != me) {
            buf_send[nsend++] = proc;
            nsend += evec->pack_exchange(i,&buf_send[nsend]);
          }
          evec->copy(nlocal-1,i,1);
          nlocal--;
        } else i++;
      }

      element->nlocal = nlocal;

      // send and recv elements from neighbor procs that touch my sub-box in dim
      // no send/recv with self
      // send size of message first
      // receiver may receive multiple messages, realloc buf_recv if needed

      nexch = nexchproc[dim];
      if (!nexch) continue;

      for (m = 0; m < nexch; m++)
        MPI_Irecv(&exchnum[dim][m],1,MPI_INT,
            exchproc[dim][m],0,world,&requests[m]);
      for (m = 0; m < nexch; m++)
        MPI_Send(&nsend,1,MPI_INT,exchproc[dim][m],0,world);
      MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);
      nrecv = 0;
      for (m = 0; m < nexch; m++) nrecv += exchnum[dim][m];
      if (nrecv > maxrecv) grow_recv(nrecv);

      offset = 0;
      for (m = 0; m < nexch; m++) {
        MPI_Irecv(&buf_recv[offset],exchnum[dim][m],
            MPI_DOUBLE,exchproc[dim][m],0,world,&requests[m]);
        offset += exchnum[dim][m];
      }
      for (m = 0; m < nexch; m++)
        MPI_Send(buf_send,nsend,MPI_DOUBLE,exchproc[dim][m],0,world);
      MPI_Waitall(nexch,requests,MPI_STATUS_IGNORE);

      // check incoming elements to see if I own it and they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      //   element may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        proc = static_cast<int> (buf_recv[m++]);
        if (proc == me) {
          value = buf_recv[m+dim+1];
          if (value >= lo && value < hi) {
            m += evec->unpack_exchange(&buf_recv[m]);
            continue;
          }
        }
        m += static_cast<int> (buf_recv[m]);
      }
    }


  }
}

/* ----------------------------------------------------------------------
   borders list nearby atoms to send to neighboring procs at every timestep
   one list is created per swap/proc that will be made
   as list is made, actually do communication
   this does equivalent of a forward_comm(), so don't need to explicitly
   call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
   ------------------------------------------------------------------------- */

void CommTiled::borders()
{

  int i,m,n,nlast,nsend,nrecv,ngroup,ncount,ncountall;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  double **x;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;

  // send/recv max one = max # of atoms/elements in single send/recv for any swap
  // send/recv max all = max # of atoms/elements in all sends/recvs within any swap

  smaxone = smaxall = 0;
  rmaxone = rmaxall = 0;

  // loop over swaps in all dimensions

  // swap atoms 
  if (atom->natoms) {

    for (int iswap = 0; iswap < nswap; iswap++) {

      // find atoms within sendboxes using >= and <
      // hi test with ">" is important b/c don't want to send an atom
      //   in lower dim (on boundary) that a proc will recv again in higher dim
      // for x-dim swaps, check owned atoms
      // for yz-dim swaps, check owned and ghost atoms
      // store sent atom indices in atom_sendlist for use in future timesteps
      // NOTE: assume SINGLE mode, add logic for MULTI mode later

      x = atom->x;

      if (iswap % 2 == 0) nlast = atom->nlocal + atom->nghost;

      ncountall = 0;

      for (m = 0; m < nsendproc[iswap]; m++) {
        bbox = sendbox[iswap][m];
        xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
        xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];

        ncount = 0;

        if (!bordergroup) {
          for (i = 0; i < nlast; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi &&
                x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == atom_maxsendlist[iswap][m]) grow_atom_list(iswap,m,ncount);
              atom_sendlist[iswap][m][ncount++] = i;
            }
          }
        } else {
          ngroup = atom->nfirst;
          for (i = 0; i < ngroup; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi &&
                x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == atom_maxsendlist[iswap][m]) grow_atom_list(iswap,m,ncount);
              atom_sendlist[iswap][m][ncount++] = i;
            }
          }
          for (i = atom->nlocal; i < nlast; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi &&
                x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == atom_maxsendlist[iswap][m]) grow_atom_list(iswap,m,ncount);
              atom_sendlist[iswap][m][ncount++] = i;
            }
          }
        }

        atom_sendnum[iswap][m] = ncount;
        smaxone = MAX(smaxone,ncount);
        ncountall += ncount;
      }
      smaxall = MAX(smaxall,ncountall);

      // send atom_sendnum counts to procs who recv from me except self
      // copy data to self if sendself is set

      nsend = nsendproc[iswap] - sendself[iswap];
      nrecv = nrecvproc[iswap] - sendself[iswap];

      if (recvother[iswap])
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&atom_recvnum[iswap][m],1,MPI_INT,
              recvproc[iswap][m],0,world,&requests[m]);
      if (sendother[iswap])
        for (m = 0; m < nsend; m++)
          MPI_Send(&atom_sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],0,world);
      if (sendself[iswap]) atom_recvnum[iswap][nrecv] = atom_sendnum[iswap][nsend];
      if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);


      // setup other per swap/proc values from atom_sendnum and atom_recvnum

      for (m = 0; m < nsendproc[iswap]; m++) {
        atom_size_reverse_recv[iswap][m] = atom_sendnum[iswap][m]*atom_size_reverse;
        if (m == 0) atom_reverse_recv_offset[iswap][0] = 0;
        else atom_reverse_recv_offset[iswap][m] =
          atom_reverse_recv_offset[iswap][m-1] + atom_sendnum[iswap][m-1];
      }

      ncountall = 0;
      for (m = 0; m < nrecvproc[iswap]; m++) {
        ncount = atom_recvnum[iswap][m];
        rmaxone = MAX(rmaxone,ncount);
        ncountall += ncount;

        atom_size_forward_recv[iswap][m] = ncount*atom_size_forward;
        atom_size_reverse_send[iswap][m] = ncount*atom_size_reverse;
        if (m == 0) {
          atom_firstrecv[iswap][0] = atom->nlocal + atom->nghost;
          atom_forward_recv_offset[iswap][0] = 0;
        } else {
          atom_firstrecv[iswap][m] = atom_firstrecv[iswap][m-1] + atom_recvnum[iswap][m-1];
          atom_forward_recv_offset[iswap][m] =
            atom_forward_recv_offset[iswap][m-1] + atom_recvnum[iswap][m-1];
        }
      }
      rmaxall = MAX(rmaxall,ncountall);

      // insure send/recv buffers are large enough for this border comm swap

      if (smaxone*atom_size_border > maxsend) grow_send(smaxone*atom_size_border,0);
      if (rmaxall*atom_size_border > maxrecv) grow_recv(rmaxall*atom_size_border);

      // swap atoms with other procs using pack_border(), unpack_border()
      // use Waitall() instead of Waitany() because calls to unpack_border()
      //   must increment per-atom arrays in ascending order

      if (ghost_velocity) {
        if (recvother[iswap]) {
          for (m = 0; m < nrecv; m++)
            MPI_Irecv(&buf_recv[atom_size_border*atom_forward_recv_offset[iswap][m]],
                atom_recvnum[iswap][m]*atom_size_border,
                MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        }
        if (sendother[iswap]) {
          for (m = 0; m < nsend; m++) {
            n = avec->pack_border_vel(atom_sendnum[iswap][m],atom_sendlist[iswap][m],
                buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][m],0,world);
          }
        }
        if (sendself[iswap]) {
          avec->pack_border_vel(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],
              pbc[iswap][nsend]);
          avec->unpack_border_vel(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
              buf_send);
        }
        if (recvother[iswap]) {
          MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
          for (m = 0; m < nrecv; m++)
            avec->unpack_border_vel(atom_recvnum[iswap][m],atom_firstrecv[iswap][m],
                &buf_recv[atom_size_border*
                atom_forward_recv_offset[iswap][m]]);
        }

      } else {
        if (recvother[iswap]) {
          for (m = 0; m < nrecv; m++)
            MPI_Irecv(&buf_recv[atom_size_border*atom_forward_recv_offset[iswap][m]],
                atom_recvnum[iswap][m]*atom_size_border,
                MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        }

        if (sendother[iswap]) {
          for (m = 0; m < nsend; m++) {
            n = avec->pack_border(atom_sendnum[iswap][m],atom_sendlist[iswap][m],
                buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][m],0,world);
          }
        }
        if (sendself[iswap]) {
          avec->pack_border(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
          avec->unpack_border(atom_recvnum[iswap][nsend],atom_firstrecv[iswap][nsend],
              buf_send);
        }

        if (recvother[iswap]) {
          MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
          for (m = 0; m < nrecv; m++)
            avec->unpack_border(atom_recvnum[iswap][m],atom_firstrecv[iswap][m],
                &buf_recv[atom_size_border*
                atom_forward_recv_offset[iswap][m]]);
        }

      }

      // increment ghost atoms

      n = nrecvproc[iswap];
      if (n)
        atom->nghost += atom_forward_recv_offset[iswap][n-1] + atom_recvnum[iswap][n-1];
    }

    // insure send/recv buffers are long enough for all forward & reverse comm
    // send buf is for one forward or reverse sends to one proc
    // recv buf is for all forward or reverse recvs in one swap

    int max = MAX(atom_maxforward*smaxone,atom_maxreverse*rmaxone);
    if (max > maxsend) grow_send(max,0);
    max = MAX(atom_maxforward*rmaxall,atom_maxreverse*smaxall);
    if (max > maxrecv) grow_recv(max);

    // reset global->local map

    if (atom_map_style) atom->map_set();
  }


  // swap elements 

  if (element->nelements) {

    for (int iswap = 0; iswap < nswap; iswap++) {

      // find elements within sendboxes using >= and <
      // hi test with ">" is important b/c don't want to send an element 
      //   in lower dim (on boundary) that a proc will recv again in higher dim
      // for x-dim swaps, check owned elements
      // for yz-dim swaps, check owned and ghost elements
      // store sent element indices in elem_sendlist for use in future timesteps
      // NOTE: assume SINGLE mode, add logic for MULTI mode later

      x = element->x;

      if (iswap % 2 == 0) nlast = element->nlocal + element->nghost;

      ncountall = 0;

      for (m = 0; m < nsendproc[iswap]; m++) {

        bbox = sendbox[iswap][m];
        xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
        xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];

        ncount = 0;

        if (!bordergroup) {
          for (i = 0; i < nlast; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi &&
                x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == elem_maxsendlist[iswap][m]) grow_elem_list(iswap,m,ncount);
              elem_sendlist[iswap][m][ncount++] = i;
            }
          }

        } else {
          error->all(FLERR,"Border group has not been considered for elements yet");
          //ngroup = elem->nfirst;
          //for (i = 0; i < ngroup; i++) {
          //  if (x[i][0] >= xlo && x[i][0] < xhi &&
          //      x[i][1] >= ylo && x[i][1] < yhi &&
          //      x[i][2] >= zlo && x[i][2] < zhi) {
          //    if (ncount == elem_maxsendlist[iswap][m]) grow_elem_list(iswap,m,ncount);
          //    elem_sendlist[iswap][m][ncount++] = i;
          //  }
          //}
          //for (i = element->nlocal; i < nlast; i++) {
          //  if (x[i][0] >= xlo && x[i][0] < xhi &&
          //      x[i][1] >= ylo && x[i][1] < yhi &&
          //      x[i][2] >= zlo && x[i][2] < zhi) {
          //    if (ncount == elem_maxsendlist[iswap][m]) grow_elem_list(iswap,m,ncount);
          //    elem_sendlist[iswap][m][ncount++] = i;
          //  }
          //}
        }

        elem_sendnum[iswap][m] = ncount;
        smaxone = MAX(smaxone,ncount);
        ncountall += ncount;

      }

      smaxall = MAX(smaxall,ncountall);

      // send elem_sendnum counts to procs who recv from me except self
      // copy data to self if sendself is set

      nsend = nsendproc[iswap] - sendself[iswap];
      nrecv = nrecvproc[iswap] - sendself[iswap];

      if (recvother[iswap])
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&elem_recvnum[iswap][m],1,MPI_INT,
              recvproc[iswap][m],0,world,&requests[m]);
      if (sendother[iswap])
        for (m = 0; m < nsend; m++)
          MPI_Send(&elem_sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],0,world);
      if (sendself[iswap]) elem_recvnum[iswap][nrecv] = elem_sendnum[iswap][nsend];
      if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

      // setup other per swap/proc values from elem_sendnum and elem_recvnum

      for (m = 0; m < nsendproc[iswap]; m++) {
        elem_size_reverse_recv[iswap][m] = elem_sendnum[iswap][m]*elem_size_reverse;
        if (m == 0) elem_reverse_recv_offset[iswap][0] = 0;
        else elem_reverse_recv_offset[iswap][m] =
          elem_reverse_recv_offset[iswap][m-1] + elem_sendnum[iswap][m-1];
      }

      ncountall = 0;
      for (m = 0; m < nrecvproc[iswap]; m++) {
        ncount = elem_recvnum[iswap][m];
        rmaxone = MAX(rmaxone,ncount);
        ncountall += ncount;

        elem_size_forward_recv[iswap][m] = ncount*elem_size_forward;
        elem_size_reverse_send[iswap][m] = ncount*elem_size_reverse;
        if (m == 0) {
          elem_firstrecv[iswap][0] = element->nlocal + element->nghost;
          elem_forward_recv_offset[iswap][0] = 0;
        } else {
          elem_firstrecv[iswap][m] = elem_firstrecv[iswap][m-1] + elem_recvnum[iswap][m-1];
          elem_forward_recv_offset[iswap][m] =
            elem_forward_recv_offset[iswap][m-1] + elem_recvnum[iswap][m-1];
        }
      }
      rmaxall = MAX(rmaxall,ncountall);

      // insure send/recv buffers are large enough for this border comm swap

      if (smaxone*elem_size_border > maxsend) grow_send(smaxone*elem_size_border,0);
      if (rmaxall*elem_size_border > maxrecv) grow_recv(rmaxall*elem_size_border);

      // swap elements with other procs using pack_border(), unpack_border()
      // use Waitall() instead of Waitany() because calls to unpack_border()
      //   must increment per-element arrays in ascending order

      if (ghost_velocity) {

        error->all(FLERR,"Ghost velocity has not been considered yet");
        //if (recvother[iswap]) {
        //  for (m = 0; m < nrecv; m++)
        //    MPI_Irecv(&buf_recv[elem_size_border*elem_forward_recv_offset[iswap][m]],
        //        elem_recvnum[iswap][m]*elem_size_border,
        //        MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        //}
        //if (sendother[iswap]) {
        //  for (m = 0; m < nsend; m++) {
        //    n = avec->pack_border_vel(elem_sendnum[iswap][m],elem_sendlist[iswap][m],
        //        buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
        //    MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][m],0,world);
        //  }
        //}
        //if (sendself[iswap]) {
        //  avec->pack_border_vel(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
        //      buf_send,pbc_flag[iswap][nsend],
        //      pbc[iswap][nsend]);
        //  avec->unpack_border_vel(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
        //      buf_send);
        //}
        //if (recvother[iswap]) {
        //  MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        //  for (m = 0; m < nrecv; m++)
        //    avec->unpack_border_vel(elem_recvnum[iswap][m],elem_firstrecv[iswap][m],
        //        &buf_recv[elem_size_border*
        //        elem_forward_recv_offset[iswap][m]]);
        //}

      } else {

        if (recvother[iswap]) {
          for (m = 0; m < nrecv; m++)
            MPI_Irecv(&buf_recv[elem_size_border*elem_forward_recv_offset[iswap][m]],
                elem_recvnum[iswap][m]*elem_size_border,
                MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        }

        if (sendother[iswap]) {
          for (m = 0; m < nsend; m++) {
            n = evec->pack_border(elem_sendnum[iswap][m],elem_sendlist[iswap][m],
                buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][m],0,world);
          }
        }
        if (sendself[iswap]) {
          evec->pack_border(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
          evec->unpack_border(elem_recvnum[iswap][nsend],elem_firstrecv[iswap][nsend],
              buf_send);
        }
        if (recvother[iswap]) {
          MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
          for (m = 0; m < nrecv; m++)
            evec->unpack_border(elem_recvnum[iswap][m],elem_firstrecv[iswap][m],
                &buf_recv[elem_size_border*
                elem_forward_recv_offset[iswap][m]]);
        }
      }

      // increment ghost elements

      n = nrecvproc[iswap];
      if (n)
        element->nghost += elem_forward_recv_offset[iswap][n-1] + elem_recvnum[iswap][n-1];


    }

    // insure send/recv buffers are long enough for all forward & reverse comm
    // send buf is for one forward or reverse sends to one proc
    // recv buf is for all forward or reverse recvs in one swap

    int max = MAX(elem_maxforward*smaxone,elem_maxreverse*rmaxone);
    if (max > maxsend) grow_send(max,0);
    max = MAX(elem_maxforward*rmaxall,elem_maxreverse*smaxall);
    if (max > maxrecv) grow_recv(max);

    // reset global->local map

    if (elem_map_style) element->map_set();
  }

}
/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_pair(Pair *pair)
{
  int i,irecv,n,nsend,nrecv;

  int nsize;

  // forward communicate atoms from pair

  if (atom->natoms) {
    nsize = pair->comm_atom_forward;

    for (int iswap = 0; iswap < nswap; iswap++) {
      nsend = nsendproc[iswap] - sendself[iswap];
      nrecv = nrecvproc[iswap] - sendself[iswap];

      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[iswap][i]],
              nsize*atom_recvnum[iswap][i],
              MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = pair->pack_atom_forward_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        pair->pack_atom_forward_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
            buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        pair->unpack_atom_forward_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          pair->unpack_atom_forward_comm(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
              &buf_recv[nsize*
              atom_forward_recv_offset[iswap][irecv]]);
        }
      }
    }
  }

  // forward communicate elements from pair

  if (element->nelements) {
    nsize = pair->comm_elem_forward;

    for (int iswap = 0; iswap < nswap; iswap++) {
      nsend = nsendproc[iswap] - sendself[iswap];
      nrecv = nrecvproc[iswap] - sendself[iswap];

      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[nsize*elem_forward_recv_offset[iswap][i]],
              nsize*elem_recvnum[iswap][i],
              MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = pair->pack_elem_forward_comm(elem_sendnum[iswap][i],elem_sendlist[iswap][i],
              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        pair->pack_elem_forward_comm(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
            buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        pair->unpack_elem_forward_comm(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
            buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          pair->unpack_elem_forward_comm(elem_recvnum[iswap][irecv],elem_firstrecv[iswap][irecv],
              &buf_recv[nsize*elem_forward_recv_offset[iswap][irecv]]);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_pair(Pair *pair)
{
  int i,irecv,n,nsend,nrecv;

  int nsize;

  // reverse communicate atoms from pair

  if (atom->natoms) {
    nsize = MAX(pair->comm_atom_reverse,pair->comm_atom_reverse_off);

    if (nsize) 
      for (int iswap = nswap-1; iswap >= 0; iswap--) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++)
            MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[iswap][i]],
                nsize*atom_sendnum[iswap][i],MPI_DOUBLE,
                sendproc[iswap][i],0,world,&requests[i]);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            n = pair->pack_atom_reverse_comm(atom_recvnum[iswap][i],atom_firstrecv[iswap][i],buf_send);
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          pair->pack_atom_reverse_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],buf_send);
          pair->unpack_atom_reverse_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],buf_send);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
            pair->unpack_atom_reverse_comm(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
                &buf_recv[nsize*atom_reverse_recv_offset[iswap][irecv]]);
          }
        }
      }
  }

  // reverse communicate elements from pair

  if (element->nelements) {
    nsize = MAX(pair->comm_elem_reverse,pair->comm_elem_reverse_off);

    if (nsize) 
      for (int iswap = nswap-1; iswap >= 0; iswap--) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++)
            MPI_Irecv(&buf_recv[nsize*elem_reverse_recv_offset[iswap][i]],
                nsize*elem_sendnum[iswap][i],MPI_DOUBLE,
                sendproc[iswap][i],0,world,&requests[i]);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            n = pair->pack_elem_reverse_comm(elem_recvnum[iswap][i],elem_firstrecv[iswap][i],buf_send);
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          pair->pack_elem_reverse_comm(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
              buf_send);
          pair->unpack_elem_reverse_comm(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],buf_send);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
            pair->unpack_elem_reverse_comm(elem_sendnum[iswap][irecv],elem_sendlist[iswap][irecv],
                &buf_recv[nsize*elem_reverse_recv_offset[iswap][irecv]]);
          }
        }
      }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_atom_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_atom_forward
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_fix(Fix *fix, int size)
{
  error->all(FLERR,"Forward comm fix not yet supported by CommTiled");
  //  int i,irecv,n,nsize,nsend,nrecv;
  //
  //  if (size) nsize = size;
  //  else nsize = fix->comm_atom_forward;
  //
  //  for (int iswap = 0; iswap < nswap; iswap++) {
  //    nsend = nsendproc[iswap] - sendself[iswap];
  //    nrecv = nrecvproc[iswap] - sendself[iswap];
  //
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[iswap][i]],
  //            nsize*atom_recvnum[iswap][i],
  //            MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
  //    }
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++) {
  //        n = fix->pack_atom_forward_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
  //            buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
  //      }
  //    }
  //    if (sendself[iswap]) {
  //      fix->pack_atom_forward_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
  //          buf_send,pbc_flag[iswap][nsend],
  //          pbc[iswap][nsend]);
  //      fix->unpack_atom_forward_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
  //          buf_send);
  //    }
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++) {
  //        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
  //        fix->unpack_atom_forward_comm(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
  //            &buf_recv[nsize*
  //            atom_forward_recv_offset[iswap][irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_atom_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_atom_forward
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix(Fix *fix, int size)
{
  error->all(FLERR,"Reverse comm fix not yet supported by CommTiled");
  //  int i,irecv,n,nsize,nsend,nrecv;
  //
  //  if (size) nsize = size;
  //  else nsize = fix->comm_atom_reverse;
  //
  //  for (int iswap = nswap-1; iswap >= 0; iswap--) {
  //    nsend = nsendproc[iswap] - sendself[iswap];
  //    nrecv = nrecvproc[iswap] - sendself[iswap];
  //
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[iswap][i]],
  //            nsize*atom_sendnum[iswap][i],MPI_DOUBLE,
  //            sendproc[iswap][i],0,world,&requests[i]);
  //    }
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++) {
  //        n = fix->pack_atom_reverse_comm(atom_recvnum[iswap][i],atom_firstrecv[iswap][i],
  //            buf_send);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
  //      }
  //    }
  //    if (sendself[iswap]) {
  //      fix->pack_atom_reverse_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
  //          buf_send);
  //      fix->unpack_atom_reverse_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
  //          buf_send);
  //    }
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++) {
  //        MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
  //        fix->unpack_atom_reverse_comm(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
  //            &buf_recv[nsize*
  //            atom_reverse_recv_offset[iswap][irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for all pack sizes to insure buf_send is big enough
   handshake sizes before irregular comm to insure buf_recv is big enough
NOTE: how to setup one big buf recv with correct offsets ??
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix_variable(Fix *fix)
{
  error->all(FLERR,"Reverse comm fix variable not yet supported by CommTiled");
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_compute(Compute *compute)
{
  int i,irecv,n,nsend,nrecv;

  int nsize;

  if (atom->natoms) {
    nsize = compute->comm_atom_forward;
    if (nsize) 
      for (int iswap = 0; iswap < nswap; iswap++) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++)
            MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[iswap][i]],
                nsize*atom_recvnum[iswap][i],
                MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            n = compute->pack_atom_forward_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
                buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          compute->pack_atom_forward_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],
              pbc[iswap][nsend]);
          compute->unpack_atom_forward_comm(atom_recvnum[iswap][nrecv],
              atom_firstrecv[iswap][nrecv],buf_send);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
            compute->unpack_atom_forward_comm(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
                &buf_recv[nsize*atom_forward_recv_offset[iswap][irecv]]);
          }
        }
      }
  }

  if (element->nelements) {
    nsize = compute->comm_elem_forward;
    if (nsize) 
      for (int iswap = 0; iswap < nswap; iswap++) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++)
            MPI_Irecv(&buf_recv[nsize*elem_forward_recv_offset[iswap][i]],
                nsize*elem_recvnum[iswap][i],
                MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            n = compute->pack_elem_forward_comm(elem_sendnum[iswap][i],elem_sendlist[iswap][i],
                buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          compute->pack_elem_forward_comm(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],
              buf_send,pbc_flag[iswap][nsend],
              pbc[iswap][nsend]);
          compute->unpack_elem_forward_comm(elem_recvnum[iswap][nrecv],
              elem_firstrecv[iswap][nrecv],buf_send);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
            compute->unpack_elem_forward_comm(elem_recvnum[iswap][irecv],elem_firstrecv[iswap][irecv],
                &buf_recv[nsize*elem_forward_recv_offset[iswap][irecv]]);
          }
        }
      }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_compute(Compute *compute)
{
  int i,irecv,n,nsend,nrecv;

  int nsize;

  // reverse communicate atoms from compute

  if (atom->natoms) {
    nsize = compute->comm_atom_reverse;

    if (nsize) 
      for (int iswap = nswap-1; iswap >= 0; iswap--) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++)
            MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[iswap][i]],
                nsize*atom_sendnum[iswap][i],MPI_DOUBLE,
                sendproc[iswap][i],0,world,&requests[i]);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            n = compute->pack_atom_reverse_comm(atom_recvnum[iswap][i],atom_firstrecv[iswap][i],buf_send);
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          compute->pack_atom_reverse_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],buf_send);
          compute->unpack_atom_reverse_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],buf_send);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
            compute->unpack_atom_reverse_comm(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
                &buf_recv[nsize*atom_reverse_recv_offset[iswap][irecv]]);
          }
        }
      }
  }

  // reverse communicate elements from compute

  if (element->nelements) {
    nsize = compute->comm_elem_reverse;

    if (nsize) 
      for (int iswap = nswap-1; iswap >= 0; iswap--) {
        nsend = nsendproc[iswap] - sendself[iswap];
        nrecv = nrecvproc[iswap] - sendself[iswap];

        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++)
            MPI_Irecv(&buf_recv[nsize*elem_reverse_recv_offset[iswap][i]],
                nsize*elem_sendnum[iswap][i],MPI_DOUBLE,
                sendproc[iswap][i],0,world,&requests[i]);
        }
        if (recvother[iswap]) {
          for (i = 0; i < nrecv; i++) {
            n = compute->pack_elem_reverse_comm(elem_recvnum[iswap][i],elem_firstrecv[iswap][i],buf_send);
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
          }
        }
        if (sendself[iswap]) {
          compute->pack_elem_reverse_comm(elem_recvnum[iswap][nrecv],elem_firstrecv[iswap][nrecv],
              buf_send);
          compute->unpack_elem_reverse_comm(elem_sendnum[iswap][nsend],elem_sendlist[iswap][nsend],buf_send);
        }
        if (sendother[iswap]) {
          for (i = 0; i < nsend; i++) {
            MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
            compute->unpack_elem_reverse_comm(elem_sendnum[iswap][irecv],elem_sendlist[iswap][irecv],
                &buf_recv[nsize*elem_reverse_recv_offset[iswap][irecv]]);
          }
        }
      }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_dump(Dump *dump)
{
  error->all(FLERR,"Forward comm dump not yet supported by CommTiled");
  //  int i,irecv,n,nsend,nrecv;
  //
  //  int nsize = dump->comm_atom_forward;
  //
  //  for (int iswap = 0; iswap < nswap; iswap++) {
  //    nsend = nsendproc[iswap] - sendself[iswap];
  //    nrecv = nrecvproc[iswap] - sendself[iswap];
  //
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[iswap][i]],
  //            nsize*atom_recvnum[iswap][i],
  //            MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
  //    }
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++) {
  //        n = dump->pack_atom_forward_comm(atom_sendnum[iswap][i],atom_sendlist[iswap][i],
  //            buf_send,pbc_flag[iswap][i],
  //            pbc[iswap][i]);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
  //      }
  //    }
  //    if (sendself[iswap]) {
  //      dump->pack_atom_forward_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
  //          buf_send,pbc_flag[iswap][nsend],
  //          pbc[iswap][nsend]);
  //      dump->unpack_atom_forward_comm(atom_recvnum[iswap][nrecv],
  //          atom_firstrecv[iswap][nrecv],buf_send);
  //    }
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++) {
  //        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
  //        dump->unpack_atom_forward_comm(atom_recvnum[iswap][irecv],atom_firstrecv[iswap][irecv],
  //            &buf_recv[nsize*
  //            atom_forward_recv_offset[iswap][irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_dump(Dump *dump)
{
  error->all(FLERR,"Forward comm dump not yet supported by CommTiled");
  //  int i,irecv,n,nsend,nrecv;
  //
  //  int nsize = dump->comm_atom_reverse;
  //
  //  for (int iswap = nswap-1; iswap >= 0; iswap--) {
  //    nsend = nsendproc[iswap] - sendself[iswap];
  //    nrecv = nrecvproc[iswap] - sendself[iswap];
  //
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[iswap][i]],
  //            nsize*atom_sendnum[iswap][i],MPI_DOUBLE,
  //            sendproc[iswap][i],0,world,&requests[i]);
  //    }
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++) {
  //        n = dump->pack_atom_reverse_comm(atom_recvnum[iswap][i],atom_firstrecv[iswap][i],
  //            buf_send);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
  //      }
  //    }
  //    if (sendself[iswap]) {
  //      dump->pack_atom_reverse_comm(atom_recvnum[iswap][nrecv],atom_firstrecv[iswap][nrecv],
  //          buf_send);
  //      dump->unpack_atom_reverse_comm(atom_sendnum[iswap][nsend],atom_sendlist[iswap][nsend],
  //          buf_send);
  //    }
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++) {
  //        MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
  //        dump->unpack_atom_reverse_comm(atom_sendnum[iswap][irecv],atom_sendlist[iswap][irecv],
  //            &buf_recv[nsize*
  //            atom_reverse_recv_offset[iswap][irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   forward communication of Nsize values in per-atom array
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_array(int nsize, double **array)
{
  error->all(FLERR,"Forward comm array not yet supported by CommTiled");
  //  int i,j,k,m,iatom,last,irecv,nsend,nrecv;
  //
  //  // insure send/recv bufs are big enough for nsize
  //  // based on smaxone/rmaxall from most recent borders() invocation
  //
  //  if (nsize > atom_maxforward) {
  //    atom_maxforward = nsize;
  //    if (atom_maxforward*smaxone > maxsend) grow_send(atom_maxforward*smaxone,0);
  //    if (atom_maxforward*rmaxall > maxrecv) grow_recv(atom_maxforward*rmaxall);
  //  }
  //
  //  for (int iswap = 0; iswap < nswap; iswap++) {
  //    nsend = nsendproc[iswap] - sendself[iswap];
  //    nrecv = nrecvproc[iswap] - sendself[iswap];
  //
  //    MPI_Barrier(world);
  //
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[iswap][i]],
  //            nsize*atom_recvnum[iswap][i],
  //            MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
  //    }
  //    if (sendother[iswap]) {
  //      for (i = 0; i < nsend; i++) {
  //        m = 0;
  //        for (iatom = 0; iatom < atom_sendnum[iswap][i]; iatom++) {
  //          j = atom_sendlist[iswap][i][iatom];
  //          for (k = 0; k < nsize; k++)
  //            buf_send[m++] = array[j][k];
  //        }
  //        MPI_Send(buf_send,nsize*atom_sendnum[iswap][i],
  //            MPI_DOUBLE,sendproc[iswap][i],0,world);
  //      }
  //    }
  //    if (sendself[iswap]) {
  //      m = 0;
  //      for (iatom = 0; iatom < atom_sendnum[iswap][nsend]; iatom++) {
  //        j = atom_sendlist[iswap][nsend][iatom];
  //        for (k = 0; k < nsize; k++)
  //          buf_send[m++] = array[j][k];
  //      }
  //      m = 0;
  //      last = atom_firstrecv[iswap][nrecv] + atom_recvnum[iswap][nrecv];
  //      for (iatom = atom_firstrecv[iswap][nrecv]; iatom < last; iatom++)
  //        for (k = 0; k < nsize; k++)
  //          array[iatom][k] = buf_send[m++];
  //    }
  //
  //    if (recvother[iswap]) {
  //      for (i = 0; i < nrecv; i++) {
  //        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
  //        m = nsize*atom_forward_recv_offset[iswap][irecv];
  //        last = atom_firstrecv[iswap][irecv] + atom_recvnum[iswap][irecv];
  //        for (iatom = atom_firstrecv[iswap][irecv]; iatom < last; iatom++)
  //          for (k = 0; k < nsize; k++)
  //            array[iatom][k] = buf_recv[m++];
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
NOTE: this method is currently not used
------------------------------------------------------------------------- */

int CommTiled::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  int nrecv = n;
  return nrecv;
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   box is owned by me and extends in dim
   ------------------------------------------------------------------------- */

void CommTiled::box_drop_brick(int idim, double *lo, double *hi, int &indexme)
{
  // NOTE: this is not triclinic compatible
  // NOTE: these error messages are internal sanity checks
  //       should not occur, can be removed at some point

  int index=-1,dir;
  if (hi[idim] == sublo[idim]) {
    index = myloc[idim] - 1;
    dir = -1;
  } else if (lo[idim] == subhi[idim]) {
    index = myloc[idim] + 1;
    dir = 1;
  } else if (hi[idim] == boxhi[idim]) {
    index = procgrid[idim] - 1;
    dir = -1;
  } else if (lo[idim] == boxlo[idim]) {
    index = 0;
    dir = 1;
  } else error->one(FLERR,"Comm tiled mis-match in box drop brick");

  int other1,other2,proc;
  double lower,upper;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  if (index < 0 || index > procgrid[idim])
    error->one(FLERR,"Comm tiled invalid index in box drop brick");

  while (1) {
    lower = boxlo[idim] + prd[idim]*split[index];
    if (index < procgrid[idim]-1)
      upper = boxlo[idim] + prd[idim]*split[index+1];
    else upper = boxhi[idim];
    if (lower >= hi[idim] || upper <= lo[idim]) break;

    if (idim == 0) proc = grid2proc[index][other1][other2];
    else if (idim == 1) proc = grid2proc[other1][index][other2];
    else proc = grid2proc[other1][other2][index];

    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap,maxoverlap,"comm:overlap");
    }

    if (proc == me) indexme = noverlap;
    overlap[noverlap++] = proc;
    index += dir;
    if (index < 0 || index >= procgrid[idim]) break;
  }
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
   ------------------------------------------------------------------------- */

void CommTiled::box_drop_tiled(int idim, double *lo, double *hi, int &indexme)
{
  box_drop_tiled_recurse(lo,hi,0,nprocs-1,indexme);
}

void CommTiled::box_drop_tiled_recurse(double *lo, double *hi,
    int proclower, int procupper,
    int &indexme)
{
  // end recursion when partition is a single proc
  // add proc to overlap list

  if (proclower == procupper) {
    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap,maxoverlap,"comm:overlap");
    }

    if (proclower == me) indexme = noverlap;
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use > and < criteria so does not include a box it only touches
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim]*rcbinfo[procmid].cutfrac;

  if (lo[idim] < cut)
    box_drop_tiled_recurse(lo,hi,proclower,procmid-1,indexme);
  if (hi[idim] > cut)
    box_drop_tiled_recurse(lo,hi,procmid,procupper,indexme);
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
   ------------------------------------------------------------------------- */

void CommTiled::box_other_brick(int idim, int idir,
    int proc, double *lo, double *hi)
{
  lo[0] = sublo[0]; lo[1] = sublo[1]; lo[2] = sublo[2];
  hi[0] = subhi[0]; hi[1] = subhi[1]; hi[2] = subhi[2];

  int other1,other2,oproc;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  int dir = -1;
  if (idir) dir = 1;
  int index = myloc[idim];
  int n = procgrid[idim];

  for (int i = 0; i < n; i++) {
    index += dir;
    if (index < 0) index = n-1;
    else if (index >= n) index = 0;

    if (idim == 0) oproc = grid2proc[index][other1][other2];
    else if (idim == 1) oproc = grid2proc[other1][index][other2];
    else oproc = grid2proc[other1][other2][index];

    if (proc == oproc) {
      lo[idim] = boxlo[idim] + prd[idim]*split[index];
      if (split[index+1] < 1.0)
        hi[idim] = boxlo[idim] + prd[idim]*split[index+1];
      else hi[idim] = boxhi[idim];
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
   ------------------------------------------------------------------------- */

void CommTiled::box_other_tiled(int idim, int idir,
    int proc, double *lo, double *hi)
{
  double (*split)[2] = rcbinfo[proc].mysplit;

  lo[0] = boxlo[0] + prd[0]*split[0][0];
  if (split[0][1] < 1.0) hi[0] = boxlo[0] + prd[0]*split[0][1];
  else hi[0] = boxhi[0];

  lo[1] = boxlo[1] + prd[1]*split[1][0];
  if (split[1][1] < 1.0) hi[1] = boxlo[1] + prd[1]*split[1][1];
  else hi[1] = boxhi[1];

  lo[2] = boxlo[2] + prd[2]*split[2][0];
  if (split[2][1] < 1.0) hi[2] = boxlo[2] + prd[2]*split[2][1];
  else hi[2] = boxhi[2];
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
   procneigh stores 6 procs that touch me
   ------------------------------------------------------------------------- */

int CommTiled::box_touch_brick(int proc, int idim, int idir)
{
  if (procneigh[idim][idir] == proc) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
   ------------------------------------------------------------------------- */

int CommTiled::box_touch_tiled(int proc, int idim, int idir)
{
  // sending to left
  // only touches if proc hi = my lo, or if proc hi = boxhi and my lo = boxlo

  if (idir == 0) {
    if (rcbinfo[proc].mysplit[idim][1] == rcbinfo[me].mysplit[idim][0])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][1] == 1.0 &&
        rcbinfo[me].mysplit[idim][0] == 0.0)
      return 1;

    // sending to right
    // only touches if proc lo = my hi, or if proc lo = boxlo and my hi = boxhi

  } else {
    if (rcbinfo[proc].mysplit[idim][0] == rcbinfo[me].mysplit[idim][1])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][0] == 0.0 &&
        rcbinfo[me].mysplit[idim][1] == 1.0)
      return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_brick(int idim, double *x)
{
  if (closer_subbox_edge(idim,x)) return procneigh[idim][1];
  return procneigh[idim][0];
}

/* ----------------------------------------------------------------------
   determine which proc owns point x via recursion thru RCB tree
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_tiled(int idim, double *x)
{
  double xnew[3];
  xnew[0] = x[0]; xnew[1] = x[1]; xnew[2] = x[2];

  if (idim == 0) {
    if (xnew[1] < sublo[1] || xnew[1] > subhi[1]) {
      if (closer_subbox_edge(1,x)) xnew[1] = subhi[1];
      else xnew[1] = sublo[1];
    }
  }
  if (idim <= 1) {
    if (xnew[2] < sublo[2] || xnew[2] > subhi[2]) {
      if (closer_subbox_edge(2,x)) xnew[2] = subhi[2];
      else xnew[2] = sublo[2];
    }
  }

  int proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
  if (proc == me) return me;

  if (idim == 0) {
    int done = 1;
    if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
      xnew[1] -= EPSILON * (subhi[1]-sublo[1]);
      done = 0;
    }
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2]-sublo[2]);
      done = 0;
    }
    if (!done) {
      proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
      done = 1;
      if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
        xnew[1] -= EPSILON * (subhi[1]-sublo[1]);
        done = 0;
      }
      if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
        xnew[2] -= EPSILON * (subhi[2]-sublo[2]);
        done = 0;
      }
      if (!done) proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
    }
  } else if (idim == 1) {
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2]-sublo[2]);
      proc = point_drop_tiled_recurse(xnew,0,nprocs-1);
    }
  }

  return proc;
}

/* ----------------------------------------------------------------------
   recursive point drop thru RCB tree
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_tiled_recurse(double *x,
    int proclower, int procupper)
{
  // end recursion when partition is a single proc
  // return proc

  if (proclower == procupper) return proclower;

  // drop point on side of cut it is on
  // use < criterion so point is not on high edge of proc sub-domain
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim]*rcbinfo[procmid].cutfrac;

  if (x[idim] < cut) return point_drop_tiled_recurse(x,proclower,procmid-1);
  else return point_drop_tiled_recurse(x,procmid,procupper);
}

/* ----------------------------------------------------------------------
   assume x[idim] is outside subbox bounds in same dim
   ------------------------------------------------------------------------- */

int CommTiled::closer_subbox_edge(int idim, double *x)
{
  double deltalo,deltahi;

  if (sublo[idim] == boxlo[idim])
    deltalo = fabs(x[idim]-prd[idim] - sublo[idim]);
  else deltalo = fabs(x[idim] - sublo[idim]);

  if (subhi[idim] == boxhi[idim])
    deltahi = fabs(x[idim]+prd[idim] - subhi[idim]);
  else deltahi = fabs(x[idim] - subhi[idim]);

  if (deltalo < deltahi) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   if RCB decomp exists and just changed, gather needed global RCB info
   ------------------------------------------------------------------------- */

void CommTiled::coord2proc_setup()
{
  if (!rcbnew) return;

  if (!rcbinfo)
    rcbinfo = (RCBinfo *)
      memory->smalloc(nprocs*sizeof(RCBinfo),"comm:rcbinfo");
  rcbnew = 0;
  RCBinfo rcbone;
  memcpy(&rcbone.mysplit[0][0],&mysplit[0][0],6*sizeof(double));
  rcbone.cutfrac = rcbcutfrac;
  rcbone.dim = rcbcutdim;
  MPI_Allgather(&rcbone,sizeof(RCBinfo),MPI_CHAR,
      rcbinfo,sizeof(RCBinfo),MPI_CHAR,world);

}

/* ----------------------------------------------------------------------
   determine which proc owns atom with coord x[3] based on current decomp
   x will be in box (orthogonal) or lamda coords (triclinic)
   if layout = UNIFORM or NONUNIFORM, invoke parent method
   if layout = TILED, use point_drop_recurse()
   return owning proc ID, ignore igx,igy,igz
   ------------------------------------------------------------------------- */

int CommTiled::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  if (layout != LAYOUT_TILED) return Comm::coord2proc(x,igx,igy,igz);
  return point_drop_tiled_recurse(x,0,nprocs-1);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
   ------------------------------------------------------------------------- */

void CommTiled::grow_send(int n, int flag)
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

void CommTiled::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap atom_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_atom_list(int iswap, int iwhich, int n)
{
  atom_maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(atom_sendlist[iswap][iwhich],atom_maxsendlist[iswap][iwhich],
      "comm:atom_sendlist[i][j]");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap atom_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_elem_list(int iswap, int iwhich, int n)
{
  elem_maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(elem_sendlist[iswap][iwhich],elem_maxsendlist[iswap][iwhich],
      "comm:elem_sendlist[i][j]");
}


/* ----------------------------------------------------------------------
   allocation of swap info
   ------------------------------------------------------------------------- */

void CommTiled::allocate_swap(int n)
{
  nsendproc = new int[n];
  nrecvproc = new int[n];
  sendother = new int[n];
  recvother = new int[n];
  sendself = new int[n];
  nprocmax = new int[n];

  sendproc = new int*[n];
  recvproc = new int*[n];
  atom_sendnum = new int*[n];
  atom_recvnum = new int*[n];
  atom_size_forward_recv = new int*[n];
  atom_firstrecv = new int*[n];
  atom_size_reverse_send = new int*[n];
  atom_size_reverse_recv = new int*[n];
  atom_forward_recv_offset = new int*[n];
  atom_reverse_recv_offset = new int*[n];

  elem_sendnum = new int*[n];
  elem_recvnum = new int*[n];
  elem_size_forward_recv = new int*[n];
  elem_firstrecv = new int*[n];
  elem_size_reverse_send = new int*[n];
  elem_size_reverse_recv = new int*[n];
  elem_forward_recv_offset = new int*[n];
  elem_reverse_recv_offset = new int*[n];

  pbc_flag = new int*[n];
  pbc = new int**[n];
  sendbox = new double**[n];
  atom_maxsendlist = new int*[n];
  atom_sendlist = new int**[n];

  elem_maxsendlist = new int*[n];
  elem_sendlist = new int**[n];

  for (int i = 0; i < n; i++) {
    sendproc[i] = recvproc[i] = NULL;
    atom_sendnum[i] = atom_recvnum[i] = NULL;
    atom_size_forward_recv[i] = atom_firstrecv[i] = NULL;
    atom_size_reverse_send[i] = atom_size_reverse_recv[i] = NULL;
    atom_forward_recv_offset[i] = atom_reverse_recv_offset[i] = NULL;
    elem_sendnum[i] = elem_recvnum[i] = NULL;
    elem_size_forward_recv[i] = elem_firstrecv[i] = NULL;
    elem_size_reverse_send[i] = elem_size_reverse_recv[i] = NULL;
    elem_forward_recv_offset[i] = elem_reverse_recv_offset[i] = NULL;


    pbc_flag[i] = NULL;
    pbc[i] = NULL;
    sendbox[i] = NULL;
    atom_maxsendlist[i] = NULL;
    atom_sendlist[i] = NULL;
    elem_maxsendlist[i] = NULL;
    elem_sendlist[i] = NULL;

  }

  maxreqstat = 0;
  requests = NULL;

  for (int i = 0; i < n; i++) {
    nprocmax[i] = DELTA_PROCS;
    grow_swap_send(i,DELTA_PROCS,0);
    grow_swap_recv(i,DELTA_PROCS);
  }

  nexchproc = new int[n/2];
  nexchprocmax = new int[n/2];
  exchproc = new int*[n/2];
  exchnum = new int*[n/2];

  for (int i = 0; i < n/2; i++) {
    nexchprocmax[i] = DELTA_PROCS;
    exchproc[i] = new int[DELTA_PROCS];
    exchnum[i] = new int[DELTA_PROCS];
  }
}

/* ----------------------------------------------------------------------
   grow info for swap I, to allow for N procs to communicate with
   ditto for complementary recv for swap I+1 or I-1, as invoked by caller
   ------------------------------------------------------------------------- */

void CommTiled::grow_swap_send(int i, int n, int nold)
{
  delete [] sendproc[i];
  sendproc[i] = new int[n];

  delete [] atom_sendnum[i];
  delete [] atom_size_reverse_recv[i];
  delete [] atom_reverse_recv_offset[i];
  delete [] elem_sendnum[i];
  delete [] elem_size_reverse_recv[i];
  delete [] elem_reverse_recv_offset[i];

  atom_sendnum[i] = new int[n];
  atom_size_reverse_recv[i] = new int[n];
  atom_reverse_recv_offset[i] = new int[n];
  elem_sendnum[i] = new int[n];
  elem_size_reverse_recv[i] = new int[n];
  elem_reverse_recv_offset[i] = new int[n];

  delete [] pbc_flag[i];
  pbc_flag[i] = new int[n];
  memory->destroy(pbc[i]);
  memory->create(pbc[i],n,6,"comm:pbc_flag");
  memory->destroy(sendbox[i]);
  memory->create(sendbox[i],n,6,"comm:sendbox");

  delete [] atom_maxsendlist[i];
  delete [] elem_maxsendlist[i];
  atom_maxsendlist[i] = new int[n];
  elem_maxsendlist[i] = new int[n];

  for (int j = 0; j < nold; j++) {
    memory->destroy(atom_sendlist[i][j]);
    memory->destroy(elem_sendlist[i][j]);
  }

  delete [] atom_sendlist[i];
  delete [] elem_sendlist[i];
  atom_sendlist[i] = new int*[n];
  elem_sendlist[i] = new int*[n];

  for (int j = 0; j < n; j++) {
    atom_maxsendlist[i][j] = BUFMIN;
    elem_maxsendlist[i][j] = BUFMIN;
    memory->create(atom_sendlist[i][j],BUFMIN,"comm:atom_sendlist[i][j]");
    memory->create(elem_sendlist[i][j],BUFMIN,"comm:atom_sendlist[i][j]");
  }

}

void CommTiled::grow_swap_recv(int i, int n)
{
  delete [] recvproc[i];
  recvproc[i] = new int[n];
  delete [] atom_recvnum[i];
  delete [] elem_recvnum[i];
  atom_recvnum[i] = new int[n];
  elem_recvnum[i] = new int[n];

  delete [] atom_size_forward_recv[i];
  delete [] elem_size_forward_recv[i];
  atom_size_forward_recv[i] = new int[n];
  elem_size_forward_recv[i] = new int[n];
  delete [] atom_firstrecv[i];
  delete [] atom_forward_recv_offset[i];
  delete [] atom_size_reverse_send[i];
  delete [] elem_firstrecv[i];
  delete [] elem_forward_recv_offset[i];
  delete [] elem_size_reverse_send[i];
  atom_firstrecv[i] = new int[n];
  atom_forward_recv_offset[i] = new int[n];
  atom_size_reverse_send[i] = new int[n];
  elem_firstrecv[i] = new int[n];
  elem_forward_recv_offset[i] = new int[n];
  elem_size_reverse_send[i] = new int[n];
}

/* ----------------------------------------------------------------------
   deallocate swap info
   ------------------------------------------------------------------------- */

void CommTiled::deallocate_swap(int n)
{
  delete [] nsendproc;
  delete [] nrecvproc;
  delete [] sendother;
  delete [] recvother;
  delete [] sendself;

  for (int i = 0; i < n; i++) {
    delete [] sendproc[i];
    delete [] recvproc[i];

    delete [] atom_sendnum[i];
    delete [] atom_recvnum[i];
    delete [] atom_size_forward_recv[i];
    delete [] atom_firstrecv[i];
    delete [] atom_size_reverse_send[i];
    delete [] atom_size_reverse_recv[i];
    delete [] atom_forward_recv_offset[i];
    delete [] atom_reverse_recv_offset[i];
    delete [] elem_sendnum[i];
    delete [] elem_recvnum[i];
    delete [] elem_size_forward_recv[i];
    delete [] elem_firstrecv[i];
    delete [] elem_size_reverse_send[i];
    delete [] elem_size_reverse_recv[i];
    delete [] elem_forward_recv_offset[i];
    delete [] elem_reverse_recv_offset[i];

    delete [] pbc_flag[i];
    memory->destroy(pbc[i]);
    memory->destroy(sendbox[i]);
    delete [] atom_maxsendlist[i];
    delete [] elem_maxsendlist[i];

    for (int j = 0; j < nprocmax[i]; j++) {
      memory->destroy(atom_sendlist[i][j]);
      memory->destroy(elem_sendlist[i][j]);
    }
    delete [] atom_sendlist[i];
    delete [] elem_sendlist[i];
  }

  delete [] sendproc;
  delete [] recvproc;

  delete [] atom_sendnum;
  delete [] atom_recvnum;
  delete [] atom_size_forward_recv;
  delete [] atom_firstrecv;
  delete [] atom_size_reverse_send;
  delete [] atom_size_reverse_recv;
  delete [] atom_forward_recv_offset;
  delete [] atom_reverse_recv_offset;

  delete [] elem_sendnum;
  delete [] elem_recvnum;
  delete [] elem_size_forward_recv;
  delete [] elem_firstrecv;
  delete [] elem_size_reverse_send;
  delete [] elem_size_reverse_recv;
  delete [] elem_forward_recv_offset;
  delete [] elem_reverse_recv_offset;

  delete [] pbc_flag;
  delete [] pbc;
  delete [] sendbox;
  delete [] atom_maxsendlist;
  delete [] atom_sendlist;
  delete [] elem_maxsendlist;
  delete [] elem_sendlist;

  delete [] requests;

  delete [] nprocmax;

  delete [] nexchproc;
  delete [] nexchprocmax;

  for (int i = 0; i < n/2; i++) {
    delete [] exchproc[i];
    delete [] exchnum[i];
  }

  delete [] exchproc;
  delete [] exchnum;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint CommTiled::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
