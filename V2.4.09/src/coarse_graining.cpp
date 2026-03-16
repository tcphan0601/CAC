#include <mpi.h>
#include <string.h>
#include "coarse_graining.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "group.h"
#include "compute.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "irregular_comm.h"
#include "universe.h"
#include "math_extra.h"

#include <map>

#define EPSILON 1e-6
#define MAXLINE 256
using namespace CAC_NS;
using namespace MathExtra;

enum{CNA,IDS};

// same as compute_cna_atom.cpp and compute_ids_atom.cpp
enum{CNA_UNKNOWN,FCC,HCP,BCC,ICOS,CNA_OTHER};
enum{IDS_OTHER,CUBIC,HEX,IDS_UNKNOWN};

/* ---------------------------------------------------------------------- */

CoarseGraining::CoarseGraining(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
}

/* ----------------------------------------------------------------------
   Called as coarse_graining command in input script.
   New atoms created will have ID starting immediately 
   following the largest atom ID existing before the
   coarse_graining command was evoked.
   Call IrregularComm to migrate atoms since atoms might move far away
   ------------------------------------------------------------------------- */

void CoarseGraining::command(int narg, char **arg)
{
  tol = 1e-3;
  if (domain->box_exist == 0)
    error->all(FLERR,"Coarse_graining command before simulation box is defined");
  if (domain->triclinic)
    error->all(FLERR,"Coarse_graining command doesn't work with triclinic box yet");
  if (atom->natoms == 0) {
    if (me == 0) error->warning(FLERR,"No atom in simulation box");
    return;
  }

  atom->init(); 
  element->init(); 

  avec = atom->avec;
  evec = element->evec;

  if (comm->nprocs != 1) error->all(FLERR,"Coarse_graining command has not been parallelized yet, use one proc only");
  if (narg < 6) error->all(FLERR,"Illegal coarse_graining command");

  if (me == 0) {
    if (screen) fprintf(screen,"Coarse graining atoms into element ...\n");
    if (logfile) fprintf(logfile,"Coarse graining atoms into element ...\n");
  }

  if (strcmp(arg[2],"Ac/fcc") == 0) 
    cut = 4.53237;
  else if (strcmp(arg[2],"Ag/fcc") == 0) 
    cut = 3.49103;
  else if (strcmp(arg[2],"Al/fcc") == 0) 
    cut = 3.45689;
  else if (strcmp(arg[2],"Ar/fcc") == 0) 
    cut = 4.48969;
  else if (strcmp(arg[2],"Au/fcc") == 0) 
    cut = 3.48250;
  else if (strcmp(arg[2],"Ba/bcc") == 0) 
    cut = 6.05968;
  else if (strcmp(arg[2],"C/dia") == 0) 
    cut = 2.74223;
  else if (strcmp(arg[2],"Ca/fcc") == 0) 
    cut = 4.76283;
  else if (strcmp(arg[2],"Ce/fcc") == 0) 
    cut = 4.40434;
  else if (strcmp(arg[2],"Cr/bcc") == 0) 
    cut = 3.47647;
  else if (strcmp(arg[2],"Cs/bcc") == 0) 
    cut = 7.30300;
  else if (strcmp(arg[2],"Cu/fcc") == 0) 
    cut = 3.08133;
  else if (strcmp(arg[2],"Eu/bcc") == 0) 
    cut = 5.56476;
  else if (strcmp(arg[2],"Fe/bcc") == 0) 
    cut = 3.46440;
  else if (strcmp(arg[2],"Ge/dia") == 0) 
    cut = 4.34762;
  else if (strcmp(arg[2],"Ir/fcc") == 0) 
    cut = 3.27764;
  else if (strcmp(arg[2],"K/bcc") == 0) 
    cut = 6.31317;
  else if (strcmp(arg[2],"Kr/fcc") == 0) 
    cut = 4.88233;
  else if (strcmp(arg[2],"Li/bcc") == 0) 
    cut = 4.21280;
  else if (strcmp(arg[2],"Mo/bcc") == 0) 
    cut = 3.80239;
  else if (strcmp(arg[2],"Na/bcc") == 0) 
    cut = 5.10606;
  else if (strcmp(arg[2],"Nb/bcc") == 0) 
    cut = 3.98345;
  else if (strcmp(arg[2],"Ne/fcc") == 0) 
    cut = 3.78124;
  else if (strcmp(arg[2],"Ni/fcc") == 0) 
    cut = 3.00451;
  else if (strcmp(arg[2],"Pb/fcc") == 0) 
    cut = 4.22509;
  else if (strcmp(arg[2],"Pd/fcc") == 0) 
    cut = 3.32032;
  else if (strcmp(arg[2],"Pt/fcc") == 0) 
    cut = 3.34593;
  else if (strcmp(arg[2],"Rb/bcc") == 0) 
    cut = 6.74773;
  else if (strcmp(arg[2],"Rh/fcc") == 0) 
    cut = 3.24350;
  else if (strcmp(arg[2],"Si/dia") == 0) 
    cut = 4.19095;
  else if (strcmp(arg[2],"Sr/fcc") == 0) 
    cut = 5.18960;
  else if (strcmp(arg[2],"Ta/bcc") == 0) 
    cut = 3.99552;
  else if (strcmp(arg[2],"Th/fcc") == 0) 
    cut = 4.33605;
  else if (strcmp(arg[2],"V/bcc") == 0) 
    cut = 3.64546;
  else if (strcmp(arg[2],"W/bcc") == 0) 
    cut = 3.81446;
  else if (strcmp(arg[2],"Xe/fcc") == 0) 
    cut = 5.29203;
  else if (strcmp(arg[2],"Yb/fcc") == 0) 
    cut = 4.68601;
  else {
    cut = universe->numeric(FLERR,arg[2]);
    if (cut < 0.0) error->all(FLERR,"Illegal coarse_graining command");
  }


  a1[0] = 1.0;  a1[1] = 0.0;  a1[2] = 0.0;
  a2[0] = 0.0;  a2[1] = 1.0;  a2[2] = 0.0;
  a3[0] = 0.0;  a3[1] = 0.0;  a3[2] = 1.0;

  if (strcmp(arg[3],"fcc") == 0) {
    structure_type = FCC;
    compute_style = CNA;
    a2[0] = 0.5; a2[1] = sqrt(3.0/4.0);
    a3[0] = 0.5; a3[1] = sqrt(1.0/12.0); a3[2] = sqrt(2.0/3.0);
  } else if (strcmp(arg[3],"bcc") == 0) {
    structure_type = BCC;
    compute_style = CNA;
    a2[0] = -1.0/3.0; a2[1] = sqrt(8.0/9.0);
    a3[0] = 1.0/3.0;  a3[1] = sqrt(2.0/9.0); a3[2] = sqrt(2.0/3.0);
  } else if (strcmp(arg[3],"cubic/diamond") == 0) {
    structure_type = CUBIC;
    compute_style = IDS;
    a2[0] = 0.5; a2[1] = sqrt(3.0/4.0);
    a3[0] = 0.5; a3[1] = sqrt(1.0/12.0); a3[2] = sqrt(2.0/3.0);
  } else error->all(FLERR,"Invalid structure type for coarse_graining command");

  ngrains = atom->ngrains;
  if (ngrains)
    grain_tag_enable = 1;
  else {
    grain_tag_enable = 0;
    ngrains = 1;
  }

  memory->create(centroid,ngrains+1,3,"coarse_graining::centroid");

  grain_etype = universe->inumeric(FLERR,arg[4]);
  grain_ctype = universe->inumeric(FLERR,arg[5]);
  if (grain_etype > element->netypes || grain_etype <= 0)
    error->all(FLERR,"Invalid etype in coarse_graining command");
  if (grain_ctype > atom->ntypes || grain_ctype <= 0)
    error->all(FLERR,"Invalid ctype in coarse_graining command");

  if (element->element_shape_ids[grain_etype] != Element::HEXAHEDRON)
    error->all(FLERR,"Coarse_graining command can only create Hexahedron element");
  ncellx = evec->ncells[grain_etype][0];
  ncelly = evec->ncells[grain_etype][1];
  ncellz = evec->ncells[grain_etype][2];
  if (ncellx != ncelly || ncellx != ncellz || ncelly != ncellz)
    error->all(FLERR,"Coarse_graining command can only create elements with uniform sizes");

  // options

  options(narg-6,&arg[6]);

  cutsq = cut*cut;


  // request a full atomic neighbor list for use by this command 
  // also sort out neighbors based on distance

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;
  neighbor->requests[irequest]->sort = 1;
  neighbor->requests[irequest]->sortnum = 14;
  neighbor->requests[irequest]->occasional = 2;
  neighbor->requests[irequest]->atomonlylist = 1;
  neighbor->requests[irequest]->atomlist = 0;
  neighbor->requests[irequest]->elemlist = 0;
  neighbor->requests[irequest]->intglist = 0;

  // use custom cutoff to save computational time

  neighbor->requests[irequest]->cut = 1;
  neighbor->requests[irequest]->cutoff = cut;


  char **newarg = new char*[4];
  if (compute_style == CNA) { 
    id_compute = (char *) "COARSE_GRAINING_CNA";
    newarg[2] = (char *) "cna/atom";
  } else if (compute_style == IDS) { 
    id_compute = (char *) "COARSE_GRAINING_IDS";
    newarg[2] = (char *) "ids/atom";
  }
  newarg[0] = id_compute;
  newarg[1] = (char *) "atom";
  newarg[3] = arg[2];
  modify->add_compute(4,newarg);
  delete [] newarg;


  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR,"Cannot create compute for coarse_graining command");
  compute = modify->compute[icompute];

  compute->preneighflag = 1;

  // init entire system since comm->borders and neighbor->build is done

  cac->init();

  // error check on cutoff

  if (force->pair == NULL)
    error->all(FLERR,"Coarse_graining requires a pair style be defined");
  if (cut > force->pair->cutforce)
    error->all(FLERR,"Coarse_graining cutoff is longer than pairwise cutoff");

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  // setup domain, communication and neighboring
  // exchange and acquire ghosts
  // build standard neighbor list to setup bins & stencils

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();
  comm->exchange();
  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost,atom->x);
    domain->lamda2x(element->nlocal+element->nghost,element->x);
    domain->lamda2nodex(element->nlocal+element->nghost,element->nodex);
  }
  neighbor->build_bins();

  //neighbor->build();

  // build neighbor list this command needs based on earlier request

  NeighList *list = neighbor->lists[irequest];
  neighbor->build_one(list,1);

  int natoms_previous = atom->natoms;
  int nelements_previous = element->nelements;
  int nlocal_previous = element->nlocal;

  int nlocal = atom->nlocal;

  // invoke compute structure and extract pattern
  // invoke forward comm to communicate ghost pattern

  compute->compute_peratom();
  comm->forward_comm_compute(compute);
  pattern = compute->vector_atom;

  int i,ii,j;
  double **x = atom->x;
  double unmap_x[3];
  imageint *image = atom->image;

  int nempty = 0;
  int *empty_list = new int[ngrains];

  // to be delete flag for each atom

  memory->create(select_flag,nlocal,"coarse_graining:select_flag");
  for (i = 0; i < nlocal; i++)
    select_flag[i] = 0;

  // copy data from neighbor list

  inum = list->ainum;
  ilist = list->ailist;
  numneigh = list->numneigha;
  firstneigh = list->firstneigha;

  grain_tag = atom->grain_tag;
  grain_seed_ids = new int[ngrains+1]; 
  for (i = 1; i <= ngrains; i++)
    grain_seed_ids[i] = -1;

  memory->create(atom_list,ncellx,ncelly,ncellz,"coarse_graining:atom_list");
  memory->create(nodecoord,8,3,"coarse_graining:nodecoord");

  if (grain_tag_enable) {

    // choose seed ids 
    // select atom closest to the grain centroid

    for (int igrain = 1; igrain <= ngrains; igrain++) {

      // calculate centroid
      // reset select_flag

      if (centroid_flag) {
        for (i = 0; i < nlocal; i++) 
          select_flag[i] = 0;

        centroid[igrain][0] = centroid[igrain][1] = centroid[igrain][2] = 0;
        centroid_count = 0;
        for (i = 0; i < nlocal; i++) {
          if (igrain == grain_tag[i] && pattern[i] == structure_type) {
            domain->unmap(x[i],image[i],unmap_x);
            centroid[igrain][0] += unmap_x[0];
            centroid[igrain][1] += unmap_x[1];
            centroid[igrain][2] += unmap_x[2];
            centroid_count++;
          }
        }

        centroid[igrain][0] /= centroid_count;  
        centroid[igrain][1] /= centroid_count;  
        centroid[igrain][2] /= centroid_count;  
      }

      domain->remap(centroid[igrain]);

      //if (igrain == 0) printf(" grain centroid = %g %g %g\n",centroid[igrain][0],centroid[0][1],centroid[0][2]);

      // find atom closest to grain centroid

      double delx,dely,delz,rsq;
      double min = 1e30;
      for (i = 0; i < nlocal; i++) {
        if (grain_tag[i] == igrain && pattern[i] == structure_type) {
          delx = centroid[igrain][0] - x[i][0];
          dely = centroid[igrain][1] - x[i][1];
          delz = centroid[igrain][2] - x[i][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < min) {
            min = rsq; 
            grain_seed_ids[igrain] = i;
          }
        }
      }
    } 


    for (i = 0; i < nlocal; i++) 
      select_flag[i] = 0;
    int tmp = element->nlocal; 
    for (int igrain = 1; igrain <= ngrains; igrain++) {
      if (grain_seed_ids[igrain] >= 0) 
        coarse_grain(grain_seed_ids[igrain],igrain,a1,a2,a3,0,0,0,1);
      if (element->nlocal == tmp)
        empty_list[nempty++] = igrain;
      else tmp = element->nlocal;
    }
  } else {

  }

  if (nempty) {
    char errstr[200];
    sprintf(errstr,"%d grains are too small for coarse "
        "graining or have no atoms\n",nempty);
    error->warning(FLERR,errstr);
  }

  // delete atoms flagged in select_flag
  // reset nlocal

  i = 0;
  while (i < nlocal) {
    if (select_flag[i]) {
      avec->copy(nlocal-1,i,1);
      select_flag[i] = select_flag[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;

  // init per-element fix/compute/variable values for created elements

  element->data_fix_compute_dump_variable(nlocal_previous,element->nlocal);

  // if compress flag set,
  // reset atom/elements tags to be contiguous
  // set all atom/element IDs to 0, call tag_extend()

  if (compress_flag) {
    tagint *tag = atom->tag;
    for (i = 0; i < atom->nlocal; i++) tag[i] = 0;
    atom->tag_extend();

    tag = element->tag;
    for (i = 0; i < element->nlocal; i++) tag[i] = 0;
    element->tag_extend();
  } else {

    // assign tags for created elements

    if (element->tag_enable) element->tag_extend();
    element->tag_check();
  }

  // remap element since it might be created outside of periodic box

  double **ex = element->x;
  double ***nodex = element->nodex;
  int *etype = element->etype;
  int *npe = element->npe;
  for (i = 0; i < element->nlocal; i++) domain->remap(ex[i],nodex[i],npe[etype[i]]);

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);

  nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&element->nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);

  bigint natoms_delete = natoms_previous - atom->natoms;
  bigint nelements_add = element->nelements - nelements_previous;

  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,"Coarse grained " BIGINT_FORMAT
          " atoms into " BIGINT_FORMAT " elements\n "
          " New total atoms = " BIGINT_FORMAT 
          ", new total elements = " BIGINT_FORMAT "\n",
          natoms_delete,nelements_add,atom->natoms,element->nelements);
    if (logfile) 
      fprintf(logfile,"Coarse grained " BIGINT_FORMAT
          " atoms into " BIGINT_FORMAT " elements\n "
          " New total atoms = " BIGINT_FORMAT 
          ", new total elements = " BIGINT_FORMAT "\n",
          natoms_delete,nelements_add,atom->natoms,element->nelements);
  }

  // clear map if needed

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (element->map_style) {
    element->nghost = 0;
    element->map_init();
    element->map_set();
  }

  // move atoms/elements back inside simulation box and to new processors
  // use irregular() in case atoms/elements moved a long distance
  // must be done after updating new global natoms/nelements

  if (migrate_flag) {
    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal,atom->x);
      domain->x2lamda(element->nlocal,element->x);
      domain->nodex2lamda(element->nlocal,element->nodex);
    }
    domain->reset_box();
    IrregularComm *irregular_comm = new IrregularComm(cac);
    irregular_comm->migrate(1);
    delete irregular_comm;
    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal,atom->x);
      domain->lamda2x(element->nlocal,element->x);
      domain->lamda2nodex(element->nlocal,element->nodex);
    }
  }

  // clean up

  memory->destroy(select_flag);
  memory->destroy(atom_list);
  memory->destroy(nodecoord);
  memory->destroy(centroid);
  delete [] grain_seed_ids;
  delete [] empty_list;

}

/* ----------------------------------------------------------------------
   process command options
   ------------------------------------------------------------------------- */

void CoarseGraining::options(int narg, char **arg)
{
  compress_flag = 1;
  migrate_flag = 1;
  wrapx = wrapy = wrapz = 1;
  centroid_flag = 1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarse_graining command");
      if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
      else error->all(FLERR,"Illegal coarse_graining command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"migrate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarse_graining command");
      if (strcmp(arg[iarg+1],"yes") == 0) migrate_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) migrate_flag = 0;
      else error->all(FLERR,"Illegal coarse_graining command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"wrap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarse_graining command");
      if (strcmp(arg[iarg+1],"yes") == 0) wrapx = wrapy = wrapz = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wrapx = wrapy = wrapz = 0;
      else if (strcmp(arg[iarg+1],"x") == 0) {
        wrapx = 1;
        wrapy = wrapz = 0;
      } else if (strcmp(arg[iarg+1],"y") == 0) {
        wrapy = 1;
        wrapx = wrapz = 0;
      } else if (strcmp(arg[iarg+1],"z") == 0) {
        wrapz = 1;
        wrapx = wrapy = 0;
      } else if (strcmp(arg[iarg+1],"xy") == 0) {
        wrapx = wrapy = 1;
        wrapz = 0;
      } else if (strcmp(arg[iarg+1],"xz") == 0) {
        wrapx = wrapz = 1;
        wrapy = 0;
      } else if (strcmp(arg[iarg+1],"yz") == 0) {
        wrapz = wrapy = 1;
        wrapx = 0;
      } else error->all(FLERR,"Illegal coarse_graining command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarse_graining command");
      tol = universe->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"a1") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal coarse_graining command");
      a1[0] = universe->numeric(FLERR,arg[iarg+1]);
      a1[1] = universe->numeric(FLERR,arg[iarg+2]);
      a1[2] = universe->numeric(FLERR,arg[iarg+3]);
      normalize3(a1,a1);
      if (a1[0] == 0 && a1[1] == 0 && a1[2] == 0) error->all(FLERR,
          "Invalid values for a1 in coarse_graining command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"a2") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal coarse_graining command");
      a2[0] = universe->numeric(FLERR,arg[iarg+1]);
      a2[1] = universe->numeric(FLERR,arg[iarg+2]);
      a2[2] = universe->numeric(FLERR,arg[iarg+3]);
      normalize3(a2,a2);
      if (a2[0] == 0 && a2[1] == 0 && a2[2] == 0) error->all(FLERR,
          "Invalid values for a2 in coarse_graining command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"a3") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal coarse_graining command");
      a3[0] = universe->numeric(FLERR,arg[iarg+1]);
      a3[1] = universe->numeric(FLERR,arg[iarg+2]);
      a3[2] = universe->numeric(FLERR,arg[iarg+3]);
      normalize3(a3,a3);
      if (a3[0] == 0 && a3[1] == 0 && a3[2] == 0) error->all(FLERR,
          "Invalid values for a3 in coarse_graining command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"centroid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal coarse_graining command");
      read_centroid(arg[iarg+1]);
      centroid_flag = 0;
      iarg += 2;
    } else error->all(FLERR,"Illegal coarse_graining command");
  } 
}
/* ----------------------------------------------------------------------
   coarse graining from root atom I to element 
   iex, iey, and iez are the # offset element in each direction from seed element
   start from 0 0 0 and proceed recursively in 6 directions
   ------------------------------------------------------------------------- */

void CoarseGraining::coarse_grain(int i, int grain_id, double *a1,
    double *a2, double *a3, int iex, int iey, int iez, int first_flag)
{
  int ii,jj,kk;
  int j,k,l,j_remap;
  for (ii = 0; ii < ncellx; ii++)
    for (jj = 0; jj < ncelly; jj++)
      for (kk = 0; kk < ncellz; kk++)
        atom_list[ii][jj][kk] = -1;
  atom_list[0][0][0] = i;
  int nlocal = atom->nlocal;
  int success;
  double pa1[3],pa2[3],pa3[3],na1[3],na2[3],na3[3];
  double **x = atom->x;

  if (first_flag) {
    j = align(i,-1,a1,0,1);
    k = align(i,j,a2,0,1);
    if (k < 0) error->one(FLERR,"Coarse_graining command can't detect crystal orientation");
    int flag = 1;
    for (int ll = 0; ll < numneigh[i]; ll++) {
      l = firstneigh[i][ll];
      if (is_neighbor(j,l) && 
          is_neighbor(k,l)) {
        flag = 0;
        break;     
      }
    }

    //int t;
    //if (grain_id == 0) {
    //  printf("%d %d %g %g %g\n",atom->tag[i],5
    //      ,atom->x[i][0]
    //      ,atom->x[i][1]
    //      ,atom->x[i][2]
    //      );
    //  for (ii = 0; ii < numneigh[i]; ii++) {
    //    if (firstneigh[i][ii] == j)
    //      t = 1;
    //    else if (firstneigh[i][ii] == k)
    //      t = 2;
    //    else if (firstneigh[i][ii] == l)
    //      t = 3;
    //    else t = 4;
    //    printf("%d %d %g %g %g\n",atom->tag[firstneigh[i][ii]],t
    //        ,atom->x[firstneigh[i][ii]][0]
    //        ,atom->x[firstneigh[i][ii]][1]
    //        ,atom->x[firstneigh[i][ii]][2]
    //        );
    //  }
    //}
    if (flag) error->one(FLERR,"Coarse_graining command can't detect crystal orientation");
    sub3(x[j],x[i],pa1);
    normalize3(pa1,pa1);
    sub3(x[k],x[i],pa2);
    normalize3(pa2,pa2);
    sub3(x[l],x[i],pa3);
    normalize3(pa3,pa3);

    success = collect_atoms(grain_id,pa1,pa2,pa3,0,0,0);
    if (!success) {
      copy3(pa1,na1); copy3(pa2,na2); copy3(pa3,na3);
      negate3(na1); negate3(na2); negate3(na3);
      if (!collect_atoms(grain_id,pa1,pa2,na3,0,0,0))
        if (!collect_atoms(grain_id,pa1,na2,pa3,0,0,0))
          if (!collect_atoms(grain_id,na1,pa2,pa3,0,0,0))
            if (!collect_atoms(grain_id,na1,na2,pa3,0,0,0))
              if (!collect_atoms(grain_id,na1,pa2,na3,0,0,0))
                if (!collect_atoms(grain_id,pa1,na2,na3,0,0,0))
                  if (!collect_atoms(grain_id,na1,na2,na3,0,0,0))
                    return;
    }
  } else {
    if (!collect_atoms(grain_id,a1,a2,a3,0,0,0)) return;
  }

  // check to make sure all atoms are found and collected

  for (ii = 0; ii < ncellx; ii++)
    for (jj = 0; jj < ncelly; jj++)
      for (kk = 0; kk < ncellz; kk++)
        if (atom_list[ii][jj][kk] < 0) return;

  // create new element

  int nodelist[8];
  nodelist[0] = atom_list[0][0][0];
  nodelist[1] = atom_list[ncellx-1][0][0];
  nodelist[2] = atom_list[ncelly-1][ncelly-1][0];
  nodelist[3] = atom_list[0][ncelly-1][0];
  nodelist[4] = atom_list[0][0][ncellz-1];
  nodelist[5] = atom_list[ncellx-1][0][ncellz-1];
  nodelist[6] = atom_list[ncellx-1][ncelly-1][ncellz-1];
  nodelist[7] = atom_list[0][ncelly-1][ncellz-1];

  for (ii = 0; ii < 8; ii++) {
    //printf("nodelist[%d] = %d nlocal = %d nghost = %d nall = %d\n",ii,nodelist[ii],atom->nlocal,atom->nghost,atom->nlocal+atom->nghost);
    copy3(x[nodelist[ii]],nodecoord[ii]);
  }

  evec->create_element(nodecoord,grain_etype,grain_ctype,0);

  // mark atoms in list for removal later
  // don't delete atoms now since it will mess up neighbor list

  for (ii = 0; ii < ncellx; ii++)
    for (jj = 0; jj < ncelly; jj++)
      for (kk = 0; kk < ncellz; kk++)
        select_flag[atom->map(atom->tag[atom_list[ii][jj][kk]])] = 1;

  // calculate new set of a vectors from new element

  for (j = 0; j < 3; j++) {
    pa1[j] = nodecoord[1][j] - nodecoord[0][j] + 
      nodecoord[2][j] - nodecoord[3][j] +
      nodecoord[5][j] - nodecoord[4][j] +
      nodecoord[6][j] - nodecoord[7][j];
    pa2[j] = nodecoord[3][j] - nodecoord[0][j] +
      nodecoord[2][j] - nodecoord[1][j] +
      nodecoord[7][j] - nodecoord[4][j] +
      nodecoord[6][j] - nodecoord[5][j];
    pa3[j] = nodecoord[4][j] - nodecoord[0][j] +
      nodecoord[7][j] - nodecoord[3][j] +
      nodecoord[5][j] - nodecoord[1][j] +
      nodecoord[6][j] - nodecoord[2][j];
    na1[j] = -pa1[j];
    na2[j] = -pa2[j];
    na3[j] = -pa3[j];
  }
  normalize3(pa1,pa1);
  normalize3(pa2,pa2);
  normalize3(pa3,pa3);
  normalize3(na1,na1);
  normalize3(na2,na2);
  normalize3(na3,na3);

  // find new seed atoms to continue coarse graining in 6 directions

  j = align(nodelist[4],-1,pa3,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,pa1,pa2,pa3,iex,iey,iez+1);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,pa1,pa2,pa3,iex,iey,iez+1);
      }
    }
  }

  j = align(nodelist[3],-1,pa2,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,pa1,pa2,pa3,iex,iey+1,iez);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,pa1,pa2,pa3,iex,iey+1,iez);
      }
    }
  }

  j = align(nodelist[1],-1,pa1,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,pa1,pa2,pa3,iex+1,iey,iez);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,pa1,pa2,pa3,iex+1,iey,iez);
      }
    }
  }

  j = align(nodelist[0],-1,na3,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,pa1,pa2,na3,iex,iey,iez-1);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,pa1,pa2,na3,iex,iey,iez-1);
      }
    }
  }

  j = align(nodelist[0],-1,na2,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,pa1,na2,pa3,iex,iey-1,iez);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,pa1,na2,pa3,iex,iey-1,iez);
      }
    }
  }

  j = align(nodelist[0],-1,na1,tol,1);
  if (j >= 0) {
    j_remap = atom->map(atom->tag[j]);
    if (!select_flag[j_remap]) {
      if (check_neighbors(j_remap,grain_id)) {
        if (j >= nlocal) {
          if ((fabs(x[j][0] - x[j_remap][0]) > EPSILON) <= wrapx &&
              (fabs(x[j][1] - x[j_remap][1]) > EPSILON) <= wrapy &&
              (fabs(x[j][2] - x[j_remap][2]) > EPSILON) <= wrapz)
            coarse_grain(j_remap,grain_id,na1,pa2,pa3,iex-1,iey,iez);
        } else if (j < nlocal)
          coarse_grain(j,grain_id,na1,pa2,pa3,iex-1,iey,iez);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   coarse graining 1 element from root atom I to element 
   iax, iay, and iaz are the # offset atom in each direction from seed atom
   start from 0 0 0 and proceed recursively in 3 directions
   ------------------------------------------------------------------------- */

int CoarseGraining::collect_atoms(int grain_id, double *old_a1, double *old_a2, 
    double *old_a3, int iax, int iay, int iaz)
{
  int i = atom_list[iax][iay][iaz]; 
  if (i < 0) return 0;

  if (pattern[i] != structure_type || grain_tag[i] != grain_id)
    return 0;
  if (select_flag[atom->map(atom->tag[i])]) return 0;

  int j,k,l,success;
  double **x = atom->x;
  double del[3];
  double new_a1[3];
  double new_a2[3];
  double new_a3[3];

  j = align(i,-1,old_a1,tol,1);
  if (j < 0) return 0;
  k = align(i,-1,old_a2,tol,1);
  if (k < 0) return 0;
  l = align(i,-1,old_a3,tol,1);
  if (l < 0) return 0;

  copy3(old_a1,new_a1);
  copy3(old_a2,new_a2);
  copy3(old_a3,new_a3);

  if (iaz+1 < ncellz) {
    atom_list[iax][iay][iaz+1] = l;
    sub3(x[l],x[i],new_a3);
    normalize3(new_a3,new_a3);
    success = collect_atoms(grain_id,new_a1,new_a2,new_a3,iax,iay,iaz+1);
    if (!success) return 0;
  } 

  if (iaz == 0 && iay+1 < ncelly) {
    atom_list[iax][iay+1][iaz] = k;
    sub3(x[k],x[i],new_a2);
    normalize3(new_a2,new_a2);
    success = collect_atoms(grain_id,new_a1,new_a2,new_a3,iax,iay+1,iaz);
    if (!success) return 0;
  } 

  if (iay == 0 && iaz == 0 && iax+1 < ncellx) {
    atom_list[iax+1][iay][iaz] = j;
    sub3(x[j],x[i],new_a1);
    normalize3(new_a1,new_a1);
    success = collect_atoms(grain_id,new_a1,new_a2,new_a3,iax+1,iay,iaz);
    if (!success) return 0;
  }


  return 1;
}

/* ----------------------------------------------------------------------
   Find neighbor of I that best align with vector a in direction dir (1/0 = same/opposit e direction)
   If I2 is given (>= 0), it has to be neighbor of I2 as well
   Return index of neighbor that align
   If positive tolerant is given, the alignment must be close within tolerant,
   otherwise return -1
   ------------------------------------------------------------------------- */

int CoarseGraining::align(int i, int i2, double *a, double tolerant, int dir)
{
  int j,jj,j_align,jnum,*jlist;
  double **x = atom->x;
  double del[3],tmp;
  int maxneigh;
  if (compute_style == CNA) {
    if (structure_type == FCC) maxneigh = 12;
    if (structure_type == BCC) maxneigh = 14;
  } else if (compute_style == IDS) maxneigh = 16;
  jnum = MIN(numneigh[i],maxneigh);
  jlist = firstneigh[i];

  j_align = -1;
  if (dir) {
    double max = -2; 
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (i2 > 0) 
        if (!is_neighbor(i2,j))
          continue;
      sub3(x[j],x[i],del);
      normalize3(del,del);
      tmp = dot3(a,del);
      if (tmp > max) {
        max = tmp;
        j_align = j;
      }
    }
    if (tolerant > 0 && 1.0-max > tolerant) return -1;
  } else {
    double min = 2; 
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (i2 > 0) 
        if (!is_neighbor(i2,j))
          continue;
      sub3(x[j],x[i],del);
      normalize3(del,del);
      tmp = dot3(a,del);
      if (tmp < min) {
        min = tmp;
        j_align = j;
      }
    }
    if (tolerant > 0 && 1.0+min > tolerant) return -1;
  }

  return j_align;
}

/* ----------------------------------------------------------------------
   check if neighbors of I are valid 
   ------------------------------------------------------------------------- */

int CoarseGraining::check_neighbors(int i, int grain_id)
{
  int j;
  int jnum = MIN(numneigh[i],12);
  int *jlist = firstneigh[i];
  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    if (grain_id != grain_tag[j] || 
        structure_type != pattern[j])
      return 0;
  }
  return 1;
}

/* ----------------------------------------------------------------------
   check if J is neighbor of I
   ------------------------------------------------------------------------- */

int CoarseGraining::is_neighbor(int i, int j)
{
  int jj,*jlist,jnum;
  jnum = numneigh[i];
  jlist = firstneigh[i];
  for (jj = 0; jj < jnum; jj++) {
    if (jlist[jj] == j) return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   read grain centroids from file
   ------------------------------------------------------------------------- */

void CoarseGraining::read_centroid(char *file)
{
  if (me == 0) {
    FILE *fp;
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open centroid file %s for coarse_graining command",file);
      error->one(FLERR,str);
    }
    char *line = new char[MAXLINE];
    char *ptr;
    for (int i = 1; i <= ngrains; i++) {
      fgets(line,MAXLINE,fp);
      ptr = strtok(line," \t\n\r\f");
      if (ptr)
        centroid[i][0] = atof(ptr);
      else error->one(FLERR,"Missing value in centroid file for coarse_graining command");
      if ((ptr = strtok(NULL," \t\n\r\f"))) 
        centroid[i][1] = atof(ptr);
      else error->one(FLERR,"Missing value in centroid file for coarse_graining command");
      if ((ptr = strtok(NULL," \t\n\r\f"))) 
        centroid[i][2] = atof(ptr);
      else error->one(FLERR,"Missing value in centroid file for coarse_graining command");
    }
    fclose(fp);
    fp = NULL;
    delete line;
  }
  MPI_Bcast(&centroid[0][0],ngrains*3+3,MPI_DOUBLE,0,world);
}


/* Debug function
   void CoarseGraining::header()
   {
   fprintf(fp,"ITEM: TIMESTEP\n");

   fprintf(fp,"%d\n",0);

   fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
   fprintf(fp,"%d\n",3*16*16);
   char boundstr[9];
   domain->boundary_string(boundstr);
   if (domain->triclinic == 0) {
   fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
   fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
   fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
   fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
   } else {
   fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
   fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
   fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
   fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
   }
   fprintf(fp,"ITEM: ATOMS id index x y z\n");
   }
   */

