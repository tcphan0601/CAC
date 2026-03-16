#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "fix_adaptive_2.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "compute.h"
#include "domain.h"
#include "memory.h"
#include "group.h"
#include "variable.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixAdaptive2::FixAdaptive2(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  id_compute(NULL), splitlist(NULL)
{
  create_attribute = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  //if (strcmp(element->element_style,"cac") != 0) 
  //  error->all(FLERR,"Fix adaptive only works for rhb element style");
  if (element->max_apc > 1) error->all(FLERR,"Fix adaptive does not work for multiple atoms per unit cell yet");

  if (narg < 6) error->all(FLERR,"Illegal fix adaptive command");

  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix adaptive command");

  int n = strlen(arg[4]) + 1; 
  id_compute = new char[n];
  strcpy(id_compute,arg[4]);
  icompute = modify->find_compute(id_compute);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix adative does not exist");
  compute = modify->compute[icompute];

  perfect_crystal = universe->inumeric(FLERR,arg[5]);

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    error->all(FLERR,"Illegal fix adaptive command");
  }
  
  maxsplit = 1024;
  grow_splitlist(); 

  pending = 0;
}

/* ---------------------------------------------------------------------- */

FixAdaptive2::~FixAdaptive2()
{
  delete [] id_compute;
  delete [] splitlist;
}

/* ---------------------------------------------------------------------- */

int FixAdaptive2::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive2::init()
{
  // need an occasional full intpl neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->atomlist = 0;
  neighbor->requests[irequest]->intglist = 0;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive2::init_list(int id, NeighList *ptr)
{
  list = ptr;
}


/* ---------------------------------------------------------------------- */

void FixAdaptive2::setup(int vflag)
{
  //next_reneighbor = (update->ntimestep/nevery)*nevery;
}

void FixAdaptive2::setup_pre_exchange()
{
  // need to exchange, acquire ghosts, build neigh list on the 
  // first step in order for compute to build neighbor list
  // same as in Verlet::setup()

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
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2nodex(element->nlocal,element->nodex);
  }
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();

  neighbor->build();
  neighbor->ncalls = 0;

  next_reneighbor = (update->ntimestep/nevery)*nevery;

  post_integrate();
  next_reneighbor = (update->ntimestep/nevery)*nevery+nevery;

  // if splitting
  // reset neighbor->lastcall to -1 
  // to force neigh list rebuild in Veret::setup()
  
  if (pending) {
    neighbor->lastcall = -1; 
    pending = 0;
  }
}

/* ----------------------------------------------------------------------
   if pending, dislocation detected at this time step,  
   set new next_reneighbor here after forcing neigh list rebuild
   ------------------------------------------------------------------------- */

void FixAdaptive2::pre_exchange()
{
  if (!pending) return;
  next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
  pending = 0;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive2::post_integrate()
{
  // return if not an adaptive timestep

  if (update->ntimestep < next_reneighbor) return;

  // update slip plane for all elements

  element->evec->update_slip_plane();

  // invoke compute peratom
  // set pre_neighbor_flag = 1;

  if (compute->invoked_peratom != update->ntimestep) {
    compute->compute_peratom();
    compute->preneighflag = 1;
  }

  // acquire ghost structure id
  
  comm->forward_comm_compute(compute);

  // build neighbor list for elements

  neighbor->build_one(list,1);

  int inum = list->einum;
  int *ilist = list->eilist;
  int *numneighe2a = list->numneighe2a;
  int *numneighe2e = list->numneighe2e;
  int **firstneighe2a = list->firstneighe2a;
  int **firstneighe2e = list->firstneighe2e;
  int i,j,k,l,ii,jj,kk,ll,jintpl,jetype;
  int jnum,knum,*jlist,*klist;
  double xtmp,ytmp,ztmp;
  double **ax = atom->x;
  double **ex = element->x;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  double *apattern = compute->vector_atom;
  double **iapattern = compute->vector_intpl_atom;
  int *nsurface = element->nsurface;
  int ***surface_intpl = element->surface_intpl;
  int **nsurface_intpl = element->nsurface_intpl;
  int cut,cutdim;
  int nsplit = 0;

  // loop through neighbors of elements to find dislocation nearby
  // add element to splitlist
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = ex[i][0];
    ytmp = ex[i][1];
    ztmp = ex[i][2];
    if (!(emask[i] & groupbit)) continue;


    // loop over atom neighbors

    //jnum = numneighe2a[i];
    //jlist = firstneighe2a[i];
    //for (jj = 0; jj < jnum; jj++) {
    //  j = jlist[jj];
    //  j &= NEIGHMASK;
    //  if (!(amask[j] & groupbit) || apattern[j] == perfect_crystal) continue;
    //}

    // loop over element neighbors

    jnum = numneighe2e[i];
    jlist = firstneighe2e[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jetype = etype[j];
      j &= NEIGHMASK;
      if (!(emask[j] & groupbit)) continue;
      for (kk = 0; kk < nsurface[jetype]; kk++) {
        knum = nsurface_intpl[jetype][kk];
        klist = surface_intpl[jetype][kk];
        for (k = 0; k < knum; k++) {
          jintpl = klist[k];
          if (iapattern[j][jintpl] == perfect_crystal) continue;
          if (element->evec->check_split_elem(i,j,kk,cut,cutdim)) {
            if (nsplit == maxsplit) grow_splitlist();
            splitlist[nsplit][0] = i;
            splitlist[nsplit][1] = cut;
            splitlist[nsplit][2] = cutdim;
            nsplit++;
          }
          break;
        }
      }
    }
  }

  int nsplit_all;
  MPI_Allreduce(&nsplit,&nsplit_all,1,MPI_INT,MPI_MAX,world);

  // return and update next reneighbor if no plane need to split
  // this will not force reneighbor at this time step

  if (nsplit_all == 0) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen,"No new dislocation found\n");
    return;
  }

  pending = 1;
  int nelocal = element->nlocal;
  int nalocal = atom->nlocal;

  // split elements

  int nesplit,nesplit_all;
  nesplit = element->evec->split_element(splitlist,nsplit);

  // extend tags to new atoms and elements

  atom->tag_extend();
  element->tag_extend();

  int nnew_atoms = atom->nlocal - nalocal;
  int nnew_elems = element->nlocal - nelocal;
  int nnew_atoms_all,nnew_elems_all;
  MPI_Allreduce(&nnew_atoms,&nnew_atoms_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nnew_elems,&nnew_elems_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nesplit,&nesplit_all,1,MPI_INT,MPI_SUM,world);

  // update number of elements and number of atoms

  atom->natoms += nnew_atoms_all;
  element->nelements += nnew_elems_all;
  if (comm->me == 0) {
    if (screen) fprintf(screen,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nesplit_all,nnew_atoms_all,nnew_elems_all+nesplit_all);
    if (logfile) fprintf(logfile,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nesplit_all,nnew_atoms_all,nnew_elems_all+nesplit_all);
  }
  
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

void FixAdaptive2::grow_splitlist()
{
  maxsplit += maxsplit/2;
  memory->grow(splitlist,maxsplit,3,"fix:splitlist");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

double FixAdaptive2::memory_usage()
{
  return (atom->nmax + element->nmax) * sizeof(double);
}

