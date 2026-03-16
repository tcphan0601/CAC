#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "fix_adaptive.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "neighbor.h"
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

FixAdaptive::FixAdaptive(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  id_compute(NULL), atom_plane_id(NULL), elem_plane_id(NULL),
  splitlist(NULL), splitflag(NULL), disorder_flag(NULL), disorder_flagall(NULL)
{
  create_attribute = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  //if (strcmp(element->element_style,"cac") != 0) 
  //  error->all(FLERR,"Fix adaptive only works for rhb element style");
  //if (element->apc > 1) error->all(FLERR,"Fix adaptive does not work for multiple atoms per unit cell yet");

  if (narg < 9) error->all(FLERR,"Illegal fix adaptive command");

  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix adaptive command");
  l = universe->numeric(FLERR,arg[4]); 
  if (strcmp(arg[5],"xy") == 0) split_dim = 2;
  else if (strcmp(arg[5],"xz") == 0) split_dim = 1;
  else if (strcmp(arg[5],"yz") == 0) split_dim = 0;
  else error->all(FLERR,"Illegal fix adaptive command");

  int iregion = domain->find_region(arg[6]);
  Region *region = domain->regions[iregion];
  if (strcmp(region->style,"block") != 0) 
    error->all(FLERR,"Must use region block for fix adaptive command");

  blo[0] = region->extent_xlo;
  bhi[0] = region->extent_xhi;
  blo[1] = region->extent_ylo;
  bhi[1] = region->extent_yhi;
  blo[2] = region->extent_zlo;
  bhi[2] = region->extent_zhi;
 
  int n = strlen(arg[7]) + 1; 
  id_compute = new char[n];
  strcpy(id_compute,arg[7]);
  icompute = modify->find_compute(id_compute);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix adative does not exist");
  compute = modify->compute[icompute];

  perfect_crystal = universe->inumeric(FLERR,arg[8]);

  // optional args

  int iarg = 9;
  while (iarg < narg) {
    error->all(FLERR,"Illegal fix adaptive command");
  }
  
  // add one plane above to make sure atoms/nodes are included
  nplanes = static_cast<int> ((bhi[split_dim]-blo[split_dim])/l+1);
  bhi[split_dim] = blo[split_dim] + l*nplanes;

  splitflag = new int[nplanes-1];
  splitlist = new int[nplanes-1];
  for (int i = 0; i < nplanes-1; i++) splitflag[i] = 0;

  disorder_flag = new int[nplanes];
  disorder_flagall = new int[nplanes];
  
  // perform initial allocation of atom-based and elem-based arrays
  // register with Atom and Element class

  grow_atom_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);
  grow_elem_arrays(element->nmax);
  element->add_callback(0);
  element->add_callback(1);
  vector_atom = atom_plane_id;
  vector_elem = elem_plane_id;

  // assign atoms/interpolated atoms
  // skip if reset from restart file

  double **ax = atom->x;
  int *amask = atom->mask;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++) 
    if (amask[i] & groupbit) 
      atom_plane_id[i] = find_plane(ax[i]);

  double ***nodex = element->nodex;
  int *emask = element->mask;
  int *etype = element->etype;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  int *npe = element->npe;

  int flag;
  for (int i = 0; i < nelocal; i++) {
    if (emask[i] & groupbit) {
      elem_plane_id[i] = 0;
      flag = 0;
      for (int j = 0; j < npe[etype[i]]; j++) {
        elem_plane_id[i] += find_plane(nodex[i][j]);
        flag = MIN(flag,find_plane(nodex[i][j]));
      }
      if (flag < 0) elem_plane_id[i] = -1;
      else elem_plane_id[i] /= npe[etype[i]];
    }
  } 

  pending = 0;
}

/* ---------------------------------------------------------------------- */

FixAdaptive::~FixAdaptive()
{
  delete [] id_compute;
  delete [] splitflag;
  delete [] splitlist;
  delete [] disorder_flag;
  delete [] disorder_flagall;
  memory->destroy(atom_plane_id);
  memory->destroy(elem_plane_id);
}

/* ---------------------------------------------------------------------- */

int FixAdaptive::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive::init()
{

}

/* ---------------------------------------------------------------------- */

void FixAdaptive::setup(int vflag)
{
  //next_reneighbor = (update->ntimestep/nevery)*nevery;
}

void FixAdaptive::setup_pre_exchange()
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

 void FixAdaptive::pre_exchange()
{
  if (!pending) return;
  next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
  pending = 0;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive::post_integrate()
{
  // return if not an adaptive timestep

  if (update->ntimestep < next_reneighbor) return;

  // invoke compute peratom
  // set pre_neighbor_flag = 1;

  if (compute->invoked_peratom != update->ntimestep) {
    compute->preneighflag = 1;
    compute->compute_peratom();
  }

  // loop through my atoms
  // check for planes that need to split
  // skip if atoms not in any planes or 
  // plane has already been marked as disordered

  // reset flag 

  for (int i = 0; i < nplanes; i++) disorder_flag[i] = 0;

  int nlocal = atom->nlocal;
  double *pattern = compute->vector_atom;
  int ipattern;
  double iplane;
  for (int i = 0; i < nlocal; i++) {
    iplane = atom_plane_id[i];
    if (iplane < 0) continue;
    int iplane_int = static_cast<int> (iplane);
    if (disorder_flag[iplane_int]) continue;
    ipattern = static_cast<int> (pattern[i]);
    if (ipattern != perfect_crystal) 
      disorder_flag[iplane_int] = 1;
  }

  MPI_Allreduce(disorder_flag,disorder_flagall,nplanes,MPI_INT,MPI_MAX,world); 

  // check for interplanes that need splitting
  // if yes, split elements
  // skip if is already splitted

  int nsplit = 0;
  for (int i = 0; i < nplanes-1; i++) {
    if (disorder_flagall[i] && disorder_flagall[i+1] && !splitflag[i]) {
      splitflag[i] = 1;
      splitlist[nsplit++] = i;
    }
  }
  // return and update next reneighbor if no plane need to split

  if (nsplit == 0) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen,"No new dislocation found\n");
    return;
  }
  nlocal = element->nlocal;

  // first pass: check for elements that need splitted and new etypes needed,
  // stored list of indices of elements need splitted

  int nlist_split_all = 0;
  int *list_split = new int[nlocal];
  int nlist_split = 0;

  for (int i = 0; i < nlocal; i++) {
    iplane = elem_plane_id[i];
    if (iplane < 0) continue;
    for (int s = 0; s < nsplit; s++) {
      if (element->evec->split_element(i,split_dim,iplane,splitlist[s],0,NULL)) {
        if (i != list_split[nlist_split])
          list_split[nlist_split++] = i;
      }
    }
  }

  MPI_Allreduce(&nlist_split,&nlist_split_all,1,MPI_INT,MPI_SUM,world);
  // continue if there is at least one element need splitting in the whole domain,
  // set pending = 1 to update next_reneighbor after neighbor->decide() in pre_exchange() 
  // otherwise, update next_reneighbor now to skip force reneighboring at this time step

  if (!nlist_split_all) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen,"Dislocation found, no need to split\n");
  } else {
    pending = 1;
    int nalocal_old = atom->nlocal;
    // add new requested etypes 

    element->evec->add_requested_etype();

    // second pass: splitting elements in list

    for (int ii = 0; ii < nlist_split; ii++) {
      int i = list_split[ii];
      iplane = elem_plane_id[i];
      for (int s = 0; s < nsplit; s++) 
        element->evec->split_element(i,split_dim,iplane,splitlist[s],1,elem_plane_id);
    }

    // extend tags to new atoms and elements

    atom->tag_extend();
    element->tag_extend();

    int nnew_atoms = atom->nlocal - nalocal_old;
    int nnew_elems = element->nlocal - nlocal;
    int nnew_atoms_all,nnew_elems_all;
    MPI_Allreduce(&nnew_atoms,&nnew_atoms_all,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nnew_elems,&nnew_elems_all,1,MPI_INT,MPI_SUM,world);

    // update number of elements and number of atoms

    atom->natoms += nnew_atoms_all;
    element->nelements += nnew_elems_all;
    if (comm->me == 0) {
      if (screen) fprintf(screen,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nlist_split_all,nnew_atoms_all,nnew_elems_all+nlist_split_all);
      if (logfile) fprintf(logfile,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nlist_split_all,nnew_atoms_all,nnew_elems_all+nlist_split_all);
    }
  } 
}

/* ---------------------------------------------------------------------- */

int FixAdaptive::find_plane(double *x)
{
  if (x[0] < blo[0] || x[0] >= bhi[0] ||
      x[1] < blo[1] || x[1] >= bhi[1] ||
      x[2] < blo[2] || x[2] >= bhi[2]) 
    return -1;
  else
    return static_cast<int> ((x[split_dim] - blo[split_dim])/l);
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
   ------------------------------------------------------------------------- */

void FixAdaptive::set_atom_arrays(int i)
{
  atom_plane_id[i] = 0;
}

/* ----------------------------------------------------------------------
   initialize one element's storage values, called when element is created
   ------------------------------------------------------------------------- */

void FixAdaptive::set_elem_arrays(int i)
{
  elem_plane_id[i] = 0;
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values passed on from another element, called when atom is created from an element
   ------------------------------------------------------------------------- */

void FixAdaptive::pass_on_atom_arrays(int i, int i_old, int iintpl)
{
  atom_plane_id[i] = 0;
}

/* ----------------------------------------------------------------------
   initialize one element's storage values passed on from another element, called when element is created from another element
   element plane id is assigned right after created
   ------------------------------------------------------------------------- */

void FixAdaptive::pass_on_elem_arrays(int i, int i_old, int *iintpl_list)
{

}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
   ------------------------------------------------------------------------- */

void FixAdaptive::copy_atom_arrays(int i, int j, int /*delflag*/)
{
  atom_plane_id[j] = atom_plane_id[i];
}
/* ----------------------------------------------------------------------
   copy values within local elem-based array
   ------------------------------------------------------------------------- */

void FixAdaptive::copy_elem_arrays(int i, int j, int /*delflag*/)
{
  elem_plane_id[j] = elem_plane_id[i];
}

/* ----------------------------------------------------------------------
   allocate atom-based array
   ------------------------------------------------------------------------- */

void FixAdaptive::grow_atom_arrays(int nmax)
{
  memory->grow(atom_plane_id,nmax,"adaptive:atom_plane_id");
}
/* ----------------------------------------------------------------------
   allocate elem-based array
   ------------------------------------------------------------------------- */

void FixAdaptive::grow_elem_arrays(int nmax)
{
  memory->grow(elem_plane_id,nmax,"adaptive:elem_plane_id");
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
   ------------------------------------------------------------------------- */

int FixAdaptive::pack_atom_exchange(int i, double *buf)
{
  buf[0] = atom_plane_id[i];
  return 1;
}

/* ----------------------------------------------------------------------
   pack values in local elem-based array for exchange with another proc
   ------------------------------------------------------------------------- */

int FixAdaptive::pack_elem_exchange(int i, double *buf)
{
  buf[0] = elem_plane_id[i];
  return 1;
}


/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
   ------------------------------------------------------------------------- */

int FixAdaptive::unpack_atom_exchange(int nlocal, double *buf)
{
  atom_plane_id[nlocal] = buf[0]; 
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local elem-based array from exchange with another proc
   ------------------------------------------------------------------------- */

int FixAdaptive::unpack_elem_exchange(int nlocal, double *buf)
{
  elem_plane_id[nlocal] = buf[0]; 
  return 1;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

double FixAdaptive::memory_usage()
{
  return (atom->nmax + element->nmax) * sizeof(double);
}
