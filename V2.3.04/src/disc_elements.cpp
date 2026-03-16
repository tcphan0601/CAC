#include <mpi.h>
#include <string.h>
#include "disc_elements.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "irregular_comm.h"

#include <map>

using namespace CAC_NS;
/* ---------------------------------------------------------------------- */

DiscElements::DiscElements(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
}

/* ----------------------------------------------------------------------
   Called as disc_elements command in input script.
   New atoms created will have ID starting immediately 
   following the largest atom ID existing before the
   disc_elements command was evoked.
   Call IrregularComm to migrate atoms since atoms might move far away
   ------------------------------------------------------------------------- */

void DiscElements::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Disc_elements command before simulation box is defined");

  if (element->nelements == 0) {
    if (me == 0) error->warning(FLERR,"No element in simulation box");
    return;
  }

  if (narg < 1) error->all(FLERR,"Illegal dis_elements command");

  // flag elements for discretizing
  
  if (strcmp(arg[0],"region") == 0) list_disc_region(narg,arg);
  else if (strcmp(arg[0],"group") == 0) list_disc_group(narg,arg);
  else error->all(FLERR,"Illegal disc_elements command");

  atom->init();
  element->init(); 

  int i = 0;
  int nnewatoms = 0;
  int ndel = 0;
  ElementVec *evec = element->evec;

  // discretize elements in list

  while (i < element->nlocal)
    if (disc_flag_list[i]) {
      ndel++;
      nnewatoms += evec->element2atom(i);
      disc_flag_list[i] = disc_flag_list[element->nlocal];
    } else i++;

  memory->destroy(disc_flag_list);

 
  // assign tag to new atoms
  // if compress flag set,
  // reset element tags to be contiguous
  // set all element IDs to 0, call tag_extend()
  
  atom->tag_extend();
  if (compress_flag) {
    tagint *tag = element->tag;
    int nlocal = element->nlocal;
    for (int i = 0; i < nlocal; i++) tag[i] = 0;
    element->tag_extend();
  }
 
  // update new global natoms and nelements

  int totalnewatoms;
  int totaldel;
  MPI_Allreduce(&nnewatoms,&totalnewatoms,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ndel,&totaldel,1,MPI_INT,MPI_SUM,world);
  atom->natoms += totalnewatoms;
  element->nelements -= totaldel;
  element->nintpls -= totalnewatoms;

  // Invoke new irregular comm, migrate atoms/elements and destroy it
  // must be done after updating new global natoms and nelements
  
  
  double **ax = atom->x;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

  double **ex = element->x;
  double ***nodex = element->nodex;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  for (i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i]);

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
  }
  domain->reset_box();


  IrregularComm *irrcomm = new IrregularComm(cac);
  irrcomm->migrate(1);
  delete irrcomm;

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
  }


  // reset element->map and atom->map if it exists
  // set nghost to 0 so old ghosts of deleted elements and new atoms won't be mapped

  if (element->map_style) {
    element->nghost = 0;
    element->map_init();
    element->map_set();
  }

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  FILE *out;
  if (me == 0) 
    for (int m = 0; m < 2; m++) {
      if (m == 0) out = logfile; 
      else out = screen;
      if (out) {
        fprintf(out,"Discretized %d elements into %d atoms\n",totaldel,totalnewatoms);
        fprintf(out,"  New total atoms = %d\n",atom->natoms);
        fprintf(out,"  New total elements = %d\n",element->nelements);
      }
    }
}

/* ----------------------------------------------------------------------
   flag all elements in group to be discretize, group will still exist
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_group(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find disc_elements group ID");

  options(narg-2,&arg[2]);

  // allocate discretization list

  int nlocal = element->nlocal;
  memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");


  int *mask = element->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) disc_flag_list[i] = 1;
    else disc_flag_list[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   flag all elements in region to discretize
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Could not find disc_elements region ID");
  domain->regions[iregion]->prematch();

  options(narg-2,&arg[2]);

  // allocate discretization list

  int nlocal = element->nlocal;
  memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");

  double **x = element->x;

  for (int i = 0; i < nlocal; i++) {
    if (domain->regions[iregion]->match(x[i])) disc_flag_list[i] = 1;
    else disc_flag_list[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   process command options
   ------------------------------------------------------------------------- */

void DiscElements::options(int narg, char **arg)
{
  compress_flag = 1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
      if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
      else error->all(FLERR,"Illegal disc_elements command");
      iarg += 2;
    } else error->all(FLERR,"Illegal disc_elements command");
  }
}

