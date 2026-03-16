#include <stdlib.h>
#include <string.h>
#include "delete_atoms.h"
#include "element.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "group.h"
#include "region.h"
#include "memory.h"
#include "irregular_comm.h"
#include "error.h"

#include <map>

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

DeleteAtoms::DeleteAtoms(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void DeleteAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Delete_atoms command before simulation box is defined");
  if (narg < 1) error->all(FLERR,"Illegal delete_atoms command");
  if (atom->natoms == 0) {
    if (comm->me == 0) 
      error->warning(FLERR,"No atoms in simulation box");
    return;
  }
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use delete_atoms unless atoms have IDs");

  // store state before delete

  bigint natoms_previous = atom->natoms;

  // flag atoms for deletion

  allflag = 0;

  if (strcmp(arg[0],"group") == 0) list_delete_group(narg,arg);
  else if (strcmp(arg[0],"region") == 0) list_delete_region(narg,arg);
  else error->all(FLERR,"Illegal delete_atoms command");

  // process option args
  
  options(narg-2,&arg[2]);

  // if allflag = 1, just reset atom->nlocal
  // else delete atoms one by one

  if (allflag) atom->nlocal = 0;
  else {

    // delete local atoms flagged in del_flag_list
    // reset nlocal

    AtomVec *avec = atom->avec;
    int nlocal = atom->nlocal;

    int i = 0;
    while (i < nlocal) {
      if (del_flag_list[i]) {
        avec->copy(nlocal-1,i,1);
        del_flag_list[i] = del_flag_list[nlocal-1];
        nlocal--;
      } else i++;
    }

    atom->nlocal = nlocal;
    memory->destroy(del_flag_list);
  }

  // if compress flag set,
  // reset atom tags to be contiguous
  // set all atom IDs to 0, call tag_extend()

  if (compress_flag) {
    tagint *tag = atom->tag;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) tag[i] = 0;
    atom->tag_extend();
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // reset box and redistribute atoms/elements since box might change a lot
  // Invoke new irregular comm, migrate atoms/elements and destroy it
  // must be done after updating new global natoms

  double **ax = atom->x;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

  double **ex = element->x;
  double ***nodex = element->nodex;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  for (int i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i]);

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

  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped


  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }


  bigint ndelete = natoms_previous - atom->natoms;

  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,"Deleted " BIGINT_FORMAT
          " atoms, new total = " BIGINT_FORMAT "\n",
          ndelete,atom->natoms);
    if (logfile) 
      fprintf(logfile,"Deleted " BIGINT_FORMAT
          " atoms, new total = " BIGINT_FORMAT "\n",
          ndelete,atom->natoms);
  }
}

/* ----------------------------------------------------------------------
   flag all atoms in group to delete, group will still exist
   ------------------------------------------------------------------------- */

void DeleteAtoms::list_delete_group(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal delete_atoms command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find delete_atoms group ID");

  // check for special case of group = all

  if ((strcmp(arg[1],"all") == 0) || (strcmp(arg[1],"atom") == 0)) {
    allflag = 1;
    return;
  }

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(del_flag_list,nlocal,"delete_atoms:del_flag_list");


  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) del_flag_list[i] = 1;
    else del_flag_list[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   flag all atoms in region to delete
   ------------------------------------------------------------------------- */

void DeleteAtoms::list_delete_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal delete_atoms command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Could not find delete_atoms region ID");
  domain->regions[iregion]->prematch();

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(del_flag_list,nlocal,"delete_atoms:del_flag_list");

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++) {
    if (domain->regions[iregion]->match(x[i])) del_flag_list[i] = 1;
    else del_flag_list[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   process command options
   ------------------------------------------------------------------------- */

void DeleteAtoms::options(int narg, char **arg)
{
  compress_flag = 1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal delete_atoms command");
      if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
      else error->all(FLERR,"Illegal delete_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal delete_atoms command");
  }
}
