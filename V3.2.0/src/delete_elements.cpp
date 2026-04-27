#include <stdlib.h>
#include <string.h>
#include "delete_elements.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "comm.h"
#include "domain.h"
#include "group.h"
#include "region.h"
#include "memory.h"
#include "irregular_comm.h"
#include "error.h"

#include <map>

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

DeleteElements::DeleteElements(CAC *cac) : Pointers(cac) 
{
  MPI_Comm_rank(world, &me);
}


/*  ----------------------------------------------------------------------  */

void DeleteElements::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Delete_elements command before simulation box is defined");
  if (narg < 1) error->all(FLERR, "Illegal delete_elements command");
  if (element->nelements == 0) {
    if (me == 0) 
      error->warning(FLERR, "No elements in simulation box");
    return;
  }
  if (element->tag_enable == 0)
    error->all(FLERR, "Cannot use delete_elements unless elements have IDs");

  // store state before delete

  bigint nelements_previous = element->nelements;

  // flag elements for deletion

  allflag = 0;

  if (strcmp(arg[0], "group") == 0) list_delete_group(narg, arg);
  else if (strcmp(arg[0], "region") == 0) list_delete_region(narg, arg);
  //else if (strcmp(arg[0], "overlap") == 0) delete_overlap(narg, arg);
  //else if (strcmp(arg[0], "porosity") == 0) delete_porosity(narg, arg);
  else error->all(FLERR, "Illegal delete_elements command");

  //if (allflag) {
  //  int igroup = group->find("all");
  //  if ((igroup >= 0) && modify->check_rigid_group_overlap(group->bitmask[igroup]))
  //    error->warning(FLERR, "Attempting to delete elements in rigid bodies");
  //} else {
  //  if (modify->check_rigid_list_overlap(del_flag_list))
  //    error->warning(FLERR, "Attempting to delete elements in rigid bodies");
  //}

  // if allflag = 1, just reset element->nlocal
  // else delete elements one by one


  if (allflag) {
    element->nlocal = 0;
    element->nucells = 0;
  } else {

    // delete local elements flagged in del_flag_list
    // reset nlocal

    int nlocal = element->nlocal;
    int *etype = element->etype;
    int *nucell = element->nucell;

    int i = 0;
    while (i < nlocal) {
      if (del_flag_list[i]) {
        nlocal--;
        element->evec->copy(nlocal, i, 1);
        del_flag_list[i] = del_flag_list[nlocal];
      } else i++;
    }

    element->nlocal = nlocal;
    memory->destroy(del_flag_list);
  }

  // if compress flag set, 
  // reset element tags to be contiguous
  // set all element IDs to 0, call tag_extend()

  if (compress_flag) {
    tagint *tag = element->tag;
    int nlocal = element->nlocal;
    for (int i = 0; i < nlocal; i++) tag[i] = 0;
    element->tag_extend();
  }


  bigint nblocal = element->nlocal;
  MPI_Allreduce(&nblocal, &element->nelements, 1, MPI_CAC_BIGINT, MPI_SUM, world);

  // recount # of nodes, node cells, vatoms, ucells
  
  element->count_nodes(1);
  element->count_node_cells(1);
  element->count_vatoms();
  element->count_ucells();

  // reset box and redistribute atoms/elements since box might change a lot
  // Invoke new irregular comm, migrate atoms/elements and destroy it
  // must be done after updating new global natoms

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal, atom->x);
    domain->x2lamda(element->nlocal, element->x);
    domain->nodex2lamda(element->nlocal, element->nodex);
  }
  domain->reset_box();
  IrregularComm *irrcomm = new IrregularComm(cac);
  irrcomm->migrate(1);
  delete irrcomm;
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal, atom->x);
    domain->lamda2x(element->nlocal, element->x);
    domain->lamda2nodex(element->nlocal, element->nodex);
  }


  // reset element->map if it exists
  // set nghost to 0 so old ghosts of deleted elements won't be mapped

  if (element->map_style) {
    element->nghost = 0;
    element->map_init();
    element->map_set();
  }

  // print before and after element and topology counts

  bigint ndelete = nelements_previous - element->nelements;

  if (me == 0) {
    if (screen) 
      fprintf(screen, "Deleted " BIGINT_FORMAT
          " elements, new total = " BIGINT_FORMAT "\n", 
          ndelete, element->nelements);
    if (logfile) 
      fprintf(logfile, "Deleted " BIGINT_FORMAT
          " elements, new total = " BIGINT_FORMAT "\n", 
          ndelete, element->nelements);
  }

}

/*  ----------------------------------------------------------------------
   flag all elements in group to be delete, group will still exist
   -------------------------------------------------------------------------  */

void DeleteElements::list_delete_group(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal delete_elements command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR, "Could not find delete_elements group ID");
  options(narg-2, &arg[2]);

  // check for special case of group = all

  if ((strcmp(arg[1], "all") == 0) || (strcmp(arg[1], "element") == 0)) {
    allflag = 1;
    return;
  }

  // allocate and initialize deletion list

  int nlocal = element->nlocal;
  memory->create(del_flag_list, nlocal, "delete_elements:del_flag_list");


  int *mask = element->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) del_flag_list[i] = 1;
    else del_flag_list[i] = 0;
  }
}

/*  ----------------------------------------------------------------------
   flag all elements in region to delete
   -------------------------------------------------------------------------  */

void DeleteElements::list_delete_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal delete_elements command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR, "Could not find delete_elements region ID");

  domain->regions[iregion]->prematch();
  options(narg-2, &arg[2]);

  // allocate and initialize deletion list

  int nlocal = element->nlocal;
  memory->create(del_flag_list, nlocal, "delete_elements:del_flag_list");

  double **x = element->x;

  for (int i = 0; i < nlocal; i++) {
    if (domain->regions[iregion]->match(x[i])) del_flag_list[i] = 1;
    else del_flag_list[i] = 0;
  }
}

/*  ----------------------------------------------------------------------
   delete elements so there are no pairs within cutoff
   which elements are deleted depends on ordering of elements within proc
   deletions can vary with processor count
   no guarantee that minimium number of elements will be deleted
   -------------------------------------------------------------------------  */
/* 
   void DeleteElements::delete_overlap(int narg, char **arg)
   {
   if (narg < 4) error->all(FLERR, "Illegal delete_elements command");

// read args

double cut = force->numeric(FLERR, arg[1]);
double cutsq = cut * cut;

int igroup1 = group->find(arg[2]);
int igroup2 = group->find(arg[3]);
if (igroup1 < 0 || igroup2 < 0)
error->all(FLERR, "Could not find delete_elements group ID");
options(narg-4, &arg[4]);

int group1bit = group->bitmask[igroup1];
int group2bit = group->bitmask[igroup2];

if (me == 0 && screen)
fprintf(screen, "System init for delete_elements ...\n");

// request a full neighbor list for use by this command

int irequest = neighbor->request(this);
neighbor->requests[irequest]->pair = 0;
neighbor->requests[irequest]->command = 1;
neighbor->requests[irequest]->half = 0;
neighbor->requests[irequest]->full = 1;
neighbor->requests[irequest]->occasional = 1;
neighbor->requests[irequest]->command_style = "delete_elements";

// init entire system since comm->borders and neighbor->build is done
// comm::init needs neighbor::init needs pair::init needs kspace::init, etc

cac->init();

// error check on cutoff
// if no pair style, neighbor list will be empty

if (force->pair == nullptr)
error->all(FLERR, "Delete_elements requires a pair style be defined");
if (cut > neighbor->cutneighmax)
error->all(FLERR, "Delete_elements cutoff > max neighbor cutoff");
if (cut > neighbor->cutneighmin && me == 0)
error->warning(FLERR, "Delete_elements cutoff > minimum neighbor cutoff");

// setup domain, communication and neighboring
// acquire ghosts and build standard neighbor lists

if (domain->triclinic) domain->x2lamda(element->nlocal);
domain->pbc();
domain->reset_box();
comm->setup();
if (neighbor->style) neighbor->setup_bins();
comm->exchange();
comm->borders();
if (domain->triclinic) domain->lamda2x(element->nlocal+element->nghost);
neighbor->build();

// build neighbor list this command needs based on earlier request

NeighList *list = neighbor->lists[irequest];
neighbor->build_one(list);

// allocate and initialize deletion list
// must be after exchange potentially changes nlocal

int nlocal = element->nlocal;
memory->create(del_flag_list, nlocal, "delete_elements:del_flag_list");
for (int i = 0; i < nlocal; i++) del_flag_list[i] = 0;

// double loop over owned elements and their full neighbor list
// at end of loop, there are no more overlaps
// only ever delete owned element I in I loop iteration, never J even if owned

tagint *tag = element->tag;
int *mask = element->mask;
double **x = element->x;
double *special_coul = force->special_coul;
double *special_lj = force->special_lj;

int i, j, ii, jj, inum, jnum;
double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
int *ilist, *jlist, *numneigh, **firstneigh;
double factor_lj, factor_coul;

inum = list->inum;
ilist = list->ilist;
numneigh = list->numneigh;
firstneigh = list->firstneigh;

for (ii = 0; ii < inum; ii++) {
  i = ilist[ii];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  jlist = firstneigh[i];
  jnum = numneigh[i];

  for (jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    factor_lj = special_lj[sbmask(j)];
    factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;

    // if both weighting factors are 0, skip this pair
    // could be 0 and still be in neigh list for long-range Coulombics
    // want consistency with non-charged pairs which wouldn't be in list

    if (factor_lj == 0.0 && factor_coul == 0.0) continue;

    // only consider deletion if I, J distance < cutoff
    // compute rsq identically on both I, J loop iterations
    // ignoring possibility that I, J tags are equal

    if (tag[i] < tag[j]) {
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
    } else {
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
    }
    rsq = delx * delx + dely * dely + delz * delz;
    if (rsq >= cutsq) continue;

    // only consider deletion if I, J are in groups 1, 2 respectively
    // true whether J is owned or ghost element

    if (!(mask[i] & group1bit)) continue;
    if (!(mask[j] & group2bit)) continue;

    // J is owned element:
    //   delete element I if element J has not already been deleted
    // J is ghost element:
    //   delete element I if J, I is not a candidate deletion pair
    //     due to being in groups 1, 2 respectively
    //   if they are candidate pair, then either:
    //      another proc owns J and could delete J
    //      J is a ghost of another of my owned elements, and I could delete J
    //   test on tags of I, J insures that only I or J is deleted

    if (j < nlocal) {
      if (del_flag_list[j]) continue;
    } else if ((mask[i] & group2bit) && (mask[j] & group1bit)) {
      if (tag[i] > tag[j]) continue;
    }

    del_flag_list[i] = 1;
    break;
  }
}
}
 */
/*  ----------------------------------------------------------------------
   create porosity by deleting elements in a specified region
   -------------------------------------------------------------------------  */
/* 
   void DeleteElements::delete_porosity(int narg, char **arg)
   {
   if (narg < 4) error->all(FLERR, "Illegal delete_elements command");

   int iregion = domain->find_region(arg[1]);
   if (iregion == -1) error->all(FLERR, "Could not find delete_elements region ID");
   domain->regions[iregion]->prematch();

   double porosity_fraction = force->numeric(FLERR, arg[2]);
   int seed = force->inumeric(FLERR, arg[3]);
   options(narg-4, &arg[4]);

   RanMars *random = new RanMars(cac, seed + me);

// allocate and initialize deletion list

int nlocal = element->nlocal;
memory->create(del_flag_list, nlocal, "delete_elements:del_flag_list");
for (int i = 0; i < nlocal; i++) del_flag_list[i] = 0;

double **x = element->x;

for (int i = 0; i < nlocal; i++)
if (domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2]))
if (random->uniform() <= porosity_fraction) del_flag_list[i] = 1;

delete random;
}

 */
/*  ----------------------------------------------------------------------
   process command options
   -------------------------------------------------------------------------  */

void DeleteElements::options(int narg, char **arg)
{
  compress_flag = 0;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal delete_elements command");
      if (strcmp(arg[iarg+1], "yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) compress_flag = 0;
      else error->all(FLERR, "Illegal delete_elements command");
      iarg += 2;
    } else error->all(FLERR, "Illegal delete_elements command");
  }
}

