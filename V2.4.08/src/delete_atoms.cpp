#include <stdlib.h>
#include <string.h>
#include "delete_atoms.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "region.h"
#include "memory.h"
#include "universe.h"
#include "irregular_comm.h"
#include "random_mars.h"
#include "error.h"

#include <map>

using namespace CAC_NS;

enum { UNKNOWN, FRACTION, COUNT };

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
  else if (strcmp(arg[0],"overlap") == 0) list_delete_overlap(narg,arg);
  else if (strcmp(arg[0],"random") == 0) list_delete_random(narg,arg);
  else error->all(FLERR,"Illegal delete_atoms command");


  //fprintf(screen,"before delete me = %d natoms = %d nlocal = %d\n",comm->me,atom->natoms,atom->nlocal);
  // if allflag = 1, just reset atom->nlocal
  // else delete atoms one by one

  if (allflag) atom->nlocal = 0;
  else {

    // delete local atoms flagged in dlist
    // reset nlocal

    AtomVec *avec = atom->avec;
    int nlocal = atom->nlocal;

    int i = 0;
    while (i < nlocal) {
      if (dlist[i]) {
        avec->copy(nlocal-1,i,1);
        dlist[i] = dlist[nlocal-1];
        nlocal--;
      } else i++;
    }

    atom->nlocal = nlocal;
    memory->destroy(dlist);
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



  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->reset_box();
  IrregularComm *irrcomm = new IrregularComm(cac);
  irrcomm->migrate(1);
  delete irrcomm;
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2nodex(element->nlocal,element->nodex);
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

  // process option args
  
  options(narg-2,&arg[2]);

  // check for special case of group = all

  if ((strcmp(arg[1],"all") == 0) || (strcmp(arg[1],"atom") == 0)) {
    allflag = 1;
    return;
  }

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");


  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) dlist[i] = 1;
    else dlist[i] = 0;
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

  // process option args
  
  options(narg-2,&arg[2]);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");

  double **x = atom->x;

  for (int i = 0; i < nlocal; i++) {
    if (domain->regions[iregion]->match(x[i])) dlist[i] = 1;
    else dlist[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   flag atoms to delete so that there are no pairs within cutoff
   which atoms are deleted depends on ordering of atoms within proc
   deletions can vary with processor count
   no guarantee that minimium number of atoms will be deleted
------------------------------------------------------------------------- */

void DeleteAtoms::list_delete_overlap(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR,"Illegal delete_atoms command");

  // read args

  double cut = universe->numeric(FLERR,arg[1]);
  double cutsq = cut*cut;

  int igroup1 = group->find(arg[2]);
  int igroup2 = group->find(arg[3]);
  if (igroup1 < 0 || igroup2 < 0)
    error->all(FLERR,"Could not find delete_atoms group ID");
  options(narg-4,&arg[4]);

  int group1bit = group->bitmask[igroup1];
  int group2bit = group->bitmask[igroup2];

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for delete_atoms ...\n");

  // request a full neighbor list for use by this command

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->elemlist = 0;
  neighbor->requests[irequest]->intglist = 0;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->command_style = "delete_atoms";

  // if no element in both groups, request atom-atom pair neighbor list only 
  
  if (group->count_elem(igroup1) == 0 || group->count_elem(igroup2) == 0) {
    neighbor->requests[irequest]->atomonlylist = 1;
    neighbor->requests[irequest]->atomlist = 0;
  }

  // init entire system since comm->borders and neighbor->build is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  cac->init();

  // error check on cutoff
  // if no pair style, neighbor list will be empty

  if (force->pair == NULL)
    error->all(FLERR,"Delete_atoms requires a pair style be defined");
  if (cut > neighbor->cutneighmax)
    error->all(FLERR,"Delete_atoms cutoff > max neighbor cutoff");
  if (cut > neighbor->cutneighmin && comm->me == 0)
    error->warning(FLERR,"Delete_atoms cutoff > minimum neighbor cutoff");

  // setup domain, communication and neighboring
  // acquire ghosts and build standard neighbor lists

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
  if (element->element_cluster_flag && element->max_apc > 1) 
    element->check_element_clusters();
  neighbor->build(1);

  // build neighbor list this command needs based on earlier request

  NeighList *list = neighbor->lists[irequest];
  neighbor->build_one(list);
 
  // allocate and initialize deletion list
  // must be after exchange potentially changes nlocal

  int nlocal = atom->nlocal;
  memory->create(dlist,nlocal,"delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // double loop over owned atoms and their full neighbor list
  // at end of loop, there are no more overlaps
  // only ever delete owned atom I in I loop iteration, never J even if owned

  tagint *tag = atom->tag;
  int *amask = atom->mask;
  double **x = atom->x;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int *emask = element->mask;
  int *etype = element->etype;
  double ***nodex = element->nodex;
  double coord[3];
   
  int i,j,ii,jj,inum,jnum,jintpl,jetype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*jlist_index,*numneigha2a,*numneigha2ia;
  int **firstneigha2a,**firstneigha2ia,**firstneigha2ia_index;
  double factor_lj,factor_coul;

  inum = list->ainum;
  ilist = list->ailist;
  numneigha2a = list->numneigha2a;
  numneigha2ia = list->numneigha2ia;
  firstneigha2a = list->firstneigha2a;
  firstneigha2ia = list->firstneigha2ia;
  firstneigha2ia_index = list->firstneigha2ia_index;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigha2a[i];
    jnum = numneigha2a[i];
    if (jnum) {
      jlist = firstneigha2a[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        factor_lj = special_lj[sbmask(j)];
        factor_coul = special_coul[sbmask(j)];
        j &= NEIGHMASK;

        // only consider deletion if I,J are in groups 1,2 respectively
        // true whether J is owned or ghost atom

        if (!(amask[i] & group1bit)) continue;
        if (!(amask[j] & group2bit)) continue;

        // if both weighting factors are 0, skip this pair
        // could be 0 and still be in neigh list for long-range Coulombics
        // want consistency with non-charged pairs which wouldn't be in list

        //if (factor_lj == 0.0 && factor_coul == 0.0) continue;

        // only consider deletion if I,J distance < cutoff
        // compute rsq identically on both I,J loop iterations
        // ignoring possibility that I,J tags are equal

        if (tag[i] < tag[j]) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
        } else {
          delx = x[j][0] - xtmp;
          dely = x[j][1] - ytmp;
          delz = x[j][2] - ztmp;
        }
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq >= cutsq) continue;

        // J is owned atom:
        //   delete atom I if atom J has not already been deleted
        // J is ghost atom:
        //   delete atom I if J,I is not a candidate deletion pair
        //     due to being in groups 1,2 respectively
        //   if they are candidate pair, then either:
        //      another proc owns J and could delete J
        //      J is a ghost of another of my owned atoms, and I could delete J
        //   test on tags of I,J insures that only I or J is deleted

        if (j < nlocal) {
          if (dlist[j]) continue;
        } else if ((amask[i] & group2bit) && (amask[j] & group1bit)) {
          if (tag[i] > tag[j]) continue;
        }

        dlist[i] = 1;
        break;
      }
    }

    if (dlist[i]) continue;

    jnum = numneigha2ia[i];
    if (jnum) {
      jlist = firstneigha2ia[i];
      jlist_index = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jlist_index[jj];
        jetype = etype[j];

        // only consider deletion if I is in groups 1 and J is in either group
       

        if (!(amask[i] & group1bit)) continue;
        if (!(emask[j] & group1bit)) continue;
        if (!(emask[j] & group2bit)) continue;

        //factor_lj = special_lj[sbmask(j)];
        //factor_coul = special_coul[sbmask(j)];
        //j &= NEIGHMASK;

        // if both weighting factors are 0, skip this pair
        // could be 0 and still be in neigh list for long-range Coulombics
        // want consistency with non-charged pairs which wouldn't be in list

        //if (factor_lj == 0.0 && factor_coul == 0.0) continue;

        // only consider deletion if I,J distance < cutoff

        delx = xtmp - element->evec->interpolate(nodex,j,jintpl,0);
        dely = ytmp - element->evec->interpolate(nodex,j,jintpl,1);
        delz = ztmp - element->evec->interpolate(nodex,j,jintpl,2);

        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq >= cutsq) continue;

        dlist[i] = 1;

        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   delete specified portion of atoms within group and/or region
------------------------------------------------------------------------- */

void DeleteAtoms::list_delete_random(int narg, char **arg)
{
  if (narg < 7) error->all(FLERR,"Illegal delete_atoms command");

  int random_style = UNKNOWN;
  int exactflag = 0;
  int errorflag = 0;
  bigint ncount = 0;
  double fraction = 0.0;

  if (strcmp(arg[1], "fraction") == 0) {
    random_style = FRACTION;
    fraction = universe->numeric(FLERR, arg[2]);
    if (strcmp(arg[3],"yes") == 0) exactflag = 1;
    else if (strcmp(arg[3],"no") == 0) exactflag = 0;

    if (fraction < 0.0 || fraction > 1.0)
      error->all(FLERR, "Delete_atoms random fraction has invalid value");
  } else if (strcmp(arg[1], "count") == 0) {
    random_style = COUNT;
    ncount = universe->bnumeric(FLERR, arg[2]);
    if (strcmp(arg[3],"yes") == 0) exactflag = 1;
    else if (strcmp(arg[3],"no") == 0) exactflag = 0;
    if (ncount < 0) error->all(FLERR, "Delete_atoms random count has invalid value");
    exactflag = true;
  } else {
    error->all(FLERR, "Unknown delete_atoms random style");
  }

  int igroup = group->find(arg[4]);
  if (igroup == -1) error->all(FLERR, "Could not find delete_atoms random group ID");

  int iregion = domain->find_region(arg[5]);
  if (iregion == -1 && (strcmp(arg[5], "NULL") != 0)) error->all(FLERR,"Could not find delete_atoms region ID");
  if (iregion >= 0)
    domain->regions[iregion]->prematch();

  int seed = universe->inumeric(FLERR, arg[6]);
  options(narg - 7, &arg[7]);

  RanMars *ranmars = new RanMars(cac, seed + comm->me);

  // allocate and initialize deletion list

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  // setup

  double **x = atom->x;
  int *mask = atom->mask;

  int groupbit = group->bitmask[igroup];

  // delete approximate fraction of atoms in both group and region

  if (random_style == FRACTION && !exactflag) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (iregion >= 0 && !domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2])) continue;
      if (ranmars->uniform() <= fraction) dlist[i] = 1;
    }

    // delete exact fraction or count of atoms in both group and region

  } else {
    double **x = atom->x;
    int *mask = atom->mask;

    // count = number of atoms this proc owns in both group and region

    int count = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (iregion >= 0 && !domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2])) continue;
      count++;
    }

    // convert specified fraction to ncount

    bigint bcount = count;
    bigint allcount;
    MPI_Allreduce(&bcount, &allcount, 1, MPI_CAC_BIGINT, MPI_SUM, world);

    if (random_style == FRACTION) {
      ncount = static_cast<bigint>(fraction * allcount);
    } else if (random_style == COUNT) {
      if (ncount > allcount) {
        if (errorflag) {
          error->all(FLERR, "Delete_atoms count exceeds number of eligible atoms");
        } else {
          if (comm->me == 0) 
            error->warning(FLERR, "Delete_atoms count exceeds number of eligible atoms");
        }
      }
    }

    // make selection

    int *flag = memory->create(flag, count, "delete_atoms:flag");
    int *work = memory->create(work, count, "delete_atoms:work");

    ranmars->select_subset(ncount, count, flag, work);

    // set dlist for atom indices in flag
    // flag vector from select_subset() is only for eligible atoms

    int j = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (iregion >= 0 && !domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2])) continue;
      if (flag[j]) dlist[i] = 1;
      j++;
    }

    memory->destroy(flag);
    memory->destroy(work);
  }

  // delete RN generator

  delete ranmars;
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
  // if more than one atom per node cluster, no compress element tags 

  if (element->element_cluster_flag) compress_flag = 0;

}
