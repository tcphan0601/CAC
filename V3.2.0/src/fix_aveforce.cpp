#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "fix_aveforce.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "force.h"
#include "pair.h"
//#include "respa.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NONE, CONSTANT, EQUAL};

/*  ----------------------------------------------------------------------  */

FixAveForce::FixAveForce(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  xstr(nullptr), ystr(nullptr), zstr(nullptr), idregion(nullptr)
{
  if (narg < 6) error->all(FLERR, "Illegal fix aveforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 1;
  //respa_level_support = 1;
  //ilevel_respa = nlevels_respa = 0;

  xstr = ystr = zstr = nullptr;

  if (strstr(arg[3], "v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr, &arg[3][2]);
  } else if (strcmp(arg[3], "NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = universe->numeric(FLERR, arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4], "v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr, &arg[4][2]);
  } else if (strcmp(arg[4], "NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = universe->numeric(FLERR, arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5], "v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr, &arg[5][2]);
  } else if (strcmp(arg[5], "NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = universe->numeric(FLERR, arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = nullptr;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix setforce command");
      if (strcmp(arg[iarg+1], "atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1], "node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1], "element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR, "Illegal fix setforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "region") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix aveforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR, "Region ID for fix aveforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion, arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR, "Illegal fix aveforce command");

  }

  foriginal_all[0] = foriginal_all[1] = 0.0;
  foriginal_all[2] = foriginal_all[3] = 0.0;
}

/*  ----------------------------------------------------------------------  */

FixAveForce::~FixAveForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
}

/*  ----------------------------------------------------------------------  */

int FixAveForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  //mask |= POST_FORCE_RESPA;
  //mask |= MIN_POST_FORCE;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixAveForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR, "Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else error->all(FLERR, "Variable for fix aveforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR, "Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else error->all(FLERR, "Variable for fix aveforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR, "Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else error->all(FLERR, "Variable for fix aveforce is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR, "Region ID for fix aveforce does not exist");
  }

  if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL) varflag = EQUAL;
  else varflag = CONSTANT;

  // check if assign force to node individually in group with atoms, might not be physically correct
  
  if (group_flag == Group::NODE)
    if (group->count_atom(igroup) && group->count_node(igroup))
      error->all(FLERR, "Cannot use fix ave on group contain both atoms and nodes for group node option");

  //if (strstr(update->integrate_style, "respa")) {
  //  nlevels_respa = ((Respa *) update->integrate)->nlevels;
  //  if (respa_level >= 0) ilevel_respa = MIN(respa_level, nlevels_respa-1);
  //  else ilevel_respa = nlevels_respa-1;
  //}
}

/*  ----------------------------------------------------------------------  */

void FixAveForce::setup(int vflag)
{
  if (strstr(update->integrate_style, "verlet"))
    post_force(vflag);
  //else
  //  for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
  //    ((Respa *) update->integrate)->copy_flevel_f(ilevel);
  //    post_force_respa(vflag, ilevel, 0);
  //    ((Respa *) update->integrate)->copy_f_flevel(ilevel);
  //  }
}

/*  ----------------------------------------------------------------------  */

void FixAveForce::min_setup(int vflag)
{
  post_force(vflag);
}

/*  ----------------------------------------------------------------------  */

void FixAveForce::post_force(int /* vflag */)
{
  // update region if necessary

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  int i, j, k;

  // sum forces on participating atoms and nodes

  double **ax = atom->x;
  double **af = atom->f;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;

  double foriginal[4];
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  for (i = 0; i < nalocal; i++)
    if (amask[i] & groupbit) {
      if (region && !region->match(ax[i][0], ax[i][1], ax[i][2])) continue;
      foriginal[0] += af[i][0];
      foriginal[1] += af[i][1];
      foriginal[2] += af[i][2];
      foriginal[3] += 1.0;
    }

  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  int *emask = element->mask;
  int ***nodemask = element->nodemask;
  int *etype = element->etype;
  int *nucell = element->nucell;
  double **ex = element->x;
  double ****nodex = element->nodex;
  double ****nodef = element->nodef;
  double *nodal_weight = element->nodal_weight;
  double fxtmp, fytmp, fztmp;

  if (group_flag == Group::ELEMENT) {
    for (i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) {
        if (region && !region->match(ex[i][0], ex[i][1], ex[i][2])) continue;
        int ietype = etype[i];
        fxtmp = fytmp = fztmp = 0.0;
        for (j = 0; j < apc[ietype]; j++) 
          for (k = 0; k < npe[ietype]; k++) {
            fxtmp += nodef[i][j][k][0];
            fytmp += nodef[i][j][k][1];
            fztmp += nodef[i][j][k][2];
          }

        foriginal[0] += fxtmp * nodal_weight[ietype];
        foriginal[1] += fytmp * nodal_weight[ietype];
        foriginal[2] += fztmp * nodal_weight[ietype];
        foriginal[3] += nucell[ietype] * apc[ietype];
      }
  } else if (group_flag == Group::NODE) {
    for (i = 0; i < nelocal; i++) {
      int ietype = etype[i];
      for (j = 0; j < apc[ietype]; j++) 
        for (k = 0; k < npe[ietype]; k++) {
          if (nodemask[i][j][k] & groupbit) {
            if (region && 
                !region->match(nodex[i][j][k][0], nodex[i][j][k][1], nodex[i][j][k][2])) 
              continue;
            foriginal[0] += nodef[i][j][k][0] * nodal_weight[ietype];
            foriginal[1] += nodef[i][j][k][1] * nodal_weight[ietype];
            foriginal[2] += nodef[i][j][k][2] * nodal_weight[ietype];
            foriginal[3] += nodal_weight[ietype];
          }
        }
    }
  }

  // average the force on participating atoms and nodes
  // add in requested amount, computed via variable evaluation if necessary
  // wrap variable evaluation with clear/add

  MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);

  int ncount = static_cast<int> (foriginal_all[3]);
  if (ncount == 0) return;

  if (varflag == EQUAL) {
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  double fave[3];
  fave[0] = foriginal_all[0]/ncount + xvalue;
  fave[1] = foriginal_all[1]/ncount + yvalue;
  fave[2] = foriginal_all[2]/ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

  for (i = 0; i < nalocal; i++) {
    if (amask[i] & groupbit) {
      if (region && !region->match(ax[i][0], ax[i][1], ax[i][2])) continue;
      if (xstyle) af[i][0] = fave[0];
      if (ystyle) af[i][1] = fave[1];
      if (zstyle) af[i][2] = fave[2];
    }
  }

  if (group_flag == Group::ELEMENT) {
    for (i = 0; i < nelocal; i++) {
      if (emask[i] & groupbit) {
        if (region && !region->match(ex[i][0], ex[i][1], ex[i][2])) continue;
        for (j = 0; j < apc[etype[i]]; j++) 
          for (k = 0; k < npe[etype[i]]; k++) {
            if (xstyle) nodef[i][j][k][0] = fave[0];
            if (ystyle) nodef[i][j][k][1] = fave[1];
            if (zstyle) nodef[i][j][k][2] = fave[2];
          }
      }
    } 
  } else if (group_flag == Group::NODE) {
    for (i = 0; i < nelocal; i++) {
      for (j = 0; j < apc[etype[i]]; j++) 
        for (k = 0; k < npe[etype[i]]; k++) {
          if (nodemask[i][j][k] & groupbit) {
            if (region && 
                !region->match(nodex[i][j][k][0], nodex[i][j][k][1], nodex[i][j][k][2])) 
              continue;
            if (xstyle) nodef[i][j][k][0] = fave[0];
            if (ystyle) nodef[i][j][k][1] = fave[1];
            if (zstyle) nodef[i][j][k][2] = fave[2];
          }
        }
    }
  }
}

/*  ----------------------------------------------------------------------  */

//void FixAveForce::post_force_respa(int vflag, int ilevel, int /* iloop */)
//{
//  // ave + extra force on selected RESPA level
//  // just ave on all other levels
//
//  if (ilevel == ilevel_respa) post_force(vflag);
//  else {
//    Region *region = nullptr;
//    if (iregion >= 0) {
//      region = domain->regions[iregion];
//      region->prematch();
//    }
//
//    double **x = atom->x;
//    double **f = atom->f;
//    int *mask = atom->mask;
//    int nlocal = atom->nlocal;
//
//    double foriginal[4];
//    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
//
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
//        foriginal[0] += f[i][0];
//        foriginal[1] += f[i][1];
//        foriginal[2] += f[i][2];
//        foriginal[3] += 1.0;
//      }
//
//    MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);
//
//    int ncount = static_cast<int> (foriginal_all[3]);
//    if (ncount == 0) return;
//
//    double fave[3];
//    fave[0] = foriginal_all[0]/ncount;
//    fave[1] = foriginal_all[1]/ncount;
//    fave[2] = foriginal_all[2]/ncount;
//
//    for (int i = 0; i < nlocal; i++)
//      if (mask[i] & groupbit) {
//        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
//        if (xstyle) f[i][0] = fave[0];
//        if (ystyle) f[i][1] = fave[1];
//        if (zstyle) f[i][2] = fave[2];
//      }
//  }
//}

/*  ----------------------------------------------------------------------  */

//void FixAveForce::min_post_force(int vflag)
//{
//  post_force(vflag);
//}

/*  ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
   -------------------------------------------------------------------------  */

double FixAveForce::compute_vector(int n)
{
  return foriginal_all[n];
}
