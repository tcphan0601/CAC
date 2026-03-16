#include <string.h>
#include <stdlib.h>
#include "fix_setforce.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "variable.h"
#include "group.h"


using namespace CAC_NS;
using namespace FixConst;

enum{NONE, CONSTANT, EQUAL, ATOM};

/* ------------------------------------------------------- */

FixSetForce::FixSetForce(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  xstr(nullptr), ystr(nullptr), zstr(nullptr), idregion(nullptr)//, sforce(nullptr)
{

  if (narg < 6) error->all(FLERR, "Illegal fix setforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  //sforce = nullptr;
  //nodesforce = nullptr;
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
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix setforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR, "Region ID for fix setforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion, arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR, "Illegal fix setforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  

  //maxatom = 1;
  //memory->create(sforce, maxatom, 3, "setforce:sforce");
  //maxelem = 1;
  //memory->create(nodesforce, maxelem, element->npe, 3, "setforce:sforce");
}

/* ---------------------------------------------------------------------------- */

FixSetForce::~FixSetForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  //memory->destroy(sforce);
}

/* -------------------------------------------------------------------------- */

int FixSetForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  //mask |= MIN_POST_FORCE;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixSetForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR, "Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR, "Variable for fix setforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR, "Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR, "Variable for fix setforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR, "Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR, "Variable for fix setforce is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR, "Region ID for fix setforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    error->all(FLERR, "Atom style variable not ready yet");
    //varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  //if (strstr(update->integrate_style, "respa")) {
  //  nlevels_respa = ((Respa *) update->integrate)->nlevels;
  //  if (respa_level >= 0) ilevel_respa = MIN(respa_level, nlevels_respa-1);
  //  else ilevel_respa = nlevels_respa-1;
  //}

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  //int flag = 0;
  //if (update->whichflag == 2) {
  //  if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
  //  if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
  //  if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
  //  if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
  //  if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
  //  if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  //}
  //if (flag)
  //  error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* --------------------------------------------------------------------------- */

void FixSetForce::setup(int vflag)
{
  if (strstr(update->integrate_style, "verlet")) post_force(vflag);
  else error->all(FLERR, "Only verlet integration style is considered");
}	

/* ---------------------------------------------------------------------------- */

void FixSetForce::post_force(int vflag)
{
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  // update region if necessary

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary
  
  //if (varflag == ATOM && atom->nmax > maxatom) {
  //  maxatom = atom->nmax;
  //  memory->destroy(sforce);
  //  memory->create(sforce, maxatom, 3, "setforce:sforce");
  //}

  //if (varflag == ATOM && element->nmax > maxelem) {
  //  maxelem = element->nmax;
  //  memory->destroy(nodesforce);
  //  memory->create(nodesforce, maxelem, element->npe, 3, "setforce:sforce");
  //}

  if (varflag == EQUAL) {
    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    //else if (xstyle == ATOM)
    //  input->variable->compute_atom(xvar, igroup, &sforce[0][0], 3, 0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    //else if (ystyle == ATOM)
    //  input->variable->compute_atom(yvar, igroup, &sforce[0][1], 3, 0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    //else if (zstyle == ATOM)
    //  input->variable->compute_atom(zvar, igroup, &sforce[0][2], 3, 0);

    modify->addstep_compute(update->ntimestep + 1);
  }


  double **ax = atom->x;
  double **af = atom->f;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;

  for (int i = 0; i < nalocal; i++) 
    if (amask[i] & groupbit) {
      if (region && !region->match(ax[i][0], ax[i][1], ax[i][2])) continue;
      foriginal[0] += af[i][0];
      foriginal[1] += af[i][1];
      foriginal[2] += af[i][2];
      if (xstyle) af[i][0] = xvalue;
      if (ystyle) af[i][1] = yvalue;
      if (zstyle) af[i][2] = zvalue;
    }

  double **x = element->x;
  double ****nodex = element->nodex;
  double ****nodef = element->nodef;
  int *emask = element->mask;
  int *etype = element->etype;
  int ***nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;

  if (group_flag == Group::NODE) {
    for (int i = 0; i < nelocal; i++) 
      for (int j = 0; j < apc[etype[i]]; j++) 
        for (int k = 0; k < npe[etype[i]]; k++) 
          if (nodemask[i][j][k] & groupbit) {
            if (region && !region->match(nodex[i][j][k][0], 
                  nodex[i][j][k][1], nodex[i][j][k][2])) continue;
            foriginal[0] += nodef[i][j][k][0];
            foriginal[1] += nodef[i][j][k][1];
            foriginal[2] += nodef[i][j][k][2];
            if (xstyle) nodef[i][j][k][0] = xvalue;
            if (ystyle) nodef[i][j][k][1] = yvalue;
            if (zstyle) nodef[i][j][k][2] = zvalue;
          }

  } else if (group_flag == Group::ELEMENT) {
    for (int i = 0; i < nelocal; i++) {
      if (!(emask[i] & groupbit)) continue;
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      for (int j = 0; j < apc[etype[i]]; j++) 
        for (int k = 0; k < npe[etype[i]]; k++) {
          foriginal[0] += nodef[i][j][k][0];
          foriginal[1] += nodef[i][j][k][1];
          foriginal[2] += nodef[i][j][k][2];
          if (xstyle) nodef[i][j][k][0] = xvalue;
          if (ystyle) nodef[i][j][k][1] = yvalue;
          if (zstyle) nodef[i][j][k][2] = zvalue;
        }
    }
  }
}

/*  ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
   -------------------------------------------------------------------------  */

double FixSetForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/*  ----------------------------------------------------------------------
   memory usage of local atom-based array
   -------------------------------------------------------------------------  */

double FixSetForce::memory_usage()
{
  double bytes = 0.0;
  //if (varflag == ATOM) bytes = (maxatom + maxelem * element->npe) * 3 * sizeof(double);
  return bytes;
}
