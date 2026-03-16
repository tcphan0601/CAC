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


using namespace CAC_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};
enum{NODE,ELEMENT};

/*-------------------------------------------------------*/

FixSetForce::FixSetForce(CAC *cac, int narg, char **arg) :
  Fix(cac,narg,arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix setforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  group_flag = NODE;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  sforce = NULL;
  xstr = ystr = zstr = NULL;

  if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = universe->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }

  if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = universe->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }

  if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = universe->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix setforce command");
      if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR, "Illegal fix setforce command");
      iarg += 2;
    } else error->all(FLERR, "Illegal fix setforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  
  // constant for now, will add variable later 
  
  varflag = CONSTANT;

  maxatom = 1;
  if (maxatom) memory->create(sforce,maxatom,3,"setforce:sforce");
}

/*----------------------------------------------------------------------------*/

FixSetForce::~FixSetForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/*--------------------------------------------------------------------------*/

int FixSetForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  //mask |= MIN_POST_FORCE;
  return mask;
}

/*---------------------------------------------------------------------------*/

void FixSetForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) post_force(vflag);
  else error->all(FLERR,"Only verlet integration style is considered");

}	

/*----------------------------------------------------------------------------*/

void FixSetForce::post_force(int vflag)
{
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (atom->nlocal) {

    double **af = atom->f;
    int *amask = atom->mask;
    int nalocal = atom->nlocal;

    if (varflag == CONSTANT) 
      for (int i = 0; i < nalocal; i++)
        if (amask[i] & groupbit) {
          foriginal[0] += af[i][0];
          foriginal[1] += af[i][1];
          foriginal[2] += af[i][2];
          if (xstyle) af[i][0] = xvalue;
          if (ystyle) af[i][1] = yvalue;
          if (zstyle) af[i][2] = zvalue;
        }
  }

  if (element->nlocal) {
    double ***nodef = element->nodef;
    int *emask = element->mask;
    int *etype = element->etype;
    int **nodemask = element->nodemask;
    int nelocal = element->nlocal;
    int npe = element->npe;
    int i,j;

    if (varflag == CONSTANT) {
      for (int i = 0; i < nelocal; i++) {
        if (group_flag == NODE) {
          for (int j = 0; j < npe; j++) 
            if (nodemask[i][j] & groupbit) {
              foriginal[0] += nodef[i][j][0];
              foriginal[1] += nodef[i][j][1];
              foriginal[2] += nodef[i][j][2];
              if (xstyle) nodef[i][j][0] = xvalue;
              if (ystyle) nodef[i][j][1] = yvalue;
              if (zstyle) nodef[i][j][2] = zvalue;
            }
        } else if (emask[i] & groupbit) {
          for (int j = 0; j < npe; j++) {
            foriginal[0] += nodef[i][j][0];
            foriginal[1] += nodef[i][j][1];
            foriginal[2] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] = xvalue;
            if (ystyle) nodef[i][j][1] = yvalue;
            if (zstyle) nodef[i][j][2] = zvalue;
          }
        }
      }
    }
  }
}
/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
   ------------------------------------------------------------------------- */

double FixSetForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double FixSetForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom * 3 * sizeof(double);
  return bytes;
}
