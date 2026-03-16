#include <string.h>
#include <stdlib.h>
#include "fix_addforce.h"
#include "atom.h"
#include "element.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};
enum{NODE,ELEMENT};

/* ---------------------------------------------------------------------- */

FixAddForce::FixAddForce(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg < 6) error->all(FLERR,"Illegal fix addforce command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  group_flag = NODE;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    xvalue = universe->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    yvalue = universe->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    zvalue = universe->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix addforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix addforce command");
      if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR, "Illegal fix addforce command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix addforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  //maxatom = 1;
  //memory->create(sforce,maxatom,3,"addforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixAddForce::~FixAddForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  //memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixAddForce::setmask()
{
  //datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  //mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }

  //if (estr) {
  //  evar = input->variable->find(estr);
  //  if (evar < 0)
  //    error->all(FLERR,"Variable name for fix addforce does not exist");
  //  if (input->variable->atomstyle(evar)) estyle = ATOM;
  //  else error->all(FLERR,"Variable for fix addforce is invalid style");
  //} else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix addforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    //varflag = ATOM;
    error->all(FLERR,"Atom style variable not ready yet");
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  //if (varflag == CONSTANT && estyle != NONE)
  //  error->all(FLERR,"Cannot use variable energy with "
  //      "constant force in fix addforce");
  //if ((varflag == EQUAL || varflag == ATOM) &&
  //    update->whichflag == 2 && estyle == NONE)
  //  error->all(FLERR,"Must use variable energy with fix addforce");

}

/* ---------------------------------------------------------------------- */

void FixAddForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) post_force(vflag);
  else error->all(FLERR,"Only verlet integration style is considered");
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  //if ((varflag == ATOM) && atom->nmax > maxatom) {
  //  maxatom = atom->nmax;
  //  memory->destroy(sforce);
  //  memory->create(sforce,maxatom,4,"addforce:sforce");
  //}

  // foriginal = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  double **ax = atom->x;
  double **af = atom->f;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodef = element->nodef;
  double **ex = element->x;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int npe = element->npe;
  int i,j;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {

    for (i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        if (region && !region->match(ax[i][0],ax[i][1],ax[i][2])) continue;

        //domain->unmap(x[i],image[i],unwrap);
        //foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
        foriginal[1] += af[i][0];
        foriginal[2] += af[i][1];
        foriginal[3] += af[i][2];
        af[i][0] += xvalue;
        af[i][1] += yvalue;
        af[i][2] += zvalue;
      }

    for (i = 0; i < nelocal; i++) 
      if (group_flag) {
        if (emask[i] & groupbit) {
          if (region && !region->match(ex[i][0],ex[i][1],ex[i][2])) continue;
          for (j = 0; j < npe; j++) {
            foriginal[1] += nodef[i][j][0];
            foriginal[2] += nodef[i][j][1];
            foriginal[3] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] += xvalue;
            if (ystyle) nodef[i][j][1] += yvalue;
            if (zstyle) nodef[i][j][2] += zvalue;
          }
        }
      } else {
        for (j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit) {
            if (region && !region->match(nodex[i][j][0],nodex[i][j][1],nodex[i][j][2])) continue;
            foriginal[1] += nodef[i][j][0];
            foriginal[2] += nodef[i][j][1];
            foriginal[3] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] += xvalue;
            if (ystyle) nodef[i][j][1] += yvalue;
            if (zstyle) nodef[i][j][2] += zvalue;
          }
      }
    

    // variable force, wrap with clear/add
    // potential energy = evar if defined, else 0.0
    // wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    //else if (xstyle == ATOM)
    //  input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    //else if (ystyle == ATOM)
    //  input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    //else if (zstyle == ATOM)
    //  input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        if (region && !region->match(ax[i][0],ax[i][1],ax[i][2])) continue;
        //if (estyle == ATOM) foriginal[0] += sforce[i][3];
        foriginal[1] += af[i][0];
        foriginal[2] += af[i][1];
        foriginal[3] += af[i][2];
        if (xstyle) af[i][0] += xvalue;
        if (ystyle) af[i][1] += yvalue;
        if (zstyle) af[i][2] += zvalue;
      }

    for (i = 0; i < nelocal; i++) 
      if (group_flag) {
        if (emask[i] & groupbit) {
          if (region && !region->match(ex[i][0],ex[i][1],ex[i][2])) continue;
          for (j = 0; j < npe; j++) {
            foriginal[1] += nodef[i][j][0];
            foriginal[2] += nodef[i][j][1];
            foriginal[3] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] += xvalue;
            if (ystyle) nodef[i][j][1] += yvalue;
            if (zstyle) nodef[i][j][2] += zvalue;
          }
        }
      } else {
        for (j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit) {
            if (region && !region->match(nodex[i][j][0],nodex[i][j][1],nodex[i][j][2])) continue;
            foriginal[1] += nodef[i][j][0];
            foriginal[2] += nodef[i][j][1];
            foriginal[3] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] += xvalue;
            if (ystyle) nodef[i][j][1] += yvalue;
            if (zstyle) nodef[i][j][2] += zvalue;
          }
      }
   
  }
}
/* ----------------------------------------------------------------------
   potential energy of added force
   ------------------------------------------------------------------------- */

double FixAddForce::compute_scalar()
{
  // only sum across procs one time

  //if (force_flag == 0) {
  //  MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
  //  force_flag = 1;
  //}
  //return foriginal_all[0];
  return 0.0;
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
   ------------------------------------------------------------------------- */

double FixAddForce::compute_vector(int n)
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

double FixAddForce::memory_usage()
{
  double bytes = 0.0;
  //if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
