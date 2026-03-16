#include <cstring>
#include <cstdlib>
#include <cmath>
#include "fix_temp_rescale.h"
#include "atom.h"
#include "element.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NOBIAS, BIAS};
enum{CONSTANT, EQUAL};

/*  ----------------------------------------------------------------------  */

FixTempRescale::FixTempRescale(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  tstr(NULL), id_temp(NULL), tflag(0)
{
  if (narg < 8) error->all(FLERR, "Illegal fix temp/rescale command");

  nevery = universe->inumeric(FLERR, arg[3]);
  if (nevery <= 0) error->all(FLERR, "Illegal fix temp/rescale command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;
  dynamic_group_allow = 1;

  // optional args

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix temp/rescale command");
      if (strcmp(arg[iarg+1], "atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1], "node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1], "element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR, "Illegal fix temp/rescale command");
      iarg += 2;
    } else error->all(FLERR, "Illegal fix temp/rescale command");
  }
  if (group_flag == Group::NODE && element->nodal_force_style == Element::CONSISTENT)
    error->all(FLERR, "Cannot use group node option with consistent mass style for fix temp/rescale command");

  tstr = NULL;
  if (strstr(arg[4], "v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    tstr = new char[n];
    strcpy(tstr, &arg[4][2]);
    tstyle = EQUAL;
  } else {
    t_start = universe->numeric(FLERR, arg[4]);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = universe->numeric(FLERR, arg[5]);
  t_window = universe->numeric(FLERR, arg[6]);
  fraction = universe->numeric(FLERR, arg[7]);

  // create a new compute temp
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp, id);
  strcat(id_temp, "_temp");

  char **newarg = new char*[5];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";
  newarg[3] = (char *) "group";
  if (group_flag == Group::NODE) 
    newarg[4] = (char *) "node";
  else if (group_flag == Group::ELEMENT) 
    newarg[4] = (char *) "element";
  else
    newarg[4] = (char *) "atom";
  modify->add_compute(5, newarg);
  delete [] newarg;
  tflag = 1;

  energy = 0.0;

}

/*  ----------------------------------------------------------------------  */

FixTempRescale::~FixTempRescale()
{
  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/*  ----------------------------------------------------------------------  */

int FixTempRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixTempRescale::init()
{
  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR, "Variable name for fix temp/rescale does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR, "Variable for fix temp/rescale is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR, "Temperature ID for fix temp/rescale does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/*  ----------------------------------------------------------------------  */

void FixTempRescale::end_of_step()
{
  double t_current = temperature->compute_scalar();

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  // protect against division by zero

  if (t_current == 0.0)
    error->all(FLERR, "Computed temperature for fix temp/rescale cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR, 
          "Fix temp/rescale variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v
  int i, j, k;
  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction * (t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    energy += (t_current-t_target) * efactor;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        v[i][0] *= factor;
        v[i][1] *= factor;
        v[i][2] *= factor;
      }

    double ****nodev = element->nodev;
    int ***nodemask = element->nodemask;
    int *emask = element->mask;
    int nelocal = element->nlocal;
    int *npe = element->npe;
    int *apc = element->apc;
    int *etype = element->etype;

    for (i = 0; i < nelocal; i++) {
      if (group_flag == Group::NODE) {
        for (j = 0; j < apc[etype[i]]; j++) 
          for (k = 0; k < npe[etype[i]]; k++) {
            if (nodemask[i][j][k] & groupbit) {
              nodev[i][j][k][0] *= factor;
              nodev[i][j][k][1] *= factor;
              nodev[i][j][k][2] *= factor;
            }  
          }
      } else if (group_flag == Group::ELEMENT) {
        if (emask[i] & groupbit) {
          for (j = 0; j < apc[etype[i]]; j++) 
            for (k = 0; k < npe[etype[i]]; k++) {
              nodev[i][j][k][0] *= factor;
              nodev[i][j][k][1] *= factor;
              nodev[i][j][k][2] *= factor;
            }  
        }
      }
    }
  }
}

/*  ----------------------------------------------------------------------  */

int FixTempRescale::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp, arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR, "Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR, 
          "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/*  ----------------------------------------------------------------------  */

void FixTempRescale::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/*  ----------------------------------------------------------------------  */

double FixTempRescale::compute_scalar()
{
  return energy;
}

/*  ----------------------------------------------------------------------
   extract thermostat properties
   -------------------------------------------------------------------------  */

void *FixTempRescale::extract(const char *str, int &dim)
{
  if (strcmp(str, "t_target") == 0) {
    dim = 0;
    return &t_target;
  }
  return NULL;
}
