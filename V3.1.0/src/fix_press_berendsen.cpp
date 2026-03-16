#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fix_press_berendsen.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "universe.h"
#include "comm.h"
#include "modify.h"
#include "group.h"
#include "compute.h"
#include "update.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NOBIAS, BIAS};
enum{NONE, XYZ, XY, YZ, XZ};
enum{ISO, ANISO};

/*  ----------------------------------------------------------------------  */

FixPressBerendsen::FixPressBerendsen(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  //id_temp(nullptr), 
  id_press(nullptr), 
  //tflag(0), 
  pflag(0)
{
  if (narg < 5) error->all(FLERR, "Illegal fix press/berendsen command");

  box_change_size = 1;

  // Berendsen barostat applied every step

  nevery = 1;

  // default values

  pcouple = NONE;
  bulkmodulus = 10.0;
  allremap = 1;

  for (int i = 0; i < 3; i++) {
    p_start[i] = p_stop[i] = p_period[i] = 0.0;
    p_flag[i] = 0;
    p_period[i] = 0.0;
  }

  // process keywords

  dimension = domain->dimension;

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "iso") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      pcouple = XYZ;
      p_start[0] = p_start[1] = p_start[2] = universe->numeric(FLERR, arg[iarg+1]);
      p_stop[0] = p_stop[1] = p_stop[2] = universe->numeric(FLERR, arg[iarg+2]);
      p_period[0] = p_period[1] = p_period[2] = universe->numeric(FLERR, arg[iarg+3]);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg], "aniso") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      pcouple = NONE;
      p_start[0] = p_start[1] = p_start[2] = universe->numeric(FLERR, arg[iarg+1]);
      p_stop[0] = p_stop[1] = p_stop[2] = universe->numeric(FLERR, arg[iarg+2]);
      p_period[0] = p_period[1] = p_period[2] = universe->numeric(FLERR, arg[iarg+3]);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg], "x") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[0] = universe->numeric(FLERR, arg[iarg+1]);
      p_stop[0] = universe->numeric(FLERR, arg[iarg+2]);
      p_period[0] = universe->numeric(FLERR, arg[iarg+3]);
      p_flag[0] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "y") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[1] = universe->numeric(FLERR, arg[iarg+1]);
      p_stop[1] = universe->numeric(FLERR, arg[iarg+2]);
      p_period[1] = universe->numeric(FLERR, arg[iarg+3]);
      p_flag[1] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "z") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[2] = universe->numeric(FLERR, arg[iarg+1]);
      p_stop[2] = universe->numeric(FLERR, arg[iarg+2]);
      p_period[2] = universe->numeric(FLERR, arg[iarg+3]);
      p_flag[2] = 1;
      iarg += 4;
      if (dimension == 2)
        error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");

    } else if (strcmp(arg[iarg], "couple") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      if (strcmp(arg[iarg+1], "xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1], "xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1], "yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1], "xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1], "none") == 0) pcouple = NONE;
      else error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;

    } else if (strcmp(arg[iarg], "modulus") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      bulkmodulus = universe->numeric(FLERR, arg[iarg+1]);
      if (bulkmodulus <= 0.0)
        error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dilate") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix press/berendsen command");
      if (strcmp(arg[iarg+1], "all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1], "partial") == 0) allremap = 0;
      else error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;
    //} else if (strcmp(arg[iarg], "group") == 0) {
    //  if (iarg+2 > narg)
    //    error->all(FLERR, "Illegal fix press/berendsen command");
    //  if (strcmp(arg[iarg+1], "atom") == 0) group_flag = Group::ATOM;
    //  else if (strcmp(arg[iarg+1], "node") == 0) group_flag = Group::NODE;
    //  else if (strcmp(arg[iarg+1], "element") == 0) group_flag = Group::ELEMENT;
    //  else error->all(FLERR, "Illegal fix press/berendsen command");
    //  iarg += 2;
    } else error->all(FLERR, "Illegal fix press/berendsen command");
  }

  if (allremap == 0) restart_pbc = 1;

  // error checks

  if (dimension == 2 && p_flag[2])
    error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR, 
               "Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR, 
               "Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR, 
               "Cannot use fix press/berendsen on a non-periodic dimension");

  if (pcouple == XYZ && dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] ||
       p_stop[0] != p_stop[1] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XYZ && dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] ||
       p_period[1] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");

  if ((p_flag[0] && p_period[0] <= 0.0) ||
      (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0))
    error->all(FLERR, "Fix press/berendsen damping parameters must be > 0.0");

  // pstyle = ISO if XYZ coupling or XY coupling in 2d -> 1 dof
  // else pstyle = ANISO -> 3 dof

  if (pcouple == XYZ || (dimension == 2 && pcouple == XY)) pstyle = ISO;
  else pstyle = ANISO;

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp, id);
  strcat(id_temp, "_temp");

  char **newarg = new char*[4];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3, newarg);
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press, id);
  strcat(id_press, "_press");

  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  modify->add_compute(4, newarg);
  delete [] newarg;
  pflag = 1;

  //nrigid = 0;
  //rfix = nullptr;
}

/*  ----------------------------------------------------------------------  */

FixPressBerendsen::~FixPressBerendsen()
{
  //delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/*  ----------------------------------------------------------------------  */

int FixPressBerendsen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixPressBerendsen::init()
{
  if (domain->triclinic)
    error->all(FLERR, "Cannot use fix press/berendsen with triclinic box");

  // insure no conflict with fix deform

  //for (int i = 0; i < modify->nfix; i++)
  //  if (strcmp(modify->fix[i]->style, "deform") == 0) {
  //    int *dimflag = ((FixDeform *) modify->fix[i])->dimflag;
  //    if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) ||
  //        (p_flag[2] && dimflag[2]))
  //      error->all(FLERR, "Cannot use fix press/berendsen and "
  //                 "fix deform on same component of stress tensor");
  //  }

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR, "Temperature ID for fix press/berendsen does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  icompute = modify->find_compute(id_press);
  if (icompute < 0)
    error->all(FLERR, "Pressure ID for fix press/berendsen does not exist");
  pressure = modify->compute[icompute];

  // detect if any rigid fixes exist so rigid bodies move when box is remapped
  // rfix[] = indices to each fix rigid

  //delete [] rfix;
  //nrigid = 0;
  //rfix = nullptr;

  //for (int i = 0; i < modify->nfix; i++)
  //  if (modify->fix[i]->rigid_flag) nrigid++;
  //if (nrigid) {
  //  rfix = new int[nrigid];
  //  nrigid = 0;
  //  for (int i = 0; i < modify->nfix; i++)
  //    if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  //}
}

/*  ----------------------------------------------------------------------
   compute T, P before integrator starts
-------------------------------------------------------------------------  */

void FixPressBerendsen::setup(int vflag)
{
  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/*  ----------------------------------------------------------------------  */

void FixPressBerendsen::end_of_step()
{
  // compute new T, P

  if (pstyle == ISO) {
    temperature->compute_scalar();
    pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      dilation[i] =
        pow(1.0 - update->dt/p_period[i] *
            (p_target[i]-p_current[i])/bulkmodulus, 1.0/3.0);
    }
  }

  // remap simulation box and atoms/elements
  // redo KSpace coeffs since volume has changed

  remap();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/*  ----------------------------------------------------------------------  */

void FixPressBerendsen::couple()
{
  double *tensor = pressure->vector;

  if (pstyle == ISO)
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  else if (pcouple == XYZ) {
    double ave = 1.0/3.0 * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/*  ----------------------------------------------------------------------
   change box size
   remap all atoms/elements or fix group atoms/elements depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
-------------------------------------------------------------------------  */

void FixPressBerendsen::remap()
{
  int i, j, k;
  double oldlo, oldhi, ctr;

  double **ax = atom->x;
  double ****nodex = element->nodex;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;

  // convert pertinent atoms/elements and rigid bodies to lamda coords

  if (allremap) {
    domain->x2lamda(nalocal, ax);
    domain->nodex2lamda(nelocal, nodex);
  } else {
    for (i = 0; i < nalocal; i++)
      if (amask[i] & groupbit)
        domain->x2lamda(ax[i], ax[i]);
    for (i = 0; i < nelocal; i++)
      if (emask[i] & groupbit)
        for (j = 0; j < apc[etype[i]]; j++)
          for (k = 0; k < npe[etype[i]]; k++)
            domain->x2lamda(nodex[i][j][k], nodex[i][j][k]);
  }

  //if (nrigid)
  //  for (i = 0; i < nrigid; i++)
  //    modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo-ctr) * dilation[i] + ctr;
      domain->boxhi[i] = (oldhi-ctr) * dilation[i] + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms/elements and rigid bodies back to box coords

  if (allremap) {
    domain->lamda2x(nalocal, ax);
    domain->lamda2nodex(nelocal, nodex);
  } else {
    for (i = 0; i < nalocal; i++)
      if (amask[i] & groupbit)
        domain->lamda2x(ax[i], ax[i]);
    for (i = 0; i < nelocal; i++)
      if (emask[i] & groupbit)
        for (j = 0; j < apc[etype[i]]; j++)
          for (k = 0; k < npe[etype[i]]; k++)
            domain->lamda2x(nodex[i][j][k], nodex[i][j][k]);
  }

  element->evec->update_center_coord();
  //if (nrigid)
  //  for (i = 0; i < nrigid; i++)
  //    modify->fix[rfix[i]]->deform(1);
}

/*  ----------------------------------------------------------------------  */

int FixPressBerendsen::modify_param(int narg, char **arg)
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

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR, "Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR, "Temperature for NPT is not for group all");

    // reset id_temp of pressure to new temperature ID

    icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR, "Pressure ID for fix press/berendsen does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0], "press") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press, arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR, "Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR, "Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}
