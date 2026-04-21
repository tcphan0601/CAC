#include "set_atoms.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <climits>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "universe.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "random_mars.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "modify.h"

using namespace CAC_NS;
using namespace MathConst;

enum{ATOM_SELECT, TYPE_SELECT, GROUP_SELECT, REGION_SELECT};

//enum{TYPE, TYPE_FRACTION, TYPE_RATIO, TYPE_SUBSET, 
//     MOLECULE, X, Y, Z, CHARGE, MASS, SHAPE, LENGTH, TRI, 
//     DIPOLE, DIPOLE_RANDOM, SPIN, SPIN_RANDOM, QUAT, QUAT_RANDOM, 
//     THETA, THETA_RANDOM, ANGMOM, OMEGA, 
//     DIAMETER, DENSITY, VOLUME, IMAGE, BOND, ANGLE, DIHEDRAL, IMPROPER, 
//     MESO_E, MESO_CV, MESO_RHO, EDPD_TEMP, EDPD_CV, CC, SMD_MASS_DENSITY, 
//     SMD_CONTACT_RADIUS, DPDTHETA, INAME, DNAME, VX, VY, VZ};
enum{TYPE, TYPE_FRACTION, TYPE_RATIO, TYPE_SUBSET, GRAIN_ID, 
     X, Y, Z, IMAGE, VX, VY, VZ};

#define BIG INT_MAX

/*  ----------------------------------------------------------------------  */

void SetAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Set_atoms command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR, "Set_atoms command with no atoms existing");
  if (narg < 3) error->all(FLERR, "Illegal set_atoms command");

  // style and ID info

  if (strcmp(arg[0], "atom") == 0) style = ATOM_SELECT;
  else if (strcmp(arg[0], "type") == 0) style = TYPE_SELECT;
  else if (strcmp(arg[0], "group") == 0) style = GROUP_SELECT;
  else if (strcmp(arg[0], "region") == 0) style = REGION_SELECT;
  else error->all(FLERR, "Illegal set_atoms command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id, arg[1]);
  select = nullptr;
  selection(atom->nlocal);

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen, "Setting atom values ...\n");

  int allcount, origarg;

  int iarg = 2;
  while (iarg < narg) {
    varflag = varflag1 = varflag2 = varflag3 = 0;
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg], "type") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      ivalue = universe->inumeric(FLERR, arg[iarg+1]);
      set(TYPE);
      iarg += 2;

    } else if (strcmp(arg[iarg], "type/fraction") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal set_atoms command");
      newtype = universe->inumeric(FLERR, arg[iarg+1]);
      fraction = universe->numeric(FLERR, arg[iarg+2]);
      ivalue = universe->inumeric(FLERR, arg[iarg+3]);
      if (newtype <= 0 || newtype > atom->ntypes)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (ivalue <= 0)
        error->all(FLERR, "Invalid random number seed in set_atoms command");
      setrandom(TYPE_FRACTION);
      iarg += 4;

    } else if (strcmp(arg[iarg], "type/ratio") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal set_atoms command");
      newtype = universe->inumeric(FLERR, arg[iarg+1]);
      fraction = universe->numeric(FLERR, arg[iarg+2]);
      ivalue = universe->inumeric(FLERR, arg[iarg+3]);
      if (newtype <= 0 || newtype > atom->ntypes)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (ivalue <= 0)
        error->all(FLERR, "Invalid random number seed in set_atoms command");
      setrandom(TYPE_RATIO);

      iarg += 4;

    } else if (strcmp(arg[iarg], "type/subset") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal set_atoms command");
      newtype = universe->inumeric(FLERR, arg[iarg+1]);
      nsubset = universe->bnumeric(FLERR, arg[iarg+2]);
      ivalue = universe->inumeric(FLERR, arg[iarg+3]);
      if (newtype <= 0 || newtype > atom->ntypes)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (nsubset < 0)
        error->all(FLERR, "Invalid value in set_atoms command");
      if (ivalue <= 0)
        error->all(FLERR, "Invalid random number seed in set_atoms command");
      setrandom(TYPE_SUBSET);
      iarg += 4;

    } else if (strcmp(arg[iarg], "grain/id") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      ivalue = universe->inumeric(FLERR, arg[iarg+1]);
      set(GRAIN_ID);
      iarg += 2;

    } else if (strcmp(arg[iarg], "x") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(X);
      iarg += 2;

    } else if (strcmp(arg[iarg], "y") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(Y);
      iarg += 2;

    } else if (strcmp(arg[iarg], "z") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(Z);
      iarg += 2;

    } else if (strcmp(arg[iarg], "vx") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(VX);
      iarg += 2;

    } else if (strcmp(arg[iarg], "vy") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(VY);
      iarg += 2;

    } else if (strcmp(arg[iarg], "vz") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
//      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//      else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(VZ);
      iarg += 2;
/* 
    } else if (strcmp(arg[iarg], "charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      else dvalue = universe->numeric(FLERR, arg[iarg+1]);
      if (!atom->q_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
      set(CHARGE);
      iarg += 2;

    } else if (strcmp(arg[iarg], "mass") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_atoms command");
      if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      else dvalue = universe->numeric(FLERR, arg[iarg+1]);
      if (!atom->rmass_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
      set(MASS);
      iarg += 2;
 */
    } else if (strcmp(arg[iarg], "image") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal set command");
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1], "NULL") != 0) {
        ximageflag = 1;
//        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
//        else 
        ximage = universe->inumeric(FLERR, arg[iarg+1]);
      }
      if (strcmp(arg[iarg+2], "NULL") != 0) {
        yimageflag = 1;
//        if (strstr(arg[iarg+2], "v_") == arg[iarg+2]) varparse(arg[iarg+2], 2);
//        else 
        yimage = universe->inumeric(FLERR, arg[iarg+2]);
      }
      if (strcmp(arg[iarg+3], "NULL") != 0) {
        zimageflag = 1;
//        if (strstr(arg[iarg+3], "v_") == arg[iarg+3]) varparse(arg[iarg+3], 3);
//        else 
        zimage = universe->inumeric(FLERR, arg[iarg+3]);
      }
      if (ximageflag && ximage && !domain->xperiodic)
        error->all(FLERR, 
                   "Cannot set_atoms non-zero image flag for non-periodic dimension");
      if (yimageflag && yimage && !domain->yperiodic)
        error->all(FLERR, 
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (zimageflag && zimage && !domain->zperiodic)
        error->all(FLERR, 
                   "Cannot set non-zero image flag for non-periodic dimension");
      set(IMAGE);
      iarg += 4;
    } else error->all(FLERR, "Illegal set_atoms command");

    // statistics
    // for CC option, include species index

    MPI_Allreduce(&count, &allcount, 1, MPI_INT, MPI_SUM, world);

    if (comm->me == 0) {

      if (screen) {
        fprintf(screen, "  %d settings made for %s\n", 
            allcount, arg[origarg]);
      }
      if (logfile) {
        fprintf(logfile, "  %d settings made for %s\n", 
            allcount, arg[origarg]);
      }
    }
  }

  // free local memory

  delete [] id;
  delete [] select;
}

/*  ----------------------------------------------------------------------
   select atoms according to ATOM, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
   -------------------------------------------------------------------------  */

void SetAtoms::selection(int n)
{
  delete [] select;
  select = new int[n];
  int nlo, nhi;

  if (style == ATOM_SELECT) {
    if (atom->tag_enable == 0)
      error->all(FLERR, "Cannot use set atom with no atom IDs defined");
    bigint nlobig, nhibig;
    universe->boundsbig(FLERR, id, MAXTAGINT, nlobig, nhibig);

    tagint *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    universe->bounds(FLERR, id, atom->ntypes, nlo, nhi);

    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR, "Could not find set_atoms group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all(FLERR, "Set_atoms region ID does not exist");
    domain->regions[iregion]->prematch();

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2]))
        select[i] = 1;
      else select[i] = 0;
  }
}

/*  ----------------------------------------------------------------------
   set owned atom properties directly
   either scalar or per-atom values from atom-style variable(s)
   -------------------------------------------------------------------------  */

void SetAtoms::set(int keyword)
{
  // evaluate atom-style variable(s) if necessary

  vec1 = vec2 = vec3 = nullptr;

  if (varflag) {
    int nlocal = atom->nlocal;
    if (varflag1) {
      memory->create(vec1, nlocal, "set_atoms:vec1");
      input->variable->compute_atom(ivar1, 0, vec1, 1, 0);
    }
    if (varflag2) {
      memory->create(vec2, nlocal, "set_atoms:vec2");
      input->variable->compute_atom(ivar2, 0, vec2, 1, 0);
    }
    if (varflag3) {
      memory->create(vec3, nlocal, "set_atoms:vec3");
      input->variable->compute_atom(ivar3, 0, vec3, 1, 0);
    }
  }

  // check if properties of atoms in rigid bodies are updated
  // that are cached as per-body data.

//  switch (keyword) {
//    case X:
//    case Y:
//    case Z:
//    case IMAGE:
//      if (modify->check_rigid_list_overlap(select))
//        error->warning(FLERR, "Changing a property of atoms in rigid bodies "
//            "that has no effect unless rigid bodies are rebuild");
//      break;
//    default: // assume no conflict for all other properties
//      break;
//  }

  // loop over selected atoms

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    // overwrite dvalue, ivalue, xyzw value if variables defined
    // else the input script scalar value remains in place

    if (varflag) {
      if (varflag1) {
        dvalue = xvalue = vec1[i];
        ivalue = static_cast<int> (dvalue);
      }
      if (varflag2) yvalue = vec2[i];
      if (varflag3) zvalue = vec3[i];
    }

    // set values in per-atom arrays
    // error check here in case atom-style variables generated bogus value

    if (keyword == TYPE) {
      if (ivalue <= 0 || ivalue > atom->ntypes)
        error->one(FLERR, "Invalid value in set_atoms command");
      atom->type[i] = ivalue;
    } else if (keyword == GRAIN_ID) {
      if (ivalue < 0 || ivalue > atom->ngrains)
        error->one(FLERR, "Invalid value in set_atoms command");
      atom->grain_tag[i] = ivalue;
    }
    else if (keyword == X) atom->x[i][0] = dvalue;
    else if (keyword == Y) atom->x[i][1] = dvalue;
    else if (keyword == Z) atom->x[i][2] = dvalue;
    else if (keyword == VX) atom->v[i][0] = dvalue;
    else if (keyword == VY) atom->v[i][1] = dvalue;
    else if (keyword == VZ) atom->v[i][2] = dvalue;


    // reset any or all of 3 image flags

    else if (keyword == IMAGE) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if (varflag1) ximage = static_cast<int>(xvalue);
      if (varflag2) yimage = static_cast<int>(yvalue);
      if (varflag3) zimage = static_cast<int>(zvalue);
      if (ximageflag) xbox = ximage;
      if (yimageflag) ybox = yimage;
      if (zimageflag) zbox = zimage;
      atom->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }

    count++;
  }

  // clear up per-atom memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
}

/*  ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom coordinates
   make atom result independent of what proc owns it
   -------------------------------------------------------------------------  */

void SetAtoms::setrandom(int keyword)
{
  int i;

  double **x = atom->x;
  int seed = ivalue;

  RanPark *ranpark = new RanPark(cac, 1);
  RanMars *ranmars = new RanMars(cac, seed + comm->me);

  // set approx fraction of atom types to newtype

  if (keyword == TYPE_FRACTION) {
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (select[i]) {
        ranpark->reset(seed, x[i]);
        if (ranpark->uniform() > fraction) continue;
        atom->type[i] = newtype;
        count++;
      }

    // set exact count of atom types to newtype
    // for TYPE_RATIO, exact = fraction out of total eligible
    // for TYPE_SUBSET, exact = nsubset out of total eligible

  } else if (keyword == TYPE_RATIO || keyword == TYPE_SUBSET) {
    int nlocal = atom->nlocal;

    // count = number of eligible atoms I own

    count = 0;
    for (i = 0; i < nlocal; i++)
      if (select[i]) count++;

    // convert specified fraction to nsubset

    bigint bcount = count;
    bigint allcount;
    MPI_Allreduce(&bcount, &allcount, 1, MPI_CAC_BIGINT, MPI_SUM, world);

    if (keyword == TYPE_RATIO) {
      nsubset = static_cast<bigint> (fraction * allcount);
    } else if (keyword == TYPE_SUBSET) {
      if (nsubset > allcount)
        error->all(FLERR, "Set_atoms type/subset value exceeds eligible atoms");
    }

    // make selection

    int *flag = memory->create(flag, count, "set_atoms:flag");
    int *work = memory->create(work, count, "set_atoms:work");

    ranmars->select_subset(nsubset, count, flag, work);

    // change types of selected atoms
    // flag vector from select_subset() is only for eligible atoms
    count = 0;
    int eligible = 0;
    for (i = 0; i < nlocal; i++) {
      if (!select[i]) continue;
      if (flag[eligible]) {
        atom->type[i] = newtype;
        count++;
      }
      eligible++;
    }

    // clean up

    memory->destroy(flag);
    memory->destroy(work);

  }

  delete ranpark;
  delete ranmars;
}


/*  ----------------------------------------------------------------------  */

void SetAtoms::varparse(char *name, int m)
{
  varflag = 1;

  name = &name[2];
  int n = strlen(name) + 1;
  char *str = new char[n];
  strcpy(str, name);

  int ivar = input->variable->find(str);
  delete [] str;

  if (ivar < 0)
    error->all(FLERR, "Variable name for set_atoms command does not exist");
  //if (!input->variable->atomstyle(ivar))
  //  error->all(FLERR, "Variable for set_atoms command is invalid style");

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  }
}
