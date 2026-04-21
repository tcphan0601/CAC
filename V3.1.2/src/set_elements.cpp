#include "set_elements.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <climits>
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "universe.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "modify.h"

using namespace CAC_NS;

enum{ELEMENT_SELECT, NODE_SELECT, ETYPE_SELECT, CTYPE_SELECT, GROUP_SELECT, GROUP_CTYPE_SELECT, REGION_SELECT};

//enum{TYPE, TYPE_FRACTION, TYPE_RATIO, TYPE_SUBSET, 
//     MOLECULE, X, Y, Z, CHARGE, MASS, SHAPE, LENGTH, TRI, 
//     DIPOLE, DIPOLE_RANDOM, SPIN, SPIN_RANDOM, QUAT, QUAT_RANDOM, 
//     THETA, THETA_RANDOM, ANGMOM, OMEGA, 
//     DIAMETER, DENSITY, VOLUME, IMAGE, BOND, ANGLE, DIHEDRAL, IMPROPER, 
//     MESO_E, MESO_CV, MESO_RHO, EDPD_TEMP, EDPD_CV, CC, SMD_MASS_DENSITY, 
//     SMD_CONTACT_RADIUS, DPDTHETA, INAME, DNAME, VX, VY, VZ};
enum{ETYPE, ETYPE_FRACTION, ETYPE_RATIO, ETYPE_SUBSET, 
     CTYPE, CTYPE_FRACTION, CTYPE_RATIO, CTYPE_SUBSET, 
     X, Y, Z, IMAGE};

#define BIG INT_MAX

/*  ----------------------------------------------------------------------  */

void SetElements::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Set_elements command before simulation box is defined");
  if (element->nelements == 0)
    error->all(FLERR, "Set_elements command with no elements existing");
  if (narg < 3) error->all(FLERR, "Illegal set_elements command");

  // style and ID info

  if (strcmp(arg[0], "element") == 0) style = ELEMENT_SELECT;
  //else if (strcmp(arg[0], "node") == 0) style = NODE_SELECT;
  else if (strcmp(arg[0], "etype") == 0) style = ETYPE_SELECT;
  else if (strcmp(arg[0], "ctype") == 0) style = CTYPE_SELECT;
  else if (strcmp(arg[0], "group") == 0) style = GROUP_SELECT;
  else if (strcmp(arg[0], "group/ctype") == 0) style = GROUP_CTYPE_SELECT;
  else if (strcmp(arg[0], "region") == 0) style = REGION_SELECT;
  else error->all(FLERR, "Illegal set_elements command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id, arg[1]);
  if (style == GROUP_CTYPE_SELECT) {
    if (narg < 4) error->all(FLERR, "Illegal set_elements command");
    n = strlen(arg[2]) + 1;
    id2 = new char[n];
    strcpy(id2, arg[2]);

  }

  select = nullptr;
  select2 = nullptr;
  selection(element->nlocal);

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen, "Setting element values ...\n");

  int allcount, origarg;

  int iarg = 2;
  if (style == GROUP_CTYPE_SELECT) iarg++;

  while (iarg < narg) {
    varflag = varflag1 = varflag2 = varflag3 = 0;
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg], "etype") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
      //        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      //        else 
      ivalue = universe->inumeric(FLERR, arg[iarg+1]);
      set(ETYPE);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ctype") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
      //        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      //        else 
      ivalue = universe->inumeric(FLERR, arg[iarg+1]);
      set(CTYPE);
      iarg += 2;

      /* 
         } else if (strcmp(arg[iarg], "x") == 0) {
         if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
      //        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      //        else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(X);
      iarg += 2;

      } else if (strcmp(arg[iarg], "y") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
      //        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      //        else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(Y);
      iarg += 2;

      } else if (strcmp(arg[iarg], "z") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
      //        if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
      //        else 
      dvalue = universe->numeric(FLERR, arg[iarg+1]);
      set(Z);
      iarg += 2;
      */
      /* 
         } else if (strcmp(arg[iarg], "charge") == 0) {
         if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
         if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
         else dvalue = universe->numeric(FLERR, arg[iarg+1]);
         if (!element->q_flag)
         error->all(FLERR, "Cannot set this attribute for this element style");
         set(CHARGE);
         iarg += 2;

         } else if (strcmp(arg[iarg], "mass") == 0) {
         if (iarg+2 > narg) error->all(FLERR, "Illegal set_elements command");
         if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
         else dvalue = universe->numeric(FLERR, arg[iarg+1]);
         if (!element->rmass_flag)
         error->all(FLERR, "Cannot set this attribute for this element style");
         set(MASS);
         iarg += 2;
         */
    } else if (strcmp(arg[iarg], "image") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal set command");
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1], "NULL") != 0) {
        ximageflag = 1;
        //          if (strstr(arg[iarg+1], "v_") == arg[iarg+1]) varparse(arg[iarg+1], 1);
        //          else 
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
            "Cannot set_elements non-zero image flag for non-periodic dimension");
      if (yimageflag && yimage && !domain->yperiodic)
        error->all(FLERR, 
            "Cannot set non-zero image flag for non-periodic dimension");
      if (zimageflag && zimage && !domain->zperiodic)
        error->all(FLERR, 
            "Cannot set non-zero image flag for non-periodic dimension");
      set(IMAGE);
      iarg += 4;
    } else error->all(FLERR, "Illegal set_elements command");

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

  element->count_ucells();
  element->count_nodes(1);

  // free local memory

  delete [] id;
  if (style == GROUP_CTYPE_SELECT) 
    delete [] id2;
  delete [] select;
  memory->destroy(select2);
}

/*  ----------------------------------------------------------------------
    select elements according to ELEMENT, TYPE, GROUP, REGION style
    n = nlocal or nlocal+nghost depending on keyword
    -------------------------------------------------------------------------  */

void SetElements::selection(int n)
{

  if (style == CTYPE_SELECT || style == GROUP_CTYPE_SELECT) {
    memory->destroy(select2);
    memory->create(select2, n, element->maxapc, "set_elements:select2");
  } else {
    delete [] select;
    select = new int[n];
  }
  int nlo, nhi;

  if (style == ELEMENT_SELECT) {
    if (element->tag_enable == 0)
      error->all(FLERR, "Cannot use set element with no element IDs defined");
    bigint nlobig, nhibig;
    universe->boundsbig(FLERR, id, MAXTAGINT, nlobig, nhibig);

    tagint *tag = element->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == ETYPE_SELECT) {
    universe->bounds(FLERR, id, element->netypes, nlo, nhi);

    int *etype = element->etype;
    for (int i = 0; i < n; i++)
      if (etype[i] >= nlo && etype[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == CTYPE_SELECT) {
    universe->bounds(FLERR, id, atom->ntypes, nlo, nhi);

    int **ctype = element->ctype;
    int *apc = element->apc;
    int *etype = element->etype;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < apc[etype[i]]; j++)
        if (ctype[i][j] >= nlo && ctype[i][j] <= nhi) select2[i][j] = 1;
        else select2[i][j] = 0;

  } else if (style == GROUP_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR, "Could not find set_elements group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = element->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_CTYPE_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR, "Could not find set_elements group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = element->mask;
    int **ctype = element->ctype;
    int *apc = element->apc;
    int *etype = element->etype;

    universe->bounds(FLERR, id2, atom->ntypes, nlo, nhi);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < apc[etype[i]]; j++) 
        if (mask[i] & groupbit && ctype[i][j] >= nlo && ctype[i][j] <= nhi) select2[i][j] = 1;
        else select2[i][j] = 0;

  } else if (style == REGION_SELECT) {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all(FLERR, "Set_elements region ID does not exist");
    domain->regions[iregion]->prematch();

    double **x = element->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0], x[i][1], x[i][2]))
        select[i] = 1;
      else select[i] = 0;
  }
}

/*  ----------------------------------------------------------------------
    set owned element properties directly
    either scalar or per-element values from element-style variable(s)
    -------------------------------------------------------------------------  */

void SetElements::set(int keyword)
{
  // evaluate element-style variable(s) if necessary

  vec1 = vec2 = vec3 = nullptr;

  //if (varflag) {
  //  int nlocal = element->nlocal;
  //  if (varflag1) {
  //    memory->create(vec1, nlocal, "set_elements:vec1");
  //    input->variable->compute_element(ivar1, 0, vec1, 1, 0);
  //  }
  //  if (varflag2) {
  //    memory->create(vec2, nlocal, "set_elements:vec2");
  //    input->variable->compute_element(ivar2, 0, vec2, 1, 0);
  //  }
  //  if (varflag3) {
  //    memory->create(vec3, nlocal, "set_elements:vec3");
  //    input->variable->compute_element(ivar3, 0, vec3, 1, 0);
  //  }
  //}

  // check if properties of elements in rigid bodies are updated
  // that are cached as per-body data.

  //  switch (keyword) {
  //    case X:
  //    case Y:
  //    case Z:
  //    case IMAGE:
  //      if (modify->check_rigid_list_overlap(select))
  //        error->warning(FLERR, "Changing a property of elements in rigid bodies "
  //            "that has no effect unless rigid bodies are rebuild");
  //      break;
  //    default: // assume no conflict for all other properties
  //      break;
  //  }

  // loop over selected elements

  int nlocal = element->nlocal;
  int *etype = element->etype;
  int *apc = element->apc;
  for (int i = 0; i < nlocal; i++) {
    if (style == CTYPE_SELECT || style == GROUP_CTYPE_SELECT) {
      if (keyword != CTYPE) error->one(FLERR, "Can only set ctype for style ctype or group/ctype");
      for (int j = 0; j < apc[etype[i]]; j++) {
        if (!select2[i][j]) continue;
        element->ctype[i][j] = ivalue;
        count++;
      }
    } else {

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

      // set values in per-element arrays
      // error check here in case element-style variables generated bogus value

      if (keyword == ETYPE) {
        if (ivalue <= 0 || ivalue > element->netypes)
          error->one(FLERR, "Invalid value in set_elements command");
        element->etype[i] = ivalue;
      } 
      else if (keyword == CTYPE) {
        for (int j = 0 ; j < apc[etype[i]]; j++) {
          if (ivalue <= 0 || ivalue > atom->ntypes)
            error->one(FLERR, "Invalid value in set_elements command");
          element->ctype[i][j] = ivalue;
        }
      } 
      //else if (keyword == X) element->x[i][0] = dvalue;
      //else if (keyword == Y) element->x[i][1] = dvalue;
      //else if (keyword == Z) element->x[i][2] = dvalue;

      // reset any or all of 3 image flags

      else if (keyword == IMAGE) {
        int xbox = (element->image[i] & IMGMASK) - IMGMAX;
        int ybox = (element->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        int zbox = (element->image[i] >> IMG2BITS) - IMGMAX;
        if (varflag1) ximage = static_cast<int>(xvalue);
        if (varflag2) yimage = static_cast<int>(yvalue);
        if (varflag3) zimage = static_cast<int>(zvalue);
        if (ximageflag) xbox = ximage;
        if (yimageflag) ybox = yimage;
        if (zimageflag) zbox = zimage;
        element->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
          (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
          (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
      }

      count++;
    }
  }

  // clear up per-element memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
}

/*  ----------------------------------------------------------------------
    set an owned element property randomly
    set seed based on element coordinates
    make element result independent of what proc owns it
    -------------------------------------------------------------------------  */
/* 
   void SetElements::setrandom(int keyword)
   {
   int i;

   double **x = element->x;
   int seed = ivalue;
   int nlocal = element->nlocal;
   int *type;

   if (keyword == ETYPE_FRACTION || 
   keyword == ETYPE_RATIO || 
   keyword == ETYPE_SUBSET) 
   type = element->etype;
   else 
   type = element->ctype;

   RanPark *ranpark = new RanPark(cac, 1);
   RanMars *ranmars = new RanMars(cac, seed + comm->me);

// set approx fraction of element types to newtype

if (keyword == ETYPE_FRACTION || keyword == CTYPE_FRACTION) {

for (i = 0; i < nlocal; i++)
if (select[i]) {
ranpark->reset(seed, x[i]);
if (ranpark->uniform() > fraction) continue;
type[i] = newtype;
count++;
}
} 

// set exact count of element types to newtype
// for TYPE_RATIO, exact = fraction out of total eligible
// for TYPE_SUBSET, exact = nsubset out of total eligible

else {

// count = number of eligible elements I own

count = 0;
for (i = 0; i < nlocal; i++)
if (select[i]) count++;

// convert specified fraction to nsubset

bigint bcount = count;
bigint allcount;
MPI_Allreduce(&bcount, &allcount, 1, MPI_CAC_BIGINT, MPI_SUM, world);

if (keyword == ETYPE_RATIO || keyword == CTYPE_RATIO) {
nsubset = static_cast<bigint> (fraction * allcount);
} else if (keyword == ETYPE_SUBSET || keyword == CTYPE_SUBSET) {
if (nsubset > allcount)
error->all(FLERR, "Set_elements type/subset value exceeds eligible elements");
}

// make selection

int *flag = memory->create(flag, count, "set_elements:flag");
int *work = memory->create(work, count, "set_elements:work");

ranmars->select_subset(nsubset, count, flag, work);

// change types of selected elements
// flag vector from select_subset() is only for eligible elements

count = 0;
int eligible = 0;
for (i = 0; i < nlocal; i++) {
if (!select[i]) continue;
if (flag[eligible]) {
  type[i] = newtype;
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
*/

/*  ----------------------------------------------------------------------  */

void SetElements::varparse(char *name, int m)
{
  error->all(FLERR, ""); 
  varflag = 1;

  name = &name[2];
  int n = strlen(name) + 1;
  char *str = new char[n];
  strcpy(str, name);

  int ivar = input->variable->find(str);
  delete [] str;

  if (ivar < 0)
    error->all(FLERR, "Variable name for set_elements command does not exist");
  //if (!input->variable->elementstyle(ivar))
  //  error->all(FLERR, "Variable for set_elements command is invalid style");

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  }
}
