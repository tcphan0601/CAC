#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "region_cylinder.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "error.h"
#include "universe.h"
#include "comm.h"

using namespace CAC_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegCylinder::RegCylinder(CAC *cac, int narg, char **arg) :
  Region(cac, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"x") && strcmp(arg[2],"y") && strcmp(arg[2],"z"))
    error->all(FLERR,"Illegal region cylinder command");
  axis = arg[2][0];

  if (axis == 'x') {
    c1 = yscale*universe->numeric(FLERR,arg[3]);
    c2 = zscale*universe->numeric(FLERR,arg[4]);
  } else if (axis == 'y') {
    c1 = xscale*universe->numeric(FLERR,arg[3]);
    c2 = zscale*universe->numeric(FLERR,arg[4]);
  } else if (axis == 'z') {
    c1 = xscale*universe->numeric(FLERR,arg[3]);
    c2 = yscale*universe->numeric(FLERR,arg[4]);
  }

  radius = universe->numeric(FLERR,arg[5]);
  if (axis == 'x') radius *= yscale;
  else radius *= xscale;

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[0];
      else lo = domain->boxlo_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[1];
      else lo = domain->boxlo_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[2];
      else lo = domain->boxlo_bound[2];
    }
  } else {
    if (axis == 'x') lo = xscale*universe->numeric(FLERR,arg[6]);
    if (axis == 'y') lo = yscale*universe->numeric(FLERR,arg[6]);
    if (axis == 'z') lo = zscale*universe->numeric(FLERR,arg[6]);
  }

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[0];
      else hi = domain->boxhi_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[1];
      else hi = domain->boxhi_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[2];
      else hi = domain->boxhi_bound[2];
    }
  } else {
    if (axis == 'x') hi = xscale*universe->numeric(FLERR,arg[7]);
    if (axis == 'y') hi = yscale*universe->numeric(FLERR,arg[7]);
    if (axis == 'z') hi = zscale*universe->numeric(FLERR,arg[7]);
  }

  // error check

  if (radius <= 0.0) error->all(FLERR,
      "Illegal radius in region cylinder command");
  if (hi <= lo) error->all(FLERR,
      "Illegal lo and hi in region cylinder command (lo shoud be < hi)");


  rsq = radius * radius;

  // extent of cylinder
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 - radius;
      extent_yhi = c1 + radius;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'y') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'z') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = c2 - radius;
      extent_yhi = c2 + radius;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else bboxflag = 0;

}

/* ---------------------------------------------------------------------- */

RegCylinder::~RegCylinder()
{
}

/* ---------------------------------------------------------------------- */

void RegCylinder::init()
{
  Region::init();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
   ------------------------------------------------------------------------- */

int RegCylinder::inside(double x, double y, double z)
{
  double del1,del2,distsq;
  int inside = 0;

  if (axis == 'x') {
    del1 = y - c1;
    del2 = z - c2;
    distsq = del1*del1 + del2*del2;
    if (distsq <= rsq && x >= lo && x <= hi) inside = 1;
  } else if (axis == 'y') {
    del1 = x - c1;
    del2 = z - c2;
    distsq = del1*del1 + del2*del2;
    if (distsq <= rsq && y >= lo && x <= hi) inside = 1;
  } else {
    del1 = x - c1;
    del2 = y - c2;
    distsq = del1*del1 + del2*del2;
    if (distsq <= rsq && z >= lo && z <= hi) inside = 1;
  }

  return inside;
}


