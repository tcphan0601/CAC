#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "region_sphere.h"
#include "update.h"
#include "input.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;

enum{CONSTANT,VARIABLE};

/*-------------------------------------------------*/

RegSphere::RegSphere(CAC *cac, int narg, char **arg) : Region(cac, narg, arg)
{
  options(narg-6,&arg[6]);

  xc = xscale*universe->numeric(FLERR,arg[2]);
  yc = yscale*universe->numeric(FLERR,arg[3]);
  zc = zscale*universe->numeric(FLERR,arg[4]);

  rstr = NULL;

  radius = xscale*universe->numeric(FLERR,arg[5]);
  rstyle = CONSTANT;

  // error check

  if (radius < 0.0) error->all(FLERR, "Illegal region sphere command");

  // extent of sphere 
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    extent_xlo = xc - radius;
    extent_xhi = xc + radius;
    extent_ylo = yc - radius;
    extent_yhi = yc + radius;
    extent_zlo = zc - radius;
    extent_zhi = zc + radius;
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax]; 
}

/*-----------------------------------------------------------*/

RegSphere::~RegSphere()
{
  delete [] rstr;
  delete [] contact;
}

/*------------------------------------------------------------*/

void RegSphere::init()
{
  Region::init();
}

/*---------------------------------------------------------------
    inside = 1 if x,y,z is inside or on surface
    inside = 0 if x,y,z is outside and not on surface
 -----------------------------------------------------------------*/

int RegSphere::inside(double x, double y, double z)
{
  double delx = x - xc;
  double dely = y - yc;
  double delz = z - zc;
  double r = sqrt(delx*delx + dely*dely + delz*delz);

  if (r <= radius) return 1;
  return 0;
}


