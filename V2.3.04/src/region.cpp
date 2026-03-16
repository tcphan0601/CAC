#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "region.h"
#include "lattice.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "error.h"
#include "universe.h"
#include "force.h"

using namespace CAC_NS;

/*---------------------------------------------------------------------*/

Region::Region(CAC *cac, int narg, char **arg) : Pointers(cac)
{
   int n = strlen(arg[0])+1;
   id = new char[n];
   strcpy(id,arg[0]);

   n = strlen(arg[1]) + 1;
   style = new char[n];
   strcpy(style,arg[1]);

   varshape = 0;
   xstr = ystr = zstr = tstr = NULL;
   dx = dy = dz = 0.0;
}

/*----------------------------------------------------------------------*/

Region::~Region()
{
  delete [] id;
  delete [] style;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] tstr;
}

/*---------------------------------------------------------------------------*/

void Region::init()
{
}

/*----------------------------------------------------------------------------
  parse optional parameters at end of region input line
------------------------------------------------------------------------------*/

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal region command");

  // option defaults

  interior = 1;
  scaleflag = 1;
  moveflag = rotateflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal region command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal region command");
      iarg +=2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;
    } else error->all(FLERR, "Illegal region command");
  }

  // setup scaling

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else xscale = yscale = zscale = 1.0;

  if (moveflag || rotateflag) dynamic = 1;
  else dynamic = 0;
}

/*------------------------------------------------------------------------------*/

void Region::prematch()
{
}

/*-------------------------------------------------------------------------------*/

int Region::match(double x, double y, double z)
{
  return !(inside(x,y,z) ^ interior);
}

/*-------------------------------------------------------------------------------*/

int Region::match(double *x)
{
  return !(inside(x[0],x[1],x[2]) ^ interior);
}
