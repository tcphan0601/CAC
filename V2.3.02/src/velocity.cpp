#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "velocity.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "universe.h"
#include "variable.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

enum{CREATE,SET,SCALE,RAMP,ZERO};
enum{ALL,LOCAL,GEOM};
enum{NONE,CONSTANT,EQUAL,ATOM};
enum{NODE,ELEMENT};

#define WARMUP 100
#define SMALL 0.001

/*----------------------------------------------------------------------------*/

Velocity::Velocity(CAC *cac) : Pointers(cac) {}

/*----------------------------------------------------------------------------*/

void Velocity::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal velocity command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Velocity command before simulation box is defined");
  if ((atom->natoms + element->nelements) == 0)
    error->all(FLERR,"Velocity command with no atoms or elements existing");

  // atoms masses must all be set

  atom->check_mass();

  //identify group

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Cound not find velocity group ID");
  groupbit = group->bitmask[igroup];

  // identify style

  if (strcmp(arg[1],"set") == 0) style = SET;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else error->all(FLERR, "Illegal velocity command");

  // set defaults

  temperature = NULL;
  dist_flag = 0;
  sum_flag = 0;
  momentum_flag = 0;
  rotation_flag = 0;
  bias_flag = 0;
  loop_flag = ALL;
  scale_flag = 1;
  rfix = -1;
  group_flag = NODE;

  // read options from end of input line
  // change defaults as options specify

  if (style == SET) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);

  // initialize velocities based on style
  // create() invoked differently, so can be called externally

  if (style == SET) set(narg-2,&arg[2]);
  else if (style == RAMP) ramp(narg-2,&arg[2]);

}

/*------------------------------------------------------------------
  parse optional parameters at end of velocity input line
  --------------------------------------------------------------------*/

void Velocity::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal velocity command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) sum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) sum_flag = 1;
      else error->all(FLERR, "Illegal velocity command");
      iarg+=2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal velocity command");
      if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR, "Illegal velocity command");
      iarg += 2;
    //} else if (strcmp(arg[iarg],"units") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR, "Illegal velocity command");
    //  if (strcmp(arg[iarg+1],"box") == 0) scale_flag = 0;
    //  else if (strcmp(arg[iarg+1],"lattice") == 0) scale_flag = 1;
    //  else error->all(FLERR, "Illegal velocity command");
    //  iarg += 2;
    } else error->all(FLERR, "Illegal velocity command");
  }


}

/*-------------------------------------------------------------------*/

void Velocity::set(int narg, char **arg)
{
  int xstyle, ystyle, zstyle, varflag;
  double vx, vy, vz;
  char *xstr, *ystr, *zstr;
  int xvar,yvar,zvar;

  xstyle = ystyle = zstyle = CONSTANT;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[0],"v_") == arg[0]) {
    int n = strlen(&arg[0][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[0][2]);
  } else if (strcmp(arg[0],"NULL") == 0) xstyle = NONE;
  else vx = universe->numeric(FLERR,arg[0]);

  if (strstr(arg[1],"v_") == arg[1]) {
    int n = strlen(&arg[1][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[1][2]);
  } else if (strcmp(arg[1],"NULL") == 0) ystyle = NONE;
  else vy = universe->numeric(FLERR, arg[1]);

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[2][2]);
  } else if (strcmp(arg[2],"NULL") == 0) zstyle = NONE;
  else vz = universe->numeric(FLERR,arg[2]);

  // set and apply scale factors

  xscale = yscale = zscale = 1.0;

  if (xstyle && !xstr) vx *= xscale;
  if (ystyle && !ystr) vy *= yscale;
  if (zstyle && !zstr) vz *= zscale;

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for velocity set does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for velocity set is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    error->all(FLERR,"ATOM style variable not ready yet");
    //varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (domain->dimension == 2) {
    if (zstyle == CONSTANT && vz != 0.0)
      error->all(FLERR,"Cannot set non-zero z velocity for 2d simulation");
    if (zstyle == EQUAL || zstyle == ATOM)
      error->all(FLERR,"Cannot set variable z velocity for 2d simulation");
  }

  // allocate vfield array if necessary

  double **vfield = NULL;
  if (varflag == ATOM) memory->create(vfield,atom->nlocal,3,"velocity:vfield");

  // add velocity to atoms/element nodes
 
  double **v = atom->v;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodev = element->nodev;
  int **nodemask = element->nodemask;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int npe = element->npe;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nalocal; i++) 
      if (amask[i] & groupbit) 
        if (sum_flag == 0) {
          if (xstyle) v[i][0] = vx;
          if (ystyle) v[i][1] = vy;
          if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle) v[i][0] += vx;
          if (ystyle) v[i][1] += vy;
          if (zstyle) v[i][2] += vz;
        }

    for (int i = 0; i < nelocal; i++) 
      if (group_flag == NODE) { 
        for (int j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit) 
            if (sum_flag == 0) {
              if (xstyle) nodev[i][j][0] = vx;
              if (ystyle) nodev[i][j][1] = vy;
              if (zstyle) nodev[i][j][2] = vz;
            } else {
              if (xstyle) nodev[i][j][0] += vx;
              if (ystyle) nodev[i][j][1] += vy;
              if (zstyle) nodev[i][j][2] += vz;
            }
      } else if (emask[i] & groupbit) {
        for (int j = 0; j < npe; j++)
          if (sum_flag == 0) {
            if (xstyle) nodev[i][j][0] = vx;
            if (ystyle) nodev[i][j][1] = vy;
            if (zstyle) nodev[i][j][2] = vz;
          } else {
            if (xstyle) nodev[i][j][0] += vx;
            if (ystyle) nodev[i][j][1] += vy;
            if (zstyle) nodev[i][j][2] += vz;
          }
      }

  // set velocities via variables
  
  } else {
    if (xstyle == EQUAL) vx = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM) {
      if (vfield) input->variable->compute_atom(xvar,igroup,&vfield[0][0],3,0);
      else input->variable->compute_atom(xvar,igroup,NULL,3,0);
    }
    if (ystyle == EQUAL) vy = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM) {
      if (vfield) input->variable->compute_atom(yvar,igroup,&vfield[0][1],3,0);
      else input->variable->compute_atom(yvar,igroup,NULL,3,0);
    }
    if (zstyle == EQUAL) vz = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM) {
      if (vfield) input->variable->compute_atom(zvar,igroup,&vfield[0][2],3,0);
      else input->variable->compute_atom(zvar,igroup,NULL,3,0);
    }

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        if (sum_flag == 0) {
          if (xstyle == ATOM) v[i][0] = vfield[i][0];
          else if (xstyle) v[i][0] = vx;
          if (ystyle == ATOM) v[i][1] = vfield[i][1];
          else if (ystyle) v[i][1] = vy;
          if (zstyle == ATOM) v[i][2] = vfield[i][2];
          else if (zstyle) v[i][2] = vz;
        } else {
          if (xstyle == ATOM) v[i][0] += vfield[i][0];
          else if (xstyle) v[i][0] += vx;
          if (ystyle == ATOM) v[i][1] += vfield[i][1];
          else if (ystyle) v[i][1] += vy;
          if (zstyle == ATOM) v[i][2] += vfield[i][2];
          else if (zstyle) v[i][2] += vz;
        }
      }

    for (int i = 0; i < nelocal; i++) 
      if (group_flag == NODE) { 
        for (int j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit) 
            if (sum_flag == 0) {
              if (xstyle) nodev[i][j][0] = vx;
              if (ystyle) nodev[i][j][1] = vy;
              if (zstyle) nodev[i][j][2] = vz;
            } else {
              if (xstyle) nodev[i][j][0] += vx;
              if (ystyle) nodev[i][j][1] += vy;
              if (zstyle) nodev[i][j][2] += vz;
            }
      } else if (emask[i] & groupbit) {
        for (int j = 0; j < npe; j++)
          if (sum_flag == 0) {
            if (xstyle) nodev[i][j][0] = vx;
            if (ystyle) nodev[i][j][1] = vy;
            if (zstyle) nodev[i][j][2] = vz;
          } else {
            if (xstyle) nodev[i][j][0] += vx;
            if (ystyle) nodev[i][j][1] += vy;
            if (zstyle) nodev[i][j][2] += vz;
          }
      }
  }

  // clean up

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
}

/*---------------------------------------------------------------------------
  apply a ramped set of velocities
  -----------------------------------------------------------------------------*/

void Velocity::ramp(int narg, char **arg)
{
  // set scale factors
  xscale = yscale = zscale = 1.0;

  int v_dim;
  if (strcmp(arg[0],"vx") == 0) v_dim = 0;
  else if (strcmp(arg[0],"vy") == 0) v_dim=1;
  else if (strcmp(arg[0],"vz") == 0) v_dim=2;
  else error->all(FLERR,"Illegal velocity command");

  if (v_dim == 2 && domain->dimension == 2)
    error->all(FLERR,"Velocity ramp in z for a 2d problem");

  double v_lo, v_hi;
  if (v_dim == 0) {
    v_lo = xscale*universe->numeric(FLERR,arg[1]);
    v_hi = xscale*universe->numeric(FLERR,arg[2]);
  } else if (v_dim == 1) {
    v_lo = yscale*universe->numeric(FLERR,arg[1]);
    v_hi = yscale*universe->numeric(FLERR,arg[2]);
  } else if (v_dim == 2) {
    v_lo = zscale*universe->numeric(FLERR,arg[1]);
    v_hi = zscale*universe->numeric(FLERR,arg[2]);
  }

  int coord_dim;
  if (strcmp(arg[3],"x") == 0) coord_dim=0;
  else if (strcmp(arg[3],"y") == 0) coord_dim=1;
  else if (strcmp(arg[3],"z") == 0) coord_dim=2;
  else error->all(FLERR,"Illegal velocity command");

  double coord_lo,coord_hi;
  if (coord_dim == 0) {
    coord_lo = xscale*universe->numeric(FLERR,arg[4]);
    coord_hi = xscale*universe->numeric(FLERR,arg[5]);
  } else if (coord_dim == 1) {
    coord_lo = yscale*universe->numeric(FLERR,arg[4]);
    coord_hi = yscale*universe->numeric(FLERR,arg[5]);
  } else if (coord_dim == 2) {
    coord_lo = zscale*universe->numeric(FLERR,arg[4]);
    coord_hi = zscale*universe->numeric(FLERR,arg[5]);
  }

  // vramp = ramped velocity component for v_dim
  // add or set based on sum_flag

  // add velocity to atoms

  double **x = atom->x;
  double **v = atom->v;
  int *amask = atom->mask;
  int nlocal = atom->nlocal;

  double fraction, vramp;

  for (int i = 0; i < nlocal; i++) {
    if (amask[i] & groupbit) {
      fraction = (x[i][coord_dim]-coord_lo) / (coord_hi - coord_lo);
      fraction = MAX(fraction,0.0);
      fraction = MIN(fraction,1.0);
      vramp = v_lo + fraction*(v_hi-v_lo);
      if (sum_flag) v[i][v_dim] += vramp;
      else v[i][v_dim] = vramp;
    }
  }

  // add velocity to element nodes

  double ***nodev = element->nodev;
  double ***nodex = element->nodex;
  int **nodemask = element->nodemask;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int npe = element->npe;

  for (int i = 0; i < nelocal; i++) {
    if (group_flag = NODE) {
      for (int j = 0; j < npe; j++) 
        if (nodemask[i][j] & groupbit) {
          fraction = (nodex[i][j][coord_dim]-coord_lo) / (coord_hi - coord_lo);
          fraction = MAX(fraction,0.0);
          fraction = MIN(fraction,1.0);
          vramp = v_lo + fraction*(v_hi-v_lo);
          if (sum_flag) nodev[i][j][v_dim] += vramp;
          else nodev[i][j][v_dim] = vramp;
        }
    } else if (emask[i] & groupbit) {
      for (int j = 0; j < npe; j++) 
        if (nodemask[i][j] & groupbit) {
          fraction = (nodex[i][j][coord_dim]-coord_lo) / (coord_hi - coord_lo);
          fraction = MAX(fraction,0.0);
          fraction = MIN(fraction,1.0);
          vramp = v_lo + fraction*(v_hi-v_lo);
          if (sum_flag) nodev[i][j][v_dim] += vramp;
          else nodev[i][j][v_dim] = vramp;
        }
    } 
  }
}
