#include <stdio.h>
#include <string.h>
#include "fix_nve.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "error.h"
#include "group.h"

using namespace CAC_NS;
using namespace FixConst;


/* ----------------------------------------------------------------- */

FixNVE::FixNVE(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix nve command");
      if (strcmp(arg[iarg+1], "atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1], "node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1], "element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR, "Illegal fix nve command");
      iarg += 2;
    } else error->all(FLERR, "Illegal fix nve command");
  }
}

/*  ----------------------------------------------------------------------  */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixNVE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/*  ----------------------------------------------------------------------
    
   -------------------------------------------------------------------------  */

void FixNVE::initial_integrate(int vflag)
{
  double dtfm;
  double **x, **v, **f;
  double ****nodex, ****nodev, ****nodef;
  double *mass;
  int *mask, ***nodemask;
  int *type;
  int *etype;
  int **ctype;
  int nlocal;
  int *npe;
  int *apc;

  // update v and x of atoms in group

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  type = atom->type;
  mask = atom->mask;
  nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }


  // update nodev, nodex, and x of elements in group

  x = element->x;
  nodex = element->nodex;
  nodev = element->nodev;
  nodef = element->nodef;
  ctype = element->ctype;
  etype = element->etype;
  mask = element->mask;
  nodemask = element->nodemask;
  nlocal = element->nlocal;
  npe = element->npe;
  apc = element->apc;

  if (group_flag == Group::NODE) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < apc[etype[i]]; j++) {
        dtfm = dtf / (mass[ctype[i][j]]);
        for (int k = 0; k < npe[etype[i]]; k++) {
          if (nodemask[i][j][k] & groupbit) {
            nodev[i][j][k][0] += dtfm * nodef[i][j][k][0];
            nodev[i][j][k][1] += dtfm * nodef[i][j][k][1];
            nodev[i][j][k][2] += dtfm * nodef[i][j][k][2];
            nodex[i][j][k][0] += dtv * nodev[i][j][k][0];
            nodex[i][j][k][1] += dtv * nodev[i][j][k][1];
            nodex[i][j][k][2] += dtv * nodev[i][j][k][2];
          }
        }
      }
    }
  } else if (group_flag == Group::ELEMENT) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue; 
      for (int j = 0; j < apc[etype[i]]; j++) {
        dtfm = dtf / (mass[ctype[i][j]]);
        for (int k = 0; k < npe[etype[i]]; k++) {
          nodev[i][j][k][0] += dtfm * nodef[i][j][k][0];
          nodev[i][j][k][1] += dtfm * nodef[i][j][k][1];
          nodev[i][j][k][2] += dtfm * nodef[i][j][k][2];
          nodex[i][j][k][0] += dtv * nodev[i][j][k][0];
          nodex[i][j][k][1] += dtv * nodev[i][j][k][1];
          nodex[i][j][k][2] += dtv * nodev[i][j][k][2];
        }
      }
    }
  }
  element->evec->update_center_coord();
}

/*  ----------------------------------------------------------------------  */

void FixNVE::final_integrate()
{
  double dtfm;
  double **v, **f;
  double ****nodev, ****nodef;
  double *mass;
  int *mask, ***nodemask;
  int *type;
  int *etype;
  int **ctype;
  int nlocal;
  int *npe;
  int *apc;

  // update v of atoms in group

  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  type = atom->type;
  mask = atom->mask;
  nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }


  // update nodev of elements in group

  nodev = element->nodev;
  nodef = element->nodef;
  mass = atom->mass;
  etype = element->etype;
  ctype = element->ctype;
  mask = element->mask;
  nodemask = element->nodemask;
  nlocal = element->nlocal;
  npe = element->npe;
  apc = element->apc;
  if (group_flag == Group::NODE) {
    for (int i = 0; i < nlocal; i++) {
      for (int j = 0; j < apc[etype[i]]; j++) {
        dtfm = dtf / (mass[ctype[i][j]]);
        for (int k = 0; k < npe[etype[i]]; k++) {
          if (nodemask[i][j][k] & groupbit) {
            nodev[i][j][k][0] += dtfm * nodef[i][j][k][0];
            nodev[i][j][k][1] += dtfm * nodef[i][j][k][1];
            nodev[i][j][k][2] += dtfm * nodef[i][j][k][2];
          }
        }
      }
    }
  } else if (group_flag == Group::ELEMENT) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      for (int j = 0; j < apc[etype[i]]; j++) {
        dtfm = dtf / (mass[ctype[i][j]]);
        for (int k = 0; k < npe[etype[i]]; k++) {
          nodev[i][j][k][0] += dtfm * nodef[i][j][k][0];
          nodev[i][j][k][1] += dtfm * nodef[i][j][k][1];
          nodev[i][j][k][2] += dtfm * nodef[i][j][k][2];
        }
      }
    }
  }
}

/*  ----------------------------------------------------------------------  */

void FixNVE::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
