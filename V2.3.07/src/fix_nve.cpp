#include <stdio.h>
#include <string.h>
#include "fix_nve.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

enum{ATOM,NODE,ELEMENT};

/*-----------------------------------------------------------------*/

FixNVE::FixNVE(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
  group_flag = NODE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve command");
      if (strcmp(arg[iarg+1],"atom") == 0) group_flag = ATOM;
      else if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR,"Illegal fix nve command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix nve command");
  }
}

/* ---------------------------------------------------------------------- */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
    
   ------------------------------------------------------------------------- */

void FixNVE::initial_integrate(int vflag)
{
  double dtfm;
  double xtmp,ytmp,ztmp;
  double **x,**v,**f;
  double ***nodex,***nodev,***nodef;
  double *mass;
  int *mask,**nodemask;
  int *type;
  int nlocal,npe;

  // update v and x of atoms in group

  if (atom->nlocal) {
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
  }

  // update nodev,nodex, and x of elements in group

  if (element->nlocal && group_flag) {
    x = element->x;
    nodex = element->nodex;
    nodev = element->nodev;
    nodef = element->nodef;
    mass = atom->mass;
    type = element->ctype;
    mask = element->mask;
    nodemask = element->nodemask;
    nlocal = element->nlocal;
    npe = element->npe;

    if (group_flag == NODE) {
      for (int i = 0; i < nlocal; i++) {
        xtmp = ytmp = ztmp = 0.0;
        dtfm = dtf / (mass[type[i]]);
        for (int j = 0; j < npe; j++) {
          if (nodemask[i][j] & groupbit) {
            nodev[i][j][0] += dtfm * nodef[i][j][0];
            nodev[i][j][1] += dtfm * nodef[i][j][1];
            nodev[i][j][2] += dtfm * nodef[i][j][2];
            nodex[i][j][0] += dtv * nodev[i][j][0];
            nodex[i][j][1] += dtv * nodev[i][j][1];
            nodex[i][j][2] += dtv * nodev[i][j][2];
          }
          xtmp += nodex[i][j][0];
          ytmp += nodex[i][j][1];
          ztmp += nodex[i][j][2];
        }
        x[i][0] = xtmp/npe;
        x[i][1] = ytmp/npe;
        x[i][2] = ztmp/npe;
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue; 
        xtmp = ytmp = ztmp = 0.0;
        dtfm = dtf / (mass[type[i]]);
        for (int j = 0; j < npe; j++) {
          nodev[i][j][0] += dtfm * nodef[i][j][0];
          nodev[i][j][1] += dtfm * nodef[i][j][1];
          nodev[i][j][2] += dtfm * nodef[i][j][2];
          nodex[i][j][0] += dtv * nodev[i][j][0];
          nodex[i][j][1] += dtv * nodev[i][j][1];
          nodex[i][j][2] += dtv * nodev[i][j][2];
          xtmp += nodex[i][j][0];
          ytmp += nodex[i][j][1];
          ztmp += nodex[i][j][2];
        }
        x[i][0] = xtmp/npe;
        x[i][1] = ytmp/npe;
        x[i][2] = ztmp/npe;
      }
    }
  }
}
/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate()
{
  double dtfm;
  double **v,**f;
  double ***nodev,***nodef;
  double *mass;
  int *mask,**nodemask;
  int *type;
  int nlocal,npe;

  // update v of atoms in group

  if (atom->nlocal) {
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
  }

  // update nodev of elements in group

  if (element->nlocal && group_flag) {
    nodev = element->nodev;
    nodef = element->nodef;
    mass = atom->mass;
    type = element->ctype;
    mask = element->mask;
    nodemask = element->nodemask;
    nlocal = element->nlocal;
    npe = element->npe;
    if (group_flag == NODE) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / (mass[type[i]]);
        for (int j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit) {
            nodev[i][j][0] += dtfm * nodef[i][j][0];
            nodev[i][j][1] += dtfm * nodef[i][j][1];
            nodev[i][j][2] += dtfm * nodef[i][j][2];
          }
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        dtfm = dtf / (mass[type[i]]);
        for (int j = 0; j < npe; j++) { 
          nodev[i][j][0] += dtfm * nodef[i][j][0];
          nodev[i][j][1] += dtfm * nodef[i][j][1];
          nodev[i][j][2] += dtfm * nodef[i][j][2];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVE::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
