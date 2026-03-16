#include <stdio.h>
#include <string.h>
#include "fix_nve.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "error.h"

#include "domain.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NODE,ELEMENT};

/*-----------------------------------------------------------------*/

FixNVE::FixNVE(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
  group_flag = NODE;
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
  double **x,**v,**f;
  double ***nodex,***nodev,***nodef;
  double *mass;
  int *mask,**nodemask;
  int *type;
  int nlocal,npe;

  // update v and x of atoms in group
  
  if (atom->nlocal){
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
  
  if (element->nlocal){
    x = element->x;
    nodex = element->nodex;
    nodev = element->nodev;
    nodef = element->nodef;
    mass = atom->mass;
    type = element->ctype;
    nodemask = element->nodemask;
    nlocal = element->nlocal;
    npe = element->npe;
    for (int i = 0; i < nlocal; i++) {
      x[i][0] = 0.0;
      x[i][1] = 0.0;
      x[i][2] = 0.0;
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
        x[i][0] += nodex[i][j][0];
        x[i][1] += nodex[i][j][1];
        x[i][2] += nodex[i][j][2];
      }
      x[i][0] /= npe;
      x[i][1] /= npe;
      x[i][2] /= npe;
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

  double **x,***nodex,coord[3];
  x = atom->x;
  nodex = element->nodex;

  // update v of atoms in group

  if (atom->nlocal){
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

  if (element->nlocal){
    nodev = element->nodev;
    nodef = element->nodef;
    mass = atom->mass;
    type = element->ctype;
    nodemask = element->nodemask;
    nlocal = element->nlocal;
    npe = element->npe;
    for (int i = 0; i < nlocal; i++) {
      dtfm = dtf / (mass[type[i]]);
      for (int j = 0; j < npe; j++) 
        if (nodemask[i][j] & groupbit) {
          nodev[i][j][0] += dtfm * nodef[i][j][0];
          nodev[i][j][1] += dtfm * nodef[i][j][1];
          nodev[i][j][2] += dtfm * nodef[i][j][2];
        }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVE::reset_dt(){
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
