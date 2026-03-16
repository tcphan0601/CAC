#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_viscous.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "error.h"
#include "group.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

/*---------------------------------------------------------------------------*/

FixViscous::FixViscous(CAC *cac, int narg, char **arg) :
  Fix(cac,narg,arg)
{

  if (narg < 4) error->all(FLERR,"Illegal fix viscous command");

  double gamma_one = universe->numeric(FLERR,arg[3]);

  int n = atom->ntypes;
  gamma = new double[n+1];
  for (int i = 1; i <= n; i++) gamma[i] = gamma_one;

  // optional args

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix viscous command");
      int itype = universe->inumeric(FLERR,arg[iarg+1]);
      double scale = universe->numeric(FLERR,arg[iarg+2]);
      if (itype <= 0 || itype > n)
        error->all(FLERR,"Illegal fix viscous command");
      gamma[itype] = gamma_one * scale;
      iarg += 3;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve command");
      if (strcmp(arg[iarg+1],"atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1],"node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR,"Illegal fix nve command");
      iarg += 2;

    } else error->all(FLERR,"Illegal fix viscous command");
  }


}

/*---------------------------------------------------------------------------*/

FixViscous::~FixViscous()
{
  delete [] gamma;
}

/*---------------------------------------------------------------------------*/

int FixViscous::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask; 
}

/*----------------------------------------------------------------------------*/

void FixViscous::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else error->all(FLERR,"only verlet integrate style for now");
}

/*---------------------------------------------------------------------------*/

void FixViscous::post_force(int vflag)
{

  // apply drag force to atoms/elements in group
  // diretion is opposed to velocity vector
  // magnitude depends on atom/element type
  
  double drag;

  double **av = atom->v;
  double **af = atom->f;
  int *amask = atom->mask;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++)
    if (amask[i] & groupbit) {
      drag = gamma[atype[i]];
      af[i][0] -= drag*av[i][0];
      af[i][1] -= drag*av[i][1];
      af[i][2] -= drag*av[i][2];
    }

  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int **nodemask = element->nodemask;
  int *emask = element->mask;
  int *ctype = element->ctype;
  int *etype = element->etype;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  if (group_flag == Group::NODE) {
    for (int i = 0; i < nelocal; i++) {
      drag = gamma[ctype[i]];
      for (int j = 0; j < npe[etype[i]]; j++) 
        if (nodemask[i][j] & groupbit) {
          nodef[i][j][0] -= drag*nodev[i][j][0];
          nodef[i][j][1] -= drag*nodev[i][j][1];
          nodef[i][j][2] -= drag*nodev[i][j][2];
        }
    }
  } else if (group_flag == Group::ELEMENT) {
    for (int i = 0; i < nelocal; i++) {
      if (!(emask[i] & groupbit)) continue;
      drag = gamma[ctype[i]];
      for (int j = 0; j < npe[etype[i]]; j++) {
        nodef[i][j][0] -= drag*nodev[i][j][0];
        nodef[i][j][1] -= drag*nodev[i][j][1];
        nodef[i][j][2] -= drag*nodev[i][j][2];
      }
    }
  }
}

