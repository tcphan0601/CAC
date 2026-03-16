#include <string.h>
#include "fix_enforce2d.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

enum{NODE,ELEMENT};

/* ---------------------------------------------------------------------- */

FixEnforce2D::FixEnforce2D(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  flist(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal fix enforce2d command");

  // optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix enforce2d command");
      if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR, "Illegal fix enforce2d command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix enforce2d command");
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

FixEnforce2D::~FixEnforce2D()
{
}

/* ---------------------------------------------------------------------- */

int FixEnforce2D::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::init()
{
  if (domain->dimension == 3)
    error->all(FLERR,"Cannot use fix enforce2d with 3d simulation");
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force(int vflag)
{
  if (atom->nlocal) {
    double **v = atom->v;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        v[i][2] = 0.0;
        f[i][2] = 0.0;
      }
  }

  if (element->nlocal) {
    double ***nodev = element->nodev;
    double ***nodef = element->nodef;
    int *mask = element->mask;
    int **nodemask = element->nodemask;
    int nlocal = element->nlocal;
    int npe = element->npe;
    for (int i = 0; i < nlocal; i++) 
      if (group_flag) {
        if (mask[i] & groupbit)
          for (int j = 0; j < npe; j++) {
            nodev[i][j][2] = 0.0;
            nodef[i][j][2] = 0.0;
          }
      } else {
        for (int j = 0; j < npe; j++) {
          if (nodemask[i][j] & groupbit) {
            nodev[i][j][2] = 0.0;
            nodef[i][j][2] = 0.0;
          }
        }
      }
  }
}
