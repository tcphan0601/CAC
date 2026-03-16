#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "scale.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "universe.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "variable.h"

using namespace CAC_NS;

#define   EPSILON 1.0e-4

/* ---------------------------------------------------------------------- */

Scale::Scale(CAC *cac) : Pointers(cac)
{
}

/* ---------------------------------------------------------------------- */

Scale::~Scale()
{
}

/* ---------------------------------------------------------------------- */

void Scale::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Scale command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal scale command");
  if (domain->triclinic) error->all(FLERR,"Scale command has not been tested with triclinic box yet");

  double scale_factor[3];
  scale_factor[0] = scale_factor[1] = scale_factor[2] = -1.0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal scale command");
      scale_factor[0] = universe->numeric(FLERR,arg[iarg+1]);
      if (scale_factor[0] <= 0.0) error->all(FLERR,"Scale factor must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal scale command");
      scale_factor[1] = universe->numeric(FLERR,arg[iarg+1]);
      if (scale_factor[1] <= 0.0) error->all(FLERR,"Scale factor must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal scale command");
      scale_factor[2] = universe->numeric(FLERR,arg[iarg+1]);
      if (scale_factor[2] <= 0.0) error->all(FLERR,"Scale factor must be positive");
      if (domain->dimension == 2) {
        if (comm->me == 0) error->warning(FLERR,"Scale factor for z direction is ignored for 2D simulation");
        scale_factor[2] = 1.0;
      }
      iarg += 2;
    }
  }

  /* no read/write restart yet
     if (modify->nfix_restart_peratom)
     error->all(FLERR,"Cannot scale after "
     "reading restart file with per-atom info");
     */

  if (comm->me == 0 && screen) fprintf(screen,"Scaling atoms/elements ...\n");
  if (comm->me == 0 && logfile) fprintf(logfile,"Scaling atoms/elements ...\n");

  // scale atoms coordinates, elements' centers and nodes coordinates
  // also scale domain and subdomain boxes 

  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;
  double **ex = element->x;
  double **ax = atom->x;
  double ***nodex = element->nodex;
  for (int dim = 0; dim < domain->dimension; dim++) {
    if (scale_factor[dim] > 0) {
      for (int i = 0; i < nalocal; i++) 
        ax[i][dim] *= scale_factor[dim];
      for (int i = 0; i < nelocal; i++) {
        ex[i][dim] *= scale_factor[dim];
        for (int j = 0; j < npe[etype[i]]; j++) {
          nodex[i][j][dim] *= scale_factor[dim];
        }
      }
      domain->boxlo[dim] *= scale_factor[dim];
      domain->boxhi[dim] *= scale_factor[dim];
      domain->sublo[dim] *= scale_factor[dim];
      domain->subhi[dim] *= scale_factor[dim];
      domain->prd[dim] *= scale_factor[dim];
      domain->prd_half[dim] *= scale_factor[dim];
      if (dim == 0) {
        domain->xprd *= scale_factor[0];
        domain->xprd_half *= scale_factor[0];
      }
      if (dim == 1) {
        domain->yprd *= scale_factor[1];
        domain->yprd_half *= scale_factor[1];
      }
      if (dim == 2) {
        domain->yprd *= scale_factor[2];
        domain->yprd_half *= scale_factor[2];
      }
    } 
  }
}

