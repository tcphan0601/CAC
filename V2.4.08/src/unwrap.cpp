#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "unwrap.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "modify.h"
#include "domain.h"
#include "group.h"
#include "comm.h"
#include "universe.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "variable.h"

using namespace CAC_NS;

#define   EPSILON 1.0e-4


/* ---------------------------------------------------------------------- */

Unwrap::Unwrap(CAC *cac) : Pointers(cac)
{
}

/* ---------------------------------------------------------------------- */

Unwrap::~Unwrap()
{
}

/* ---------------------------------------------------------------------- */

void Unwrap::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Unwrap command before simulation box is defined");
  if (narg < 3) error->all(FLERR,"Illegal unwrap command");
  if (domain->triclinic) error->all(FLERR,"Unwrap command has not been tested with triclinic box yet");

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find unwrap group ID");
  int groupbit = group->bitmask[igroup];

  if (comm->me == 0 && screen) fprintf(screen,"Unwrapping atoms/elements ...\n");
  if (comm->me == 0 && logfile) fprintf(logfile,"Unwrapping atoms/elements ...\n");

  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;
  imageint *aimage = atom->image;
  imageint *eimage = element->image;
  double **ex = element->x;
  double **ax = atom->x;
  double ***nodex = element->nodex;

  int pbc_flag[3];
  pbc_flag[0] = pbc_flag[1] = pbc_flag[2] = 0;
  double boxlo[3],boxhi[3];
  double *prd = domain->prd;

  // unwrap from atom/element periodic image values

  if (strcmp(arg[1],"image") == 0) {
    error->all(FLERR,"Image unwrap has not been finished yet");
  } 

  // unwrap base on bound box
  
  else if (strcmp(arg[1],"box") == 0) {
    int iarg = 2;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"x") == 0) {
        pbc_flag[0] = 1;
        if (iarg+3 > narg) error->all(FLERR,"Illegal unwrap command");
        boxlo[0] = universe->numeric(FLERR,arg[iarg+1]);
        boxhi[0] = universe->numeric(FLERR,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg],"y") == 0) {
        pbc_flag[1] = 1;
        if (iarg+3 > narg) error->all(FLERR,"Illegal unwrap command");
        boxlo[1] = universe->numeric(FLERR,arg[iarg+1]);
        boxhi[1] = universe->numeric(FLERR,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg],"z") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal unwrap command");
        if (domain->dimension == 2) {
          if (comm->me == 0) error->warning(FLERR,"Unwrap factor for z direction is ignored for 2D simulation");
        } else {
          pbc_flag[2] = 1;
          boxlo[2] = universe->numeric(FLERR,arg[iarg+1]);
          boxhi[2] = universe->numeric(FLERR,arg[iarg+2]);
        }
        iarg += 3;
      } else error->all(FLERR,"Invalid unwrap direction");
    }

    imageint idim,otherdims;
    for (int i = 0; i < nalocal; i++) {
      if (pbc_flag[0]) {
        if (ax[i][0] < boxlo[0]) {
          ax[i][0] += prd[0];
          idim = aimage[i] & IMGMASK;
          otherdims = aimage[i] ^ idim;
          idim--;
          idim &= IMGMASK;
          aimage[i] = otherdims | idim;
        }
        if (ax[i][0] > boxhi[0]) {
          ax[i][0] -= prd[0];
          idim = aimage[i] & IMGMASK;
          otherdims = aimage[i] ^ idim;
          idim++;
          idim &= IMGMASK;
          aimage[i] = otherdims | idim;
        }

        if (ax[i][0] < boxlo[0] || ax[i][0] > boxhi[0]) {
          printf("i = %d x = %g lo/hi = %g %g\n",i,ax[i][0],boxlo[0],boxhi[0]);
          error->one(FLERR,"Bound box does not cover all atoms");
        }
      }

      if (pbc_flag[1]) {
        if (ax[i][1] < boxlo[1]) {
          ax[i][1] += prd[1];
          idim = (aimage[i] >> IMGBITS) & IMGMASK;
          otherdims = aimage[i] ^ (idim << IMGBITS);
          idim--;
          idim &= IMGMASK;
          aimage[i] = otherdims | (idim << IMGBITS);
        }
        if (ax[i][1] > boxhi[1]) {
          ax[i][1] -= prd[1];
          idim = (aimage[i] >> IMGBITS) & IMGMASK;
          otherdims = aimage[i] ^ (idim << IMGBITS);
          idim++;
          idim &= IMGMASK;
          aimage[i] = otherdims | (idim << IMGBITS);
        }

        if (ax[i][1] < boxlo[1] || ax[i][1] > boxhi[1]) 
          error->one(FLERR,"Bound box does not cover all atoms");
      }

      if (pbc_flag[2]) {
        if (ax[i][2] < boxlo[2]) {
          ax[i][2] += prd[2];
          idim = aimage[i] >> IMG2BITS;
          otherdims = aimage[i] ^ (idim << IMG2BITS);
          idim--;
          idim &= IMGMASK;
          aimage[i] = otherdims | (idim << IMG2BITS);
        }
        if (ax[i][2] > boxhi[2]) {
          ax[i][2] -= prd[2];
          idim = aimage[i] >> IMG2BITS;
          otherdims = aimage[i] ^ (idim << IMG2BITS);
          idim++;
          idim &= IMGMASK;
          aimage[i] = otherdims | (idim << IMG2BITS);
        }

        if (ax[i][2] < boxlo[2] || ax[i][2] > boxhi[2]) 
          error->one(FLERR,"Bound box does not cover all atoms");
      }

    }

    for (int i = 0; i < nelocal; i++) {
      if (pbc_flag[0]) {
        if (ex[i][0] < boxlo[0]) {
          ex[i][0] += prd[0];
          element->evec->update_node_coord(i,0);
          idim = eimage[i] & IMGMASK;
          otherdims = eimage[i] ^ idim;
          idim--;
          idim &= IMGMASK;
          eimage[i] = otherdims | idim;
        }
        if (ex[i][0] > boxhi[0]) {
          ex[i][0] -= prd[0];
          element->evec->update_node_coord(i,0);
          idim = eimage[i] & IMGMASK;
          otherdims = eimage[i] ^ idim;
          idim++;
          idim &= IMGMASK;
          eimage[i] = otherdims | idim;
        }

        if (ex[i][0] < boxlo[0] || ex[i][0] > boxhi[0]) 
          error->one(FLERR,"Bound box does not cover all atoms");
      }

      if (pbc_flag[1]) {
        if (ex[i][1] < boxlo[1]) {
          ex[i][1] += prd[1];
          element->evec->update_node_coord(i,1);
          idim = (eimage[i] >> IMGBITS) & IMGMASK;
          otherdims = eimage[i] ^ (idim << IMGBITS);
          idim--;
          idim &= IMGMASK;
          eimage[i] = otherdims | (idim << IMGBITS);
        }
        if (ex[i][1] > boxhi[1]) {
          ex[i][1] -= prd[1];
          element->evec->update_node_coord(i,1);
          idim = (eimage[i] >> IMGBITS) & IMGMASK;
          otherdims = eimage[i] ^ (idim << IMGBITS);
          idim++;
          idim &= IMGMASK;
          eimage[i] = otherdims | (idim << IMGBITS);
        }

        if (ex[i][1] < boxlo[1] || ex[i][1] > boxhi[1]) 
          error->one(FLERR,"Bound box does not cover all atoms");
      }

      if (pbc_flag[2]) {
        if (ex[i][2] < boxlo[2]) {
          ex[i][2] += prd[2];
          element->evec->update_node_coord(i,2);
          idim = eimage[i] >> IMG2BITS;
          otherdims = eimage[i] ^ (idim << IMG2BITS);
          idim--;
          idim &= IMGMASK;
          eimage[i] = otherdims | (idim << IMG2BITS);
        }
        if (ex[i][2] > boxhi[2]) {
          ex[i][2] -= prd[2];
          element->evec->update_node_coord(i,2);
          idim = eimage[i] >> IMG2BITS;
          otherdims = eimage[i] ^ (idim << IMG2BITS);
          idim++;
          idim &= IMGMASK;
          eimage[i] = otherdims | (idim << IMG2BITS);
        }

        if (ex[i][2] < boxlo[2] || ex[i][2] > boxhi[2]) 
          error->one(FLERR,"Bound box does not cover all atoms");
      }
    }
  } else 
    error->all(FLERR,"Invalid unwrap style");
}
