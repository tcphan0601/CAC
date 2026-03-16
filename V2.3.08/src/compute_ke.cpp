#include <mpi.h>
#include <string.h>
#include "compute_ke.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeKE::ComputeKE(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute ke command");
  if (igroup) error->all(FLERR,"Compute ke must use group all");
  
  scalar_flag = 1;
  extscalar = 1;
  
  if (narg == 3) {
    pairflag = 1;
    fixflag = 1;
	sumflag = 1;
  } else {
	pairflag = 1;
    fixflag = 1;
	int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"simple") == 0) sumflag = 1;
      else if (strcmp(arg[iarg],"accurate") == 0) sumflag = 2;
      else error->all(FLERR,"Illegal compute ke command");
      iarg++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeKE::init()
{
  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeKE::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  
  double **v = atom->v;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int *etype  = element->etype;
  int *ectype = element->ctype;
  int nlocal = atom->nlocal;
  int i,j;
  
  double ke = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      ke += mass[type[i]] *
        (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
		
  double ***nodev = element->nodev;
  int **nodemask = element->nodemask;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int npe = element->npe;
  double weight;
  int *nintpl = element->nintpl;
  
  for (i = 0; i < nelocal; i++) {
	if (sumflag == 1) {
      weight =  double(nintpl[etype[i]]) / double(npe);
	  for (j = 0; j < npe; j++) {
		if (nodemask[i][j] & groupbit) {
	      ke += mass[ectype[i]] * weight *
		    (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]);
	    }
	  } 
	}else {
	    error->all(FLERR,"accurate ke has not been implemented yet");
	}	
  }

  MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}
