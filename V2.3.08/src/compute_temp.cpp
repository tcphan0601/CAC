#include <mpi.h>
#include <cstring>
#include "compute_temp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  
  if (narg == 3) sumflag = 1;
  else {
	int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"simple") == 0) sumflag = 1;
      else if (strcmp(arg[iarg],"accurate") == 0) sumflag = 2;
      else error->all(FLERR,"Illegal compute temp command");
      iarg++;
    }
  }

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTemp::~ComputeTemp()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::init()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::dof_compute()
{
  //adjust_dof_fix();
  natoms_temp = group->count_atom(igroup) + group->count_node(igroup);
  dof = domain->dimension * natoms_temp;
  //dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTemp::compute_scalar()
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
  
  double t = 0.0;
  
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      t += mass[type[i]] *
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
	      t += mass[ectype[i]] * weight *
		    (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]);
	    }
	  } 
	}else {
	  error->all(FLERR,"accurate ke has not been implemented yet");
	}	
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::compute_vector()
{

  int i;

  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;

}
