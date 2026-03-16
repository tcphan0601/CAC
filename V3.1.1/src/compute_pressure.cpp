#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "compute_pressure.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "error.h"

#include "element.h"
#include "atom.h"
#include "comm.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

ComputePressure::ComputePressure(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg), 
  vptr(nullptr), id_temp(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal compute pressure command");
  if (igroup) error->all(FLERR, "Compute pressure must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  if (strcmp(arg[3], "NULL") == 0) id_temp = nullptr;
  else {
    int n = strlen(arg[3]) + 1;
    id_temp = new char[n];
    strcpy(id_temp, arg[3]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR, "Could not find compute pressure temperature ID");
    if (modify->compute[icompute]->tempflag == 0)
      error->all(FLERR, "Compute pressure temperature ID does not "
                 "compute temperature");
  }

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg], "pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg], "fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg], "virial") == 0) {
        pairflag = 1;
        fixflag = 1;
      } else error->all(FLERR, "Illegal compute pressure command");
      iarg++;
    }
  }

  // error check

  if (keflag && id_temp == nullptr)
    error->all(FLERR, "Compute pressure requires temperature ID "
	       "to include kinetic energy");

  vector = new double[6];
  nvirial = 0;
  vptr = nullptr;
}

/*  ----------------------------------------------------------------------  */

ComputePressure::~ComputePressure()
{
  delete [] id_temp;
  delete [] vector;
  delete [] vptr;
}

/*  ----------------------------------------------------------------------  */

void ComputePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (keflag) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR, "Could not find compute pressure temperature ID");
    temperature = modify->compute[icompute];
  }

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = nullptr;

  if (pairflag && force->pair) nvirial++;
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->virial_flag) nvirial++;

  if (nvirial) {
    vptr = new double * [nvirial];
    nvirial = 0;
    if (pairflag && force->pair) vptr[nvirial++] = force->pair->virial;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->virial_flag)
          vptr[nvirial++] = modify->fix[i]->virial;
  }
}

/*  ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
-------------------------------------------------------------------------  */

double ComputePressure::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR, "Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t = temperature->compute_scalar();
    else t = temperature->scalar;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3, 3);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
          virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2, 2);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
          virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}

/*  ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
   -------------------------------------------------------------------------  */

void ComputePressure::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR, "Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double *ke_tensor;
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    ke_tensor = temperature->vector;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6, 3);
    if (keflag) 
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4, 2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}

/*  ----------------------------------------------------------------------  */

void ComputePressure::virial_compute(int n, int ndiag)
{
  int i, j;
  double v[6], *vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v, virial, n, MPI_DOUBLE, MPI_SUM, world);

  // LJ long-range tail correction, only if pair contributions are included

  if (force->pair && pairflag && force->pair->tail_flag)
    for (i = 0; i < ndiag; i++) virial[i] += force->pair->ptail * inv_volume;
}

/*  ----------------------------------------------------------------------  */

void ComputePressure::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  int n = strlen(id_new) + 1;
  id_temp = new char[n];
  strcpy(id_temp, id_new);
}

