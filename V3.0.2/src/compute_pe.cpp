#include <mpi.h>
#include <string.h>
#include "compute_pe.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

ComputePE::ComputePE(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal compute pe command");
  if (igroup) error->all(FLERR, "Compute pe must use group all");

  scalar_flag = 1;
  extscalar = 1;
  peflag = 1;
  timeflag = 1;

  if (narg == 3) {
    pairflag = 1;
    fixflag = 1;
  } else {
    pairflag = 0;
    fixflag = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg], "fix") == 0) fixflag = 1;
      else error->all(FLERR, "Illegal compute pe command");
      iarg++;
    }
  }
}

/*  ----------------------------------------------------------------------  */

double ComputePE::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR, "Energy was not tallied on needed timestep");

  double one = 0.0;
  if (pairflag && force->pair)
    one += force->pair->eng_vdwl + force->pair->eng_coul;

  MPI_Allreduce(&one, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  if (pairflag && force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    scalar += force->pair->etail / volume;
  }

//  if (fixflag && modify->n_thermo_energy) scalar += modify->thermo_energy();

  return scalar;
}
