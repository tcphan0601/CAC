#include "min_sd.h"
#include <cmath>
#include "error.h"
#include "update.h"
#include "output.h"
#include "timer.h"

using namespace CAC_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

/* ---------------------------------------------------------------------- */

MinSD::MinSD(CAC *cac) : MinLineSearch(cac) {}

/* ----------------------------------------------------------------------
   minimization via steepest descent
------------------------------------------------------------------------- */

int MinSD::iterate(int maxiter)
{
  int i, m, n, fail, ntimestep;
  double fdotf;
  double *fatom, *hatom;

  // initialize working vectors

  for (i = 0; i < nvec; i++) h[i] = fvec[i];
  /*
  if (nextra_atom)
    for (m = 0; m < nextra_atom; m++) {
      fatom = fextra_atom[m];
      hatom = hextra_atom[m];
      n = extra_nlen[m];
      for (i = 0; i < n; i++) hatom[i] = fatom[i];
    }
  if (nextra_global)
    for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];
  */
  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // line minimization along h from current position x
    // h = downhill gradient direction

    eprevious = ecurrent;
    fail = (this->*linemin)(ecurrent, alpha_final);
    if (fail) return fail;

    // function evaluation criterion

    if (neval >= update->max_eval) return MAXEVAL;

    // energy tolerance criterion

    if (fabs(ecurrent-eprevious) <
        update->etol * 0.5 * (fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
      return ETOL;

    // force tolerance criterion

    fdotf = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotf = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf = fnorm_inf();   // infinite force norm
      else if (normstyle == TWO) fdotf = fnorm_sqr();   // Euclidean force 2-norm
      else error->all(FLERR, "Illegal min_modify command");
      if (fdotf < update->ftol*update->ftol) return FTOL;
    }

    // set new search direction h to f = -Grad(x)

    for (i = 0; i < nvec; i++) h[i] = fvec[i];
    /*
    if (nextra_atom)
      for (m = 0; m < nextra_atom; m++) {
        fatom = fextra_atom[m];
        hatom = hextra_atom[m];
        n = extra_nlen[m];
        for (i = 0; i < n; i++) hatom[i] = fatom[i];
      }
    if (nextra_global)
      for (i = 0; i < nextra_global; i++) hextra[i] = fextra[i];
    */
    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
