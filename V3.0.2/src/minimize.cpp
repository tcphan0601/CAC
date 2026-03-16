#include "minimize.h"
#include "domain.h"
#include "update.h"
#include "min.h"
#include "finish.h"
#include "timer.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

Minimize::Minimize(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void Minimize::command(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Illegal minimize command");

  if (domain->box_exist == 0)
    error->all(FLERR, "Minimize command before simulation box is defined");

  // ignore minimize command, if walltime limit was already reached
  if (timer->is_timeout()) return;

  update->etol = universe->numeric(FLERR, arg[0]);
  update->ftol = universe->numeric(FLERR, arg[1]);
  update->nsteps = universe->inumeric(FLERR, arg[2]);
  update->max_eval = universe->inumeric(FLERR, arg[3]);

  if (update->etol < 0.0 || update->ftol < 0.0)
    error->all(FLERR, "Illegal minimize command");

  update->whichflag = 2;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  if (update->laststep < 0)
    error->all(FLERR, "Too many iterations");

  cac->init();
  timer->init_timeout();
  update->minimize->setup();

  timer->init();
  timer->barrier_start();
  update->minimize->run(update->nsteps);
  timer->barrier_stop();

  update->minimize->cleanup();

  Finish finish(cac);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}
