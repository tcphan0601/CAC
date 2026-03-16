#include <stdlib.h>
#include <string.h>
#include "run.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "output.h"
#include "finish.h"
#include "universe.h"
#include "input.h"
#include "timer.h"
#include "error.h"
#include "comm.h"

using namespace CAC_NS;

#define MAXLINE 2048

/* ---------------------------------------------------------------------- */

Run::Run(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void Run::command(int narg, char **arg)
{

  if (narg < 1) error->all(FLERR,"Illegal run command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Run command before simulation box is defined");

  bigint nsteps_input = universe->bnumeric(FLERR,arg[0]);

  // parse optional args

  int uptoflag = 0;
  int startflag = 0;
  int stopflag = 0;
  bigint start,stop;
  int preflag = 1;
  int postflag = 1;
  int nevery = 0;
  int ncommands = 0;
  int first,last;

  // set nsteps as integer, using upto value if specified

  int nsteps;
  if (!uptoflag) {
    if (nsteps_input < 0 || nsteps_input > MAXSMALLINT)
      error->all(FLERR,"Invalid run command N value");
    nsteps = static_cast<int> (nsteps_input);
  }
 
  // if nevery, make copies of arg strings that are commands
  // required because re-parsing commands via input->one() will wipe out args 

  char **commands = NULL;
  if (nevery && ncommands > 0) {
    commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands],arg[i]);
      ncommands++;
    }
  }

  // perform a single run
  // use start/stop to set begin/end step
  // if pre or 1st run, do System init/setup,
  // else just init timer and setup output
  // if post, do full Finish, else just print time

  update->whichflag = 1;
  timer->init_timeout();

  if (nevery == 0) {
    update->nsteps = nsteps;
    update->firststep = update->ntimestep;
    update->laststep = update->ntimestep + nsteps;

    if (update->laststep < 0 || update->laststep > MAXBIGINT)
      error->all(FLERR,"too many time step"); 

    if (startflag) update->beginstep = start;
    else update->beginstep = update->firststep;
    if (stopflag) update->endstep = stop;
    else update->endstep = update->laststep;

    if (preflag || update->first_update == 0) {
        cac->init();
        update->integrate->setup();
    }

    timer->init();
    timer->barrier_start();
    update->integrate->run(nsteps);
    timer->barrier_stop();
    update->integrate->cleanup();
    Finish finish(cac);
    finish.end(postflag);  
  }

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  if (commands) {
    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
  }
}
