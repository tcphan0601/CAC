#include <mpi.h>
#include <string.h>
#include <ctype.h>
#include "cac.h"
#include "style_atom.h"
#include "style_command.h"
#include "style_pair.h"
#include "style_fix.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "input.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "comm_brick.h"
#include "group.h"
#include "output.h"
#include "domain.h"
#include "force.h"
#include "timer.h"

using namespace CAC_NS;

/* ----------------------------------------------------------------------
   start up CAC
   allocate fundamental classes (memory, error, universe, input)
   parse input switches
   initialize communicators, screen & logfile output
   input is allocated at end after MPI info is setup
   ------------------------------------------------------------------------- */

CAC::CAC(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this, communicator);

  screen = NULL;
  logfile = NULL;
  infile = NULL;

  initclock = MPI_Wtime();

  // parse input switches 

  int inflag = 0;
  int screenflag = 0;
  int logflag = 0;

  nclass = 0;
  suffix = suffix2 = NULL;
  suffix_enable = 0;


  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-echo") == 0 ||
               strcmp(arg[iarg],"-e") == 0) {

      // only check that sufficient arguments exist, actual argument will be checked in input.cpp
      
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      iarg += 2;

    } else if (strcmp(arg[iarg],"-in") == 0 ||
        strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      inflag = iarg + 1;
      iarg +=2;
    } else if (strcmp(arg[iarg],"-log") == 0 ||
               strcmp(arg[iarg],"-l") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      logflag = iarg + 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"-screen") == 0 ||
               strcmp(arg[iarg],"-sc") == 0) {
      if (iarg+2 > narg)
        error->universe_all(FLERR,"Invalid command-line argument");
      screenflag = iarg + 1;
      iarg += 2;
    }

  }

  // set universe screen and logfile

  if (universe->me == 0) {
    if (screenflag == 0)
      universe->uscreen = stdout;
    else if (strcmp(arg[screenflag],"none") == 0)
      universe->uscreen = NULL;
    else {
      universe->uscreen = fopen(arg[screenflag],"w");
      if (universe->uscreen == NULL)
        error->universe_one(FLERR,"Cannot open universe screen file");
    }
    if (logflag == 0) {
      //if (helpflag == 0) {
      universe->ulogfile = fopen("log.cac","w");
      if (universe->ulogfile == NULL)
        error->universe_warn(FLERR,"Cannot open log.cac for writing");
      //}
    } else if (strcmp(arg[logflag],"none") == 0)
      universe->ulogfile = NULL;
    else {
      universe->ulogfile = fopen(arg[logflag],"w");
      if (universe->ulogfile == NULL)
        error->universe_one(FLERR,"Cannot open universe log file");
    }
  }

  if (universe->me > 0) {
    if (screenflag == 0) universe->uscreen = stdout;
    else universe->uscreen = NULL;
    universe->ulogfile = NULL;
  }

  // make universe and single world the same, since no partition switch
  // world inherits settings from universe
  // set world screen, logfile, communicator, infile
  // open input script if from file
 
  screen = universe->uscreen;
  logfile = universe->ulogfile;
  world = universe->uworld;

  if (universe->me == 0) {
    infile = fopen(arg[inflag],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s", arg[inflag]);
      error->one(FLERR,str);
    }

    if (screen) fprintf(screen,"Simulation running on CAC (%s)\n",universe->version);
    if (logfile) fprintf(logfile,"Simulation running on CAC (%s)\n",universe->version);
  }

  // check consistency of datatype settings in cactype.h
  
  if (sizeof(smallint) != sizeof(int))
    error->all(FLERR,"Smallint setting in cactype.h is invalid");
  if (sizeof(imageint) < sizeof(smallint))
    error->all(FLERR,"Imageint setting in cactype.h is invalid");
  if (sizeof(tagint) < sizeof(smallint))
    error->all(FLERR,"Tagint setting in cactype.h is invalid");
  if (sizeof(bigint) < sizeof(imageint) || sizeof(bigint) < sizeof(tagint))
    error->all(FLERR,"Bigint setting in cactype.h is invalid");

  int mpisize;
  MPI_Type_size(MPI_CAC_TAGINT,&mpisize);
  if (mpisize != sizeof(tagint))
    error->all(FLERR,"MPI_CAC_TAGINT and tagint in "
        "cactype.h are not compatible");
  MPI_Type_size(MPI_CAC_BIGINT,&mpisize);
  if (mpisize != sizeof(bigint))
    error->all(FLERR,"MPI_CAC_BIGINT and bigint in "
        "cactype.h are not compatible");

#ifdef CAC_SMALLBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || 
      sizeof(tagint) != 4 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef CAC_BIGBIG
  if (sizeof(smallint) != 4 || sizeof(imageint) != 8 || 
      sizeof(tagint) != 8 || sizeof(bigint) != 8)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif
#ifdef CAC_SMALLSMALL
  if (sizeof(smallint) != 4 || sizeof(imageint) != 4 || 
      sizeof(tagint) != 4 || sizeof(bigint) != 4)
    error->all(FLERR,"Small to big integers are not sized correctly");
#endif

  input = new Input(this,narg,arg);

  create();

}

/* ----------------------------------------------------------------------
   shutdown CAC
   delete top-level classes
   close screen and log files in world and universe
   output files were already closed in destroy()
   delete fundamental classes
------------------------------------------------------------------------- */

CAC::~CAC()
{
  const int me = comm->me;

  destroy();

  double totalclock = MPI_Wtime() - initclock;
  if ((me == 0) && (screen || logfile)) {
    char outtime[128];
    int seconds = fmod(totalclock,60.0);
    totalclock  = (totalclock - seconds) / 60.0;
    int minutes = fmod(totalclock,60.0);
    int hours = (totalclock - minutes) / 60.0;
    sprintf(outtime,"Total wall time: "
            "%d:%02d:%02d\n", hours, minutes, seconds);
    if (screen) fputs(outtime,screen);
    if (logfile) fputs(outtime,logfile);
  }

  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);
  logfile = NULL;
  if (screen != stdout) screen = NULL;
  if (infile && infile != stdin) fclose(infile);
  if (world != universe->uworld) MPI_Comm_free(&world);

  delete input;
  delete universe;
  delete error;
  delete memory;
}

void CAC::create()
{
  comm = new CommBrick(this);
  neighbor = new Neighbor(this);
  domain = new Domain(this);
  atom = new Atom(this);
  atom->create_avec("atomic",0,NULL);
  element = new Element(this); 
  
  char **args = new char*[2];
  args[0] = (char *) "8";
  args[1] = (char *) "1";
  element->create_evec("cac",2,args);
  delete args;
  group = new Group(this);
  force = new Force(this);
  modify = new Modify(this);
  output = new Output(this);
  update = new Update(this);
  timer = new Timer(this); 
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
   ------------------------------------------------------------------------- */

void CAC::destroy()
{

  delete update;
  update = NULL;

  delete neighbor;
  neighbor = NULL;

  delete comm;
  comm = NULL;

  delete force;
  force = NULL;

  delete group;
  group = NULL;

  delete output;
  output = NULL;

  delete modify;          // modify must come after output, force, update
                          //   since they delete fixes
  modify = NULL;

  delete domain;          // domain must come after modify
                          //   since fix destructors access domain
  domain = NULL;

  delete atom;        
  atom = NULL;

  delete element;        
  element = NULL;

  delete timer;
  timer = NULL;
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
   ------------------------------------------------------------------------- */

void CAC::init()
{
  update->init();
  force->init();
  domain->init();
  if (atom->natoms) atom->init();
  if (element->nelements) element->init(); // must before neighbor for element size
  modify->init();
  neighbor->init();
  comm->init();
  output->init();
}
