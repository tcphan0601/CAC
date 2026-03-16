#include <string.h>
#include <stdlib.h>
#include "fix_balance.h"
#include "balance.h"
#include "update.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "irregular_comm.h"
#include "universe.h"
#include "modify.h"
#include "rcb.h"
#include "timer.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

enum{SHIFT,BISECTION};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), balance(NULL), irregular_comm(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix balance command");

  box_change_domain = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  screenflag = 0;
  // parse required arguments

  int dimension = domain->dimension;

  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix balance command");
  thresh = universe->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"shift") == 0) lbstyle = SHIFT;
  else if (strcmp(arg[5],"rcb") == 0) lbstyle = BISECTION;
  else error->all(FLERR,"Illegal fix balance command");

  int iarg = 5;
  if (lbstyle == SHIFT) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[iarg+1]) > 3)
      error->all(FLERR,"Illegal fix balance command");
    strcpy(bstr,arg[iarg+1]);
    nitermax = universe->inumeric(FLERR,arg[iarg+2]);
    if (nitermax <= 0) error->all(FLERR,"Illegal fix balance command");
    stopthresh = universe->numeric(FLERR,arg[iarg+3]);
    if (stopthresh < 1.0) error->all(FLERR,"Illegal fix balance command");
    iarg += 4;
  } else if (lbstyle == BISECTION) {
    iarg++;
  }

  // error checks

  if (lbstyle == SHIFT) {
    int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Fix balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Fix balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Fix balance shift string is invalid");
    }
  }

  if (lbstyle == BISECTION && comm->style == 0)
    error->all(FLERR,"Fix balance rcb cannot be used with comm_style brick");

  // create instance of Balance class
  // if SHIFT, initialize it with params
  // process remaining optional args via Balance

  balance = new Balance(cac);
  if (lbstyle == SHIFT) balance->shift_setup(bstr,nitermax,thresh);
  screenflag = balance->options(iarg,narg,arg);

  // create instance of IrregularComm class

  irregular_comm = new IrregularComm(cac);

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;
  lastbalance = -1;
  
  // compute initial outputs

  itercount = 0;
  pending = 0;
  imbfinal = imbprev = maxloadperproc = minloadperproc = 0.0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete balance;
  delete irregular_comm;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::post_constructor()
{
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{

}

/* ---------------------------------------------------------------------- */

void FixBalance::setup(int vflag)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // has to be be done before rebalance() invokes IrregularComm::migrate()
  //   since it requires atoms be inside simulation box
  //   even though pbc() will be done again in Verlet::run()
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->x2lamda(element->nlocal,element->npe,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2x(element->nlocal,element->npe,element->nodex);
  }

  // perform a rebalance if threshold exceeded

  imbnow = balance->imbalance_factor(maxloadperproc,minloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
  }
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
  }


  // perform a rebalance if threshold exceeded

  imbnow = balance->imbalance_factor(maxloadperproc,minloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occured
   ------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_factor(maxloadperproc,minloadperproc);
  pending = 0;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
   ------------------------------------------------------------------------- */

void FixBalance::rebalance()
{

  if (screenflag && comm->me == 0) fprintf(screen,"Rebalancing...\n  timestep = " BIGINT_FORMAT "\n",update->ntimestep);

  imbprev = imbnow;

  // invoke balancer and reset comm->uniform flag

  int *sendproc;
  if (lbstyle == SHIFT) {
    itercount = balance->shift();
    comm->layout = LAYOUT_NONUNIFORM;
  } else if (lbstyle == BISECTION) {
    sendproc = balance->bisection();
    comm->layout = LAYOUT_TILED;
  }

  // output of new decomposition

  if (balance->outflag) balance->dumpout(update->ntimestep);

  // reset proc sub-domains
  // check and warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  domain->subbox_too_small_check(neighbor->skin);

  // move atoms to new processors via irregular_comm()
  // only needed if migrate_check() says an atom moves to far
  // else allow caller's comm->exchange() to do it

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
  }
  if (lbstyle == BISECTION) irregular_comm->migrate(0,1,sendproc);
  else if (irregular_comm->migrate_check()) irregular_comm->migrate();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
  }

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in comm->exchange()

  pending = 1;

  double imbnew = balance->imbalance_factor(maxloadperproc,minloadperproc);
  if (screenflag && comm->me == 0) fprintf(screen,"  old/new imbalance factor = %g/%g\n",imbprev,imbnew);
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
   ------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
   ------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return maxloadperproc;
  if (i == 1) return minloadperproc;
  if (i == 2) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = irregular_comm->memory_usage();
  if (balance->rcb) bytes += balance->rcb->memory_usage();
  return bytes;
}
