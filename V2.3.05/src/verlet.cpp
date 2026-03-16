#include <string.h>
#include "verlet.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "element_vec.h"
#include "force.h"
#include "pair.h"
#include "output.h"
#include "dump.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "balance.h"

using namespace CAC_NS;

Verlet::Verlet(CAC *cac, int narg, char **arg) : 
  Integrate(cac, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
   ------------------------------------------------------------------------- */

void Verlet::init()
{

  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_fdotr_compute via sum over ghosts

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;

}

/* ----------------------------------------------------------------------
   setup before run
   ------------------------------------------------------------------------- */

void Verlet::setup(int flag)
{

  FILE *out;
  if (comm->me == 0) {
    for (int m = 0; m < 2; m++) {
      if (m == 0) out = screen;
      else out = logfile;

      fprintf(out,"Setting up Verlet run ...\n");
      if (flag) {
        fprintf(out,"  Unit style       : %s\n", update->unit_style);
        fprintf(out,"  Current step     : " BIGINT_FORMAT "\n", update->ntimestep);
        fprintf(out,"  Time step        : %g\n", update->dt);
      }
    }
  }

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  element->setup();
  modify->setup_pre_exchange();

  // print out imbalance info


  Balance *balance = new Balance(cac);
  double maxcost,mincost;
  double imbalance = balance->imbalance_factor(maxcost,mincost);

  if (comm->me == 0) 
    for (int m = 0; m < 2; m++) {
      if (m == 0) out = screen;
      else out = logfile;
      if (flag) {
        fprintf(out,"  Imbalance factor : %g\n",imbalance);
        fprintf(out,"  Max/min load/proc: %g/%g\n",maxcost,mincost);
        timer->print_timeout(out);
      }
    }
  
  delete balance;

  if (triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->x2lamda(element->nlocal,element->npe,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins(); 
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (triclinic) {
    domain->lamda2x(atom->nlocal + atom->nghost,atom->x);
    domain->lamda2x(element->nlocal + element->nghost,element->x);
    domain->lamda2x(element->nlocal + element->nghost,element->npe,element->nodex);
  }

  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  force->setup();
  ev_set(update->ntimestep);
  force_clear(); 
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  output->setup(flag);
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms/element nodes
   clear other arrays as needed
   ------------------------------------------------------------------------- */

void Verlet::force_clear()
{
  int i;
  size_t nbytes;

  if (external_force_clear) return;

  // clear force on all particles/element nodes
  // when using threads always clear all forces.

  nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;
  if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

  nbytes = sizeof(double) * element->nlocal;
  if (force->newton) nbytes += sizeof(double) * element->nghost;
  if (nbytes) memset(&element->nodef[0][0][0],0,3*element->npe*nbytes);
}


/* ----------------------------------------------------------------------
   run for N steps
   ------------------------------------------------------------------------- */

void Verlet::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;


  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {
    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // check and update element size 

    if (element->nelements) element->evec->check_element_size();

    // initial time integration

    timer->stamp();
    modify->initial_integrate(vflag);    
    if (n_post_integrate) modify->post_integrate();
    timer->stamp(Timer::MODIFY);

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(Timer::COMM);
    } else {
      if (n_pre_exchange) {
        timer->stamp();
        modify->pre_exchange();
        timer->stamp(Timer::MODIFY);
      }

      if (triclinic) {
        domain->x2lamda(atom->nlocal,atom->x);
        domain->x2lamda(element->nlocal,element->x);
        domain->x2lamda(element->nlocal,element->npe,element->nodex);
      }
      domain->pbc();

      // update cutoff based on max element size

      neighbor->update_cut();

      if (domain->box_change) {
        domain->reset_box();
        comm->setup();
        if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();
      comm->exchange();
      if (sortflag && ntimestep >= atom->nextsort) atom->sort();
      comm->borders();
      if (triclinic) {
        domain->lamda2x(atom->nlocal+atom->nghost,atom->x);
        domain->lamda2x(element->nlocal+element->nghost,element->x);
        domain->lamda2x(element->nlocal+element->nghost,element->npe,element->nodex);
      }
      timer->stamp(Timer::COMM);
      if (n_pre_neighbor) {
        modify->pre_neighbor();
        timer->stamp(Timer::MODIFY);
      }
      neighbor->build();
      timer->stamp(Timer::NEIGH);
    }

    // force computations

    force_clear();
    timer->stamp();

    if (n_pre_force) {
      modify->pre_force(vflag);
      timer->stamp(Timer::MODIFY);
    }

    if (pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }

    if (n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);
    }

    // reverse communication of forces

    if (force->newton) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
    }

    if (n_post_force) modify->post_force(vflag);
    modify->final_integrate();

    if (n_end_of_step) modify->end_of_step();
    timer->stamp(Timer::MODIFY);

    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Verlet::cleanup()
{
  modify->post_run();
  domain->box_too_small_check();
  update->update_time();
}
