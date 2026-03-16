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
#include "fix_balance.h"
#include "input.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "balance.h"

using namespace CAC_NS;

Verlet::Verlet(CAC *cac, int narg, char **arg) : 
  Integrate(cac, narg, arg) {}

/*  ----------------------------------------------------------------------
   initialization before run
   -------------------------------------------------------------------------  */

void Verlet::init()
{

  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR, "No fixes defined, atoms won't move");

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

/*  ----------------------------------------------------------------------
   setup before run
   -------------------------------------------------------------------------  */

void Verlet::setup(int forceflag, int flag)
{

  if (element->nelements + atom->natoms < comm->nprocs && comm->style) 
    error->all(FLERR, "Total number of procs can't be more than total number of (atoms + elements)");
  FILE *out;
  if (comm->me == 0) {
    for (int m = 0; m < 2; m++) {
      if (m == 0) out = screen;
      else out = logfile;

      if (out) {
        fprintf(out, "Setting up Verlet run ...\n");
        if (flag) {
          fprintf(out, "  Unit style       : %s\n", update->unit_style);
          fprintf(out, "  Current step     : " BIGINT_FORMAT "\n", update->ntimestep);
          fprintf(out, "  Time step        : %g\n", update->dt);
        }
      }
    }
  }

  MPI_Barrier(world);

  update->setupflag = 1;
  neighbor->lastcall = -1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  element->setup();
  modify->setup_pre_exchange();

  if (triclinic) {
    domain->x2lamda(atom->nlocal, atom->x);
    domain->x2lamda(element->nlocal, element->x);
    domain->nodex2lamda(element->nlocal, element->nodex);
  }

  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();

  comm->exchange();

  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();

  if (triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost, atom->x);
    domain->lamda2x(element->nlocal+element->nghost, element->x);
    domain->lamda2nodex(element->nlocal+element->nghost, element->nodex);
  }

  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();

  // neigh list first build might occur already
  // in setup_pre_exchange from fix_adaptive

  if (forceflag) 
    if (neighbor->lastcall < 0) neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  force->setup();
  ev_set(update->ntimestep);
  force_clear(); 

  modify->setup_pre_force(vflag);

  if (forceflag) {
    if (pair_compute_flag) force->pair->compute(eflag, vflag);
    else if (force->pair) force->pair->compute_dummy(eflag, vflag);
  }

  modify->setup_pre_reverse(eflag, vflag);

  if (forceflag) {
    if (force->newton) comm->reverse_comm();
    if (element->nelements) element->compute_nodef();
  }
  modify->setup(vflag);

  // print out imbalance info

  double imbalance;
  double maxload, minload;
  Fix *fix_balance = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "balance") == 0)
      fix_balance = modify->fix[i];
  if (fix_balance) {
    imbalance = fix_balance->compute_scalar();
    maxload = fix_balance->compute_vector(0);
    minload = fix_balance->compute_vector(1);
  } else {
    Balance *balance = new Balance(cac);
    imbalance = balance->imbalance_factor(maxload, minload);
    delete balance;
  }
  if (comm->me == 0) 
    for (int m = 0; m < 2; m++) {
      if (m == 0) out = screen;
      else out = logfile;
      if (out && flag) {
        fprintf(out, "  Imbalance factor : %g\n", imbalance);
        fprintf(out, "  Max/min load/proc: %g/%g\n", maxload, minload);
        timer->print_timeout(out);
      }
    }
  output->setup(flag);
  update->setupflag = 0;
}

/*  ----------------------------------------------------------------------
   clear force on own & ghost atoms/element nodes
   clear other arrays as needed
   -------------------------------------------------------------------------  */

void Verlet::force_clear()
{
  int i;
  size_t nbytes;

  if (external_force_clear) return;

  // clear force on all particles/element nodes
  // when using threads always clear all forces.

  nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;
  if (nbytes) memset(&atom->f[0][0], 0, 3 * nbytes);

  nbytes = sizeof(double) * element->nlocal;
  if (force->newton) nbytes += sizeof(double) * element->nghost;
  if (nbytes) {
    memset(&element->nodef[0][0][0][0], 0, 3 * element->maxapc * element->maxnpe * nbytes);
    memset(&element->gaussf[0][0][0][0], 0, 3 * element->max_ngcell * element->maxapc * nbytes);
  }
}

/*  ----------------------------------------------------------------------
   run for N steps
   -------------------------------------------------------------------------  */

void Verlet::run(int n, int forceflag)
{
  bigint ntimestep;
  int nflag, sortflag;

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

    // initial time integration

    timer->stamp();
    modify->initial_integrate(vflag);    
    if (n_post_integrate) modify->post_integrate();
    timer->stamp(Timer::MODIFY);

    // regular communication vs neighbor list rebuild

    if (forceflag) {
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

        timer->stamp();
        if (triclinic) {
          domain->x2lamda(atom->nlocal, atom->x);
          domain->x2lamda(element->nlocal, element->x);
          domain->nodex2lamda(element->nlocal, element->nodex);
        }
        domain->pbc();

        // need to exchange first before setup borders cutghost

        if (domain->box_change) {
          comm->setup_exchange();
          domain->reset_box();
        }

        comm->exchange();

        // if there are elements, force comm->setup_borders() and neighbor->setup_bins() to 
        // update cutghost from element->local_element_bound_box

        if (element->nelements || domain->box_change) {
          comm->setup_borders();
          if (neighbor->style) neighbor->setup_bins();
        }

        if (sortflag && ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) {
          domain->lamda2x(atom->nlocal+atom->nghost, atom->x);
          domain->lamda2x(element->nlocal+element->nghost, element->x);
          domain->lamda2nodex(element->nlocal+element->nghost, element->nodex);
        }
        timer->stamp(Timer::COMM);
        if (forceflag) {
          if (n_pre_neighbor) {
            modify->pre_neighbor();
            timer->stamp(Timer::MODIFY);
          }
          neighbor->build();
        }
        timer->stamp(Timer::NEIGH);
      }

      // force computations

      force_clear();
    }
    timer->stamp();

    if (n_pre_force) {
      modify->pre_force(vflag);
      timer->stamp(Timer::MODIFY);
    }

    if (forceflag)
      if (pair_compute_flag) {
        force->pair->compute(eflag, vflag);
        timer->stamp(Timer::PAIR);
      }

    if (n_pre_reverse) {
      modify->pre_reverse(eflag, vflag);
      timer->stamp(Timer::MODIFY);
    }

    // reverse communication of forces

    if (forceflag)
    if (force->newton) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
      if (element->nelements) {
        element->compute_nodef();
        timer->stamp(Timer::PAIR);
      }
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

    // reset skip_exchange 
    if (comm->skip_exchange) comm->skip_exchange = 0;
  }
}

/*  ----------------------------------------------------------------------  */

void Verlet::cleanup()
{
  modify->post_run();
  domain->box_too_small_check();
  update->update_time();

  // reset neighbor->lastcall

  neighbor->lastcall = -1;
}
