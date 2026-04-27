#include <mpi.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair.h"
#include "atom.h"
#include "element.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "math_const.h"

using namespace CAC_NS;
using namespace MathConst;

//enum{NONE, RLINEAR, RSQ, BMP};

// allocate space for static class instance variable and initialize it
int Pair::instance_total = 0;

/*  ----------------------------------------------------------------------  */

Pair::Pair(CAC *cac) : Pointers(cac)
{
  instance_me = instance_total++;
  
  eng_vdwl = eng_coul = 0.0;
  nodal_energy_weight = -1;

  comm_atom_forward = comm_elem_forward = comm_atom_reverse = comm_elem_reverse = comm_atom_reverse_off = comm_elem_reverse_off = 0;

  single_enable = 1;
  restartinfo = 1;
  respa_enable = 0;
  one_coeff = 0;
  no_virial_fdotr_compute = 1;
  writedata = 0;
  ghostneigh = 0;

  nextra = 0;
  pvector = nullptr;
  single_extra = 0;
  svector = nullptr;

  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = dipoleflag = 0;
  reinitflag = 1;
  centroidstressflag = CENTROID_SAME;

  // pair_modify settingsx
  
  compute_flag = 1;
  manybody_flag = 0;
  threebody_flag = 0;
  offset_flag = 0;
  mix_flag = GEOMETRIC;
  tail_flag = 0;
  etail = ptail = etail_ij = ptail_ij = 0.0;
  ncoultablebits = 12;
  ndisptablebits = 12;
  tabinner = sqrt(2.0);
  tabinner_disp = sqrt(2.0);

  allocated = 0;
//  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = maxcvatom = 0;
  maxeelem = maxvelem = maxcvelem = 0;
  eatom = nullptr;
  vatom = nullptr;
  cvatom = nullptr;
  enode = nullptr;
  vnode = nullptr;
  cvnode = nullptr;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

//  execution_space = Host;
//  datamask_read = ALL_MASK;
//  datamask_modify = ALL_MASK;

  copymode = 0;
  dvalue = 0;

  // debug variables

  debug_mode = 0;
}

Pair::~Pair()
{
  if (copymode) return;
  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(cvatom);
  memory->destroy(enode);
  memory->destroy(vnode);
  memory->destroy(cvnode);
}

/*  ----------------------------------------------------------------------  */

void Pair::init()
{
  int i, j;

  if (offset_flag && tail_flag)
    error->all(FLERR, "Cannot have both pair_modify shift and tail set to yes");
  if (tail_flag && domain->dimension == 2)
    error->all(FLERR, "Cannot use pair tail corrections with 2d simulations");
  if (tail_flag && domain->nonperiodic && comm->me == 0)
    error->warning(FLERR, "Using pair tail corrections with nonperiodic system");
  if (!compute_flag && tail_flag)
    error->warning(FLERR, "Using pair tail corrections with compute set to no");
  if (!compute_flag && offset_flag)
    error->warning(FLERR, "Using pair potential shift with compute set to no");

  // for manybody potentials
  // check if bonded exclusions could invalidate the neighbor list

  // I, I coeffs must be set\
  // init_one() will check if I, J is set explicitly or inferred by mixing

  if (!allocated) error->all(FLERR, "All pair coeffs are not set");

  for (i = 1; i <= atom->ntypes; i++)
    if (setflag[i][i] == 0) error->all(FLERR, "All pair coeffs are not set");

  // style-specific initialization

  init_style();
 
  // call init_one() for each I, J
  // set cutsq for each I, J, used to neighbor
  // cutforce = max of all I, J cutoffs

  cutforce = 0.0;
  etail = ptail = 0.0;
  double cut;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      cut = init_one(i, j);
      cutsq[i][j] = cutsq[j][i] = cut * cut;
      cutforce = MAX(cutforce, cut);
      if (tail_flag) {
        etail += etail_ij;
        ptail += ptail_ij;
        if (i != j) {
          etail += etail_ij;
          ptail += ptail_ij;
        }
      }
  }


}

/*  ----------------------------------------------------------------------
   init specific to a pair style
   specific pair style can override this function
     if needs its own error checks
     if needs another kind of neighbor list
   request default neighbor list = half list
-------------------------------------------------------------------------  */

void Pair::init_style()
{
  neighbor->request(this, instance_me);
}

/*  ----------------------------------------------------------------------
    mixing of pair potential prefactors (epsilon)
-------------------------------------------------------------------------  */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
  if (mix_flag == GEOMETRIC)
    return sqrt(eps1 * eps2);
  else if (mix_flag == ARITHMETIC)
    return sqrt(eps1 * eps2);
  else if (mix_flag == SIXTHPOWER)
    return (2.0 * sqrt(eps1 * eps2) *
      pow(sig1, 3.0) * pow(sig2, 3.0) / (pow(sig1, 6.0) + pow(sig2, 6.0)));
  else return 0.0;
}

/*  ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
-------------------------------------------------------------------------  */

double Pair::mix_distance(double sig1, double sig2)
{
  if (mix_flag == GEOMETRIC)
    return sqrt(sig1 * sig2);
  else if (mix_flag == ARITHMETIC)
    return (0.5 * (sig1+sig2));
  else if (mix_flag == SIXTHPOWER)
    return pow((0.5 * (pow(sig1, 6.0) + pow(sig2, 6.0))), 1.0/6.0);
  else return 0.0;
}

/*  ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   specific pair style can override this function
-------------------------------------------------------------------------  */

void Pair::init_list(int which, NeighList *ptr)
{
  list = ptr;
}

/*  ----------------------------------------------------------------------  */

void Pair::compute_dummy(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/*  ----------------------------------------------------------------------  */

double Pair::memory_usage()
{
  int maxnpe = element->maxnpe;
  int n = comm->nthreads;
  double bytes = n * maxeatom * sizeof(double);
  bytes += n * maxvatom * 6 * sizeof(double);
  bytes += n * maxeelem * maxnpe * sizeof(double);
  bytes += n * maxvelem * maxnpe * 6 * sizeof(double);
  return bytes;
}

/*  ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for bitwise settings of eflag/vflag
   set the following flags, values are otherwise set to 0:
     eflag_global != 0 if ENERGY_GLOBAL bit of eflag set
     eflag_atom   != 0 if ENERGY_ATOM bit of eflag set
     eflag_either != 0 if eflag_global or eflag_atom is set
     vflag_global != 0 if VIRIAL_PAIR bit of vflag set, OR
                       if VIRIAL_FDOTR bit of vflag is set but no_virial_fdotr = 1
     vflag_fdotr  != 0 if VIRIAL_FDOTR bit of vflag set and no_virial_fdotr = 0
     vflag_atom   != 0 if VIRIAL_ATOM bit of vflag set, OR
                       if VIRIAL_CENTROID bit of vflag set
                       and centroidstressflag != CENTROID_AVAIL
     cvflag_atom  != 0 if VIRIAL_CENTROID bit of vflag set
                       and centroidstressflag = CENTROID_AVAIL
     vflag_either != 0 if any of vflag_global, vflag_atom, cvflag_atom is set
     evflag       != 0 if eflag_either or vflag_either is set
   centroidstressflag is set by the pair style to one of these values:
     CENTROID_SAME = same as two-body stress
     CENTROID_AVAIL = different and implemented
     CENTROID_NOTAVAIL = different but not yet implemented
-------------------------------------------------------------------------  */

void Pair::ev_setup(int eflag, int vflag, int alloc)
{
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;
  int i, j, k, n;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_global = vflag & VIRIAL_PAIR;
  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 1) vflag_global = 1;
  vflag_fdotr = 0;
  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 0) vflag_fdotr = 1;
  vflag_atom = vflag & VIRIAL_ATOM;
  if (vflag & VIRIAL_CENTROID && centroidstressflag != CENTROID_AVAIL) vflag_atom = 1;
  cvflag_atom = 0;
  if (vflag & VIRIAL_CENTROID && centroidstressflag == CENTROID_AVAIL) cvflag_atom = 1;
  vflag_either = vflag_global || vflag_atom || cvflag_atom;

  evflag = eflag_either || vflag_either;

  // reallocate per-atom/per-node arrays if necessary

  if (eflag_atom) {
    if (atom->nmax > maxeatom) {
      maxeatom = atom->nmax;
      if (alloc) {
        memory->destroy(eatom);
        memory->create(eatom, comm->nthreads * maxeatom, "pair:eatom");
      }
    }
    if (element->nmax > maxeelem) {
      maxeelem = element->nmax;
      if (alloc) {
        memory->destroy(enode);
        memory->create(enode, comm->nthreads * maxeelem, maxapc, maxnpe, "pair:enode");
      }
    }
  }

  if (vflag_atom) {
    if (atom->nmax > maxvatom) {
      maxvatom = atom->nmax;
      if (alloc) {
        memory->destroy(vatom);
        memory->create(vatom, comm->nthreads * maxvatom, 6, "pair:vatom");
      }
    }
    if (element->nmax > maxvelem) {
      maxvelem = element->nmax;
      if (alloc) {
        memory->destroy(vnode);
        memory->create(vnode, comm->nthreads * maxvelem, maxapc, maxnpe, 6, "pair:vnode");
      }
    }
  }

  if (cvflag_atom) {
    if (atom->nmax > maxcvatom) {
      maxcvatom = atom->nmax;
      if (alloc) {
        memory->destroy(cvatom);
        memory->create(cvatom, comm->nthreads * maxcvatom, 9, "pair:cvatom");
      }
    }
    if (element->nmax > maxcvelem) {
      maxcvelem = element->nmax;
      if (alloc) {
        memory->destroy(cvnode);
        memory->create(cvnode, comm->nthreads * maxcvelem, maxapc, maxnpe, 9, "pair:cvnode");
      }
    }
  }

  // zero accumulators
  // use force->newton instead of newton_pair
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info

  if (eflag_global) eng_vdwl = eng_coul = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) 
          enode[i][j][k] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) {
          vnode[i][j][k][0] = 0.0;
          vnode[i][j][k][1] = 0.0;
          vnode[i][j][k][2] = 0.0;
          vnode[i][j][k][3] = 0.0;
          vnode[i][j][k][4] = 0.0;
          vnode[i][j][k][5] = 0.0;
        }
  }
  if (cvflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
    }
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) {
          cvnode[i][j][k][0] = 0.0;
          cvnode[i][j][k][1] = 0.0;
          cvnode[i][j][k][2] = 0.0;
          cvnode[i][j][k][3] = 0.0;
          cvnode[i][j][k][4] = 0.0;
          cvnode[i][j][k][5] = 0.0;
          cvnode[i][j][k][6] = 0.0;
          cvnode[i][j][k][7] = 0.0;
          cvnode[i][j][k][8] = 0.0;
        }
  }

  // if vflag_global = 2 and no element in simulation and pair::compute() calls virial_fdotr_compute() 
  // compute global virial via (F dot r) instead of via pairwise summation
  // unset other flags as appropriate

  if (vflag_global == 2 && no_virial_fdotr_compute == 0 && element->nelements == 0) {
    vflag_fdotr = 1;
  } else vflag_fdotr = 0;
}

/*  ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators 
   i and j can be either atoms or nodes
   check if jeptr/jvptr == nullptr
   -------------------------------------------------------------------------  */

void Pair::ev_tally(double escale, double vscale, 
    double *ieptr, double *jeptr, double *ivptr, double *jvptr, 
    double evdwl, double ecoul, double fpair, 
    double delx, double dely, double delz)
{
  double epairhalf, v[6];

  if (eflag_either) {
    if (eflag_global) {
      eng_vdwl += evdwl * escale;
      eng_coul += ecoul * escale;
    }
    if (eflag_atom) {
      epairhalf = 0.5 * (evdwl + ecoul);
      if (ieptr != nullptr) *ieptr += epairhalf;
      if (jeptr != nullptr) *jeptr += epairhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx * delx * fpair;
    v[1] = dely * dely * fpair;
    v[2] = delz * delz * fpair;
    v[3] = delx * dely * fpair;
    v[4] = delx * delz * fpair;
    v[5] = dely * delz * fpair;

    if (vflag_global) {
      virial[0] += v[0] * vscale;
      virial[1] += v[1] * vscale;
      virial[2] += v[2] * vscale;
      virial[3] += v[3] * vscale;
      virial[4] += v[4] * vscale;
      virial[5] += v[5] * vscale;
    }
    if (vflag_atom) {
      if (ivptr != nullptr) {
        ivptr[0] += 0.5 * v[0]; ivptr[1] += 0.5 * v[1];
        ivptr[2] += 0.5 * v[2]; ivptr[3] += 0.5 * v[3];
        ivptr[4] += 0.5 * v[4]; ivptr[5] += 0.5 * v[5];
      }
      if (jvptr != nullptr) {
        jvptr[0] += 0.5 * v[0]; jvptr[1] += 0.5 * v[1];
        jvptr[2] += 0.5 * v[2]; jvptr[3] += 0.5 * v[3];
        jvptr[4] += 0.5 * v[4]; jvptr[5] += 0.5 * v[5];
      }
    }
  }
}

/*  ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji * fj + drki * fk
 -------------------------------------------------------------------------  */

void Pair::ev_tally3(double escale, double vscale, double *ieptr, double *jeptr, 
    double *ivptr, double *jvptr, double *keptr, double *kvptr, 
    double evdwl, double ecoul, double *fj, double *fk, double *drji, double *drki)
{
  double epairthird, v[6];

  if (eflag_either) {
    if (eflag_global) {
      eng_vdwl += evdwl * escale;
      eng_coul += ecoul * escale;
    }
    if (eflag_atom) {
      epairthird = THIRD * (evdwl + ecoul);
      if (ieptr != nullptr) *ieptr += epairthird;
      if (jeptr != nullptr) *jeptr += epairthird;
      if (keptr != nullptr) *keptr += epairthird;
    }
  }

  if (vflag_either) {
    v[0] = drji[0] * fj[0] + drki[0] * fk[0];
    v[1] = drji[1] * fj[1] + drki[1] * fk[1];
    v[2] = drji[2] * fj[2] + drki[2] * fk[2];
    v[3] = drji[0] * fj[1] + drki[0] * fk[1];
    v[4] = drji[0] * fj[2] + drki[0] * fk[2];
    v[5] = drji[1] * fj[2] + drki[1] * fk[2];

    if (vflag_global) {
      virial[0] += v[0] * vscale;
      virial[1] += v[1] * vscale;
      virial[2] += v[2] * vscale;
      virial[3] += v[3] * vscale;
      virial[4] += v[4] * vscale;
      virial[5] += v[5] * vscale;
    }

    if (vflag_atom) {
      if (ivptr != nullptr) {
        ivptr[0] += THIRD * v[0]; ivptr[1] += THIRD * v[1];
        ivptr[2] += THIRD * v[2]; ivptr[3] += THIRD * v[3];
        ivptr[4] += THIRD * v[4]; ivptr[5] += THIRD * v[5];
      }
      if (jvptr != nullptr) {
        jvptr[0] += THIRD * v[0]; jvptr[1] += THIRD * v[1];
        jvptr[2] += THIRD * v[2]; jvptr[3] += THIRD * v[3];
        jvptr[4] += THIRD * v[4]; jvptr[5] += THIRD * v[5];
      }
      if (kvptr != nullptr) {
        kvptr[0] += THIRD * v[0]; kvptr[1] += THIRD * v[1];
        kvptr[2] += THIRD * v[2]; kvptr[3] += THIRD * v[3];
        kvptr[4] += THIRD * v[4]; kvptr[5] += THIRD * v[5];
      }
    }
  }
}


/*  ----------------------------------------------------------------------
   compute global pair virial via summing F dot r over own & ghost atoms
   at this point, only pairwise forces have been accumulated in atom->f
   only called when no element in simulation
   -------------------------------------------------------------------------  */

void Pair::virial_fdotr_compute()
{
  double **x = atom->x;
  double **f = atom->f;

  // sum over force on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++) {
    virial[0] += f[i][0] * x[i][0];
    virial[1] += f[i][1] * x[i][1];
    virial[2] += f[i][2] * x[i][2];
    virial[3] += f[i][1] * x[i][0];
    virial[4] += f[i][2] * x[i][0];
    virial[5] += f[i][2] * x[i][1];
  }

  // prevent multiple calls to update the virial
  // when a hybrid pair style uses both a gpu and non-gpu pair style
  // or when respa is used with gpu pair styles

  vflag_fdotr = 0;
}

/*  ----------------------------------------------------------------------
   set all flags to zero for energy, virial computation
   called by some complicated many-body potentials that use individual flags
   to insure no holdover of flags from previous timestep
   -------------------------------------------------------------------------  */

void Pair::ev_unset()
{
  evflag = 0;

  eflag_either = 0;
  eflag_global = 0;
  eflag_atom = 0;

  vflag_either = 0;
  vflag_global = 0;
  vflag_atom = 0;
  vflag_fdotr = 0;
}

NeighList* Pair::get_list()
{
  return list;
}
