#include "atom.h"
#include "comm.h"
#include "element.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair_lj_mdf.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathConst;

PairLJMDF::PairLJMDF(CAC *cac) : Pair(cac) {
  // rRESPA not developed yet
  // respa_enable = 1;
  writedata = 1;

  nemax = 0;
}

PairLJMDF::~PairLJMDF() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(aparm);
    memory->destroy(bparm);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairLJMDF::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut_inner, n + 1, n + 1, "pair:cut_inner");
  memory->create(cut_inner_sq, n + 1, n + 1, "pair:cut_inner_sq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(aparm, n + 1, n + 1, "pair:aparm");
  memory->create(bparm, n + 1, n + 1, "pair:bparm");
  memory->create(lj1, n + 1, n + 1, "pair:lj1");
  memory->create(lj2, n + 1, n + 1, "pair:lj2");
  memory->create(lj3, n + 1, n + 1, "pair:lj3");
  memory->create(lj4, n + 1, n + 1, "pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairLJMDF::settings(int narg, char **arg) {
  if (narg != 2)
    error->all(FLERR, "Illegal pair_style command");

  cut_inner_global = universe->numeric(FLERR, arg[0]);
  cut_global = universe->numeric(FLERR, arg[1]);

  // reset cutoffs that have been explicitly set
  //
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i + 1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairLJMDF::coeff(int narg, char **arg) {

  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double aparm_one = universe->numeric(FLERR, arg[2]);
  double bparm_one = universe->numeric(FLERR, arg[3]);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 6) {
    cut_inner_one = universe->numeric(FLERR, arg[4]);
    cut_one = universe->numeric(FLERR, arg[5]);
  }
  if (cut_inner_one <= 0.0 || cut_inner_one > cut_one)
    error->all(FLERR, "Illegal pair_coeff command");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      aparm[i][j] = aparm_one;
      bparm[i][j] = bparm_one;
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
   ------------------------------------------------------------------------- */

void PairLJMDF::write_data_all(FILE *fp) {
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, aparm[i][j], bparm[i][j],
              cut[i][j]);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairLJMDF::init_one(int i, int j) {
  if (setflag[i][j] == 0)
    error->all(FLERR, "All pair coeffs are not set");

  cut_inner_sq[i][j] = cut_inner[i][j] * cut_inner[i][j];

  lj1[i][j] = 12.0 * aparm[i][j];
  lj2[i][j] = 6.0 * bparm[i][j];
  lj3[i][j] = aparm[i][j];
  lj4[i][j] = bparm[i][j];

  cut[j][i] = cut[i][j];
  cut_inner[j][i] = cut_inner[i][j];
  cut_inner_sq[j][i] = cut_inner_sq[i][j];

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairLJMDF::compute(int eflag, int vflag) {
  int jnum;
  double r2inv, r6inv, forcelj, fpair;

  int i, j, ii, jj;
  int ictype, jctype, ietype,
      jetype; // note: ictype and jctype are atomic (chemical) type, same as
              // itype and jtype in LAMMPS

  int node, inode, jintpl, iintpl, iintg, iintg_local;

  double evdwl = 0.0;

  int *jlist, *jindexlist;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  double fx, fy, fz;

  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;

  int newton_pair = force->newton_pair;
  // double *special_lj = force->special_lj;

  int *npe = element->npe;

  double ***nodex = element->nodex;

  int *nintg = element->nintg;

  int *etype = element->etype;
  int *ctype = element->ctype;

  int **i2ia = element->i2ia;
  int **i2n = element->i2n;
  double ***shape_array = element->shape_array;
  double ***weighted_shape_array = element->weighted_shape_array;

  int ainum = list->ainum;
  int einum = list->einum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;

  int *e2ilist = list->e2ilist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighi2a = list->numneighi2a;
  int *numneighi2ia = list->numneighi2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighi2a = list->firstneighi2a;
  int **firstneighi2ia = list->firstneighi2ia;
  int **firstneighi2ia_index = list->firstneighi2ia_index;

  // compute forces on each ATOM
  // loop over neighbors of my atoms

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];
    ictype = atype[i];

    // loop over neighboring atoms, skip half if newton is on or j is owned

    jnum = numneigha2a[i];
    if (jnum) {
      jlist = firstneigha2a[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jctype = atype[j];

        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (rsq < cutsq[ictype][jctype]) {
          r2inv = 1.0 / rsq;
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
          fpair = forcelj * r2inv;
          // fpair = factor_lj*forcelj*r2inv;

          af[i][0] += delx * fpair;
          af[i][1] += dely * fpair;
          af[i][2] += delz * fpair;

          if (newton_pair || j < nalocal) {
            af[j][0] -= delx * fpair;
            af[j][1] -= dely * fpair;
            af[j][2] -= delz * fpair;
          }

          if (eflag)
            evdwl = r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);

          if (evflag)
            atom_ev_tally(i, j, nalocal, newton_pair, evdwl, 0.0, fpair, delx,
                          dely, delz);
        }
      }
    }

    // loop over neighboring interpolated atoms

    jnum = numneigha2ia[i];
    if (jnum) {
      jlist = firstneigha2ia[i];
      jindexlist = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jetype = etype[j];
        jintpl = jindexlist[jj];

        // interpolate positions at site J

        delx = xtmp;
        dely = ytmp;
        delz = ztmp;
        for (node = 0; node < npe[jetype]; node++) {
          delx -= shape_array[jetype][jintpl][node] * nodex[j][node][0];
          dely -= shape_array[jetype][jintpl][node] * nodex[j][node][1];
          delz -= shape_array[jetype][jintpl][node] * nodex[j][node][2];
        }
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq[ictype][jctype]) {
          r2inv = 1.0 / rsq;
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
          fpair = forcelj * r2inv;
          // fpair = factor_lj*forcelj*r2inv;

          af[i][0] += delx * fpair;
          af[i][1] += dely * fpair;
          af[i][2] += delz * fpair;

          if (eflag)
            evdwl = r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);

          if (evflag)
            atom_ev_tally_full(i, evdwl, 0.0, fpair, delx, dely, delz);
        }
      }
    }
  }

  // compute forces on each NODE
  // compute forces on each integration points and distribute to nodes
  // loop over my integration points iing in element i
  // fx fy fz = force at current integration point.
  // use e2ilist to get the local index of the first integration point
  //   of element i
  // 1. Loop over neighboring atoms. Add force contribution to fx fy fz.
  // 2. Loop over neighboring interpolated atoms. Add force contribution to fx
  // fy fz
  // 3. After getting total forces at all integration points in element i,
  //    distribute forces at integration points to nodef in element i.

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    ictype = ctype[i];

    reset_force_columns();

    for (iintg = 0; iintg < nintg[ietype]; iintg++) {
      fx = fy = fz = 0.0;

      // interpolate postions of integration point from node values

      iintpl = i2ia[ietype][iintg];
      inode = i2n[ietype][iintg];
      xtmp = ytmp = ztmp = 0.0;
      for (node = 0; node < npe[ietype]; node++) {
        xtmp += shape_array[ietype][iintpl][node] * nodex[i][node][0];
        ytmp += shape_array[ietype][iintpl][node] * nodex[i][node][1];
        ztmp += shape_array[ietype][iintpl][node] * nodex[i][node][2];
      }

      // loop over neighboring atoms and add force contribution

      iintg_local = e2ilist[ii] + iintg;
      jnum = numneighi2a[iintg_local];
      if (jnum) {
        jlist = firstneighi2a[iintg_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          jctype = atype[j];
          delx = xtmp - ax[j][0];
          dely = ytmp - ax[j][1];
          delz = ztmp - ax[j][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq[ictype][jctype]) {
            r2inv = 1.0 / rsq;
            r6inv = r2inv * r2inv * r2inv;
            forcelj =
                r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
            fpair = forcelj * r2inv;

            fx += delx * fpair;
            fy += dely * fpair;
            fz += delz * fpair;

            // tally eng_dwl and virial if integration point is a node

            if (eflag)
              evdwl =
                  r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);
            if (evflag)
              node_ev_tally(i, inode, evdwl, 0.0, fpair, delx, dely, delz);
          }
        }
      }

      // loop over neighboring interpolated atoms and accumulate force
      // contribution

      jnum = numneighi2ia[iintg_local];
      if (jnum) {
        jlist = firstneighi2ia[iintg_local];
        jindexlist = firstneighi2ia_index[iintg_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          jintpl = jindexlist[jj];
          jetype = etype[j];
          jctype = ctype[j];

          // interpolating positions at site J from nodex

          delx = xtmp;
          dely = ytmp;
          delz = ztmp;
          for (node = 0; node < npe[jetype]; node++) {
            delx -= shape_array[jetype][jintpl][node] * nodex[j][node][0];
            dely -= shape_array[jetype][jintpl][node] * nodex[j][node][1];
            delz -= shape_array[jetype][jintpl][node] * nodex[j][node][2];
          }
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq[ictype][jctype]) {
            r2inv = 1.0 / rsq;
            r6inv = r2inv * r2inv * r2inv;
            forcelj =
                r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
            fpair = forcelj * r2inv;

            fx += delx * fpair;
            fy += dely * fpair;
            fz += delz * fpair;

            // tally eng_dwl and virial if integration point is a node

            if (eflag)
              evdwl =
                  r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);
            if (evflag)
              node_ev_tally(i, inode, evdwl, 0.0, fpair, delx, dely, delz);
          }
        }
      }

      // distribute forces at integration point to nodal force columns

      for (node = 0; node < npe[ietype]; node++) {
        force_columns[node][0] +=
            weighted_shape_array[ietype][iintg][node] * fx;
        force_columns[node][1] +=
            weighted_shape_array[ietype][iintg][node] * fy;
        force_columns[node][2] +=
            weighted_shape_array[ietype][iintg][node] * fz;
      }
    }
    compute_nodef(i);
  }

  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   compute force from atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairLJMDF::pair_a2a(int i, int j, double *f) {
  double delx, dely, delz;
  double rsq, r2inv, r6inv, forcelj, fpair;
  int ictype = atom->type[i];
  int jctype = atom->type[j];
  double **x = atom->x;

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  rsq = delx * delx + dely * dely + delz * delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force

  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
  fpair = forcelj * r2inv;

  f[0] = delx * fpair;
  f[1] = dely * fpair;
  f[2] = delz * fpair;
}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairLJMDF::pair_a2ia(int i, int j, int jintpl, double *jx, double *f) {
  double delx, dely, delz;
  double rsq, r2inv, r6inv, forcelj, fpair;
  int ictype = atom->type[i];
  int jctype = element->ctype[j];

  double **x = atom->x;

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = x[i][0] - jx[0];
  dely = x[i][1] - jx[1];
  delz = x[i][2] - jx[2];
  rsq = delx * delx + dely * dely + delz * delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force

  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
  fpair = forcelj * r2inv;

  f[0] = delx * fpair;
  f[1] = dely * fpair;
  f[2] = delz * fpair;
}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on interpolated atom I
   ------------------------------------------------------------------------- */

void PairLJMDF::pair_ia2ia(int i, int iintpl, double *xtmp, int j, int jintpl,
                           double *jx, double *f) {
  double delx, dely, delz;
  double rsq, r2inv, r6inv, forcelj, fpair;
  int ictype = element->ctype[i];
  int jctype = element->ctype[j];

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = xtmp[0] - jx[0];
  dely = xtmp[1] - jx[1];
  delz = xtmp[2] - jx[2];
  rsq = delx * delx + dely * dely + delz * delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force

  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
  fpair = forcelj * r2inv;

  f[0] = delx * fpair;
  f[1] = dely * fpair;
  f[2] = delz * fpair;
}
