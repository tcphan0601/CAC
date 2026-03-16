#include "pair_lj_mdf.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

PairLJMDF::PairLJMDF(CAC *cac) : Pair(cac) {
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJMDF::~PairLJMDF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
    memory->destroy(cut_inner_sq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJMDF::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, ictype, jctype, ietype, jetype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
  double *fi, *fj;
  double *ieptr, *jeptr, *ivptr, *jvptr;
  double iescale, ivscale, jescale, jvscale;
  int iindex, jindex;

  int *ilist, *iindexlist, *jlist, *jindexlist, *numneigh, **firstneigh, **firstneighindex;
  int inode, jnode, node, jucell, iucell, igcell, jgcell, inpe, iapc;
  int japc, jnpe, jbasis, ibasis;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;

  int *npe = element->npe;
  int *apc = element->apc;
  double **ex = element->x;
  double ****nodex = element->nodex;
  double ****gaussf = element->gaussf;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int nelocal = element->nlocal;
  int **g2u = element->g2u;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  double ***shape_array = element->shape_array;
  double ***weighted_shape_array = element->weighted_shape_array;
  double *nodal_weight = element->nodal_weight;

  //double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double rr, d, dd, tt, dt, dp, philj;

  inum = list->inum;
  ilist = list->ilist;
  iindexlist = list->iindexlist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firstneighindex = list->firstneighindex;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];

    ieptr = ivptr = NULL;
    iescale = ivscale = 0.0;

    // i is atom

    if (iindex < 0) {
      ictype = atype[i];
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
      fi = af[i];
      if (eflag_global) iescale = 0.5;
      if (vflag_global) ivscale = 0.5;
      if (eflag_atom) ieptr = &eatom[i];
      if (vflag_atom) ivptr = vatom[i];
    } 

    // i is gauss point

    else {
      ietype = etype[i];
      iapc = apc[ietype];
      inpe = npe[ietype];
      igcell = iindex / iapc;
      ibasis = iindex % iapc;
      ictype = ctype[i][ibasis];
      inode = g2n[ietype][igcell];
      fi = gaussf[i][ibasis][igcell];
      if (inode < 0) {
        iucell = g2u[ietype][igcell];
        xtmp = ytmp = ztmp = 0.0;
        for (node = 0; node < inpe; node++) {
          xtmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][0];
          ytmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][1];
          ztmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][2];
        }
      } else {
        xtmp = nodex[i][ibasis][inode][0];
        ytmp = nodex[i][ibasis][inode][1];
        ztmp = nodex[i][ibasis][inode][2];
        if (eflag_global) iescale = nodal_weight[ietype] / 2.0;
        if (vflag_global) ivscale = nodal_weight[ietype] / 2.0;
        if (eflag_atom) ieptr = &enode[i][ibasis][inode];
        if (vflag_atom) ivptr = vnode[i][ibasis][inode];
      }

    }

    // loop over neighbors

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      //factor_lj = special_lj[sbmask(j)];
      //j &= NEIGHMASK;

      delx = xtmp; dely = ytmp; delz = ztmp;
      jescale = jvscale = 0.0;
      fj = jeptr = jvptr = NULL;

      // j is atom

      if (jindex < 0) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
        jctype = atype[j];
        if (newton_pair || j < nalocal) {
          fj = af[j];
          if (eflag_atom) jeptr = &eatom[j];
          if (vflag_atom) jvptr = vatom[j];
          if (eflag_global) jescale = 0.5;
          if (vflag_global) jvscale = 0.5;
        }
      } 

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jnpe = npe[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = ctype[j][jbasis];
        for (node = 0; node < jnpe; node++) {
          delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
          dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
          delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
        }
        if (newton_pair || j < nelocal) {
          jgcell = u2g[jetype][jucell];
          if (jgcell >= 0) {
            jnode = g2n[jetype][jgcell];
            fj = gaussf[j][jbasis][jgcell];
            if (jnode >= 0) {
              if (eflag_atom) jeptr = &enode[j][jbasis][jnode];
              if (vflag_atom) jvptr = vnode[j][jbasis][jnode];
              if (eflag_global) jescale = nodal_weight[jetype]/2.0;
              if (vflag_global) jvscale = nodal_weight[jetype]/2.0;
            }
          }
        }
      }

      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[ictype][jctype]) {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;
        forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);

        if (rsq > cut_inner_sq[ictype][jctype]) {
          philj = r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);

          rr = sqrt(rsq);
          dp = (cut[ictype][jctype] - cut_inner[ictype][jctype]);
          d = (rr - cut_inner[ictype][jctype]) / dp;
          dd = 1.0 - d;
          // taperig function - mdf style
          tt = (1.0 + 3.0 * d + 6.0 * d * d) * dd * dd * dd;
          // minus the derivative of the tapering function
          dt = 30.0 * d * d * dd * dd * rr / dp;

          forcelj = forcelj * tt + philj * dt;
        } else {
          tt = 1;
        }
        //fpair = factor_lj * forcelj * r2inv;
        fpair = forcelj * r2inv;

        fi[0] += delx * fpair;
        fi[1] += dely * fpair;
        fi[2] += delz * fpair;
        if (fj != NULL) {
          fj[0] -= delx * fpair;
          fj[1] -= dely * fpair;
          fj[2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = r6inv * (lj3[ictype][jctype] * r6inv - lj4[ictype][jctype]);
          if (rsq > cut_inner_sq[ictype][jctype]) evdwl *= tt;

          //evdwl *= factor_lj;

          if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
              ieptr, jeptr, ivptr, jvptr, 
              evdwl, 0.0, fpair, delx, dely, delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairLJMDF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  memory->create(cut, n+1, n+1, "pair:cut");
  memory->create(cut_inner, n+1, n+1, "pair:cut_inner");
  memory->create(cut_inner_sq, n+1, n+1, "pair:cut_inner_sq");
  memory->create(epsilon, n+1, n+1, "pair:epsilon");
  memory->create(sigma, n+1, n+1, "pair:sigma");
  memory->create(lj1, n+1, n+1, "pair:lj1");
  memory->create(lj2, n+1, n+1, "pair:lj2");
  memory->create(lj3, n+1, n+1, "pair:lj3");
  memory->create(lj4, n+1, n+1, "pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairLJMDF::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  cut_inner_global = universe->numeric(FLERR, arg[0]);
  cut_global = universe->numeric(FLERR, arg[1]);

  if (cut_inner_global <= 0.0 || cut_inner_global > cut_global)
    error->all(FLERR, "Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairLJMDF::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 6)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double epsilon_one = universe->numeric(FLERR, arg[2]);
  double sigma_one = universe->numeric(FLERR, arg[3]);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 6) {
    cut_inner_one = universe->numeric(FLERR, arg[4]);
    cut_one = universe->numeric(FLERR, arg[5]);
  }
  if (cut_inner_global <= 0.0 || cut_inner_global > cut_global)
    error->all(FLERR, "Illegal pair_style command");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
   ------------------------------------------------------------------------- */

double PairLJMDF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], 
        sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut_inner[i][j] = mix_distance(cut_inner[i][i], cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  cut_inner_sq[i][j] = cut_inner[i][j]*cut_inner[i][j];
  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  cut[j][i] = cut[i][j];  // BUG FIX
  cut_inner[j][i] = cut_inner[i][j];
  cut_inner_sq[j][i] = cut_inner_sq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::write_restart(FILE *fp)
   {
   write_restart_settings(fp);

   int i, j;
   for (i = 1; i <= atom->ntypes; i++)
   for (j = i; j <= atom->ntypes; j++) {
   fwrite(&setflag[i][j], sizeof(int), 1, fp);
   if (setflag[i][j]) {
   fwrite(&epsilon[i][j], sizeof(double), 1, fp);
   fwrite(&sigma[i][j], sizeof(double), 1, fp);
   fwrite(&cut_inner[i][j], sizeof(double), 1, fp);
   fwrite(&cut[i][j], sizeof(double), 1, fp);
   }
   }
   }
   */
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::read_restart(FILE *fp)
   {
   read_restart_settings(fp);
   allocate();

   int i, j;
   int me = comm->me;
   for (i = 1; i <= atom->ntypes; i++)
   for (j = i; j <= atom->ntypes; j++) {
   if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, NULL, error);
   MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
   if (setflag[i][j]) {
   if (me == 0) {
   utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, NULL, error);
   utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, NULL, error);
   utils::sfread(FLERR, &cut_inner[i][j], sizeof(double), 1, fp, NULL, error);
   utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, NULL, error);
   }
   MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
   MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
   MPI_Bcast(&cut_inner[i][j], 1, MPI_DOUBLE, 0, world);
   MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
   }
   }
   }
   */
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::write_restart_settings(FILE *fp)
   {
   fwrite(&mix_flag, sizeof(int), 1, fp);
   }
   */
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::read_restart_settings(FILE *fp)
   {
   int me = comm->me;
   if (me == 0) {
   utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, NULL, error);
   }
   MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
   }
   */
/* ----------------------------------------------------------------------
   proc 0 writes to data file
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::write_data(FILE *fp)
   {
   for (int i = 1; i <= atom->ntypes; i++)
   fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
   }
   */
/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
   ------------------------------------------------------------------------- */
/*
   void PairLJMDF::write_data_all(FILE *fp)
   {
   for (int i = 1; i <= atom->ntypes; i++)
   for (int j = i; j <= atom->ntypes; j++)
   fprintf(fp, "%d %d %g %g %g %g\n", 
   i, j, epsilon[i][j], sigma[i][j], 
   cut_inner[i][j], cut[i][j]);
   }
   */
/* ---------------------------------------------------------------------- */

//double PairLJMDF::single(int /*i*/, int /*j*/, int ictype, int jctype, 
//                             double rsq, 
//                             double /*factor_coul*/, double factor_lj, 
//                             double &fforce)
/*{
  double r2inv, r6inv, forcelj, philj;
  double rr, dp, d, tt, dt, dd;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;

  philj = r6inv*(lj3[ictype][jctype]*r6inv-lj4[ictype][jctype]);
  forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);

  if (rsq > cut_inner_sq[ictype][jctype]) {

  rr = sqrt(rsq);
  dp = (cut[ictype][jctype] - cut_inner[ictype][jctype]);
  d = (rr - cut_inner[ictype][jctype]) / dp;
  dd = 1-d;
  tt = (1. + 3.*d + 6.*d*d)* dd*dd*dd;
  dt = 30.* d*d * dd*dd * rr / dp;

  forcelj = forcelj*tt + philj*dt;
  philj *= tt;
  }

  fforce = factor_lj*forcelj*r2inv;

  return factor_lj*philj;
  }
  */
/* ---------------------------------------------------------------------- */
/*
   void *PairLJMDF::extract(const char *str, int &dim)
   {
   dim = 2;
   if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
   if (strcmp(str, "sigma") == 0) return (void *) sigma;
   return NULL;
   }
   */
