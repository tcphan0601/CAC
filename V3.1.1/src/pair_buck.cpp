#include "pair_buck.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathConst;

//#define TEST 48
//#define TESTUCELL 2743
#define TEST 46
#define TESTUCELL 13
#define TESTBASIS 1
#define DEBUG  0
/* ---------------------------------------------------------------------- */

PairBuck::PairBuck(CAC *cac) : Pair(cac)
{
  writedata = 1;
  //centroidstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairBuck::~PairBuck()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(a);
    memory->destroy(rho);
    memory->destroy(c);
    memory->destroy(rhoinv);
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBuck::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, ictype, jctype, ietype, jetype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcebuck, factor_lj;
  double r, rexp;
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

    ieptr = ivptr = nullptr;
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
      iucell = g2u[ietype][igcell];
      if (inode < 0) {
        //iucell = g2u[ietype][igcell];
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
      fj = jeptr = jvptr = nullptr;

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
        r = sqrt(rsq);
        rexp = exp(-r * rhoinv[ictype][jctype]);
        forcebuck = buck1[ictype][jctype] * r * rexp - buck2[ictype][jctype] * r6inv;
        //fpair = factor_lj*forcebuck*r2inv;
        fpair = forcebuck * r2inv;

        fi[0] += delx * fpair;
        fi[1] += dely * fpair;
        fi[2] += delz * fpair;
        if (fj != nullptr) {
          fj[0] -= delx * fpair;
          fj[1] -= dely * fpair;
          fj[2] -= delz * fpair;
        }
        if (DEBUG) {
/*
        if (atom->tag[i] == TEST)
          printf("pair force buck i = %d fpair = %g f = %g %g %g j = %d del = %g %g %g\n",atom->tag[i],fpair,delx*fpair,dely*fpair,delz*fpair,atom->tag[j],delx,dely,delz);
        if (atom->tag[j] == TEST)
          printf("pair force buck j = %d fpair = %g f = %g %g %g i = %d del = %g %g %g\n",atom->tag[j],fpair,-delx*fpair,-dely*fpair,-delz*fpair,atom->tag[i],delx,dely,delz);

*/
        if (element->tag[i] == TEST && iucell == TESTUCELL && ibasis == TESTBASIS)
          //printf("pair force buck i = %d %d %d r = %g j = %d %d %d fpair = %g f = %g %g %g jx = %g %g %g \n"
          //printf("%d %d %d %g %d %d %d %g %-1.16e %-1.16e %-1.16e %g %g %g \n"
          printf("%-1.16e %-1.16e %-1.16e\n"
              //,element->tag[i],iucell,ibasis,sqrt(rsq)
              //,element->tag[j],jucell,jbasis
              //,fpair
              ,delx*fpair,dely*fpair,delz*fpair
              //,xtmp-delx,ytmp-dely,ztmp-delz
              );
        if (element->tag[j] == TEST && jucell == TESTUCELL && jbasis == TESTBASIS)
          //printf("pair force buck i = %d %d %d r = %g j = %d %d %d fpair = %g f = %g %g %g jx = %g %g %g \n"
          //printf("%d %d %d %g %d %d %d %g %-1.16e %-1.16e %-1.16e %g %g %g \n"
          printf("%-1.16e %-1.16e %-1.16e\n"
              //,element->tag[j],jucell,jbasis,sqrt(rsq)
              //,element->tag[i],iucell,ibasis
              //,fpair
              ,-delx*fpair,-dely*fpair,-delz*fpair
              //,xtmp,ytmp,ztmp
              );
        }
        if (eflag) {
          evdwl = a[ictype][jctype] * rexp - c[ictype][jctype] * r6inv -
            offset[ictype][jctype];
          //evdwl *= factor_lj;
        }

        if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
            ieptr, jeptr, ivptr, jvptr, 
            evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBuck::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  memory->create(cut, n+1, n+1, "pair:cut_lj");
  memory->create(a, n+1, n+1, "pair:a");
  memory->create(rho, n+1, n+1, "pair:rho");
  memory->create(c, n+1, n+1, "pair:c");
  memory->create(rhoinv, n+1, n+1, "pair:rhoinv");
  memory->create(buck1, n+1, n+1, "pair:buck1");
  memory->create(buck2, n+1, n+1, "pair:buck2");
  memory->create(offset, n+1, n+1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBuck::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = universe->numeric(FLERR, arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuck::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double a_one = universe->numeric(FLERR, arg[2]);
  double rho_one = universe->numeric(FLERR, arg[3]);
  if (rho_one <= 0) error->all(FLERR, "Incorrect args for pair coefficients");
  double c_one = universe->numeric(FLERR, arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = universe->numeric(FLERR, arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      c[i][j] = c_one;
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

double PairBuck::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  rhoinv[i][j] = 1.0 / rho[i][j];
  buck1[i][j] = a[i][j] / rho[i][j];
  buck2[i][j] = 6.0 * c[i][j];

  if (offset_flag && (cut[i][j] > 0.0)) {
    double rexp = exp(-cut[i][j] / rho[i][j]);
    offset[i][j] = a[i][j] * rexp - c[i][j] / pow(cut[i][j], 6.0);
  } else offset[i][j] = 0.0;

  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  // compute I, J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double rho1 = rho[i][j];
    double rho2 = rho1 * rho1;
    double rho3 = rho2 * rho1;
    double rc = cut[i][j];
    double rc2 = rc * rc;
    double rc3 = rc2 * rc;
    etail_ij = 2.0 * MY_PI * all[0] * all[1] *
      (a[i][j] * exp(-rc/rho1) * rho1 * (rc2 + 2.0 * rho1 * rc + 2.0 * rho2) -
       c[i][j] / (3.0 * rc3));
    ptail_ij = (-1 / 3.0) * 2.0 * MY_PI * all[0] * all[1] *
      (-a[i][j] * exp(-rc / rho1) *
       (rc3 + 3.0 * rho1 * rc2 + 6.0 * rho2 * rc + 6.0 * rho3) + 2.0 * c[i][j] / rc3);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */
/*
void PairBuck::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j], sizeof(double), 1, fp);
        fwrite(&rho[i][j], sizeof(double), 1, fp);
        fwrite(&c[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}
*/
/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */
/*
void PairBuck::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &a[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &rho[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &c[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&a[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&rho[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&c[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}
*/
/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */
/*
void PairBuck::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}
*/
/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */
/*
void PairBuck::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}
*/
/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */
/*
void PairBuck::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, a[i][i], rho[i][i], c[i][i]);
}
*/
/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */
/*
void PairBuck::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, 
              a[i][j], rho[i][j], c[i][j], cut[i][j]);
}
*/
/* ---------------------------------------------------------------------- */

//double PairBuck::single(int /*i*/, int /*j*/, int ictype, int jctype, 
//                        double rsq, double /*factor_coul*/, double factor_lj, 
//                        double &fforce)
/*
{
  double r2inv, r6inv, r, rexp, forcebuck, phibuck;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  r = sqrt(rsq);
  rexp = exp(-r*rhoinv[ictype][jctype]);
  forcebuck = buck1[ictype][jctype]*r*rexp - buck2[ictype][jctype]*r6inv;
  fforce = factor_lj*forcebuck*r2inv;

  phibuck = a[ictype][jctype]*rexp - c[ictype][jctype]*r6inv -
    offset[ictype][jctype];
  return factor_lj*phibuck;
}
*/
/* ---------------------------------------------------------------------- */
/*
void *PairBuck::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "a") == 0) return (void *) a;
  if (strcmp(str, "c") == 0) return (void *) c;
  return nullptr;
}
*/
