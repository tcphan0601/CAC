#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut.h"
#include "atom.h"
#include "force.h"
#include "element.h"
#include "comm.h"
#include "neighbor.h"
#include "universe.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "math_const.h"

#include "update.h"

using namespace CAC_NS;
using namespace MathConst;

PairLJCut::PairLJCut(CAC *cac) : Pair(cac)
{
  // rRESPA not developed yet
  //respa_enable = 1;
  writedata = 1;
}


PairLJCut::~PairLJCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}

/*  ----------------------------------------------------------------------
   allocate all arrays
   -------------------------------------------------------------------------  */

void PairLJCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(cut, n+1, n+1, "pair:cut");
  memory->create(epsilon, n+1, n+1, "pair:epsilon");
  memory->create(sigma, n+1, n+1, "pair:sigma");
  memory->create(lj1, n+1, n+1, "pair:lj1");
  memory->create(lj2, n+1, n+1, "pair:lj2");
  memory->create(lj3, n+1, n+1, "pair:lj3");
  memory->create(lj4, n+1, n+1, "pair:lj4");
  memory->create(offset, n+1, n+1, "pair:offset");
}

/*  ----------------------------------------------------------------------
   global settings
   -------------------------------------------------------------------------  */

void PairLJCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = universe->numeric(FLERR, arg[0]);

  // reset cutoffs that have been explicitly set
  //
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

}

/*  ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   -------------------------------------------------------------------------  */

void PairLJCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double epsilon_one = universe->numeric(FLERR, arg[2]);
  double sigma_one = universe->numeric(FLERR, arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = universe->numeric(FLERR, arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*  ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
   -------------------------------------------------------------------------  */

void PairLJCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/*  ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
   -------------------------------------------------------------------------  */

double PairLJCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], 
        sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR, "Pair cutoff < Respa interior cutoff");

  // compute I, J contribution to long-range tail correction
  // count total # of atoms & virtual atoms of type I and J via Allreduce

  if (tail_flag) {
    int *atype = atom->type;
    int nalocal = atom->nlocal;
    int **ctype = element->ctype;
    int *etype = element->etype;
    int *apc = element->apc;
    int nelocal = element->nlocal;
    int *nucell = element->nucell;


    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nalocal; k++) {
      if (atype[k] == i) count[0] += 1.0;
      if (atype[k] == j) count[1] += 1.0;
    }
    for (int k = 0; k < nelocal; k++) 
      for (int l = 0; l < apc[etype[k]]; l++) {
        if (ctype[k][l] == i) count[0] += nucell[etype[k]];
        if (ctype[k][l] == j) count[1] += nucell[etype[k]];
      }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;

    etail_ij = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] *
      sig6 * (sig6 - 3.0 * rc6) / (9.0 * rc9);
    ptail_ij = 16.0 * MY_PI * all[0] * all[1] * epsilon[i][j] *
      sig6 * (2.0 * sig6 - 3.0 * rc6) / (9.0 * rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairLJCut::init_style()
{
  int irequest;

  irequest = neighbor->request(this,instance_me);
  cut_respa = NULL;
}

/*  ----------------------------------------------------------------------  */

void PairLJCut::compute(int eflag, int vflag)
{
  double r2inv, r6inv, forcelj, factor_lj;
  int i, j, ii, jj, m, n, jnum;
  int ictype, jctype, ietype, jetype; // note: ictype and jctype are atomic (chemical) type, same as itype and jtype in LAMMPS

  int inode, jnode, node, jucell, iucell, igcell, jgcell, inpe, iapc;
  int japc, jnpe, jbasis, ibasis;

  double evdwl = 0.0;
  double iescale, ivscale, jescale, jvscale;
  int *jlist, *jindexlist;
  int iindex, jindex;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, fpair;
  double *fi, *fj;
  double *ieptr, *jeptr, *ivptr, *jvptr;

  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  int newton_pair = force->newton_pair;

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


  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  // compute forces on each ATOM
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];
    ieptr = ivptr = NULL;
    iescale = ivscale = 0.0;

    // i is atom

    if (iindex < 0) {
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
      ictype = atype[i];
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
        if (eflag_global) iescale = nodal_weight[ietype]/2.0;
        if (vflag_global) ivscale = nodal_weight[ietype]/2.0;
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
        r2inv = 1.0/rsq;
        r6inv = r2inv * r2inv * r2inv;
        forcelj = r6inv * (lj1[ictype][jctype] * r6inv - lj2[ictype][jctype]);
        fpair = forcelj * r2inv;
        //fpair = factor_lj * forcelj * r2inv;

        fi[0] += delx * fpair;
        fi[1] += dely * fpair;
        fi[2] += delz * fpair;

        if (fj) {
          fj[0] -= delx * fpair;
          fj[1] -= dely * fpair;
          fj[2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = r6inv * (lj3[ictype][jctype] * r6inv-lj4[ictype][jctype]) -
            offset[ictype][jctype];
          //evdwl *= factor_lj;
        }
        if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
            ieptr, jeptr, ivptr, jvptr, 
            evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  // debug
  if (debug_write_flag) {
    fclose(fp);
    debug_write_flag = 0;
  }
  // end debug

  if (vflag_fdotr) virial_fdotr_compute();
}

