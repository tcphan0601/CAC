#include "pair_gauss_cut.h"
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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathConst;

PairGaussCut::PairGaussCut(CAC *cac) : Pair(cac) 
{
  // rRESPA not developed yet
  // respa_enable = 1;
  writedata = 1;
}

PairGaussCut::~PairGaussCut() 
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(hgauss);
    memory->destroy(sigmah);
    memory->destroy(rmh);
    memory->destroy(pgauss);
    memory->destroy(offset);
  }
}

/*  ----------------------------------------------------------------------  */

void PairGaussCut::compute(int eflag, int vflag) 
{
  int i, j, ii, jj, m, n, jnum;
  int ictype, jctype, ietype, jetype; // note: ictype and jctype are atomic (chemical) type, same as itype and jtype in LAMMPS

  int inode, jnode, node, jucell, iucell, igcell, jgcell, inpe, iapc;
  int japc, jnpe, jbasis, ibasis;

  double evdwl = 0.0;
  double iescale, ivscale, jescale, jvscale;
  int *jlist, *jindexlist;
  int iindex, jindex;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, fpair;
  double r, rexp, ugauss;
  double *fi, *fj;
  double *ieptr, *jeptr, *ivptr, *jvptr;

  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  //double *special_lj = force->special_lj;

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
        r = sqrt(rsq);
        rexp = (r-rmh[ictype][jctype])/sigmah[ictype][jctype];
        ugauss = pgauss[ictype][jctype] * exp(-0.5 * rexp * rexp);
        fpair = rexp/r * ugauss/sigmah[ictype][jctype];

        fi[0] += delx * fpair;
        fi[1] += dely * fpair;
        fi[2] += delz * fpair;

        if (fj != NULL) {
          fj[0] -= delx * fpair;
          fj[1] -= dely * fpair;
          fj[2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = ugauss - offset[ictype][jctype];
        }
        if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
            ieptr, jeptr, ivptr, jvptr, 
            evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/*  ----------------------------------------------------------------------
    allocate all arrays
    -------------------------------------------------------------------------  */

void PairGaussCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  memory->create(cut, n+1, n+1, "pair:cut");
  memory->create(hgauss, n+1, n+1, "pair:hgauss");
  memory->create(sigmah, n+1, n+1, "pair:sigmah");
  memory->create(rmh, n+1, n+1, "pair:rmh");
  memory->create(pgauss, n+1, n+1, "pair:pgauss");
  memory->create(offset, n+1, n+1, "pair:offset");
}

/*  ----------------------------------------------------------------------
    global settings
    -------------------------------------------------------------------------  */

void PairGaussCut::settings(int narg, char **arg)
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

/*  ----------------------------------------------------------------------
    set coeffs for one or more type pairs
    -------------------------------------------------------------------------  */

void PairGaussCut::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double hgauss_one = universe->numeric(FLERR, arg[2]);
  double rmh_one = universe->numeric(FLERR, arg[3]);
  double sigmah_one = universe->numeric(FLERR, arg[4]);
  if (sigmah_one <= 0.0)
    error->all(FLERR, "Incorrect args for pair coefficients");


  double cut_one = cut_global;
  if (narg == 6) cut_one = universe->numeric(FLERR, arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      hgauss[i][j] = hgauss_one;
      sigmah[i][j] = sigmah_one;
      rmh[i][j] = rmh_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*  ----------------------------------------------------------------------
    init for one type pair i, j and corresponding j, i
    -------------------------------------------------------------------------  */

double PairGaussCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    hgauss[i][j] = mix_energy(fabs(hgauss[i][i]), fabs(hgauss[j][j]), 
        fabs(sigmah[i][i]), fabs(sigmah[j][j]));

    // If either of the particles is repulsive (ie, if hgauss > 0), 
    // then the interaction between both is repulsive.
    double sign_hi = (hgauss[i][i] >= 0.0) ? 1.0 : -1.0;
    double sign_hj = (hgauss[j][j] >= 0.0) ? 1.0 : -1.0;
    hgauss[i][j] *= MAX(sign_hi, sign_hj);

    sigmah[i][j] = mix_distance(sigmah[i][i], sigmah[j][j]);
    rmh[i][j] = mix_distance(rmh[i][i], rmh[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  pgauss[i][j] = hgauss[i][j] / sqrt(MY_2PI) / sigmah[i][j];

  if (offset_flag) {
    double rexp = (cut[i][j]-rmh[i][j])/sigmah[i][j];
    offset[i][j] = pgauss[i][j] * exp(-0.5 * rexp * rexp);
  } else offset[i][j] = 0.0;

  hgauss[j][i] = hgauss[i][j];
  sigmah[j][i] = sigmah[i][j];
  rmh[j][i] = rmh[i][j];
  pgauss[j][i] = pgauss[i][j];
  offset[j][i] = offset[i][j];
  cut[j][i] = cut[i][j];

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
  }

  return cut[i][j];
}

/*  ----------------------------------------------------------------------
    proc 0 writes to data file
    -------------------------------------------------------------------------  */

void PairGaussCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, hgauss[i][i], rmh[i][i], sigmah[i][i]);
}

/*  ----------------------------------------------------------------------
    proc 0 writes all pairs to data file
    -------------------------------------------------------------------------  */

void PairGaussCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, hgauss[i][j], rmh[i][j], sigmah[i][j], cut[i][j]);
}

/*  ----------------------------------------------------------------------  */
double PairGaussCut::memory_usage()
{
  const int n = atom->ntypes;

  double bytes = Pair::memory_usage();

  bytes += 7 * ((n+1) * (n+1) * sizeof(double) + (n+1) * sizeof(double *));
  bytes += 1 * ((n+1) * (n+1) * sizeof(int) + (n+1) * sizeof(int *));

  return bytes;
}
