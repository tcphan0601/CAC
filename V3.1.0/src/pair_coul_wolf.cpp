#include "pair_coul_wolf.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "force.h"
#include "universe.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace MathConst;

//#define TEST 48
//#define TESTUCELL 2743
#define TEST 46
#define TESTUCELL 13
#define TESTBASIS 1
#define DEBUG 0
/* ---------------------------------------------------------------------- */

PairCoulWolf::PairCoulWolf(CAC *cac) : Pair(cac)
{
  single_enable = 0;        // NOTE: single() method below is not yet correct
  //centroidstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairCoulWolf::~PairCoulWolf()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulWolf::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double *fi, *fj;
  double *ieptr, *jeptr, *ivptr, *jvptr;
  double iescale, ivscale, jescale, jvscale;
  int iindex, jindex;
  double rsq, forcecoul, factor_coul;
  double prefactor;
  double r;
  int *ilist, *iindexlist, *jlist, *jindexlist, *numneigh, **firstneigh, **firstneighindex;
  double erfcc, erfcd, v_sh, dvdrr, e_self, e_shift, f_shift, qisq;

  int ictype, jctype, ietype, jetype; // note: ictype and jctype are atomic (chemical) type, same as itype and jtype in LAMMPS

  int inode, jnode, node, jucell, iucell, igcell, jgcell, inpe, iapc;
  int japc, jnpe, jbasis, ibasis;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  double *q = atom->q;  // note: this q array is per-type, not per-atom as in LAMMPS

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

  //double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  // self and shifted coulombic energy

  e_self = v_sh = 0.0;
  e_shift = erfc(alf * cut_coul) / cut_coul;
  f_shift = -(e_shift + 2.0 * alf / MY_PIS * exp(-alf * alf * cut_coul * cut_coul)) /
    cut_coul;

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
      qtmp = q[ictype];
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
      qtmp = q[ictype];
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

    qisq = qtmp * qtmp;
    e_self = -(e_shift / 2.0 + alf / MY_PIS) * qisq * qqrd2e;
    if (evflag) ev_tally(iescale * 2, ivscale * 2, 
        ieptr, ieptr, ivptr, ivptr, 
        0.0, e_self, 0.0, 0.0, 0.0, 0.0);


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];

      // special bonds not implemented yet
      //factor_coul = special_coul[sbmask(j)];
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

      if (rsq < cut_coulsq) {
        r = sqrt(rsq);
        prefactor = qqrd2e * qtmp * q[jctype] / r;
        erfcc = erfc(alf * r);
        erfcd = exp(-alf * alf * rsq);
        v_sh = (erfcc - e_shift * r) * prefactor;
        dvdrr = (erfcc / rsq + 2.0 * alf / MY_PIS * erfcd / r) + f_shift;
        forcecoul = dvdrr * rsq * prefactor;
        //if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
        fpair = forcecoul / rsq;

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
          printf("pair force coul i = %d fpair = %g f = %g %g %g j = %d del = %g %g %g\n",atom->tag[i],fpair,delx*fpair,dely*fpair,delz*fpair,atom->tag[j],delx,dely,delz);
        if (atom->tag[j] == TEST)
          printf("pair force coul j = %d fpair = %g f = %g %g %g i = %d del = %g %g %g\n",atom->tag[j],fpair,-delx*fpair,-dely*fpair,-delz*fpair,atom->tag[i],delx,dely,delz);
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
          ecoul = v_sh;
          //if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        } else ecoul = 0.0;

        if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
            ieptr, jeptr, ivptr, jvptr, 
            0.0, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairCoulWolf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
   unlike other pair styles, 
   there are no individual pair settings that these override
   ------------------------------------------------------------------------- */

void PairCoulWolf::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  alf = universe->numeric(FLERR, arg[0]);
  cut_coul = universe->numeric(FLERR, arg[1]);
}

/* ----------------------------------------------------------------------
   set cutoffs for one or more type pairs, optional
   ------------------------------------------------------------------------- */

void PairCoulWolf::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairCoulWolf::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair coul/wolf requires per-type attribute q through charge command");
  atom->check_charge(FLERR);

  neighbor->request(this, instance_me);

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
   ------------------------------------------------------------------------- */

double PairCoulWolf::init_one(int /*i*/, int /*j*/)
{
  return cut_coul;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
   ------------------------------------------------------------------------- */
/*
   void PairCoulWolf::write_restart(FILE *fp)
   {
   write_restart_settings(fp);

   int i, j;
   for (i = 1; i <= atom->ntypes; i++)
   for (j = i; j <= atom->ntypes; j++)
   fwrite(&setflag[i][j], sizeof(int), 1, fp);
   }
   */
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */
/*
   void PairCoulWolf::read_restart(FILE *fp)
   {
   read_restart_settings(fp);
   allocate();

   int i, j;
   int me = comm->me;
   for (i = 1; i <= atom->ntypes; i++)
   for (j = i; j <= atom->ntypes; j++) {
   if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
   MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
   }
   }
   */
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
   ------------------------------------------------------------------------- */
/*
   void PairCoulWolf::write_restart_settings(FILE *fp)
   {
   fwrite(&alf, sizeof(double), 1, fp);
   fwrite(&cut_coul, sizeof(double), 1, fp);
   fwrite(&offset_flag, sizeof(int), 1, fp);
   fwrite(&mix_flag, sizeof(int), 1, fp);
   }
   */
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */
/*
   void PairCoulWolf::read_restart_settings(FILE *fp)
   {
   if (comm->me == 0) {
   utils::sfread(FLERR, &alf, sizeof(double), 1, fp, nullptr, error);
   utils::sfread(FLERR, &cut_coul, sizeof(double), 1, fp, nullptr, error);
   utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
   utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
   }
   MPI_Bcast(&alf, 1, MPI_DOUBLE, 0, world);
   MPI_Bcast(&cut_coul, 1, MPI_DOUBLE, 0, world);
   MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
   MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
   }
   */
/* ----------------------------------------------------------------------
   only the pair part is calculated here
   ------------------------------------------------------------------------- */

//double PairCoulWolf::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq, 
//                           double factor_coul, double /*factor_lj*/, 
//                            double &fforce)
/*
   {
   double r, prefactor;
   double forcecoul, phicoul;
   double e_shift, f_shift, dvdrr, erfcc, erfcd;

   e_shift = erfc(alf*cut_coul) / cut_coul;
   f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
   cut_coul;

   if (rsq < cut_coulsq) {
   r = sqrt(rsq);
   prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
   erfcc = erfc(alf*r);
   erfcd = exp(-alf*alf*r*r);
   dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
   forcecoul = dvdrr*rsq*prefactor;
   if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
   } else forcecoul = 0.0;
   fforce = forcecoul / rsq;

   double eng = 0.0;
   if (rsq < cut_coulsq) {
   phicoul = prefactor * (erfcc-e_shift*r);
   if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
   eng += phicoul;
   }
   return eng;
   }
   */
