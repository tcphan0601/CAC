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

using namespace CAC_NS;
using namespace MathConst;

PairLJCut::PairLJCut(CAC *cac) : Pair(cac)
{
  // rRESPA not developed yet
  //respa_enable = 1;
  writedata = 1;

  nemax = 0;
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

/* ----------------------------------------------------------------------
   allocate all arrays
   ------------------------------------------------------------------------- */

void PairLJCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairLJCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = universe->numeric(FLERR,arg[0]);

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairLJCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  universe->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  universe->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = universe->numeric(FLERR,arg[2]);
  double sigma_one = universe->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = universe->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
   ------------------------------------------------------------------------- */

void PairLJCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairLJCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
        sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms & interpolated atoms of type I and J via Allreduce

  if (tail_flag) {
    int *atype = atom->type;
    int nalocal = atom->nlocal;
    int *ctype = element->ctype;
    int *etype = element->etype;
    int nelocal = element->nlocal;
    int *nintpl = element->nintpl;


    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nalocal; k++) {
      if (atype[k] == i) count[0] += 1.0;
      if (atype[k] == j) count[1] += 1.0;
    }
    for (int k = 0; k < nelocal; k++) {
      if (ctype[k] == i) count[0] += nintpl[etype[k]];
      if (ctype[k] == j) count[1] += nintpl[etype[k]];
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;

    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
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



/* ---------------------------------------------------------------------- */

void PairLJCut::compute(int eflag, int vflag)
{
  int jnum;
  double xtmp,ytmp,ztmp;
  double r2inv,r6inv,forcelj,factor_lj,fpair;

  int i,j,k,ii,jj,m,n;
  int ictype,jctype,ietype,jetype; // note: ictype and jctype are atomic (chemical) type, same as itype and jtype in LAMMPS

  int node,inode,jintpl,iintpl,iintg,jintg,iintpl_local,iintg_local;

  double evdwl = 0.0;
  double p,*coeff;
  int *jlist,*jindexlist;
  double ix,iy,iz,jx,jy,jz,delx,dely,delz,rsq;
  double fx,fy,fz;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;
  int newton_pair = force->newton_pair;
  //double *special_lj = force->special_lj;

  double **ex = element->x;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  double **weight = element->weight;
  int *nintg = element->nintg;
  int *nintpl = element->nintpl;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  int **i2ia = element->i2ia;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int max_nintg = element->max_nintg;
  double ***shape_array = element->shape_array;

  tagint *atag = atom->tag;
  tagint itag,jtag;

  int ainum = list->ainum;
  int einum = list->einum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int **n2ilist = list->n2ilist;
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
    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];
    itag = atag[i];
    ictype = atype[i];

    // loop over neighboring atoms, skip half if newton is on or j is owned

    jnum = numneigha2a[i];
    if(jnum) {
      jlist = firstneigha2a[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jctype = atype[j];
        jtag = atag[j];
        //factor_lj = special_lj[sbmask(j)];

        // check for unique pair only if newton is on or j is owned atom

        if (newton_pair || j < nalocal) {

          // i and j are distint atoms, only count as one pair
          // odd and even sum of tag condition is to have 
          // the pairs distributed among procs more evenly

          if (itag > jtag) {
            if ((itag+jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag+jtag) % 2 == 1) continue;
          } else {

            // if i and j have the same tag, atom j is a ghost image of atom i
            // only consider the pair from upper side in all directions since
            // force contribution from the othersize will be added in reverse_comm
            // see LAMMPS documents

            if (ax[j][2] < iz) continue;
            if (ax[j][2] == iz && ax[j][1] < iy) continue;
            if (ax[j][2] == iz && ax[j][1] == iy && ax[j][0] < ix) continue;
          }
        }

        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < cutsq[ictype][jctype]) {
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
          fpair = forcelj*r2inv;
          //fpair = factor_lj*forcelj*r2inv;

          af[i][0] += delx*fpair;
          af[i][1] += dely*fpair;
          af[i][2] += delz*fpair;

          if (newton_pair || j < nalocal) {
            af[j][0] -= delx*fpair;
            af[j][1] -= dely*fpair;
            af[j][2] -= delz*fpair;
          }

          if (eflag) {
            evdwl = r6inv*(lj3[ictype][jctype]*r6inv-lj4[ictype][jctype]) -
              offset[ictype][jctype];
            //evdwl *= factor_lj;
          }
          if (evflag) atom_ev_tally(i,j,nalocal,newton_pair,
              evdwl,0.0,fpair,delx,dely,delz);
        }
      }
    }

    // loop over neighboring interpolated atoms

    jnum = numneigha2ia[i];
    if(jnum) {
      jlist = firstneigha2ia[i];
      jindexlist = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jetype = etype[j];
        jintpl = jindexlist[jj];
        jintg = ia2i[jetype][jintpl];

        // interpolate positions at site J

        jx = jy = jz = 0.0;
        for (node = 0; node < npe; node++) {
          jx += shape_array[jetype][jintpl][node]*nodex[j][node][0];
          jy += shape_array[jetype][jintpl][node]*nodex[j][node][1];
          jz += shape_array[jetype][jintpl][node]*nodex[j][node][2];
        }

        delx = ix - jx;
        dely = iy - jy;
        delz = iz - jz;
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq[ictype][jctype]) {
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
          fpair = forcelj*r2inv;
          //fpair = factor_lj*forcelj*r2inv;

          af[i][0] += delx*fpair;
          af[i][1] += dely*fpair;
          af[i][2] += delz*fpair;

          if (eflag) {
            evdwl = r6inv*(lj3[ictype][jctype]*r6inv-lj4[ictype][jctype]) -
              offset[ictype][jctype];
            //evdwl *= factor_lj;
          }

          if (evflag) atom_ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);

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
  // 2. Loop over neighboring interpolated atoms. Add force contribution to fx fy fz
  // 3. After getting total forces at all integration points in element i, 
  //    distribute forces at integration points to nodef in element i.
  // 4. Scale down forces on nodes 

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    ictype = ctype[i];

    for (iintg = 0; iintg < nintg[ietype]; iintg++) {
      fx = fy = fz = 0.0;

      // interpolate postions of integration point from node values

      iintpl = i2ia[ietype][iintg];
      inode = i2n[ietype][iintg];
      ix = iy = iz = 0.0;
      for (node = 0; node < npe; node++) {
        ix += shape_array[ietype][iintpl][node]*nodex[i][node][0];
        iy += shape_array[ietype][iintpl][node]*nodex[i][node][1];
        iz += shape_array[ietype][iintpl][node]*nodex[i][node][2];
      }

      // loop over neighboring atoms and add force contribution

      iintg_local = e2ilist[ii] + iintg;
      jnum = numneighi2a[iintg_local];
      if (jnum) {
        jlist = firstneighi2a[iintg_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          jctype = atype[j];
          delx = ix - ax[j][0];
          dely = iy - ax[j][1];
          delz = iz - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq[ictype][jctype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
            fpair = forcelj*r2inv;

            fx += delx*fpair;
            fy += dely*fpair;
            fz += delz*fpair;

            // tally eng_dwl and virial if integration point is a node

            if (eflag) {
              evdwl = r6inv*(lj3[ictype][jctype]*r6inv-lj4[ictype][jctype]) -
                offset[ictype][jctype];
              //evdwl *= factor_lj;
            }
            if (evflag) node_ev_tally(i,inode,evdwl,0.0,fpair,delx,dely,delz);
          }
        }
      }

      // loop over neighboring interpolated atoms and accumulate force contribution

      jnum = numneighi2ia[iintg_local];
      if (jnum) {
        jlist = firstneighi2ia[iintg_local];
        jindexlist = firstneighi2ia_index[iintg_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          jintpl = jindexlist[jj];
          jetype = etype[j];
          jctype = ctype[j];
          jintg = ia2i[jetype][jintpl];

          // interpolating positions at site J from nodex

          jx = jy = jz = 0.0;
          for (node = 0; node < npe; node++) {
            jx += shape_array[jetype][jintpl][node]*nodex[j][node][0];
            jy += shape_array[jetype][jintpl][node]*nodex[j][node][1];
            jz += shape_array[jetype][jintpl][node]*nodex[j][node][2];
          }
          delx = ix - jx;
          dely = iy - jy;
          delz = iz - jz;
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq[ictype][jctype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
            fpair = forcelj*r2inv;

            fx += delx*fpair;
            fy += dely*fpair;
            fz += delz*fpair;

            // tally eng_dwl and virial if integration point is a node

            if (eflag) {
              evdwl = r6inv*(lj3[ictype][jctype]*r6inv-lj4[ictype][jctype]) -
                offset[ictype][jctype];
            }
            if (evflag) node_ev_tally(i,inode,evdwl,0.0,fpair,delx,dely,delz);
          }
        }
      }

      // distribute forces at integration point to forces at nodes 

      for (node = 0; node < npe; node++) {
        nodef[i][node][0] += weight[ietype][iintg]*shape_array[ietype][iintpl][node]*fx;
        nodef[i][node][1] += weight[ietype][iintg]*shape_array[ietype][iintpl][node]*fy;
        nodef[i][node][2] += weight[ietype][iintg]*shape_array[ietype][iintpl][node]*fz;
      }
    }

    // scale back forces on node

    for (node = 0; node < npe; node++) {
      nodef[i][node][0] *= elemscale[ietype];
      nodef[i][node][1] *= elemscale[ietype];
      nodef[i][node][2] *= elemscale[ietype];
    }
  }
}

/* ----------------------------------------------------------------------
   compute force from atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairLJCut::pair_a2a(int i, int j, double *f)
{
  double delx,dely,delz;
  double rsq,r2inv,r6inv,forcelj,fpair;
  int ictype = atom->type[i];
  int jctype = atom->type[j];
  double **x = atom->x;

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force 

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
  fpair = forcelj*r2inv;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;

}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairLJCut::pair_a2ia(int i, int j, int jintpl, double *jx, double *f)
{
  double delx,dely,delz;
  double rsq,r2inv,r6inv,forcelj,fpair;
  int ictype = atom->type[i];
  int jctype = element->ctype[j];
  int jetype = element->etype[j];
  double **x = atom->x;

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = x[i][0] - jx[0];
  dely = x[i][1] - jx[1];
  delz = x[i][2] - jx[2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force 

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
  fpair = forcelj*r2inv;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;
}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on interpolated atom I
   ------------------------------------------------------------------------- */

void PairLJCut::pair_ia2ia(int i, int iintpl, double *ix, int j, int jintpl, double *jx, double *f)
{
  double delx,dely,delz;
  double rsq,r2inv,r6inv,forcelj,fpair;
  int ictype = element->ctype[i];
  int ietype = element->etype[i];
  int jctype = element->ctype[j];
  int jetype = element->etype[j];

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = ix[0] - jx[0];
  dely = ix[1] - jx[1];
  delz = ix[2] - jx[2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutsq[ictype][jctype]) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force 

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[ictype][jctype]*r6inv - lj2[ictype][jctype]);
  fpair = forcelj*r2inv;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;
}

/* ----------------------------------------------------------------------
   memory usage of local elem-based arrays
   ------------------------------------------------------------------------- */

double PairLJCut::memory_usage()
{
  double bytes = nemax * element->max_nintg * 3 * sizeof(double);
  return bytes;
}
