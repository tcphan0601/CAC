#include "meam.h"
#include <cmath>
#include "memory.h"
#include "math_special.h"
#include "atom.h"
#include "element.h"
#include "error.h"

using namespace CAC_NS;
#define TEST 12
#define TESTINTG 5

void
MEAM::meam_dens_setup(int atom_nmax, int elem_nmax, int atom_nall, int elem_nall, int n_neigh)
{
  int i, j, k;

  // grow local arrays if necessary

  if (atom_nmax > namax) {
    memory->destroy(atomrho);
    memory->destroy(atomrho0);
    memory->destroy(atomrho1);
    memory->destroy(atomrho2);
    memory->destroy(atomrho3);
    memory->destroy(atomfrhop);
    memory->destroy(atomgamma);
    memory->destroy(atomdgamma1);
    memory->destroy(atomdgamma2);
    memory->destroy(atomdgamma3);
    memory->destroy(atomarho2b);
    memory->destroy(atomarho1);
    memory->destroy(atomarho2);
    memory->destroy(atomarho3);
    memory->destroy(atomarho3b);
    memory->destroy(atomt_ave);
    memory->destroy(atomtsq_ave);

    namax = atom_nmax;

    memory->create(atomrho, namax, "pair:atomrho");
    memory->create(atomrho0, namax, "pair:atomrho0");
    memory->create(atomrho1, namax, "pair:atomrho1");
    memory->create(atomrho2, namax, "pair:atomrho2");
    memory->create(atomrho3, namax, "pair:atomrho3");
    memory->create(atomfrhop, namax, "pair:atomfrhop");
    memory->create(atomgamma, namax, "pair:atomgamma");
    memory->create(atomdgamma1, namax, "pair:atomdgamma1");
    memory->create(atomdgamma2, namax, "pair:atomdgamma2");
    memory->create(atomdgamma3, namax, "pair:atomdgamma3");
    memory->create(atomarho2b, namax, "pair:atomarho2b");
    memory->create(atomarho1, namax, 3, "pair:atomarho1");
    memory->create(atomarho2, namax, 6, "pair:atomarho2");
    memory->create(atomarho3, namax, 10, "pair:atomarho3");
    memory->create(atomarho3b, namax, 3, "pair:atomarho3b");
    memory->create(atomt_ave, namax, 3, "pair:atomt_ave");
    memory->create(atomtsq_ave, namax, 3, "pair:atomtsq_ave");
  }
  int max_npe = element->max_npe;
  if (elem_nmax > nemax) {
    memory->destroy(noderho);
    memory->destroy(noderho0);
    memory->destroy(noderho1);
    memory->destroy(noderho2);
    memory->destroy(noderho3);
    memory->destroy(nodefrhop);
    memory->destroy(nodegamma);
    memory->destroy(nodedgamma1);
    memory->destroy(nodedgamma2);
    memory->destroy(nodedgamma3);
    memory->destroy(nodearho2b);
    memory->destroy(nodearho1);
    memory->destroy(nodearho2);
    memory->destroy(nodearho3);
    memory->destroy(nodearho3b);
    memory->destroy(nodet_ave);
    memory->destroy(nodetsq_ave);

    nemax = elem_nmax;

    memory->create(noderho, nemax, max_npe, "pair:noderho");
    memory->create(noderho0, nemax, max_npe, "pair:noderho0");
    memory->create(noderho1, nemax, max_npe, "pair:noderho1");
    memory->create(noderho2, nemax, max_npe, "pair:noderho2");
    memory->create(noderho3, nemax, max_npe, "pair:noderho3");
    memory->create(nodefrhop, nemax, max_npe, "pair:nodefrhop");
    memory->create(nodegamma, nemax, max_npe, "pair:nodegamma");
    memory->create(nodedgamma1, nemax, max_npe, "pair:nodedgamma1");
    memory->create(nodedgamma2, nemax, max_npe, "pair:nodedgamma2");
    memory->create(nodedgamma3, nemax, max_npe, "pair:nodedgamma3");
    memory->create(nodearho2b, nemax, max_npe, "pair:nodearho2b");
    memory->create(nodearho1, nemax, max_npe, 3, "pair:nodearho1");
    memory->create(nodearho2, nemax, max_npe, 6, "pair:nodearho2");
    memory->create(nodearho3, nemax, max_npe, 10, "pair:nodearho3");
    memory->create(nodearho3b, nemax, max_npe, 3, "pair:nodearho3b");
    memory->create(nodet_ave, nemax, max_npe, 3, "pair:nodet_ave");
    memory->create(nodetsq_ave, nemax, max_npe, 3, "pair:nodetsq_ave");
  }


  if (n_neigh > maxneigh) {
    memory->destroy(scrfcn);
    memory->destroy(dscrfcn);
    memory->destroy(fcpair);
    maxneigh = n_neigh;
    memory->create(scrfcn, maxneigh, "pair:scrfcn");
    memory->create(dscrfcn, maxneigh, "pair:dscrfcn");
    memory->create(fcpair, maxneigh, "pair:fcpair");
  }

  // zero out local arrays

  for (i = 0; i < atom_nall; i++) {
    atomrho0[i] = 0.0;
    atomarho2b[i] = 0.0;
    atomarho1[i][0] = atomarho1[i][1] = atomarho1[i][2] = 0.0;
    for (j = 0; j < 6; j++)
      atomarho2[i][j] = 0.0;
    for (j = 0; j < 10; j++)
      atomarho3[i][j] = 0.0;
    atomarho3b[i][0] = atomarho3b[i][1] = atomarho3b[i][2] = 0.0;
    atomt_ave[i][0] = atomt_ave[i][1] = atomt_ave[i][2] = 0.0;
    atomtsq_ave[i][0] = atomtsq_ave[i][1] = atomtsq_ave[i][2] = 0.0;
  }

  for (i = 0; i < elem_nall; i++) {
    for (j = 0; j < max_npe; j++) {
      noderho0[i][j] = 0.0;
      nodearho2b[i][j] = 0.0;
      nodearho1[i][j][0] = nodearho1[i][j][1] = nodearho1[i][j][2] = 0.0;
      for (k = 0; k < 6; k++)
        nodearho2[i][j][k] = 0.0;
      for (k = 0; k < 10; k++)
        nodearho3[i][j][k] = 0.0;
      nodearho3b[i][j][0] = nodearho3b[i][j][1] = nodearho3b[i][j][2] = 0.0;
      nodet_ave[i][j][0] = nodet_ave[i][j][1] = nodet_ave[i][j][2] = 0.0;
      nodetsq_ave[i][j][0] = nodetsq_ave[i][j][1] = nodetsq_ave[i][j][2] = 0.0;
    }
  }

}

  
void MEAM::meam_dens_init(int i, int iintg, int inode, int *fmap,
    int numneigha_half, int *firstneigha_half, 
    int numneigha_full, int *firstneigha_full, 
    int numneighia_half, int *firstneighia_half, int *firstneighia_index_half, 
    int numneighia_full, int *firstneighia_full, int *firstneighia_index_full, 
    int fnoffset)
{

  //     Compute screening function and derivatives (scrfcn, dscrfcn, and fcpair)
  getscreen(i, iintg, inode, &scrfcn[fnoffset], &dscrfcn[fnoffset], &fcpair[fnoffset], fmap, 
      numneigha_half, firstneigha_half, numneigha_full, firstneigha_full,
      numneighia_half, firstneighia_half, firstneighia_index_half, 
      numneighia_full, firstneighia_full, firstneighia_index_full);

  //     Calculate intermediate density terms 
  //     Only values at nodes are actually tallied
    calc_rho1(i, iintg, inode, fmap, &scrfcn[fnoffset], &fcpair[fnoffset],
        numneigha_half, firstneigha_half, numneigha_full, firstneigha_full,
        numneighia_half, firstneighia_half, firstneighia_index_half, 
        numneighia_full, firstneighia_full, firstneighia_index_full);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


void MEAM::getscreen(int i, int iintg, int inode, double *scrfcn, double *dscrfcn, double *fcpair, int *fmap,
    int numneigha_half, int *firstneigha_half, int numneigha_full, int *firstneigha_full, 
    int numneighia_half, int *firstneighia_half, int *firstneighia_index_half, 
    int numneighia_full, int *firstneighia_full, int *firstneighia_index_full)
{
  int jn, j, kn, k, iintpl, jintpl, kintpl, node;
  int elti, eltj, eltk; // elt = element, index in meam map
  double delxij, delyij, delzij, rij2, rij;
  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp, delxik, delyik, delzik, rik2;
  double xktmp, yktmp, zktmp, delxjk, delyjk, delzjk, rjk2;
  double xik, xjk, sij, fcij, sfcij, dfikj, dfcij, sikj, cikj;
  double Cmin, Cmax, delc, a, coef1, coef2;
  double dCikj;
  double rnorm, dfc, drinv;

  drinv = 1.0 / this->delr_meam;

  tagint *tag = atom->tag;
  tagint itag,jtag;
  double **ax = atom->x;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;
  double ***shape_array = element->shape_array;
  double ***nodex = element->nodex;

  if (iintg < 0) {
    itag = tag[i];
    xitmp = ax[i][0];
    yitmp = ax[i][1];
    zitmp = ax[i][2];
    elti = fmap[atype[i]];
  } else if (inode >= 0){
    xitmp = nodex[i][inode][0];
    yitmp = nodex[i][inode][1];
    zitmp = nodex[i][inode][2];
    elti = fmap[ctype[i]];
  } else {
    int ietype = etype[i];
    iintpl = element->i2ia[ietype][iintg];
    xitmp = yitmp = zitmp = 0.0;
    for (int node = 0; node < npe[ietype]; node++) {
      xitmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
      yitmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
      zitmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
    }
    elti = fmap[ctype[i]];

  }

  // atom neighbor list

    int counter=0;
  for (jn = 0; jn < numneigha_half + numneighia_half; jn++) {
    if (jn < numneigha_half) {
      j = firstneigha_half[jn];
      eltj = fmap[atype[j]];
      if (eltj < 0) continue;

      // check i j pair to count only half

      jtag = tag[j]; 

      xjtmp = ax[j][0];
      yjtmp = ax[j][1];
      zjtmp = ax[j][2];

    } else {
      j = firstneighia_half[jn - numneigha_half];
      jintpl = firstneighia_index_half[jn - numneigha_half];
      eltj = fmap[ctype[j]];
      if (eltj < 0) continue;

      xjtmp = yjtmp = zjtmp = 0.0;
      for (node = 0; node < npe[etype[j]]; node++) {
        xjtmp += shape_array[etype[j]][jintpl][node] * nodex[j][node][0];
        yjtmp += shape_array[etype[j]][jintpl][node] * nodex[j][node][1];
        zjtmp += shape_array[etype[j]][jintpl][node] * nodex[j][node][2];
      }
    }

    //     First compute screening function itself, sij

    delxij = xjtmp - xitmp;
    delyij = yjtmp - yitmp;
    delzij = zjtmp - zitmp;
    rij2 = delxij * delxij + delyij * delyij + delzij * delzij;

    if (rij2 >= this->cutforcesq) {
      dscrfcn[jn] = 0.0;
      scrfcn[jn] = 0.0;
      fcpair[jn] = 0.0;
      continue;
    }

    const double rbound = this->ebound_meam[elti][eltj] * rij2;
    rij = sqrt(rij2);
    rnorm = (this->cutforce - rij) * drinv;
    sij = 1.0;
    dscrfcn[jn] = 0.0;

    //     if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
    //     merge sij and dscrfcn calculation into same loop to save computing as compared to original code
    for (kn = 0; kn < numneigha_full + numneighia_full; kn++) {
      if (kn < numneigha_full) {
        k = firstneigha_full[kn];
        if (k == j && jn < numneigha_half) continue;
        eltk = fmap[atype[k]];
        if (eltk < 0) continue;
        xktmp = ax[k][0];
        yktmp = ax[k][1];
        zktmp = ax[k][2];

      } else {
        k = firstneighia_full[kn - numneigha_full];
        kintpl = firstneighia_index_full[kn - numneigha_full];
        if (jn >= numneigha_half && k == j && kintpl == jintpl) continue; 
        eltk = fmap[ctype[k]];
        if (eltk < 0) continue;
        xktmp = yktmp = zktmp = 0.0;
        for (node = 0; node < npe[etype[k]]; node++) {
          xktmp += shape_array[etype[k]][kintpl][node]*nodex[k][node][0];
          yktmp += shape_array[etype[k]][kintpl][node]*nodex[k][node][1];
          zktmp += shape_array[etype[k]][kintpl][node]*nodex[k][node][2];
        }
      }
counter++;
      delxjk = xktmp - xjtmp;
      delyjk = yktmp - yjtmp;
      delzjk = zktmp - zjtmp;
      rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
      if (rjk2 > rbound) continue;

      delxik = xktmp - xitmp;
      delyik = yktmp - yitmp;
      delzik = zktmp - zitmp;
      rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
      if (rik2 > rbound) continue;

      xik = rik2 / rij2;
      xjk = rjk2 / rij2;
      a = 1 - (xik - xjk) * (xik - xjk);
      //     if a <= 0, then ellipse equation doesn't describe this case and
      //     atom k can't possibly screen i-j
      if (a <= 0.0) continue;

      cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
      Cmax = this->Cmax_meam[elti][eltj][eltk];
      Cmin = this->Cmin_meam[elti][eltj][eltk];
      if (cikj >= Cmax) continue;
      //     note that cikj may be slightly negative (within numerical
      //     tolerance) if atoms are colinear, so don't reject that case here
      //     (other negative cikj cases were handled by the test on "a" above)
      else if (cikj <= Cmin) {
        sij = 0;
        dscrfcn[jn] = 0.0;
        break;
      } 

      delc = Cmax - Cmin;
      cikj = (cikj - Cmin) / delc;
      sikj = dfcut(cikj, dfikj);
      coef1 = dfikj / (delc * sikj);
      dCikj = dCfunc(rij2, rik2, rjk2);
      sij *= fcut(cikj);
      dscrfcn[jn] += coef1 * dCikj;

      //if (debug) error->all(FLERR,"TEST"); 

/*
      int flag = compute_screen(xitmp,yitmp,zitmp,xjtmp,yjtmp,zjtmp,xktmp,yktmp,zktmp,elti,eltj,eltk,rbound,rij2,sikj,coef1,dCikj);
      if (flag == 1) {

        //        if (iintg>=0)
        //          if ((element->tag[i]-1)*element->max_nintpl+iintpl+1 == 4563) 
        //            if ((element->tag[j]-1)*element->max_nintpl+jintpl+1== 5651) {
        //              printf("sikj = %-1.16e k = %d\n",sikj,(element->tag[k]-1)*element->max_nintpl+kintpl+1);
        //            }

        sij *= sikj;
        dscrfcn[jn] += coef1 * dCikj;
      } else if (flag == 2) {
        sij = 0;
        dscrfcn[jn] = 0.0;
        break;
      }
      */
    }
    fcij = dfcut(rnorm, dfc);
    scrfcn[jn] = sij;
    fcpair[jn] = fcij;

    sfcij = sij * fcij;
    if (!iszero(sfcij) && !isone(sfcij)) {
      coef1 = sfcij;
      dfcij = dfc * drinv;
      coef2 = sij * dfcij / rij;
      dscrfcn[jn] = dscrfcn[jn] * coef1 - coef2;
    } else dscrfcn[jn] = 0.0;

    //if (!iszero(dscrfcn[jn])) printf("%d %g %g %g\n",i,rij,dscrfcn[jn],sij);

    //if ((atom->tag[i] == 180 && atom->tag[j] == 120) ||
    //    (atom->tag[j] == 180 && atom->tag[i] == 120))
    //printf("i = %d j = %d sij = %g fcpair = %g dscrfcn = %g n = %d\n",std::max(itag,jtag),std::min(itag,jtag),sij,fcij,dscrfcn[n],n);
  }


    //printf("%d\n",counter);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


void MEAM::calc_rho1(int i, int iintg, int inode, int *fmap, double *scrfcn, double *fcpair,
    int numneigha_half, int *firstneigha_half, int numneigha_full, int *firstneigha_full, 
    int numneighia_half, int *firstneighia_half, int *firstneighia_index_half, 
    int numneighia_full, int *firstneighia_full, int *firstneighia_index_full)
{
  int ietype, iintpl;
  int jn, j, jnode, jintg, jintpl, jetype, m, n, p, elti, eltj, node;
  int nv2, nv3;
  double delij[3], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;

  // double G,Gbar,gam,shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;
  double xitmp,yitmp,zitmp,xjtmp,yjtmp,zjtmp;

  double *rho0i,*rho0j,*t_avei,*t_avej,*tsq_avei,*tsq_avej,*arho2bi,*arho2bj;
  double *arho1i,*arho1j,*arho3bi,*arho3bj,*arho2i,*arho2j,*arho3i,*arho3j;

  tagint *tag = atom->tag;
  tagint itag,jtag;
  int iflag,jflag;

  double **ax = atom->x;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int **i2ia = element->i2ia;
  double ***shape_array = element->shape_array;
  double ***nodex = element->nodex;

  // atom neighbor list 

  if (iintg < 0) {
    xitmp = ax[i][0];
    yitmp = ax[i][1];
    zitmp = ax[i][2];
    elti = fmap[atype[i]];
    itag = atom->tag[i];
    rho0i = &atomrho0[i];
    t_avei = atomt_ave[i];
    tsq_avei = atomtsq_ave[i];
    arho2bi = &atomarho2b[i];
    arho1i = atomarho1[i];
    arho3bi = atomarho3b[i];
    arho2i = atomarho2[i];
    arho3i = atomarho3[i];
    iflag = 1;
  } else if (inode < 0) {
    ietype = etype[i];
    iintpl = i2ia[ietype][iintg];

    xitmp = yitmp = zitmp = 0.0;
    for (node = 0; node < npe[etype[i]]; node++) {
      xitmp += shape_array[ietype][iintpl][node] * nodex[i][node][0];
      yitmp += shape_array[ietype][iintpl][node] * nodex[i][node][1];
      zitmp += shape_array[ietype][iintpl][node] * nodex[i][node][2];
    }
    elti = fmap[ctype[i]];
    iflag = 0;
  } else {
    xitmp = nodex[i][inode][0];
    yitmp = nodex[i][inode][1];
    zitmp = nodex[i][inode][2];
    elti = fmap[ctype[i]];
    rho0i = &noderho0[i][inode];
    t_avei = nodet_ave[i][inode];
    tsq_avei = nodetsq_ave[i][inode];
    arho2bi = &nodearho2b[i][inode];
    arho1i = nodearho1[i][inode];
    arho3bi = nodearho3b[i][inode];
    arho2i = nodearho2[i][inode];
    arho3i = nodearho3[i][inode];
    iflag = 1;
  }

  for (jn = 0; jn < numneigha_half + numneighia_half; jn++) {
    if (jn < numneigha_half) {
      j = firstneigha_half[jn];

      xjtmp = ax[j][0];
      yjtmp = ax[j][1];
      zjtmp = ax[j][2];

      eltj = fmap[atype[j]];
      if (eltj < 0) continue;
      rho0j = &atomrho0[j];
      t_avej = atomt_ave[j];
      tsq_avej = atomtsq_ave[j];
      arho2bj = &atomarho2b[j];
      arho1j = atomarho1[j];
      arho3bj = atomarho3b[j];
      arho2j = atomarho2[j];
      arho3j = atomarho3[j];
      jflag = 1;
    } else {

      j = firstneighia_half[jn - numneigha_half];
      jintpl = firstneighia_index_half[jn - numneigha_half];
      jetype = etype[j];
      jintg = ia2i[jetype][jintpl];
      if (jintg >= 0) jnode = i2n[jetype][jintg];
      else jnode = -1;

      if (jnode < 0) {
        xjtmp = yjtmp = zjtmp = 0.0;
        for (node = 0; node < npe[jetype]; node++) {
          xjtmp += shape_array[jetype][jintpl][node] * nodex[j][node][0];
          yjtmp += shape_array[jetype][jintpl][node] * nodex[j][node][1];
          zjtmp += shape_array[jetype][jintpl][node] * nodex[j][node][2];
        }
        jflag = 0;
      } else {
        xjtmp = nodex[j][jnode][0];
        yjtmp = nodex[j][jnode][1];
        zjtmp = nodex[j][jnode][2];
        rho0j = &noderho0[j][jnode];
        t_avej = nodet_ave[j][jnode];
        tsq_avej = nodetsq_ave[j][jnode];
        arho2bj = &nodearho2b[j][jnode];
        arho1j = nodearho1[j][jnode];
        arho3bj = nodearho3b[j][jnode];
        arho2j = nodearho2[j][jnode];
        arho3j = nodearho3[j][jnode];
        jflag = 1;
      }
      eltj = fmap[ctype[j]];
      if (eltj < 0) continue;
    }

    if (!iflag && !jflag) continue;
    if (iszero(scrfcn[jn])) continue;

    sij = scrfcn[jn] * fcpair[jn];
    delij[0] = xjtmp - xitmp;
    delij[1] = yjtmp - yitmp;
    delij[2] = zjtmp - zitmp;
    rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
    if (rij2 < this->cutforcesq) {
      rij = sqrt(rij2);
      if (iflag) {
        aj = rij / this->re_meam[eltj][eltj] - 1.0;
        ro0j = this->rho0_meam[eltj];
        rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj) * sij;
        rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj) * sij;
        rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj) * sij;
        rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj) * sij;
      }

      if (jflag) {
        ai = rij / this->re_meam[elti][elti] - 1.0;

        ro0i = this->rho0_meam[elti];

        rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai) * sij;
        rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai) * sij;
        rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai) * sij;
        rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai) * sij;
      }

      if (this->ialloy == 1) {
        if (iflag) {
          rhoa1j *= this->t1_meam[eltj];
          rhoa2j *= this->t2_meam[eltj];
          rhoa3j *= this->t3_meam[eltj];
        }
        if (jflag) {
          rhoa1i *= this->t1_meam[elti];
          rhoa2i *= this->t2_meam[elti];
          rhoa3i *= this->t3_meam[elti];
        }
      }

      // For ialloy = 2, use single-element value (not average)
      if (this->ialloy != 2) {
        if (iflag) {
          t_avei[0] += this->t1_meam[eltj] * rhoa0j;
          t_avei[1] += this->t2_meam[eltj] * rhoa0j;
          t_avei[2] += this->t3_meam[eltj] * rhoa0j;
        }
        if (jflag) {
          t_avej[0] += this->t1_meam[elti] * rhoa0i;
          t_avej[1] += this->t2_meam[elti] * rhoa0i;
          t_avej[2] += this->t3_meam[elti] * rhoa0i;
        }
      }
      if (this->ialloy == 1) {
        if (iflag) {
          tsq_avei[0] += this->t1_meam[eltj] * this->t1_meam[eltj] * rhoa0j;
          tsq_avei[1] += this->t2_meam[eltj] * this->t2_meam[eltj] * rhoa0j;
          tsq_avei[2] += this->t3_meam[eltj] * this->t3_meam[eltj] * rhoa0j;
        }
        if (jflag) {
          tsq_avej[0] += this->t1_meam[elti] * this->t1_meam[elti] * rhoa0i;
          tsq_avej[1] += this->t2_meam[elti] * this->t2_meam[elti] * rhoa0i;
          tsq_avej[2] += this->t3_meam[elti] * this->t3_meam[elti] * rhoa0i;
        }
      }

      if (iflag) {
        *rho0i += rhoa0j;
        *arho2bi += rhoa2j;
        A1j = rhoa1j / rij;
        A2j = rhoa2j / rij2;
        A3j = rhoa3j / (rij2 * rij);
      }
      if (jflag) {
        *rho0j += rhoa0i;
        *arho2bj += rhoa2i;
        A1i = rhoa1i / rij;
        A2i = rhoa2i / rij2;
        A3i = rhoa3i / (rij2 * rij);
      }

      nv2 = 0;
      nv3 = 0;
      for (m = 0; m < 3; m++) {
        if (iflag) {
          arho1i[m] += A1j * delij[m];
          arho3bi[m] += rhoa3j * delij[m] / rij;
        }
        if (jflag) {
          arho1j[m] -= A1i * delij[m];
          arho3bj[m] -= rhoa3i * delij[m] / rij;
        }
        for (n = m; n < 3; n++) {
          if (iflag) arho2i[nv2] += A2j * delij[m] * delij[n];
          if (jflag) arho2j[nv2] += A2i * delij[m] * delij[n];
          nv2++; 
          for (p = n; p < 3; p++) {
            if (iflag) arho3i[nv3] += A3j * delij[m] * delij[n] * delij[p];
            if (jflag) arho3j[nv3] -= A3i * delij[m] * delij[n] * delij[p];
            nv3++;
          }
        }
      }
    }
  }
}

int MEAM::compute_screen(double xi, double yi, double zi, double xj, double yj, double zj, 
    double xk, double yk, double zk, int elti, int eltj, int eltk, double rbound, double rij2, double &sikj, double &coef1, double &dCikj)
{
  //if (debug)
  //  printf("xi = %-1.16e %-1.16e %-1.16e xj = %-1.16e %-1.16e %-1.16e xk = %-1.16e %-1.16e %-1.16e rij2 = %d\n",xi,yi,zi,xj,yj,zj,xk,yk,zk,rij2);
  double delc, cikj, dfikj;
  double delxjk = xk - xj;
  double delyjk = yk - yj;
  double delzjk = zk - zj;
  double rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
  if (rjk2 > rbound) return 0;

  double delxik = xk - xi;
  double delyik = yk - yi;
  double delzik = zk - zi;
  double rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
  if (rik2 > rbound) return 0;

  double xik = rik2 / rij2;
  double xjk = rjk2 / rij2;
  double a = 1 - (xik - xjk) * (xik - xjk);
  //     if a <= 0, then ellipse equation doesn't describe this case and
  //     atom k can't possibly screen i-j
  if (a <= 0.0) return 0;

  cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
  double Cmax = this->Cmax_meam[elti][eltj][eltk];
  double Cmin = this->Cmin_meam[elti][eltj][eltk];
  if (cikj >= Cmax) return 0;
  //     note that cikj may be slightly negative (within numerical
  //     tolerance) if atoms are colinear, so don't reject that case here
  //     (other negative cikj cases were handled by the test on "a" above)
  else if (cikj <= Cmin) return 2; 

  delc = Cmax - Cmin;
  cikj = (cikj - Cmin) / delc;
  sikj = dfcut(cikj, dfikj);
  coef1 = dfikj / (delc * sikj);
  dCikj = dCfunc(rij2, rik2, rjk2);
  sikj = fcut(cikj);
  //if (debug) error->all(FLERR,"TEST"); 
  return 1;
}

