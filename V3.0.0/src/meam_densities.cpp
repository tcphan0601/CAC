#include "meam.h"
#include <cmath>
#include "memory.h"
#include "math_special.h"
#include "atom.h"
#include "element.h"
#include "error.h"

using namespace CAC_NS;

void
MEAM::meam_dens_setup(int atom_nmax, int elem_nmax, int atom_nall, int elem_nall, int n_neigh)
{
  int i, j, k, l;

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
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;
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

    memory->create(noderho, nemax, maxapc, maxnpe, "pair:noderho");
    memory->create(noderho0, nemax, maxapc, maxnpe, "pair:noderho0");
    memory->create(noderho1, nemax, maxapc, maxnpe, "pair:noderho1");
    memory->create(noderho2, nemax, maxapc, maxnpe, "pair:noderho2");
    memory->create(noderho3, nemax, maxapc, maxnpe, "pair:noderho3");
    memory->create(nodefrhop, nemax, maxapc, maxnpe, "pair:nodefrhop");
    memory->create(nodegamma, nemax, maxapc, maxnpe, "pair:nodegamma");
    memory->create(nodedgamma1, nemax, maxapc, maxnpe, "pair:nodedgamma1");
    memory->create(nodedgamma2, nemax, maxapc, maxnpe, "pair:nodedgamma2");
    memory->create(nodedgamma3, nemax, maxapc, maxnpe, "pair:nodedgamma3");
    memory->create(nodearho2b, nemax, maxapc, maxnpe, "pair:nodearho2b");
    memory->create(nodearho1, nemax, maxapc, maxnpe, 3, "pair:nodearho1");
    memory->create(nodearho2, nemax, maxapc, maxnpe, 6, "pair:nodearho2");
    memory->create(nodearho3, nemax, maxapc, maxnpe, 10, "pair:nodearho3");
    memory->create(nodearho3b, nemax, maxapc, maxnpe, 3, "pair:nodearho3b");
    memory->create(nodet_ave, nemax, maxapc, maxnpe, 3, "pair:nodet_ave");
    memory->create(nodetsq_ave, nemax, maxapc, maxnpe, 3, "pair:nodetsq_ave");
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

  for (i = 0; i < elem_nall; i++) 
    for (j = 0; j < maxapc; j++) 
      for (k = 0; k < maxnpe; k++) {
        noderho0[i][j][k] = 0.0;
        nodearho2b[i][j][k] = 0.0;
        for (l = 0; l < 3; l++) {
          nodearho1[i][j][k][l] = 0.0;
          nodearho3b[i][j][k][l] = 0.0;
          nodet_ave[i][j][k][l] = 0.0;
          nodetsq_ave[i][j][k][l] = 0.0;
        }
        for (l = 0; l < 6; l++)
          nodearho2[i][j][k][l] = 0.0;
        for (l = 0; l < 10; l++)
          nodearho3[i][j][k][l] = 0.0;
      }
}


void MEAM::meam_dens_init(int i, int iindex, int *fmap, int fnoffset, 
    int numneigh, int *firstneigh, int *firstneighindex, 
    int numneigh_full, int *firstneigh_full, int *firstneighindex_full)
{

  //     Compute screening function and derivatives (scrfcn, dscrfcn, and fcpair)
  getscreen(i, iindex, &scrfcn[fnoffset], &dscrfcn[fnoffset], &fcpair[fnoffset], fmap, 
      numneigh, firstneigh, firstneighindex, 
      numneigh_full, firstneigh_full, firstneighindex_full);

  //     Calculate intermediate density terms 
  //     Only values at nodes are actually calculated
  calc_rho1(i, iindex, fmap, &scrfcn[fnoffset], &fcpair[fnoffset], 
      numneigh, firstneigh, firstneighindex);
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


void MEAM::getscreen(int i, int iindex, double *scrfcn, double *dscrfcn, double *fcpair, int *fmap, 
    int numneigh, int *firstneigh, int *firstneighindex, 
    int numneigh_full, int *firstneigh_full, int *firstneighindex_full)
{
  int jn, j, kn, k, node;
  int jindex, jucell, jetype, japc, jnpe, jbasis;
  int kindex, kucell, ketype, kapc, knpe, kbasis;
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

  double **ax = atom->x;
  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;
  int **g2u = element->g2u;
  int **g2n = element->g2n;
  double ***shape_array = element->shape_array;
  double ****nodex = element->nodex;

  // i is atom

  if (iindex < 0) {
    xitmp = ax[i][0];
    yitmp = ax[i][1];
    zitmp = ax[i][2];
    elti = fmap[atype[i]];
  } 

  // i is gauss point

  else {
    int ietype = etype[i];
    int iapc = apc[ietype];
    int inpe = npe[ietype];
    int igcell = iindex / iapc;
    int ibasis = iindex % iapc;
    elti = fmap[ctype[i][ibasis]];
    int inode = g2n[ietype][igcell];
    if (inode >= 0){
      xitmp = nodex[i][ibasis][inode][0];
      yitmp = nodex[i][ibasis][inode][1];
      zitmp = nodex[i][ibasis][inode][2];
    } else {
      int iucell = element->g2u[ietype][igcell];
      xitmp = yitmp = zitmp = 0.0;
      for (node = 0; node < inpe; node++) {
        xitmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][0];
        yitmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][1];
        zitmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][2];
      }
    }
  }

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    jindex = firstneighindex[jn];

    // j is atom

    if (jindex < 0) {
      eltj = fmap[atype[j]];
      if (eltj < 0) continue;
      xjtmp = ax[j][0];
      yjtmp = ax[j][1];
      zjtmp = ax[j][2];
    } 

    // j is virtual atom

    else {
      jetype = etype[j];
      japc = apc[jetype];
      jnpe = npe[jetype];
      jucell = jindex / japc;
      jbasis = jindex % japc;
      eltj = fmap[ctype[j][jbasis]];
      if (eltj < 0) continue;

      xjtmp = yjtmp = zjtmp = 0.0;
      for (node = 0; node < jnpe; node++) {
        xjtmp += shape_array[jetype][jucell][node] * 
          nodex[j][jbasis][node][0];
        yjtmp += shape_array[jetype][jucell][node] * 
          nodex[j][jbasis][node][1];
        zjtmp += shape_array[jetype][jucell][node] * 
          nodex[j][jbasis][node][2];
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

    //     if rjk2 > ebound * rijsq, atom k is definitely outside the ellipse
    //     merge sij and dscrfcn calculation into same loop to save computing as compared to original code

    for (kn = 0; kn < numneigh_full; kn++) {
      k = firstneigh_full[kn];
      kindex = firstneighindex_full[kn];

      // k is atom

      if (kindex < 0) {
        eltk = fmap[atype[k]];
        if (eltk < 0) continue;
        xktmp = ax[k][0];
        yktmp = ax[k][1];
        zktmp = ax[k][2];
      } 

      // k is virtual atom

      else {
        ketype = etype[k];
        kapc = apc[ketype];
        knpe = npe[ketype];
        kucell = kindex / kapc;
        kbasis = kindex % kapc;
        eltk = fmap[ctype[k][kbasis]];
        if (eltk < 0) continue;
        xktmp = yktmp = zktmp = 0.0;
        for (node = 0; node < knpe; node++) {
          xktmp += shape_array[ketype][kucell][node] * 
            nodex[k][kbasis][node][0];
          yktmp += shape_array[ketype][kucell][node] * 
            nodex[k][kbasis][node][1];
          zktmp += shape_array[ketype][kucell][node] * 
            nodex[k][kbasis][node][2];
        }
      }
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

    //if (!iszero(dscrfcn[jn])) printf("%d %g %g %g\n", i, rij, dscrfcn[jn], sij);

    //if ((atom->tag[i] == 180 && atom->tag[j] == 120) ||
    //    (atom->tag[j] == 180 && atom->tag[i] == 120))
    //printf("i = %d j = %d sij = %g fcpair = %g dscrfcn = %g n = %d\n", std::max(itag, jtag), std::min(itag, jtag), sij, fcij, dscrfcn[n], n);
  }


}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


void MEAM::calc_rho1(int i, int iindex, int *fmap, double *scrfcn, double *fcpair, 
    int numneigh, int *firstneigh, int *firstneighindex)
{
  int ietype, iucell, inode, igcell, iapc, inpe, ibasis;
  int jetype, jucell, jnode, jgcell, japc, jnpe, jbasis, jindex;
  int jn, j, m, n, p, elti, eltj, node;
  int nv2, nv3;
  double delij[3], rij2, rij, sij;
  double ai, aj, rhoa0j, rhoa1j, rhoa2j, rhoa3j, A1j, A2j, A3j;

  // double G, Gbar, gam, shp[3+1];
  double ro0i, ro0j;
  double rhoa0i, rhoa1i, rhoa2i, rhoa3i, A1i, A2i, A3i;
  double xtmp, ytmp, ztmp;

  double *rho0i, *rho0j, *t_avei, *t_avej, *tsq_avei, *tsq_avej, *arho2bi, *arho2bj;
  double *arho1i, *arho1j, *arho3bi, *arho3bj, *arho2i, *arho2j, *arho3i, *arho3j;

  int iflag, jflag;

  double **ax = atom->x;
  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  int **g2u = element->g2u;
  double ***shape_array = element->shape_array;
  double ****nodex = element->nodex;

  // i is atom

  if (iindex < 0) {
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];
    elti = fmap[atype[i]];
    rho0i = &atomrho0[i];
    t_avei = atomt_ave[i];
    tsq_avei = atomtsq_ave[i];
    arho2bi = &atomarho2b[i];
    arho1i = atomarho1[i];
    arho3bi = atomarho3b[i];
    arho2i = atomarho2[i];
    arho3i = atomarho3[i];
    iflag = 1;
  } 

  // i is gauss point

  else {
    ietype = etype[i];
    iapc = apc[ietype];
    inpe = npe[ietype];
    igcell = iindex / iapc;
    ibasis = iindex % iapc;
    elti = fmap[ctype[i][ibasis]];
    inode = g2n[ietype][igcell];
    if (inode < 0) {
      iucell = element->g2u[ietype][igcell];
      xtmp = ytmp = ztmp = 0.0;
      for (node = 0; node < inpe; node++) {
        xtmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][0];
        ytmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][1];
        ztmp += shape_array[ietype][iucell][node] *
          nodex[i][ibasis][node][2];
      }
      iflag = 0;
    } else {
      xtmp = nodex[i][ibasis][inode][0];
      ytmp = nodex[i][ibasis][inode][1];
      ztmp = nodex[i][ibasis][inode][2];
      rho0i = &noderho0[i][ibasis][inode];
      t_avei = nodet_ave[i][ibasis][inode];
      tsq_avei = nodetsq_ave[i][ibasis][inode];
      arho2bi = &nodearho2b[i][ibasis][inode];
      arho1i = nodearho1[i][ibasis][inode];
      arho3bi = nodearho3b[i][ibasis][inode];
      arho2i = nodearho2[i][ibasis][inode];
      arho3i = nodearho3[i][ibasis][inode];
      iflag = 1;
    }
  }

  for (jn = 0; jn < numneigh; jn++) {
    j = firstneigh[jn];
    jindex = firstneighindex[jn];
    delij[0] = -xtmp;
    delij[1] = -ytmp;
    delij[2] = -ztmp;

    // j is atom

    if (jindex < 0) {
      delij[0] += ax[j][0];
      delij[1] += ax[j][1];
      delij[2] += ax[j][2];

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
    } 

    // j is virtual atom

    else {
      jetype = etype[j];
      japc = apc[jetype];
      jnpe = npe[jetype];
      jucell = jindex / japc;
      jbasis = jindex % japc;
      eltj = fmap[ctype[j][jbasis]];
      if (eltj < 0) continue;

      jgcell = u2g[jetype][jucell];
      if (jgcell >= 0) jnode = g2n[jetype][jgcell];
      else jnode = -1;

      if (jnode < 0) {
        for (node = 0; node < jnpe; node++) {
          delij[0] += shape_array[jetype][jucell][node] * 
            nodex[j][jbasis][node][0];
          delij[1] += shape_array[jetype][jucell][node] * 
            nodex[j][jbasis][node][1];
          delij[2] += shape_array[jetype][jucell][node] * 
            nodex[j][jbasis][node][2];
        }
        jflag = 0;
      } else {
        delij[0] += nodex[j][jbasis][jnode][0];
        delij[1] += nodex[j][jbasis][jnode][1];
        delij[2] += nodex[j][jbasis][jnode][2];
        rho0j = &noderho0[j][jbasis][jnode];
        t_avej = nodet_ave[j][jbasis][jnode];
        tsq_avej = nodetsq_ave[j][jbasis][jnode];
        arho2bj = &nodearho2b[j][jbasis][jnode];
        arho1j = nodearho1[j][jbasis][jnode];
        arho3bj = nodearho3b[j][jbasis][jnode];
        arho2j = nodearho2[j][jbasis][jnode];
        arho3j = nodearho3[j][jbasis][jnode];
        jflag = 1;
      }
    }

    if (!iflag && !jflag) continue;
    if (iszero(scrfcn[jn])) continue;

    sij = scrfcn[jn] * fcpair[jn];
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

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int MEAM::compute_screen(double xi, double yi, double zi, double xj, double yj, double zj, 
    double xk, double yk, double zk, int elti, int eltj, int eltk, double rbound, double rij2, double &sikj, double &coef1, double &dCikj)
{
  //if (debug)
  //  printf("xi = %-1.16e %-1.16e %-1.16e xj = %-1.16e %-1.16e %-1.16e xk = %-1.16e %-1.16e %-1.16e rij2 = %d\n", xi, yi, zi, xj, yj, zj, xk, yk, zk, rij2);
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
  //if (debug) error->all(FLERR, "TEST"); 
  return 1;
}

// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

void MEAM::meam_dens_final(int eflag_either, int eflag_global, int eflag_atom, double *eng_vdwl, 
    double *eatom, double ***enode, int *fmap, double **scale, int &errorflag)
{
  int i, elti;
  int m;
  double rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  double denom, rho_bkgd, Fl;
  double scaleii;
  int inode, ibasis, ietype, ictype;
  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;

  //     Complete the calculation of density for atoms

  for (i = 0; i < atom->nlocal; i++) {
    elti = fmap[atype[i]];
    if (elti >= 0) {
      scaleii = scale[atype[i]][atype[i]];
      atomrho1[i] = 0.0;
      atomrho2[i] = -1.0 / 3.0 * atomarho2b[i] * atomarho2b[i];
      atomrho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        atomrho1[i] += atomarho1[i][m] * atomarho1[i][m];
        atomrho3[i] -= 3.0 / 5.0 * atomarho3b[i][m] * atomarho3b[i][m];
      }
      for (m = 0; m < 6; m++) {
        atomrho2[i] += this->v2D[m] * atomarho2[i][m] * atomarho2[i][m];
      }
      for (m = 0; m < 10; m++) {
        atomrho3[i] += this->v3D[m] * atomarho3[i][m] * atomarho3[i][m];
      }

      if (atomrho0[i] > 0.0) {
        if (this->ialloy == 1) {
          atomt_ave[i][0] = fdiv_zero(atomt_ave[i][0], atomtsq_ave[i][0]);
          atomt_ave[i][1] = fdiv_zero(atomt_ave[i][1], atomtsq_ave[i][1]);
          atomt_ave[i][2] = fdiv_zero(atomt_ave[i][2], atomtsq_ave[i][2]);
        } else if (this->ialloy == 2) {
          atomt_ave[i][0] = this->t1_meam[elti];
          atomt_ave[i][1] = this->t2_meam[elti];
          atomt_ave[i][2] = this->t3_meam[elti];
        } else {
          atomt_ave[i][0] /= atomrho0[i];
          atomt_ave[i][1] /= atomrho0[i];
          atomt_ave[i][2] /= atomrho0[i];
        }
      }

      atomgamma[i] = atomt_ave[i][0] * atomrho1[i] + atomt_ave[i][1] * atomrho2[i] + atomt_ave[i][2] * atomrho3[i];

      if (atomrho0[i] > 0.0) {
        atomgamma[i] /= (atomrho0[i] * atomrho0[i]);
      }

      Z = get_Zij(this->lattce_meam[elti][elti]);

      G = G_gam(atomgamma[i], this->ibar_meam[elti], errorflag);
      if (errorflag != 0)
        return;

      get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);

      if (this->ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        if (this->mix_ref_t == 1) {
          gam = (atomt_ave[i][0] * shp[0] + atomt_ave[i][1] * shp[1] + atomt_ave[i][2] * shp[2]) / (Z * Z);
        } else {
          gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
            (Z * Z);
        }
        Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
      }
      atomrho[i] = atomrho0[i] * G;

      if (this->mix_ref_t == 1) {
        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          gam = (atomt_ave[i][0] * shp[0] + atomt_ave[i][1] * shp[1] + atomt_ave[i][2] * shp[2]) / (Z * Z);
          Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
        }
        rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
      } else {
        if (this->bkgd_dyn == 1) {
          rho_bkgd = this->rho0_meam[elti] * Z;
        } else {
          rho_bkgd = this->rho_ref_meam[elti];
        }
      }
      rhob = atomrho[i] / rho_bkgd;
      denom = 1.0 / rho_bkgd;

      G = dG_gam(atomgamma[i], this->ibar_meam[elti], dG);

      atomdgamma1[i] = (G - 2 * dG * atomgamma[i]) * denom;

      if (!iszero(atomrho0[i])) {
        atomdgamma2[i] = (dG / atomrho0[i]) * denom;
      } else {
        atomdgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (this->mix_ref_t == 1) {
        atomdgamma3[i] = atomrho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        atomdgamma3[i] = 0.0;
      }

      Fl = embedding(this->A_meam[elti], this->Ec_meam[elti][elti], rhob, atomfrhop[i]);

      if (eflag_either != 0) {
        Fl *= scaleii;
        if (eflag_global != 0) {
          *eng_vdwl = *eng_vdwl + Fl;
        }
        if (eflag_atom != 0) {
          eatom[i] = eatom[i] + Fl;
        }
      }
    }
  }

  //     Complete the calculation of density for nodes

  for (i = 0; i < element->nlocal; i++) {
    ietype = etype[i];
    for (ibasis = 0; ibasis < apc[ietype]; ibasis++) {
      ictype = ctype[i][ibasis];
      elti = fmap[ictype];
      if (elti < 0) continue;
      for (inode = 0; inode < npe[ietype]; inode++) {
        scaleii = scale[ictype][ictype];
        noderho1[i][ibasis][inode] = 0.0;
        noderho2[i][ibasis][inode] = -1.0 / 3.0 * nodearho2b[i][ibasis][inode] * nodearho2b[i][ibasis][inode];
        noderho3[i][ibasis][inode] = 0.0;
        for (m = 0; m < 3; m++) {
          noderho1[i][ibasis][inode] += nodearho1[i][ibasis][inode][m] * nodearho1[i][ibasis][inode][m];
          noderho3[i][ibasis][inode] -= 3.0 / 5.0 * nodearho3b[i][ibasis][inode][m] * nodearho3b[i][ibasis][inode][m];
        }
        for (m = 0; m < 6; m++) {
          noderho2[i][ibasis][inode] += this->v2D[m] * nodearho2[i][ibasis][inode][m] * nodearho2[i][ibasis][inode][m];
        }
        for (m = 0; m < 10; m++) {
          noderho3[i][ibasis][inode] += this->v3D[m] * nodearho3[i][ibasis][inode][m] * nodearho3[i][ibasis][inode][m];
        }

        if (noderho0[i][ibasis][inode] > 0.0) {
          if (this->ialloy == 1) {
            nodet_ave[i][ibasis][inode][0] = fdiv_zero(nodet_ave[i][ibasis][inode][0], nodetsq_ave[i][ibasis][inode][0]);
            nodet_ave[i][ibasis][inode][1] = fdiv_zero(nodet_ave[i][ibasis][inode][1], nodetsq_ave[i][ibasis][inode][1]);
            nodet_ave[i][ibasis][inode][2] = fdiv_zero(nodet_ave[i][ibasis][inode][2], nodetsq_ave[i][ibasis][inode][2]);
          } else if (this->ialloy == 2) {
            nodet_ave[i][ibasis][inode][0] = this->t1_meam[elti];
            nodet_ave[i][ibasis][inode][1] = this->t2_meam[elti];
            nodet_ave[i][ibasis][inode][2] = this->t3_meam[elti];
          } else {
            nodet_ave[i][ibasis][inode][0] /= noderho0[i][ibasis][inode];
            nodet_ave[i][ibasis][inode][1] /= noderho0[i][ibasis][inode];
            nodet_ave[i][ibasis][inode][2] /= noderho0[i][ibasis][inode];
          }
        }

        nodegamma[i][ibasis][inode] = nodet_ave[i][ibasis][inode][0] * noderho1[i][ibasis][inode] 
          + nodet_ave[i][ibasis][inode][1] * noderho2[i][ibasis][inode] 
          + nodet_ave[i][ibasis][inode][2] * noderho3[i][ibasis][inode];

        if (noderho0[i][ibasis][inode] > 0.0) {
          nodegamma[i][ibasis][inode] /= (noderho0[i][ibasis][inode] * noderho0[i][ibasis][inode]);
        }

        Z = get_Zij(this->lattce_meam[elti][elti]);

        G = G_gam(nodegamma[i][ibasis][inode], this->ibar_meam[elti], errorflag);
        if (errorflag != 0)
          return;

        get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);

        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          if (this->mix_ref_t == 1) {
            gam = (nodet_ave[i][ibasis][inode][0] * shp[0] + nodet_ave[i][ibasis][inode][1] * shp[1] + nodet_ave[i][ibasis][inode][2] * shp[2]) / (Z * Z);
          } else {
            gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
              (Z * Z);
          }
          Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
        }
        noderho[i][ibasis][inode] = noderho0[i][ibasis][inode] * G;

        if (this->mix_ref_t == 1) {
          if (this->ibar_meam[elti] <= 0) {
            Gbar = 1.0;
            dGbar = 0.0;
          } else {
            gam = (nodet_ave[i][ibasis][inode][0] * shp[0] 
                + nodet_ave[i][ibasis][inode][1] * shp[1] 
                + nodet_ave[i][ibasis][inode][2] * shp[2]) / (Z * Z);
            Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
          }
          rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
        } else {
          if (this->bkgd_dyn == 1) {
            rho_bkgd = this->rho0_meam[elti] * Z;
          } else {
            rho_bkgd = this->rho_ref_meam[elti];
          }
        }
        rhob = noderho[i][ibasis][inode] / rho_bkgd;
        denom = 1.0 / rho_bkgd;

        G = dG_gam(nodegamma[i][ibasis][inode], this->ibar_meam[elti], dG);

        nodedgamma1[i][ibasis][inode] = (G - 2 * dG * nodegamma[i][ibasis][inode]) * denom;

        if (!iszero(noderho0[i][ibasis][inode])) {
          nodedgamma2[i][ibasis][inode] = (dG / noderho0[i][ibasis][inode]) * denom;
        } else {
          nodedgamma2[i][ibasis][inode] = 0.0;
        }

        //     dgamma3 is nonzero only if we are using the "mixed" rule for
        //     computing t in the reference system (which is not correct, but
        //     included for backward compatibility
        if (this->mix_ref_t == 1) {
          nodedgamma3[i][ibasis][inode] = noderho0[i][ibasis][inode] * G * dGbar / (Gbar * Z * Z) * denom;
        } else {
          nodedgamma3[i][ibasis][inode] = 0.0;
        }

        Fl = embedding(this->A_meam[elti], this->Ec_meam[elti][elti], rhob, nodefrhop[i][ibasis][inode]);

        if (eflag_either != 0) {
          Fl *= scaleii;
          if (eflag_global != 0) {
            *eng_vdwl = *eng_vdwl + Fl;
          }
          if (eflag_atom != 0) {
            enode[i][ibasis][inode] = enode[i][ibasis][inode] + Fl;
          }
        }
      }
    }
  }
}

