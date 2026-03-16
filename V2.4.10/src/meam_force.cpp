#include "meam.h"
#include <cmath>
#include <algorithm>
#include "math_special.h"
#include "atom.h"
#include "element.h"
#include "error.h"

using namespace CAC_NS;
#define NDENS 37

void MEAM::meam_force(int i, int iintg, int inode, int eflag_either, int eflag_global, int eflag_atom, 
    int vflag_either, int vflag_global, int vflag_atom, double *virial,
    double *eng_vdwl, double *eatom, double **enode, int *fmap, double **scale, 
    int numneigha_half, int *firstneigha_half, int numneigha_full, int *firstneigha_full, 
    int numneighia_half, int *firstneighia_half, int *firstneighia_index_half, 
    int numneighia_full, int *firstneighia_full, int *firstneighia_index_full, 
    int *numneighnia2ia, int **firstneighnia2ia, int **firstneighnia2ia_index, 
    int *numneighnia2a, int **firstneighnia2a, int **ia2nialist,
    int fnoffset, double **vatom, double ***vnode, double ***intgf)
{
  int kk, node, m, n, p, q, ind;
  int ictype, ietype, iintpl;
  int j, jn, jintg, jintpl, jnode, jetype; 
  int k, kn, kintg, kintpl, knode, ketype; 
  int l, ln, lintg, lintpl, lnode, letype; 
  int elti, eltj, eltk, eltl;

  double tmp;
  double xitmp, yitmp, zitmp;
  double xjtmp, yjtmp, zjtmp;
  double xktmp, yktmp, zktmp;
  double xltmp, yltmp, zltmp;
  double delij[3], deljk[3], rij2, rij;
  double xij, xik, xjk, cikj, sikj, dfc, a;
  double dCikj1, dCikj2;
  double delc, rik2, rjk2;
  double v[6], fij[3], fik[3], fjk[3];
  double third, sixth;
  double dUdsij, dUdrijm[3], force;
  double dsij1, dsij2, force1, force2;
  double pp, phi, phip, etmp;
  double sij;
  double *fi,*fj,*fk;
  double *vir_i,*vir_j,*vir_k; // virial pointers
  double ivirial_scale,jvirial_scale,kvirial_scale;

  double *idensities = new double[NDENS];
  double *jdensities = new double[NDENS];
  double *kdensities = new double[NDENS];

  double scaleij;

  double drinv = 1.0 / this->delr_meam;
  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  int counter = 0;
  int counter2 = 0;
  double **ax = atom->x;
  double **atomf = atom->f;

  double ***nodex = element->nodex;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;
  double ***shape_array = element->shape_array;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int **i2ia = element->i2ia;
  tagint *etag = element->tag;
  int max_nintpl = element->max_nintpl;
  double *nodal_weight = element->nodal_weight;
  vir_i = NULL;
  ivirial_scale = 0.0;
  if (iintg < 0) {
    ictype = atype[i];
    elti = fmap[ictype];
    if (elti < 0) return;
    if (vflag_atom) vir_i = vatom[i];
    if (vflag_global) ivirial_scale = 1.0;
    iintpl = -1;
    xitmp = ax[i][0];
    yitmp = ax[i][1];
    zitmp = ax[i][2];
    copy_densities(idensities,i,-1,-1);
    fi = atomf[i];
  } else {
    ictype = ctype[i];
    elti = fmap[ictype];
    if (elti < 0) return;
    ietype = etype[i];
    iintpl = i2ia[ietype][iintg];
    fi = intgf[i][iintg];

    if (inode >= 0) {
      if (vflag_atom) vir_i = vnode[i][inode];
      if (vflag_global) ivirial_scale = nodal_weight[ietype];
      xitmp = nodex[i][inode][0];
      yitmp = nodex[i][inode][1];
      zitmp = nodex[i][inode][2];
    } else {
      xitmp = yitmp = zitmp = 0.0;
      for (node = 0; node < npe[ietype]; node++) {
        xitmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
        yitmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
        zitmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
      }
    }
    copy_densities(idensities,i,inode,iintpl);
  }

  for (jn = 0; jn < numneigha_half + numneighia_half; jn++) {

    counter++;
    if (iszero(scrfcn[fnoffset + jn])) continue;

    fj = NULL;
    vir_j = NULL;
    jvirial_scale = 0.0;
    jnode = -1;
    jintg = -1;
    jintpl = -1;

    // neighboring atoms

    if (jn < numneigha_half) {

      j = firstneigha_half[jn];
      eltj = fmap[atype[j]];
      if (eltj < 0) continue;
      xjtmp = ax[j][0];
      yjtmp = ax[j][1];
      zjtmp = ax[j][2];
      scaleij = scale[ictype][atype[j]];
      fj = atomf[j];
      if (vflag_atom) vir_j = vatom[j];
      if (vflag_global) jvirial_scale = 1.0;
    }

    // neighboring interpolated atoms

    else {
      j = firstneighia_half[jn-numneigha_half];
      eltj = fmap[ctype[j]];
      if (eltj < 0) continue;
      jintpl = firstneighia_index_half[jn-numneigha_half];
      jetype = etype[j];
      jintg = ia2i[jetype][jintpl];
      if (jintg >= 0)
        jnode = i2n[jetype][jintg];
      if (jnode >= 0) {
        xjtmp = nodex[j][jnode][0];
        yjtmp = nodex[j][jnode][1];
        zjtmp = nodex[j][jnode][2];
        if (vflag_atom) vir_j = vnode[j][jnode];
        if (vflag_global) jvirial_scale = nodal_weight[jetype];
      } else {
        xjtmp = yjtmp = zjtmp = 0.0;
        for (node = 0; node < npe[jetype]; node++) {
          xjtmp += shape_array[jetype][jintpl][node]*nodex[j][node][0];
          yjtmp += shape_array[jetype][jintpl][node]*nodex[j][node][1];
          zjtmp += shape_array[jetype][jintpl][node]*nodex[j][node][2];
        }
      }
      scaleij = scale[ictype][ctype[j]];
      if (jintg >= 0) fj = intgf[j][jintg];
    }

    sij = scrfcn[fnoffset + jn] * fcpair[fnoffset + jn];
    delij[0] = xjtmp - xitmp;
    delij[1] = yjtmp - yitmp;
    delij[2] = zjtmp - zitmp;
    rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
    rij = sqrt(rij2);
    ind = this->eltind[elti][eltj];
    pp = rij * this->rdrar;
    kk = (int)pp;
    kk = std::min(kk, this->nrar - 2);
    pp = pp - kk;
    pp = std::min(pp, 1.0);
    phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + this->phirar1[ind][kk]) * pp + this->phirar[ind][kk];
    phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + this->phirar4[ind][kk];

    if (eflag_either != 0) {
      double phi_sc = phi * scaleij;
      if (eflag_global != 0)
        *eng_vdwl = *eng_vdwl + phi_sc * sij;
      if (eflag_atom != 0) {
        etmp = 0.5 * phi_sc * sij;

        if (iintg < 0) eatom[i] += etmp;
        else if (inode >= 0) enode[i][inode] += etmp;

        if (jn < numneigha_half) eatom[j] += etmp;
        else if (jnode >= 0) enode[j][jnode] += etmp;
      }
    }


    // Compute pair force
    copy_densities(jdensities,j,jnode,jintpl);

    compute_pair_force(dUdrijm,dUdsij,force,dscrfcn[fnoffset+jn],elti,eltj,scaleij,
        idensities,jdensities,delij,rij,rij2,phi,phip,sij);
    for (m = 0; m < 3; m++) {
      fij[m] = delij[m] * force + dUdrijm[m];
      fi[m] += fij[m];
      if (fj) fj[m] -= fij[m];
    }

    //     Tabulate per-atom virial as symmetrized stress tensor

    if (vflag_either && (iintpl < 0 || inode >= 0 || jintpl < 0 || jnode >= 0)) {
      v[0] = -0.5 * (delij[0] * fij[0]);
      v[1] = -0.5 * (delij[1] * fij[1]);
      v[2] = -0.5 * (delij[2] * fij[2]);
      v[3] = -0.25 * (delij[0] * fij[1] + delij[1] * fij[0]);
      v[4] = -0.25 * (delij[0] * fij[2] + delij[2] * fij[0]);
      v[5] = -0.25 * (delij[1] * fij[2] + delij[2] * fij[1]);

      for (m = 0; m < 6; m++) {
        if (vir_i) vir_i[m] += v[m];
        if (vir_j) vir_j[m] += v[m];
        virial[m] += v[m] * (ivirial_scale + jvirial_scale);
      }
    }

    //     Now compute forces on other atoms k due to change in sij
    //     2 cases here, k is atom or interpolated atom


    double dxik(0), dyik(0), dzik(0);


    const double rbound = rij2 * this->ebound_meam[elti][eltj];
    for (kn = 0; kn < numneigha_full + numneighia_full; kn++) {
      vir_k = NULL;
      kvirial_scale = 0;
      fk = NULL;
      kintpl = -1;
      kintg = -1;
      knode = -1;

      if (kn < numneigha_full) {
        k = firstneigha_full[kn];
        if ((k == j && jn < numneigha_half)) continue;
        eltk = fmap[atype[k]];
        if (eltk < 0) continue;
        xktmp = ax[k][0];
        yktmp = ax[k][1];
        zktmp = ax[k][2];
        if (vflag_atom) vir_k = vatom[k];
        if (vflag_global) kvirial_scale = 1.0;
        fk = atomf[k];
      } else {
        k = firstneighia_full[kn-numneigha_full];
        eltk = fmap[ctype[k]];
        if (eltk < 0) continue;
        kintpl = firstneighia_index_full[kn-numneigha_full];
        if (jn >= numneigha_half && k == j && kintpl == jintpl) continue;
        ketype = etype[k];
        kintg = ia2i[ketype][kintpl];
        knode = -1;
        if (kintg >= 0)
          knode = i2n[ketype][kintg];
        xktmp = yktmp = zktmp = 0.0;
        if (knode >= 0) {
          xktmp = nodex[k][knode][0];
          yktmp = nodex[k][knode][1];
          zktmp = nodex[k][knode][2];
          if (vflag_atom) vir_k = vnode[k][knode];
          if (vflag_global) kvirial_scale = nodal_weight[ketype];
        } else {
          xktmp = yktmp = zktmp = 0.0;
          for (node = 0; node < npe[ketype]; node++) {
            xktmp += shape_array[ketype][kintpl][node]*nodex[k][node][0];
            yktmp += shape_array[ketype][kintpl][node]*nodex[k][node][1];
            zktmp += shape_array[ketype][kintpl][node]*nodex[k][node][2];
          }
        }
        if (kintg >= 0) fk = intgf[k][kintg];
      }

      double Cmax = this->Cmax_meam[elti][eltj][eltk];
      double Cmin = this->Cmin_meam[elti][eltj][eltk];

      dsij1 = 0.0;
      dsij2 = 0.0;
      rik2 = 0.0;
      delc = Cmax - Cmin;
      deljk[0] = xktmp - xjtmp;
      deljk[1] = yktmp - yjtmp;
      deljk[2] = zktmp - zjtmp;
      rjk2 = deljk[0] * deljk[0] + deljk[1] * deljk[1] + deljk[2] * deljk[2];

      if (rjk2 <= rbound) {
        dxik = xktmp - xitmp;
        dyik = yktmp - yitmp;
        dzik = zktmp - zitmp;
        rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
        if (rik2 <= rbound) {
          if (!iszero(sij) && !isone(sij)) {
            xik = rik2 / rij2;
            xjk = rjk2 / rij2;
            a = 1 - (xik - xjk) * (xik - xjk);
            if (!iszero(a)) {
              cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
              if (cikj >= Cmin && cikj <= Cmax) {
                cikj = (cikj - Cmin) / delc;
                sikj = dfcut(cikj, dfc);
                dCfunc2(rij2, rik2, rjk2, dCikj1, dCikj2);
                a = sij / delc * dfc / sikj;
                dsij1 = a * dCikj1;
                dsij2 = a * dCikj2;
              }
            }
          }
        }
      }
      if (!iszero(dsij1) || !iszero(dsij2)) {
        force1 = dUdsij * dsij1;
        force2 = dUdsij * dsij2;
        fik[0] = force1 * dxik;
        fik[1] = force1 * dyik;
        fik[2] = force1 * dzik;
        fjk[0] = force2 * deljk[0];
        fjk[1] = force2 * deljk[1];
        fjk[2] = force2 * deljk[2];
        for (m = 0; m < 3; m++) {
          fi[m] += fik[m];
          if (fj) fj[m] += fjk[m];
          if (fk) fk[m] -= fik[m] + fjk[m];
        }

        //     Tabulate per-atom virial as symmetrized stress tensor

        if (vflag_either && (iintpl < 0 || inode >= 0 || 
                             jintpl < 0 || jnode >= 0 || 
                             kintpl < 0 || knode >= 0)) {
          v[0] = -third * (dxik * fik[0] + deljk[0] * fjk[0]);
          v[1] = -third * (dyik * fik[1] + deljk[1] * fjk[1]);
          v[2] = -third * (dzik * fik[2] + deljk[2] * fjk[2]);
          v[3] = -sixth * (dxik * fik[1] + deljk[0] * fjk[1] + dyik * fik[0] + deljk[1] * fjk[0]);
          v[4] = -sixth * (dxik * fik[2] + deljk[0] * fjk[2] + dzik * fik[0] + deljk[2] * fjk[0]);
          v[5] = -sixth * (dyik * fik[2] + deljk[1] * fjk[2] + dzik * fik[1] + deljk[2] * fjk[1]);

          for (m = 0; m < 6; m++) {
            if (vir_i) vir_i[m] += v[m];
            if (vir_j) vir_j[m] += v[m];
            if (vir_k) vir_k[m] += v[m];
            virial[m] += v[m] * (ivirial_scale + jvirial_scale + kvirial_scale);
          }
        }
      }

      // special case where j and k are both interpolated atoms and not integration points 
      // since j-k pair will not be in the loop but also have threebody component to i
      // skip half since j-k show up twice

      if (kintpl < 0 || kintg >= 0) continue;
      if (jintpl < 0 || jintg >= 0) continue;

      if ((etag[k]-1)*max_nintpl + kintpl < 
          (etag[j]-1)*max_nintpl + jintpl)
        continue;

      if (rjk2 >= this->cutforcesq) continue;

      
      double rjkbound = rjk2 * this->ebound_meam[eltj][eltk];

      if (rij2 > rjkbound) continue;
      if (iszero(rik2)) {
        dxik = xktmp - xitmp;
        dyik = yktmp - yitmp;
        dzik = zktmp - zitmp;
        rik2 = dxik * dxik + dyik * dyik + dzik * dzik;
      }
      if (rik2 > rjkbound) continue;

      double sjk = 1.0;
      double rjk = sqrt(rjk2);
      double rnorm = (this->cutforce - rjk) * drinv;
      double fcjk,sfcjk,dscrfcnjk;
      double sjlk,dCjlk,coef1;

      // need to calculate screening functions for j-k pair, 
      // use neighbor list of j to screen j-k pair
      
      int jnia = ia2nialist[j][jintpl];
      int *lalist = firstneighnia2a[jnia];
      int *lialist = firstneighnia2ia[jnia];
      int *lialist_index = firstneighnia2ia_index[jnia];
      int lanum = numneighnia2a[jnia];
      int lianum = numneighnia2ia[jnia];

      for (ln = 0; ln < lanum + lianum; ln++) {
        if (ln < lanum) {
          l = lalist[ln];
          eltl = fmap[atype[l]];
          if (eltl < 0) continue;
          xltmp = ax[l][0];
          yltmp = ax[l][1];
          zltmp = ax[l][2];
        } else {
          l = lialist[ln-lanum];
          lintpl = lialist_index[ln-lanum];
          if (l == k && lintpl == kintpl) continue;
          eltl = fmap[ctype[l]];
          if (eltl < 0) continue;
          letype = etype[l];
          lintg = ia2i[letype][lintpl];
          if (lintg >= 0) lnode = i2n[letype][lintg];
          else lnode = -1;
          if (lnode < 0) {
            xltmp = yltmp = zltmp = 0.0;
            for (node = 0; node < npe[letype]; node++) {
              xltmp += shape_array[letype][lintpl][node]*nodex[l][node][0];
              yltmp += shape_array[letype][lintpl][node]*nodex[l][node][1];
              zltmp += shape_array[letype][lintpl][node]*nodex[l][node][2];
            }
          } else {
            xltmp = nodex[l][lnode][0];
            yltmp = nodex[l][lnode][1];
            zltmp = nodex[l][lnode][2];
          }
        }
        int flag = compute_screen(xjtmp,yjtmp,zjtmp,xktmp,yktmp,zktmp,xltmp,yltmp,zltmp,eltj,eltk,eltl,rjkbound,rjk2,sjlk,coef1,dCjlk);
        if (flag == 1) {
          sjk *= sjlk;
          dscrfcnjk += coef1 * dCjlk;
        } else if (flag == 2) {
          sjk = 0;
          dscrfcnjk = 0.0;
          break;
        }

      } // end of l loop

      if (iszero(sjk) || isone(sjk)) continue;


      fcjk = dfcut(rnorm, dfc);  
      dscrfcnjk = dscrfcnjk * sjk * fcjk - sjk * dfc * drinv / rjk;

      sjk *= fcjk;

      // densities for j and k will be interpolated.

      ind = this->eltind[eltj][eltk];
      pp = rjk * this->rdrar;
      kk = (int)pp;
      kk = std::min(kk, this->nrar - 2);
      pp = pp - kk;
      pp = std::min(pp, 1.0);

      phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + this->phirar1[ind][kk]) * pp + this->phirar[ind][kk];
      phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + this->phirar4[ind][kk];
      double dUdrjkm[3],dUdsjk;

      copy_densities(kdensities,k,-1,kintpl);

      compute_pair_force(dUdrjkm,dUdsjk,force,dscrfcnjk,eltj,eltk,scale[ctype[j]][ctype[k]],
          jdensities,kdensities,deljk,rjk,rjk2,phi,phip,sjk);
      xij = rij2 / rjk2;
      xik = rik2 / rjk2;
      a = 1 - (xij - xik) * (xij - xik);
      if (iszero(a)) continue;

      double cjik;
      cjik = (2.0 * (xij + xik) + a - 2.0) / a;
      const double Cmax2 = this->Cmax_meam[eltj][eltk][elti];
      const double Cmin2 = this->Cmin_meam[eltj][eltk][elti];
      if (cjik < Cmin2 || cjik > Cmax2) continue;

      double dsjk1,dsjk2,dCjik1,dCjik2,sjik;
      delc = Cmax2 - Cmin2;
      cjik = (cjik - Cmin) / delc;
      sjik = dfcut(cjik, dfc);
      dCfunc2(rjk2, rij2, rik2, dCjik1, dCjik2);
      a = sjk / delc * dfc / sjik;
      dsjk1 = a * dCjik1;
      dsjk2 = a * dCjik2;

      force1 = dUdsjk * dsjk1;
      force2 = dUdsjk * dsjk2;

      // reverse the sign for f the sign of delij and dxik are from i to j and k 
      // virials scale with r^2 so not reversed

      fij[0] = force1 * delij[0];
      fij[1] = force1 * delij[1];
      fij[2] = force1 * delij[2];
      fik[0] = force2 * dxik;
      fik[1] = force2 * dyik;
      fik[2] = force2 * dzik;
      fi[0] += fij[0] + fik[0];
      fi[1] += fij[1] + fik[1];
      fi[2] += fij[2] + fik[2];

      if (vflag_either && (iintpl < 0 || inode >= 0)) {
        v[0] = -third * (delij[0] * fij[0] + dxik * fik[0]);
        v[1] = -third * (delij[1] * fij[1] + dyik * fik[1]);
        v[2] = -third * (delij[2] * fij[2] + dzik * fik[2]);
        v[3] = -sixth * (delij[0] * fij[1] + dxik * fik[1] + delij[1] * fij[0] + dyik * fik[0]);
        v[4] = -sixth * (delij[0] * fij[2] + dxik * fik[2] + delij[2] * fij[0] + dzik * fik[0]);
        v[5] = -sixth * (delij[1] * fij[2] + dyik * fik[2] + delij[2] * fij[1] + dzik * fik[1]);
        for (m = 0; m < 6; m++) {
          if (vir_i) vir_i[m] += v[m];
          virial[m] += v[m] * ivirial_scale;
        }
      }
    } //     end of k loop
  } //     end of j loop

  delete [] idensities;
  delete [] jdensities;
  delete [] kdensities;
}

/**************************************
Output: dUdrijm[3], dUdsij, and force
 ***************************************/

void MEAM::compute_pair_force(double *dUdrijm, double &dUdsij, double &force,
    double pair_dscrfcn, int elti, int eltj, double scaleij, double *idensities, double *jdensities, 
    double *delij, double rij, double rij2, double phi, double phip, double sij) 
{
  int n, m, p, q, nv2, nv3;
  double shpi[3], shpj[3];
  double dUdrij;
  double a1, a1i, a1j, a2, a2i, a2j, a3, a3i, a3j, rij3;
  double ai, aj, ro0i, ro0j, invrei, invrej;
  double rhoa0j, drhoa0j, rhoa0i, drhoa0i;
  double rhoa1j, drhoa1j, rhoa1i, drhoa1i;
  double rhoa2j, drhoa2j, rhoa2i, drhoa2i;
  double a3a, rhoa3j, drhoa3j, rhoa3i, drhoa3i;
  double drho0dr1, drho0dr2, drho0ds1, drho0ds2;
  double drho1dr1, drho1dr2, drho1ds1, drho1ds2;
  double drho1drm1[3], drho1drm2[3];
  double drho2dr1, drho2dr2, drho2ds1, drho2ds2;
  double drho2drm1[3], drho2drm2[3];
  double drho3dr1, drho3dr2, drho3ds1, drho3ds2;
  double drho3drm1[3], drho3drm2[3];
  double dt1dr1, dt1dr2, dt1ds1, dt1ds2;
  double dt2dr1, dt2dr2, dt2ds1, dt2ds2;
  double dt3dr1, dt3dr2, dt3ds1, dt3ds2;
  double drhodr1, drhodr2, drhods1, drhods2, drhodrm1[3], drhodrm2[3];
  double arg, arg1i1, arg1j1, arg1i2, arg1j2, arg1i3, arg1j3, arg3i3, arg3j3;
  double t1i, t2i, t3i, t1j, t2j, t3j;

  double rho0i,rho1i,rho2i,rho3i,frhopi,dgamma1i,dgamma2i,dgamma3i,arho2bi; 
  double *arho1i,*arho2i,*arho3i,*arho3bi,*t_avei,*tsq_avei;
  double rho0j,rho1j,rho2j,rho3j,frhopj,dgamma1j,dgamma2j,dgamma3j,arho2bj; 
  double *arho1j,*arho2j,*arho3j,*arho3bj,*t_avej,*tsq_avej;
  rho0i = idensities[0];
  rho1i = idensities[1];
  rho2i = idensities[2];
  rho3i = idensities[3];
  frhopi = idensities[4];
  dgamma1i = idensities[5];
  dgamma2i = idensities[6];
  dgamma3i = idensities[7];
  arho1i = &idensities[8];
  arho2i = &idensities[11];
  arho2bi = idensities[17];
  arho3i = &idensities[18];
  arho3bi = &idensities[28];
  t_avei = &idensities[31];
  tsq_avei = &idensities[34];

  rho0j = jdensities[0];
  rho1j = jdensities[1];
  rho2j = jdensities[2];
  rho3j = jdensities[3];
  frhopj = jdensities[4];
  dgamma1j = jdensities[5];
  dgamma2j = jdensities[6];
  dgamma3j = jdensities[7];
  arho1j = &jdensities[8];
  arho2j = &jdensities[11];
  arho2bj = jdensities[17];
  arho3j = &jdensities[18];
  arho3bj = &jdensities[28];
  t_avej = &jdensities[31];
  tsq_avej = &jdensities[34];

  invrei = 1.0 / this->re_meam[elti][elti];
  ai = rij * invrei - 1.0;
  ro0i = this->rho0_meam[elti];
  rhoa0i = ro0i * MathSpecial::fm_exp(-this->beta0_meam[elti] * ai);
  drhoa0i = -this->beta0_meam[elti] * invrei * rhoa0i;
  rhoa1i = ro0i * MathSpecial::fm_exp(-this->beta1_meam[elti] * ai);
  drhoa1i = -this->beta1_meam[elti] * invrei * rhoa1i;
  rhoa2i = ro0i * MathSpecial::fm_exp(-this->beta2_meam[elti] * ai);
  drhoa2i = -this->beta2_meam[elti] * invrei * rhoa2i;
  rhoa3i = ro0i * MathSpecial::fm_exp(-this->beta3_meam[elti] * ai);
  drhoa3i = -this->beta3_meam[elti] * invrei * rhoa3i;

  if (elti != eltj) {
    invrej = 1.0 / this->re_meam[eltj][eltj];
    aj = rij * invrej - 1.0;
    ro0j = this->rho0_meam[eltj];
    rhoa0j = ro0j * MathSpecial::fm_exp(-this->beta0_meam[eltj] * aj);
    drhoa0j = -this->beta0_meam[eltj] * invrej * rhoa0j;
    rhoa1j = ro0j * MathSpecial::fm_exp(-this->beta1_meam[eltj] * aj);
    drhoa1j = -this->beta1_meam[eltj] * invrej * rhoa1j;
    rhoa2j = ro0j * MathSpecial::fm_exp(-this->beta2_meam[eltj] * aj);
    drhoa2j = -this->beta2_meam[eltj] * invrej * rhoa2j;
    rhoa3j = ro0j * MathSpecial::fm_exp(-this->beta3_meam[eltj] * aj);
    drhoa3j = -this->beta3_meam[eltj] * invrej * rhoa3j;
  } else {
    rhoa0j = rhoa0i;
    drhoa0j = drhoa0i;
    rhoa1j = rhoa1i;
    drhoa1j = drhoa1i;
    rhoa2j = rhoa2i;
    drhoa2j = drhoa2i;
    rhoa3j = rhoa3i;
    drhoa3j = drhoa3i;
  }

  const double t1mi = this->t1_meam[elti];
  const double t2mi = this->t2_meam[elti];
  const double t3mi = this->t3_meam[elti];
  const double t1mj = this->t1_meam[eltj];
  const double t2mj = this->t2_meam[eltj];
  const double t3mj = this->t3_meam[eltj];

  if (this->ialloy == 1) {
    rhoa1j  *= t1mj;
    rhoa2j  *= t2mj;
    rhoa3j  *= t3mj;
    rhoa1i  *= t1mi;
    rhoa2i  *= t2mi;
    rhoa3i  *= t3mi;
    drhoa1j *= t1mj;
    drhoa2j *= t2mj;
    drhoa3j *= t3mj;
    drhoa1i *= t1mi;
    drhoa2i *= t2mi;
    drhoa3i *= t3mi;
  }

  nv2 = 0;
  nv3 = 0;
  arg1i1 = 0.0;
  arg1j1 = 0.0;
  arg1i2 = 0.0;
  arg1j2 = 0.0;
  arg1i3 = 0.0;
  arg1j3 = 0.0;
  arg3i3 = 0.0;
  arg3j3 = 0.0;
  for (n = 0; n < 3; n++) {
    for (p = n; p < 3; p++) {
      for (q = p; q < 3; q++) {
        arg = delij[n] * delij[p] * delij[q] * this->v3D[nv3];
        arg1i3 = arg1i3 + arho3i[nv3] * arg;
        arg1j3 = arg1j3 - arho3j[nv3] * arg;
        nv3 = nv3 + 1;
      }
      arg = delij[n] * delij[p] * this->v2D[nv2];
      arg1i2 = arg1i2 + arho2i[nv2] * arg;
      arg1j2 = arg1j2 + arho2j[nv2] * arg;
      nv2 = nv2 + 1;
    }
    arg1i1 = arg1i1 + arho1i[n] * delij[n];
    arg1j1 = arg1j1 - arho1j[n] * delij[n];
    arg3i3 = arg3i3 + arho3bi[n] * delij[n];
    arg3j3 = arg3j3 - arho3bj[n] * delij[n];
  }

  //     rho0 terms
  drho0dr1 = drhoa0j * sij;
  drho0dr2 = drhoa0i * sij;

  //     rho1 terms
  a1 = 2 * sij / rij;
  drho1dr1 = a1 * (drhoa1j - rhoa1j / rij) * arg1i1;
  drho1dr2 = a1 * (drhoa1i - rhoa1i / rij) * arg1j1;
  a1 = 2.0 * sij / rij;
  for (m = 0; m < 3; m++) {
    drho1drm1[m] = a1 * rhoa1j * arho1i[m];
    drho1drm2[m] = -a1 * rhoa1i * arho1j[m];
  }

  //     rho2 terms
  a2 = 2 * sij / rij2;
  drho2dr1 = a2 * (drhoa2j - 2 * rhoa2j / rij) * arg1i2 - 2.0 / 3.0 * arho2bi * drhoa2j * sij;
  drho2dr2 = a2 * (drhoa2i - 2 * rhoa2i / rij) * arg1j2 - 2.0 / 3.0 * arho2bj * drhoa2i * sij;
  a2 = 4 * sij / rij2;
  for (m = 0; m < 3; m++) {
    drho2drm1[m] = 0.0;
    drho2drm2[m] = 0.0;
    for (n = 0; n < 3; n++) {
      drho2drm1[m] = drho2drm1[m] + arho2i[this->vind2D[m][n]] * delij[n];
      drho2drm2[m] = drho2drm2[m] - arho2j[this->vind2D[m][n]] * delij[n];
    }
    drho2drm1[m] = a2 * rhoa2j * drho2drm1[m];
    drho2drm2[m] = -a2 * rhoa2i * drho2drm2[m];
  }

  //     rho3 terms
  rij3 = rij * rij2;
  a3 = 2 * sij / rij3;
  a3a = 6.0 / 5.0 * sij / rij;
  drho3dr1 = a3 * (drhoa3j - 3 * rhoa3j / rij) * arg1i3 - a3a * (drhoa3j - rhoa3j / rij) * arg3i3;
  drho3dr2 = a3 * (drhoa3i - 3 * rhoa3i / rij) * arg1j3 - a3a * (drhoa3i - rhoa3i / rij) * arg3j3;
  a3 = 6 * sij / rij3;
  a3a = 6 * sij / (5 * rij);
  for (m = 0; m < 3; m++) {
    drho3drm1[m] = 0.0;
    drho3drm2[m] = 0.0;
    nv2 = 0;
    for (n = 0; n < 3; n++) {
      for (p = n; p < 3; p++) {
        arg = delij[n] * delij[p] * this->v2D[nv2];
        drho3drm1[m] = drho3drm1[m] + arho3i[this->vind3D[m][n][p]] * arg;
        drho3drm2[m] = drho3drm2[m] + arho3j[this->vind3D[m][n][p]] * arg;
        nv2 = nv2 + 1;
      }
    }
    drho3drm1[m] = (a3 * drho3drm1[m] - a3a * arho3bi[m]) * rhoa3j;
    drho3drm2[m] = (-a3 * drho3drm2[m] + a3a * arho3bj[m]) * rhoa3i;
  }

  //     Compute derivatives of weighting functions t wrt rij
  t1i = t_avei[0];
  t2i = t_avei[1];
  t3i = t_avei[2];
  t1j = t_avej[0];
  t2j = t_avej[1];
  t3j = t_avej[2];

  if (this->ialloy == 1) {

    a1i = fdiv_zero(drhoa0j * sij, tsq_avei[0]);
    a1j = fdiv_zero(drhoa0i * sij, tsq_avej[0]);
    a2i = fdiv_zero(drhoa0j * sij, tsq_avei[1]);
    a2j = fdiv_zero(drhoa0i * sij, tsq_avej[1]);
    a3i = fdiv_zero(drhoa0j * sij, tsq_avei[2]);
    a3j = fdiv_zero(drhoa0i * sij, tsq_avej[2]);

    dt1dr1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
    dt1dr2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
    dt2dr1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
    dt2dr2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
    dt3dr1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
    dt3dr2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

  } else if (this->ialloy == 2) {

    dt1dr1 = 0.0;
    dt1dr2 = 0.0;
    dt2dr1 = 0.0;
    dt2dr2 = 0.0;
    dt3dr1 = 0.0;
    dt3dr2 = 0.0;

  } else {

    ai = 0.0;
    if (!iszero(rho0i))
      ai = drhoa0j * sij / rho0i;
    aj = 0.0;
    if (!iszero(rho0j))
      aj = drhoa0i * sij / rho0j;

    dt1dr1 = ai * (t1mj - t1i);
    dt1dr2 = aj * (t1mi - t1j);
    dt2dr1 = ai * (t2mj - t2i);
    dt2dr2 = aj * (t2mi - t2j);
    dt3dr1 = ai * (t3mj - t3i);
    dt3dr2 = aj * (t3mi - t3j);
  }

  //     Compute derivatives of total density wrt rij, sij and rij(3)
  get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpi);
  get_shpfcn(this->lattce_meam[eltj][eltj], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shpj);

  drhodr1 = dgamma1i * drho0dr1 +
    dgamma2i * (dt1dr1 * rho1i + t1i * drho1dr1 + dt2dr1 * rho2i + t2i * drho2dr1 +
        dt3dr1 * rho3i + t3i * drho3dr1) -
    dgamma3i * (shpi[0] * dt1dr1 + shpi[1] * dt2dr1 + shpi[2] * dt3dr1);
  drhodr2 = dgamma1j * drho0dr2 +
    dgamma2j * (dt1dr2 * rho1j + t1j * drho1dr2 + dt2dr2 * rho2j + t2j * drho2dr2 +
        dt3dr2 * rho3j + t3j * drho3dr2) -
    dgamma3j * (shpj[0] * dt1dr2 + shpj[1] * dt2dr2 + shpj[2] * dt3dr2);
  for (m = 0; m < 3; m++) {
    drhodrm1[m] = dgamma2i * (t1i * drho1drm1[m] + t2i * drho2drm1[m] + t3i * drho3drm1[m]);
    drhodrm2[m] = dgamma2j * (t1j * drho1drm2[m] + t2j * drho2drm2[m] + t3j * drho3drm2[m]);
  }

  //     Compute derivatives wrt sij, but only if necessary
  if (!iszero(pair_dscrfcn)) {
    drho0ds1 = rhoa0j;
    drho0ds2 = rhoa0i;
    a1 = 2.0 / rij;
    drho1ds1 = a1 * rhoa1j * arg1i1;
    drho1ds2 = a1 * rhoa1i * arg1j1;
    a2 = 2.0 / rij2;
    drho2ds1 = a2 * rhoa2j * arg1i2 - 2.0 / 3.0 * arho2bi * rhoa2j;
    drho2ds2 = a2 * rhoa2i * arg1j2 - 2.0 / 3.0 * arho2bj * rhoa2i;
    a3 = 2.0 / rij3;
    a3a = 6.0 / (5.0 * rij);
    drho3ds1 = a3 * rhoa3j * arg1i3 - a3a * rhoa3j * arg3i3;
    drho3ds2 = a3 * rhoa3i * arg1j3 - a3a * rhoa3i * arg3j3;

    if (this->ialloy == 1) {
      a1i = fdiv_zero(rhoa0j, tsq_avei[0]);
      a1j = fdiv_zero(rhoa0i, tsq_avej[0]);
      a2i = fdiv_zero(rhoa0j, tsq_avei[1]);
      a2j = fdiv_zero(rhoa0i, tsq_avej[1]);
      a3i = fdiv_zero(rhoa0j, tsq_avei[2]);
      a3j = fdiv_zero(rhoa0i, tsq_avej[2]);

      dt1ds1 = a1i * (t1mj - t1i * MathSpecial::square(t1mj));
      dt1ds2 = a1j * (t1mi - t1j * MathSpecial::square(t1mi));
      dt2ds1 = a2i * (t2mj - t2i * MathSpecial::square(t2mj));
      dt2ds2 = a2j * (t2mi - t2j * MathSpecial::square(t2mi));
      dt3ds1 = a3i * (t3mj - t3i * MathSpecial::square(t3mj));
      dt3ds2 = a3j * (t3mi - t3j * MathSpecial::square(t3mi));

    } else if (this->ialloy == 2) {

      dt1ds1 = 0.0;
      dt1ds2 = 0.0;
      dt2ds1 = 0.0;
      dt2ds2 = 0.0;
      dt3ds1 = 0.0;
      dt3ds2 = 0.0;

    } else {

      ai = 0.0;
      if (!iszero(rho0i))
        ai = rhoa0j / rho0i;
      aj = 0.0;
      if (!iszero(rho0j))
        aj = rhoa0i / rho0j;

      dt1ds1 = ai * (t1mj - t1i);
      dt1ds2 = aj * (t1mi - t1j);
      dt2ds1 = ai * (t2mj - t2i);
      dt2ds2 = aj * (t2mi - t2j);
      dt3ds1 = ai * (t3mj - t3i);
      dt3ds2 = aj * (t3mi - t3j);
    }

    drhods1 = dgamma1i * drho0ds1 +
      dgamma2i * (dt1ds1 * rho1i + t1i * drho1ds1 + dt2ds1 * rho2i + t2i * drho2ds1 +
          dt3ds1 * rho3i + t3i * drho3ds1) -
      dgamma3i * (shpi[0] * dt1ds1 + shpi[1] * dt2ds1 + shpi[2] * dt3ds1);
    drhods2 = dgamma1j * drho0ds2 +
      dgamma2j * (dt1ds2 * rho1j + t1j * drho1ds2 + dt2ds2 * rho2j + t2j * drho2ds2 +
          dt3ds2 * rho3j + t3j * drho3ds2) -
      dgamma3j * (shpj[0] * dt1ds2 + shpj[1] * dt2ds2 + shpj[2] * dt3ds2);
  }

  //     Compute derivatives of energy wrt rij, sij and rij[3]
  dUdrij = phip * sij + frhopi * drhodr1 + frhopj * drhodr2;
  dUdsij = 0.0;
  if (!iszero(pair_dscrfcn)) {
    dUdsij = phi + frhopi * drhods1 + frhopj * drhods2;
  }
  for (m = 0; m < 3; m++) {
    dUdrijm[m] = frhopi * drhodrm1[m] + frhopj * drhodrm2[m];
  }
  if (!isone(scaleij)) {
    dUdrij *= scaleij;
    dUdsij *= scaleij;
    dUdrijm[0] *= scaleij;
    dUdrijm[1] *= scaleij;
    dUdrijm[2] *= scaleij;
  }
  //     Return the part of the force due to dUdrij and dUdsij

  force = dUdrij / rij + dUdsij * pair_dscrfcn;
}


void MEAM::copy_densities(double *densities, int i, int inode, int iintpl)
{
  int m = 0;
  double tmp;

  // atoms
  if (iintpl < 0) {
    densities[m++] = atomrho0[i];
    densities[m++] = atomrho1[i];
    densities[m++] = atomrho2[i];
    densities[m++] = atomrho3[i];
    densities[m++] = atomfrhop[i];
    densities[m++] = atomdgamma1[i];
    densities[m++] = atomdgamma2[i];
    densities[m++] = atomdgamma3[i];
    for (int j = 0; j < 3; j++)
      densities[m++] = atomarho1[i][j];
    for (int j = 0; j < 6; j++)
      densities[m++] = atomarho2[i][j];
    densities[m++] = atomarho2b[i];
    for (int j = 0; j < 10; j++)
      densities[m++] = atomarho3[i][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = atomarho3b[i][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = atomt_ave[i][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = atomtsq_ave[i][j];

  } 

  // nodes
  else if (inode >= 0) {
    densities[m++] = noderho0[i][inode];
    densities[m++] = noderho1[i][inode];
    densities[m++] = noderho2[i][inode];
    densities[m++] = noderho3[i][inode];
    densities[m++] = nodefrhop[i][inode];
    densities[m++] = nodedgamma1[i][inode];
    densities[m++] = nodedgamma2[i][inode];
    densities[m++] = nodedgamma3[i][inode];
    for (int j = 0; j < 3; j++)
      densities[m++] = nodearho1[i][inode][j];
    for (int j = 0; j < 6; j++)
      densities[m++] = nodearho2[i][inode][j];
    densities[m++] = nodearho2b[i][inode];
    for (int j = 0; j < 10; j++)
      densities[m++] = nodearho3[i][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodearho3b[i][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodet_ave[i][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodetsq_ave[i][inode][j];

  } 

  // interpolated atoms
  else {
    for (int j = 0; j < NDENS; j++)
      densities[j] = 0.0;
    int *npe = element->npe;
    int ietype = element->etype[i];
    double ***shape_array = element->shape_array;

    for (int node = 0; node < npe[ietype]; node++) {
      tmp = shape_array[ietype][iintpl][node];
      if (iszero(tmp)) continue; 
      m = 0;
      densities[m++] += tmp*noderho0[i][node];
      densities[m++] += tmp*noderho1[i][node];
      densities[m++] += tmp*noderho2[i][node];
      densities[m++] += tmp*noderho3[i][node];
      densities[m++] += tmp*nodefrhop[i][node];
      densities[m++] += tmp*nodedgamma1[i][node];
      densities[m++] += tmp*nodedgamma2[i][node];
      densities[m++] += tmp*nodedgamma3[i][node];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp*nodearho1[i][node][j];
      for (int j = 0; j < 6; j++)
        densities[m++] += tmp*nodearho2[i][node][j];
      densities[m++] += tmp*nodearho2b[i][node];
      for (int j = 0; j < 10; j++)
        densities[m++] += tmp*nodearho3[i][node][j];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp*nodearho3b[i][node][j];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp*nodet_ave[i][node][j];
      for (int j = 0; j < 3; j++) 
        densities[m++] += tmp*nodetsq_ave[i][node][j];
    }
  }
}
