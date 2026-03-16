#include "meam.h"
#include <cmath>
#include <algorithm>
#include "math_special.h"
#include "atom.h"
#include "element.h"
#include "error.h"

using namespace CAC_NS;
#define NDENS 37
//#define TEST 453172
//#define TEST 436793
#define TEST 311
#define TESTGCELL 12
#define TESTBASIS 3

void MEAM::meam_force(int i, int iindex, int eflag_either, int eflag_global, int eflag_atom, 
    int vflag_either, int vflag_global, int vflag_atom, double *virial, 
    double *eng_vdwl, double *eatom, double ***enode, int *fmap, double **scale, 
    int numneigh, int *firstneigh, int *firstneighindex, 
    int numneigh_full, int *firstneigh_full, int *firstneighindex_full, 
    int *nvanumneigh, int **nvafirstneigh, int **nvafirstneighindex, 
    int ***va2nvalist, int fnoffset, double **vatom, double ****vnode)
{
  int kk, node, m, n, p, q, ind;
  int ictype, ietype, igcell, iucell, ibasis, iapc, inpe, inode;
  int j, jn, jindex, jgcell, jucell, jbasis, jnode, jetype, jctype, japc, jnpe; 
  int k, kn, kindex, kgcell, kucell, kbasis, knode, ketype, kctype, kapc, knpe;
  int l, ln, lindex, lgcell, lucell, lbasis, lnode, letype, lapc, lnpe;
  int elti, eltj, eltk, eltl;

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
  double *fi, *fj, *fk;
  double *ivptr, *jvptr, *kvptr; // virial pointers
  double ivscale, jvscale, kvscale;
  double iescale, jescale;

  double *idensities = new double[NDENS];
  double *jdensities = new double[NDENS];
  double *kdensities = new double[NDENS];

  tagint *tag = element->tag;
  double scaleij;

  double drinv = 1.0 / this->delr_meam;
  third = 1.0 / 3.0;
  sixth = 1.0 / 6.0;

  int counter = 0;
  int counter2 = 0;
  double **ax = atom->x;
  double **atomf = atom->f;

  double ****nodex = element->nodex;
  double ****gaussf = element->gaussf;
  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;
  double ***shape_array = element->shape_array;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  int **g2u = element->g2u;
  tagint *etag = element->tag;
  int max_nucell = element->max_nucell;
  int maxapc = element->maxapc;
  double *nodal_weight = element->nodal_weight;
  ivptr = NULL;
  ivscale = iescale = 0.0;
 
  // i is atom
  
  if (iindex < 0) {
    ictype = atype[i];
    elti = fmap[ictype];
    if (elti < 0) return;
    if (vflag_atom) ivptr = vatom[i];
    if (vflag_global) ivscale = 1.0;
    if (eflag_global) iescale = 1.0;
    iucell = -1;
    xitmp = ax[i][0];
    yitmp = ax[i][1];
    zitmp = ax[i][2];
    copy_densities(idensities, i, -1, -1, -1);
    fi = atomf[i];
  } 

  // i is gauss point
  
  else {
    ietype = etype[i];
    iapc = apc[ietype];
    inpe = npe[ietype];
    igcell = iindex / iapc;
    ibasis = iindex % iapc;
    ictype = ctype[i][ibasis];
    elti = fmap[ictype];
    if (elti < 0) return;
    inode = g2n[ietype][igcell];
    iucell = g2u[ietype][igcell];
    fi = gaussf[i][ibasis][igcell];
    if (inode >= 0) {
      if (vflag_atom) ivptr = vnode[i][ibasis][inode];
      if (vflag_global) ivscale = nodal_weight[ietype];
      //if (eflag_global) iescale = nodal_weight[ietype];
      if (eflag_global) iescale = 1.0;
      xitmp = nodex[i][ibasis][inode][0];
      yitmp = nodex[i][ibasis][inode][1];
      zitmp = nodex[i][ibasis][inode][2];
    } else {
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
    copy_densities(idensities, i, inode, iucell, ibasis);
  }

  for (jn = 0; jn < numneigh; jn++) {
    if (iszero(scrfcn[fnoffset + jn])) continue;
    j = firstneigh[jn];
    jindex = firstneighindex[jn];
    fj = NULL;
    jvptr = NULL;
    jvscale = jescale = 0.0;
    jnode = jgcell = jucell = jbasis = -1;

    // j is atom

    if (jindex < 0) {
      eltj = fmap[atype[j]];
      if (eltj < 0) continue;
      xjtmp = ax[j][0];
      yjtmp = ax[j][1];
      zjtmp = ax[j][2];
      jctype = atype[j];
      scaleij = scale[ictype][jctype];
      fj = atomf[j];
      if (vflag_atom) jvptr = vatom[j];
      if (vflag_global) jvscale = 1.0;
      if (eflag_global) jescale = 1.0;
    }

    // j is virtual atom

    else {
      jetype = etype[j];
      japc = apc[jetype];
      jnpe = npe[jetype];
      jucell = jindex / japc;
      jbasis = jindex % japc;
      jctype = ctype[j][jbasis];
      eltj = fmap[jctype];
      if (eltj < 0) continue;
      jgcell = u2g[jetype][jucell];
      if (jgcell >= 0) {
        jnode = g2n[jetype][jgcell];
        fj = gaussf[j][jbasis][jgcell];
      }
      if (jnode >= 0) {
        xjtmp = nodex[j][jbasis][jnode][0];
        yjtmp = nodex[j][jbasis][jnode][1];
        zjtmp = nodex[j][jbasis][jnode][2];
        if (vflag_atom) jvptr = vnode[j][jbasis][jnode];
        if (vflag_global) jvscale = nodal_weight[jetype];
        //if (eflag_global) jescale = nodal_weight[jetype];
        if (eflag_global) jescale = 1.0;
      } else {
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
      scaleij = scale[ictype][jctype];
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
      double phi_sc = 0.5 * phi * scaleij * sij;
      if (!std::isfinite(phi_sc)) {
        printf("rij = %g i = %d %d %d j = %d %d %d i = %d iindex = %d jn = %d j = %d jindex = %d\n"
            ,rij,element->tag[i],ibasis,iucell,element->tag[j],jbasis,jucell,i,iindex,jn,j,jindex);
        printf("xitmp = %g %g %g xjtmp = %g %g %g rij2 = %g\n",xitmp,yitmp,zitmp,xjtmp,yjtmp,zjtmp,rij2);
        for (inode = 0; inode < 8; inode++)
          printf("nodex[%d][%d][%d] = %g %g %g\n",i,ibasis,inode
              ,nodex[i][ibasis][inode][0]
              ,nodex[i][ibasis][inode][1]
              ,nodex[i][ibasis][inode][2]);
        error->all(FLERR,"TEST"); 
      }
      if (eflag_global != 0)
        *eng_vdwl += phi_sc * (iescale + jescale);

//      printf("eflag_global=%d eng_vdwl = %g\n",eflag_global,*eng_vdwl);
      if (eflag_atom != 0) {
        if (iindex < 0) eatom[i] += phi_sc;
        else if (inode >= 0) enode[i][ibasis][inode] += phi_sc;
        if (jindex < 0) eatom[j] += phi_sc;
        else if (jnode >= 0) enode[j][jbasis][jnode] += phi_sc;
      }
    }

    // Compute pair force

    copy_densities(jdensities, j, jnode, jucell, jbasis);
    compute_pair_force(dUdrijm, dUdsij, force, dscrfcn[fnoffset + jn], elti, eltj, scaleij, 
        idensities, jdensities, delij, rij, rij2, phi, phip, sij);

    for (m = 0; m < 3; m++) {
      fij[m] = delij[m] * force + dUdrijm[m];
      fi[m] += fij[m];
      if (fj) fj[m] -= fij[m];
    }
    /*
       if (iindex >= 0)
       if (element->tag[i] == TEST && ibasis == TESTBASIS && igcell == TESTGCELL) 
       printf("%g %g %g\n",fij[0],fij[1],fij[2]);
       if (jindex >= 0)
       if (element->tag[j] == TEST && jbasis == TESTBASIS && jgcell == TESTGCELL) 
       printf("%g %g %g\n",-fij[0],-fij[1],-fij[2]);
       */

    //     Tabulate per-atom virial as symmetrized stress tensor

    if (vflag_either && (iucell < 0 || inode >= 0 || jucell < 0 || jnode >= 0)) {
      v[0] = -0.5 * (delij[0] * fij[0]);
      v[1] = -0.5 * (delij[1] * fij[1]);
      v[2] = -0.5 * (delij[2] * fij[2]);
      v[3] = -0.25 * (delij[0] * fij[1] + delij[1] * fij[0]);
      v[4] = -0.25 * (delij[0] * fij[2] + delij[2] * fij[0]);
      v[5] = -0.25 * (delij[1] * fij[2] + delij[2] * fij[1]);

      for (m = 0; m < 6; m++) {
        if (ivptr) ivptr[m] += v[m];
        if (jvptr) jvptr[m] += v[m];
        virial[m] += v[m] * (ivscale + jvscale);
      }
    }

    //     Now compute forces on other atoms k due to change in sij
    //     2 cases here, k is atom or virtual atom


    double dxik(0), dyik(0), dzik(0);


    const double rbound = rij2 * this->ebound_meam[elti][eltj];

    for (kn = 0; kn < numneigh_full; kn++) {

      k = firstneigh_full[kn];
      kindex = firstneighindex_full[kn];

      // skip if j and k are the same atom/virtual atom

      if (k == j && kindex == jindex) continue;

      kvptr = NULL;
      kvscale = 0;
      fk = NULL;
      knode = kgcell = kucell = kbasis = -1;

      // k is atom

      if (kindex < 0) {
        eltk = fmap[atype[k]];
        if (eltk < 0) continue;
        xktmp = ax[k][0];
        yktmp = ax[k][1];
        zktmp = ax[k][2];
        if (vflag_atom) kvptr = vatom[k];
        if (vflag_global) kvscale = 1.0;
        fk = atomf[k];
      } 

      // k is virtual atom

      else {
        ketype = etype[k];
        kapc = apc[ketype];
        knpe = npe[ketype];
        kucell = kindex / kapc;
        kbasis = kindex % kapc;
        kctype = ctype[k][kbasis];
        eltk = fmap[kctype];
        if (eltk < 0) continue;
        kgcell = u2g[ketype][kucell];
        if (kgcell >= 0) {
          knode = g2n[ketype][kgcell];
          fk = gaussf[k][kbasis][kgcell];
        }
        xktmp = yktmp = zktmp = 0.0;
        if (knode >= 0) {
          xktmp = nodex[k][kbasis][knode][0];
          yktmp = nodex[k][kbasis][knode][1];
          zktmp = nodex[k][kbasis][knode][2];
          if (vflag_atom) kvptr = vnode[k][kbasis][knode];
          if (vflag_global) kvscale = nodal_weight[ketype];
        } else {
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
        /*
           if (iindex >= 0)
           if (element->tag[i] == TEST && ibasis == TESTBASIS && igcell == TESTGCELL) 
           printf("%g %g %g\n",fik[0],fik[1],fik[2]);

           if (jindex >= 0)
           if (element->tag[j] == TEST && jbasis == TESTBASIS && jgcell == TESTGCELL) 
           printf("%g %g %g\n",fjk[0],fjk[1],fjk[2]);

           if (kindex >= 0)
           if (element->tag[i] == TEST && kbasis == TESTBASIS && kgcell == TESTGCELL) 
           printf("%g %g %g\n",-fik[0]-fjk[0],-fik[1]-fjk[1],-fik[2]-fjk[2]);
           */


        //     Tabulate per-atom virial as symmetrized stress tensor

        if (vflag_either && 
            (iucell < 0 || inode >= 0 || 
             jucell < 0 || jnode >= 0 || 
             kucell < 0 || knode >= 0)) {
          v[0] = -third * (dxik * fik[0] + deljk[0] * fjk[0]);
          v[1] = -third * (dyik * fik[1] + deljk[1] * fjk[1]);
          v[2] = -third * (dzik * fik[2] + deljk[2] * fjk[2]);
          v[3] = -sixth * (dxik * fik[1] + deljk[0] * fjk[1] + dyik * fik[0] + deljk[1] * fjk[0]);
          v[4] = -sixth * (dxik * fik[2] + deljk[0] * fjk[2] + dzik * fik[0] + deljk[2] * fjk[0]);
          v[5] = -sixth * (dyik * fik[2] + deljk[1] * fjk[2] + dzik * fik[1] + deljk[2] * fjk[1]);

          for (m = 0; m < 6; m++) {
            if (ivptr) ivptr[m] += v[m];
            if (jvptr) jvptr[m] += v[m];
            if (kvptr) kvptr[m] += v[m];
            virial[m] += v[m] * (ivscale + jvscale + kvscale);
          }
        }
      }

      // special case where j and k are both virtual atoms and NOT gaussian points
      // since j-k pair will not be in the loop but also have threebody component to i
      // skip half since j-k show up twice

      if (kindex < 0 || kgcell >= 0) continue;
      if (jindex < 0 || jgcell >= 0) continue;

      if ((etag[k]-1) * max_nucell * maxapc + kindex < 
          (etag[j]-1) * max_nucell * maxapc + jindex)
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
      double fcjk, sfcjk, dscrfcnjk;
      double sjlk, dCjlk, coef1;

      // need to calculate screening functions for j-k pair, 
      // use neighbor list of j to screen j-k pair

      int jnva = va2nvalist[j][jbasis][jucell];
      int *llist = nvafirstneigh[jnva];
      int *lindexlist = nvafirstneighindex[jnva];
      int lnum = nvanumneigh[jnva];

      for (ln = 0; ln < lnum; ln++) {
        l = llist[ln];
        lindex = lindexlist[ln];

        // skip if k and l are the same virtual atom

        if (l == k && lindex == kindex) continue;

        // l is atom

        if (lindex < 0) {
          eltl = fmap[atype[l]];
          if (eltl < 0) continue;
          xltmp = ax[l][0];
          yltmp = ax[l][1];
          zltmp = ax[l][2];
        } 

        // l is virtual atom

        else {
          letype = etype[l];
          lapc = apc[letype];
          lnpe = npe[letype];
          lucell = lindex / lapc;
          lbasis = lindex % lapc;
          eltl = fmap[ctype[l][lbasis]];
          if (eltl < 0) continue;
          lgcell = u2g[letype][lucell];
          if (lgcell >= 0) lnode = g2n[letype][lgcell];
          else lnode = -1;
          if (lnode < 0) {
            xltmp = yltmp = zltmp = 0.0;
            for (node = 0; node < lnpe; node++) {
              xltmp += shape_array[letype][lucell][node] *
                nodex[l][lbasis][node][0];
              yltmp += shape_array[letype][lucell][node] *
                nodex[l][lbasis][node][1];
              zltmp += shape_array[letype][lucell][node] *
                nodex[l][lbasis][node][2];
            }
          } else {
            xltmp = nodex[l][lbasis][lnode][0];
            yltmp = nodex[l][lbasis][lnode][1];
            zltmp = nodex[l][lbasis][lnode][2];
          }
        }
        int flag = compute_screen(xjtmp, yjtmp, zjtmp, xktmp, yktmp, zktmp, xltmp, yltmp, zltmp, 
            eltj, eltk, eltl, rjkbound, rjk2, sjlk, coef1, dCjlk);
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

      double sjktmp = sjk;
      sjk *= fcjk;

      // densities for j and k will be interpolated.

      ind = this->eltind[eltj][eltk];
      pp = rjk * this->rdrar;
      kk = (int)pp;
      kk = std::min(kk, this->nrar - 2);
      pp = pp - kk;
      pp = std::min(pp, 1.0);

      phi = ((this->phirar3[ind][kk] * pp + this->phirar2[ind][kk]) * pp + 
          this->phirar1[ind][kk]) * pp + this->phirar[ind][kk];
      phip = (this->phirar6[ind][kk] * pp + this->phirar5[ind][kk]) * pp + 
        this->phirar4[ind][kk];
      double dUdrjkm[3], dUdsjk;

      copy_densities(kdensities, k, knode, kucell, kbasis);

      compute_pair_force(dUdrjkm, dUdsjk, force, dscrfcnjk, 
          eltj, eltk, scale[jctype][kctype], jdensities, kdensities,
          deljk, rjk, rjk2, phi, phip, sjk);
      xij = rij2 / rjk2;
      xik = rik2 / rjk2;
      a = 1 - (xij - xik) * (xij - xik);
      if (iszero(a)) continue;

      double cjik;
      cjik = (2.0 * (xij + xik) + a - 2.0) / a;
      const double Cmax2 = this->Cmax_meam[eltj][eltk][elti];
      const double Cmin2 = this->Cmin_meam[eltj][eltk][elti];
      if (cjik < Cmin2 || cjik > Cmax2) continue;

      double dsjk1, dsjk2, dCjik1, dCjik2, sjik;
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

      /*
         if (iindex >= 0)
         if (element->tag[i] == TEST && ibasis == TESTBASIS && igcell == TESTGCELL) 
         printf("%g %g %g\n",fij[0]+fik[0],fij[1]+fik[1],fij[2]+fik[2]);
         */
      /*
         if (iindex < 0) {
         if (atom->tag[i] == TEST) {
         printf("Threebody force k nva f = %g %g %g lnum = %d j = %d %d %d nlocal = %d jtag = %d\n",fij[0]+fik[0],fij[1]+fik[1],fij[2]+fik[2],lnum,j,jbasis,jucell,element->nlocal,element->tag[j]);
         char file[100];
         sprintf(file,"dump_neigh_%d_%d_%d.atom",j,jbasis,jucell);
         FILE *fp = fopen(file,"w");
         fprintf(fp,"ITEM: TIMESTEP\n");
         fprintf(fp,"0\n");
         fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
         fprintf(fp,"%d\n",lnum+1);
         fprintf(fp, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
         fprintf(fp, "-4.6805000000000000e+01 4.8805000000000000e+01 0.0000000000000000e+00\n");
         fprintf(fp, "-1.1267903340238979e+03 1.1287903340238979e+03 0.0000000000000000e+00\n");
         fprintf(fp, "-1.0558860256821749e+03 1.0578860256821749e+03 1.3735085464227600e+02\n");

         fprintf(fp, "ITEM: ATOMS id type x y z aid eid ibasis iucell\n" );
         int count=1;
         fprintf(fp, "%d 1 %g %g %g %d %d %d %d\n",count++,xjtmp,yjtmp,zjtmp,0,element->tag[j],jbasis,jucell);
         int laid,leid;
         laid = leid = -1;
         for (ln = 0; ln < lnum; ln++) {
         l = llist[ln];
         lindex = lindexlist[ln];
         if (lindex < 0) {
         xltmp = ax[l][0];
         yltmp = ax[l][1];
         zltmp = ax[l][2];
         lbasis = lucell = -1;
         laid = atom->tag[l];
         } else {
         leid = element->tag[l];
         letype = etype[l];
         lapc = apc[letype];
         lnpe = npe[letype];
         lucell = lindex / lapc;
         lbasis = lindex % lapc;
         lgcell = u2g[letype][lucell];
         if (lgcell >= 0) lnode = g2n[letype][lgcell];
         else lnode = -1;
         if (lnode < 0) {
         xltmp = yltmp = zltmp = 0.0;
         for (node = 0; node < lnpe; node++) {
         xltmp += shape_array[letype][lucell][node] *
         nodex[l][lbasis][node][0];
         yltmp += shape_array[letype][lucell][node] *
         nodex[l][lbasis][node][1];
         zltmp += shape_array[letype][lucell][node] *
         nodex[l][lbasis][node][2];
         }
         } else {
         xltmp = nodex[l][lbasis][lnode][0];
         yltmp = nodex[l][lbasis][lnode][1];
         zltmp = nodex[l][lbasis][lnode][2];
         }
         }

         fprintf(fp, "%d 2 %g %g %g %d %d %d %d\n",count++,xltmp,yltmp,zltmp,laid,leid,lbasis,lucell);
         }

         fclose(fp);
         }
         }
         */

      int ll,kk,llindex,kkindex;
      if (l < k) {
        kk = l;
        ll = k;
        llindex = kindex;
        kkindex = lindex;

      } else if (l > k){
        ll = l;
        kk = k;
        llindex = lindex;
        kkindex = kindex;
      } else {
        ll = kk = l;
        if (kindex < lindex) {
          kkindex = kindex;
          llindex = lindex;

        } else {
          kkindex = lindex;
          llindex = kindex;
        }
      }

      //      printf("%d %d %d %d\n",kk,kindex,ll,lindex);
      if (vflag_either && (iucell < 0 || inode >= 0)) {
        v[0] = -third * (delij[0] * fij[0] + dxik * fik[0]);
        v[1] = -third * (delij[1] * fij[1] + dyik * fik[1]);
        v[2] = -third * (delij[2] * fij[2] + dzik * fik[2]);
        v[3] = -sixth * (delij[0] * fij[1] + dxik * fik[1] + delij[1] * fij[0] + dyik * fik[0]);
        v[4] = -sixth * (delij[0] * fij[2] + dxik * fik[2] + delij[2] * fij[0] + dzik * fik[0]);
        v[5] = -sixth * (delij[1] * fij[2] + dyik * fik[2] + delij[2] * fij[1] + dzik * fik[1]);
        for (m = 0; m < 6; m++) {
          if (ivptr) ivptr[m] += v[m];
          virial[m] += v[m] * ivscale;
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
 ************************************** */

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

  double rho0i, rho1i, rho2i, rho3i, frhopi, dgamma1i, dgamma2i, dgamma3i, arho2bi; 
  double *arho1i, *arho2i, *arho3i, *arho3bi, *t_avei, *tsq_avei;
  double rho0j, rho1j, rho2j, rho3j, frhopj, dgamma1j, dgamma2j, dgamma3j, arho2bj; 
  double *arho1j, *arho2j, *arho3j, *arho3bj, *t_avej, *tsq_avej;
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


void MEAM::copy_densities(double *densities, int i, int inode, int iucell, int ibasis)
{
  int m = 0;
  double tmp;

  // atoms
  if (iucell < 0) {
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
    densities[m++] = noderho0[i][ibasis][inode];
    densities[m++] = noderho1[i][ibasis][inode];
    densities[m++] = noderho2[i][ibasis][inode];
    densities[m++] = noderho3[i][ibasis][inode];
    densities[m++] = nodefrhop[i][ibasis][inode];
    densities[m++] = nodedgamma1[i][ibasis][inode];
    densities[m++] = nodedgamma2[i][ibasis][inode];
    densities[m++] = nodedgamma3[i][ibasis][inode];
    for (int j = 0; j < 3; j++)
      densities[m++] = nodearho1[i][ibasis][inode][j];
    for (int j = 0; j < 6; j++)
      densities[m++] = nodearho2[i][ibasis][inode][j];
    densities[m++] = nodearho2b[i][ibasis][inode];
    for (int j = 0; j < 10; j++)
      densities[m++] = nodearho3[i][ibasis][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodearho3b[i][ibasis][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodet_ave[i][ibasis][inode][j];
    for (int j = 0; j < 3; j++) 
      densities[m++] = nodetsq_ave[i][ibasis][inode][j];

  } 

  // virtual atoms
  else {
    for (int j = 0; j < NDENS; j++)
      densities[j] = 0.0;
    int *npe = element->npe;
    int ietype = element->etype[i];
    double ***shape_array = element->shape_array;

    for (int node = 0; node < npe[ietype]; node++) {
      tmp = shape_array[ietype][iucell][node];
      if (iszero(tmp)) continue; 
      m = 0;
      densities[m++] += tmp * noderho0[i][ibasis][node];
      densities[m++] += tmp * noderho1[i][ibasis][node];
      densities[m++] += tmp * noderho2[i][ibasis][node];
      densities[m++] += tmp * noderho3[i][ibasis][node];
      densities[m++] += tmp * nodefrhop[i][ibasis][node];
      densities[m++] += tmp * nodedgamma1[i][ibasis][node];
      densities[m++] += tmp * nodedgamma2[i][ibasis][node];
      densities[m++] += tmp * nodedgamma3[i][ibasis][node];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp * nodearho1[i][ibasis][node][j];
      for (int j = 0; j < 6; j++)
        densities[m++] += tmp * nodearho2[i][ibasis][node][j];
      densities[m++] += tmp * nodearho2b[i][ibasis][node];
      for (int j = 0; j < 10; j++)
        densities[m++] += tmp * nodearho3[i][ibasis][node][j];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp * nodearho3b[i][ibasis][node][j];
      for (int j = 0; j < 3; j++)
        densities[m++] += tmp * nodet_ave[i][ibasis][node][j];
      for (int j = 0; j < 3; j++) 
        densities[m++] += tmp * nodetsq_ave[i][ibasis][node][j];
    }
  }
}
