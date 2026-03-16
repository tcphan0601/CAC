#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_sw.h"
#include "atom.h"
#include "element.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

#include "update.h"

using namespace CAC_NS;

#define MAXLINE 1024
#define DELTA 4

/*  ----------------------------------------------------------------------  */

PairSW::PairSW(CAC *cac) : Pair(cac)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  threebody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  maxshort = 10;
  neighshort = neighindexshort = NULL;
  rsqshort = NULL;
  delshort = NULL;

  nnvamax = 0;
  shortnva_list = shortnva_flag = NULL;

}

/*  ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
   -------------------------------------------------------------------------  */

PairSW::~PairSW()
{
  if (copymode) return;
  memory->destroy(shortnva_flag); 
  memory->destroy(shortnva_list);

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
    memory->destroy(neighindexshort);
    memory->destroy(rsqshort);
    memory->destroy(delshort);
    delete [] map;
  }
}

/*  ----------------------------------------------------------------------  */

void PairSW::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, kk, inode, jnode, knode, node, inva;
  int jnum, *jlist, *jindexlist;
  int ictype, jctype, kctype, ietype, jetype, ketype;
  int iucell, jucell, kucell, igcell, jgcell, kgcell;
  int ibasis, jbasis, kbasis, iapc, japc, kapc;
  int iindex, jindex, kindex;
  int inna, inpe, jnpe;
  int ijparam, ikparam, ijkparam;
  int numshort;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, evdwl, fpair;
  double iescale, ivscale, jescale, jvscale, kescale, kvscale;
  double fix, fiy, fiz, fjx, fjy, fjz; // force accumulations for i and j
  double rsq, rsq1, rsq2;
  double *delr1, *delr2;
  double fj3[3], fk3[3]; // force component from threebody term
  double *fi, *fj, *fk; // force pointers for i, j, and k
  double *ieptr, *jeptr, *keptr, *ivptr, *jvptr, *kvptr;
  double third = 1.0/3.0;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;
  double ****nodex = element->nodex;
  double ****gaussf = element->gaussf;
  int *ngcell = element->ngcell;
  int *nucell = element->nucell;
  int **g2u = element->g2u;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  int max_ngcell = element->max_ngcell;
  double ***shape_array = element->shape_array;
  double ***weighted_shape_array = element->weighted_shape_array;
  double *nodal_weight = element->nodal_weight;

  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  int maxapc = element->maxapc;
  tagint *atag = atom->tag;
  tagint *etag = element->tag;


  // pointers and parameters from neighbor list

  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  int nvainum = list->nvainum;
  int ***va2nvalist = list->va2nvalist;
  int *nvailist = list->nvailist;
  int *nvaindexlist = list->nvaindexlist;
  int *nvanumneigh = listfull->nvanumneigh;
  int **nvafirstneigh = listfull->nvafirstneigh;
  int **nvafirstneighindex = listfull->nvafirstneighindex;

  // grow shortnva arrays if necessary
  // need to be at least nvainum in length
  // grow to nvainum * 1.5 to reduce growth frequency

  if (element->nelements && nnvamax < nvainum) {
    memory->destroy(shortnva_list);
    memory->destroy(shortnva_flag);
    nnvamax = static_cast<int> (nvainum * 1.5);
    memory->create(shortnva_list, nnvamax, "pair:shortnva_list");
    memory->create(shortnva_flag, nnvamax, "pair:shortnva_flag");
  }

  // initialztmpe shortnva_flag to 0

  int numshortnva = 0;
  int nbytes = sizeof(int) * nnvamax;
  if (nbytes) memset(shortnva_flag, 0, nbytes);

  //-----------------------------------------
  // Outline:
  // 1. Loop through neighbor list to calculate pair component, 
  //    extract short neighborlist, mark virtual atoms in short 
  //    neighbor list that are not gaussian points as short 
  //    nva (neighboring virtual atoms, including ghost), then loop through 
  //    short neighbor list to calculate threebody terms.
  // 2. Loop through neighbor of short nvas from neighbor list, 
  //    extract short neighborlist, skip if there is an atom or gauss point in the short neighborlist, 
  //    loop through short neighborlist to calculate threebody terms and 
  //    accumulate forces on its neighboring atoms/gauss points
  //    NOTE::no tallying forces on ghost
  // ---------------------------------------------


  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];
    igcell = -1;
    ieptr = ivptr = NULL;
    iescale = ivscale = 0.0;
    fix = fiy = fiz = 0.0;

    // i is atom
    
    if (iindex < 0) {
      itag = atag[i];
      ictype = map[atype[i]];
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
      if (eflag_global) iescale = 0.5;
      if (vflag_global) ivscale = 0.5;
      if (eflag_atom) ieptr = &eatom[i];
      if (vflag_atom) ivptr = vatom[i];
      fi = af[i];
    }

    // i is gauss point
    
    else {
      itag = etag[i];
      ietype = etype[i];
      iapc = apc[ietype];
      inpe = npe[ietype];
      igcell = iindex / iapc;
      ibasis = iindex % iapc;
      ictype = map[ctype[i][ibasis]];
      inode = g2n[ietype][igcell];
      fi = gaussf[i][ibasis][igcell];
      if (inode < 0) {
        iucell = g2u[ietype][igcell];
        xtmp = ytmp = ztmp = 0.0;
        for (node = 0; node < npe[ietype]; node++) {
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

    //**************************************
    // two-body interactions
    // skip half
    //**************************************

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    numshort = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      delx = xtmp; dely = ytmp; delz = ztmp;
      jescale = jvscale = 0.0;
      fj = jeptr = jvptr = NULL;

      if (jindex < 0) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
        jtag = atag[j];
        jctype = map[atype[j]];
        fj = af[j];
        if (eflag_atom) jeptr = &eatom[j];
        if (vflag_atom) jvptr = vatom[j];
        if (eflag_global) jescale = 0.5;
        if (vflag_global) jvscale = 0.5;
      } 

      // j is virtual atom

      else {
        jtag = etag[j];
        jetype = etype[j];
        japc = apc[jetype];
        jnpe = npe[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = map[ctype[j][jbasis]];
        for (node = 0; node < jnpe; node++) {
          delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
          dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
          delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
        }
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

      rsq = delx * delx + dely * dely + delz * delz;
      ijparam = elem2param[ictype][jctype][jctype];
      if (rsq >= params[ijparam].cutsq) continue;

      // mark virtual atom j that are not gauss point as short nia for
      // later calculation

      if (jindex >= 0 && jgcell < 0) {
        inna = va2nvalist[j][jbasis][jucell];  
        if (!shortnva_flag[inna]) {
          shortnva_list[numshortnva++] = inva;
          shortnva_flag[inva] = 1;
        }
      }

      if (numshort >= maxshort) {
        maxshort += maxshort/2;
        memory->grow(neighshort, maxshort, "pair:neighshort");
        memory->grow(neighindexshort, maxshort, "pair:neighindexshort");
        memory->grow(rsqshort, maxshort, "pair:rsqshort");
        memory->grow(delshort, maxshort, 3, "pair:delshort");
      }

      neighshort[numshort] = j;
      neighindexshort[numshort] = jindex;
      rsqshort[numshort] = rsq;
      delshort[numshort][0] = -delx;
      delshort[numshort][1] = -dely;
      delshort[numshort][2] = -delz;
      numshort++;

      if (iindex < 0 && jindex < 0) {
        // i and j are distint atoms, only count as one pair
        // odd and even sum of tag condition is to have 
        // the pairs distributed among procs more evenly

        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        }

        // if i and j have the same tag, atom j is a ghost image of atom i
        // only consider the pair from upper right side since
        // force contribution from the otherside will be added in reverse_comm
        // see LAMMPS documents

        else {
          if (delz > 0) continue;
          if (delz == 0 && dely > 0) continue;
          if (delz == 0 && dely == 0 && delx > 0) continue;
        }
      } else if (igcell >= 0 && jgcell >= 0) {
        if ((itag * max_ngcell + igcell) * maxapc + ibasis > 
            (jtag * max_ngcell + jgcell) * maxapc + jbasis) {
          if (((itag * max_ngcell + igcell) * maxapc + ibasis +
               (jtag * max_ngcell + jgcell) * maxapc + jbasis) % 2 == 0) 
            continue;
        } else if ((itag * max_ngcell + igcell) * maxapc + ibasis < 
                   (jtag * max_ngcell + jgcell) * maxapc + jbasis) {
          if (((itag * max_ngcell + igcell) * maxapc + ibasis +
               (jtag * max_ngcell + jgcell) * maxapc + jbasis) % 2 == 1) 
            continue;
        } else {
          if (delz > 0) continue;
          if (delz == 0 && dely > 0) continue;
          if (delz == 0 && dely == 0 && delx > 0) continue;
        }
      } 

      // skip if i is atom and j is gauss point

      else if (iindex < 0 && jgcell >= 0) continue;

      twobody(&params[ijparam], rsq, fpair, eflag, evdwl);

      fix += delx * fpair;
      fiy += dely * fpair;
      fiz += delz * fpair;
      if (fj != NULL) {
        fj[0] -= delx * fpair;
        fj[1] -= dely * fpair;
        fj[2] -= delz * fpair;
      }
      if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
          ieptr, jeptr, ivptr, jvptr, evdwl, 0.0, fpair, delx, dely, delz);
    }

    iescale *= 2.0 * third;
    ivscale *= 2.0 * third;

    //************************************** 
    // three-body interaction 
    //************************************** 

    // loop over pairs j-k in short neighbor list

    for (jj = 0; jj < numshort-1; jj++) {
      j = neighshort[jj];
      jindex = neighindexshort[jj];
      jescale = jvscale = 0.0;
      fj = jeptr = jvptr = NULL;

      // j is atom

      if (jindex < 0) {
        jctype = map[atype[j]];
        fj = af[j];
        if (eflag_atom) jeptr = &eatom[j];
        if (vflag_atom) jvptr = vatom[j];
        if (eflag_global) jescale = third;
        if (vflag_global) jvscale = third;
      }

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = map[ctype[j][jbasis]];
        jgcell = u2g[jetype][jucell];
        if (jgcell >= 0) {
          fj = gaussf[j][jbasis][jgcell];
          jnode = g2n[jetype][jgcell];
          if (jnode >= 0) {
            if (eflag_atom) jeptr = &enode[j][jbasis][jnode];
            if (vflag_atom) jvptr = vnode[j][jbasis][jnode];
            if (eflag_global) jescale = nodal_weight[jetype] * third;
            if (vflag_global) jvscale = nodal_weight[jetype] * third;
          }
        }
      }

      ijparam = elem2param[ictype][jctype][jctype];
      delr1 = delshort[jj];
      rsq1 = rsqshort[jj];
      fjx = fjy = fjz = 0.0;

      // loop over short neighbors

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        kindex = neighindexshort[kk];
        kescale = kvscale = 0.0;
        fk = keptr = kvptr = NULL;
        // k is atom

        if (kindex < 0) {
          kctype = map[atype[k]];
          fk = af[k];
          if (eflag_atom) keptr = &eatom[k];
          if (vflag_atom) kvptr = vatom[k];
          if (eflag_global) kescale = third;
          if (vflag_global) kvscale = third;
        }

        // k is virtual atom

        else {
          ketype = etype[k];
          kapc = apc[ketype];
          kucell = kindex / kapc;
          kbasis = kindex % kapc;
          kctype = map[ctype[k][kbasis]];
          kgcell = u2g[ketype][kucell];
          if (kgcell >= 0) {
            fk = gaussf[k][kbasis][kgcell];
            knode = g2n[ketype][kgcell];
            if (knode >= 0) {
              if (eflag_atom) keptr = &enode[k][kbasis][knode];
              if (vflag_atom) kvptr = vnode[k][kbasis][knode];
              if (eflag_global) kescale = nodal_weight[ketype] * third;
              if (vflag_global) kvscale = nodal_weight[ketype] * third;
            }
          }
        }

        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        delr2 = delshort[kk];
        rsq2 = rsqshort[kk];

        threebody(&params[ijparam], &params[ikparam], &params[ijkparam], 
            rsq1, rsq2, delr1, delr2, fj3, fk3, eflag, evdwl);

        fix -= fj3[0] + fk3[0];
        fiy -= fj3[1] + fk3[1];
        fiz -= fj3[2] + fk3[2];
        if (fj != NULL) {
          fjx += fj3[0];
          fjy += fj3[1];
          fjz += fj3[2];
        }
        if (fk != NULL) {
          fk[0] += fk3[0];
          fk[1] += fk3[1];
          fk[2] += fk3[2];
        }
        if (evflag) ev_tally3(iescale + jescale + kescale, ivscale + jvscale + kvscale, 
            ieptr, jeptr, keptr, ivptr, jvptr, kvptr, evdwl, 0.0, fj3, fk3, delr1, delr2);
      }

      if (fj != NULL) {
        fj[0] += fjx;
        fj[1] += fjy;
        fj[2] += fjz;
      }
    }
  }

  //--------------------------
  // 3. THREEBODY COMPONENTS FROM SHORT NVA
  //--------------------------
  // no force tally for ghost to avoid double count

  // loop over shortnva

  for (ii = 0; ii < numshortnva; ii++) {
    inva = shortnva_list[ii];
    i = nvailist[inva];
    iindex = nvaindexlist[inva];
    ietype = etype[i];
    inpe = npe[ietype];
    iapc = apc[ietype];
    ibasis = iindex % iapc;
    iucell = iindex / iapc;
    ictype = map[ctype[i][ibasis]];

    xtmp = ytmp = ztmp = 0.0;
    for (node = 0; node < inpe; node++) {
      xtmp += shape_array[ietype][iucell][node] * 
        nodex[i][ibasis][node][0];
      ytmp += shape_array[ietype][iucell][node] *
        nodex[i][ibasis][node][1];
      ztmp += shape_array[ietype][iucell][node] *
        nodex[i][ibasis][node][2];
    }

    // loop over its neighbors and extract short neighbor list

    // int continue_flag = 0;

    // loop over neighboring atoms

    jnum = nvanumneigh[inva];
    jlist = nvafirstneigh[inva];
    jindexlist = nvafirstneighindex[inva];
    numshort = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      delx = xtmp; dely = ytmp; delz = ztmp;
      if (jindex < 0) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
        jctype = map[atype[j]];
      } 

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = map[ctype[j][jbasis]];
        for (node = 0; node < jnpe; node++) {
          delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
          dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
          delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
        }
      }

      rsq = delx * delx + dely * dely + delz * delz;
      ijparam = elem2param[ictype][jctype][jctype];
      if (rsq >= params[ijparam].cutsq) continue;

      if (numshort >= maxshort) {
        maxshort += maxshort/2;
        memory->grow(neighshort, maxshort, "pair:neighshort");
        memory->grow(neighindexshort, maxshort, "pair:neighindexshort");
        memory->grow(rsqshort, maxshort, "pair:rsqshort");
        memory->grow(delshort, maxshort, 3, "pair:delshort");
      }

      neighshort[numshort] = j;
      neighindexshort[numshort] = jindex;
      rsqshort[numshort] = rsq;
      delshort[numshort][0] = -delx;
      delshort[numshort][1] = -dely;
      delshort[numshort][2] = -delz;
      numshort++;
    }

    // loop over short neighbor list pair j-k to calculate threebody term
    // only tally force for OWN atoms/gauss points 

    for (jj = 0; jj < numshort-1; jj++) {
      j = neighshort[jj];
      jindex = neighindexshort[jj];
      jescale = jvscale = 0.0;
      fj = jeptr = jvptr = NULL;

      // j is atom

      if (jindex < 0) {
        jctype = map[atype[j]];
        if (j < nalocal) {
          fj = af[j];
          if (eflag_atom) jeptr = &eatom[j];
          if (vflag_atom) jvptr = vatom[j];
          if (eflag_global) jescale = third;
          if (vflag_global) jvscale = third;
        }
      }

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = map[ctype[j][jbasis]];
        if (j < nelocal) {
          jgcell = u2g[jetype][jucell];
          if (jgcell >= 0) {
            fj = gaussf[j][jbasis][jgcell];
            jnode = g2n[jetype][jgcell];
            if (jnode >= 0) {
              if (eflag_atom) jeptr = &enode[j][jbasis][jnode];
              if (vflag_atom) jvptr = vnode[j][jbasis][jnode];
              if (eflag_global) jescale = nodal_weight[jetype] * third;
              if (vflag_global) jvscale = nodal_weight[jetype] * third;
            }
          }
        }
      }

      ijparam = elem2param[ictype][jctype][jctype];
      delr1 = delshort[jj];
      rsq1 = rsqshort[jj];
      fjx = fjy = fjz = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        kindex = neighindexshort[kk];
        kescale = kvscale = 0.0;
        fk = keptr = kvptr = NULL;

        // k is atom

        if (kindex < 0) {
          kctype = map[atype[k]];
          if (k < nalocal) {
            fk = af[k];
            if (eflag_atom) keptr = &eatom[k];
            if (vflag_atom) kvptr = vatom[k];
            if (eflag_global) kescale = third;
            if (vflag_global) kvscale = third;
          }
        }

        // k is virtual atom

        else {
          ketype = etype[k];
          kapc = apc[ketype];
          kucell = kindex / kapc;
          kbasis = kindex % kapc;
          kctype = map[ctype[k][kbasis]];
          if (k < nelocal) {
            kgcell = u2g[ketype][kucell];
            if (kgcell >= 0) {
              fk = gaussf[k][kbasis][kgcell];
              knode = g2n[ketype][kgcell];
              if (knode >= 0) {
                if (eflag_atom) keptr = &enode[k][kbasis][knode];
                if (vflag_atom) kvptr = vnode[k][kbasis][knode];
                if (eflag_global) kescale = nodal_weight[ketype] * third;
                if (vflag_global) kvscale = nodal_weight[ketype] * third;
              }
            }
          }
        }

        if (fj == NULL && fk == NULL) continue;

        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];
        delr2 = delshort[kk];
        rsq2 = rsqshort[kk];
        threebody(&params[ijparam], &params[ikparam], &params[ijkparam], 
            rsq1, rsq2, delr1, delr2, fj3, fk3, eflag, evdwl);

        if (fj != NULL) {
          fjx += fj3[0];
          fjy += fj3[1];
          fjz += fj3[2];
        }
        if (fk != NULL) {
          fk[0] += fk3[0];
          fk[1] += fk3[1];
          fk[2] += fk3[2];
        }
        if (evflag) ev_tally3(iescale + jescale + kescale, ivscale + jvscale + kvscale, 
            ieptr, jeptr, keptr, ivptr, jvptr, kvptr, evdwl, 0.0, fj3, fk3, delr1, delr2);
      }

      if (fj != NULL) {
        fj[0] += fjx;
        fj[1] += fjy;
        fj[2] += fjz;
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/*  ----------------------------------------------------------------------  */

void PairSW::allocate()
{

  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(neighshort, maxshort, "pair:neighshort");
  memory->create(neighindexshort, maxshort, "pair:neighindexshort");
  memory->create(rsqshort, maxshort, "pair:rsqshort");
  memory->create(delshort, maxshort, 3, "pair:delshort");
  map = new int[n+1];
}

/*  ----------------------------------------------------------------------
    global settings
    -------------------------------------------------------------------------  */

void PairSW::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/*  ----------------------------------------------------------------------
    set coeffs for one or more type pairs
    -------------------------------------------------------------------------  */

void PairSW::coeff(int narg, char **arg)
{
  int i, j, n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // insure I, J args are **

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i], elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j], arg[i]);
      nelements++;
    }
  }

  // read potential file and initialztmpe potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I, J = **

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i, j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*  ----------------------------------------------------------------------
    init specific to this pair style
    -------------------------------------------------------------------------  */

void PairSW::init_style()
{
  if (atom->tag_enable == 0 || element->tag_enable == 0)
    error->all(FLERR, "Pair style Stillinger-Weber requires atom/element IDs");
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor list and threebody 

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->threebody = 1;
}

/*  ----------------------------------------------------------------------
    init for one type pair i, j and corresponding j, i
    -------------------------------------------------------------------------  */

double PairSW::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  return cutmax;
}

/*  ----------------------------------------------------------------------  */

void PairSW::read_file(char *file)
{
  int params_per_line = 14;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str, "Cannot open Stillinger-Weber potential file %s", file);
      error->one(FLERR, str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n, nwords, ielement, jelement, kelement;
  char line[MAXLINE], *ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    nwords = universe->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n], MAXLINE-n, fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof, 1, MPI_INT, 0, world);
      if (eof) break;
      MPI_Bcast(&n, 1, MPI_INT, 0, world);
      MPI_Bcast(line, n, MPI_CHAR, 0, world);
      if ((ptr = strchr(line, '#'))) *ptr = '\0';
      nwords = universe->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR, "Incorrect format in Stillinger-Weber potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line, " \t\n\r\f");
    while ((words[nwords++] = strtok(NULL, " \t\n\r\f"))) continue;

    // ielement, jelement, kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0], elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1], elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2], elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), 
          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].epsilon = atof(words[3]);
    params[nparams].sigma = atof(words[4]);
    params[nparams].littlea = atof(words[5]);
    params[nparams].lambda = atof(words[6]);
    params[nparams].gamma = atof(words[7]);
    params[nparams].costheta = atof(words[8]);
    params[nparams].biga = atof(words[9]);
    params[nparams].bigb = atof(words[10]);
    params[nparams].powerp = atof(words[11]);
    params[nparams].powerq = atof(words[12]);
    params[nparams].tol = atof(words[13]);

    if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
        params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
        params[nparams].gamma < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
        params[nparams].powerq < 0.0 || params[nparams].tol < 0.0)
      error->all(FLERR, "Illegal Stillinger-Weber parameter");

    nparams++;
  }

  delete [] words;
}

/*  ----------------------------------------------------------------------  */

void PairSW::setup_params()
{
  int i, j, k, m, n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param, nelements, nelements, nelements, "pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR, "Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR, "Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a * sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].sigma * params[m].littlea;

    rtmp = params[m].cut;
    if (params[m].tol > 0.0) {
      if (params[m].tol > 0.01) params[m].tol = 0.01;
      if (params[m].gamma < 1.0)
        rtmp = rtmp +
          params[m].gamma * params[m].sigma / log(params[m].tol);
      else rtmp = rtmp +
        params[m].sigma / log(params[m].tol);
    }
    params[m].cutsq = rtmp * rtmp;

    params[m].sigma_gamma = params[m].sigma * params[m].gamma;
    params[m].lambda_epsilon = params[m].lambda * params[m].epsilon;
    params[m].lambda_epsilon2 = 2.0 * params[m].lambda * params[m].epsilon;
    params[m].c1 = params[m].biga * params[m].epsilon *
      params[m].powerp * params[m].bigb *
      pow(params[m].sigma, params[m].powerp);
    params[m].c2 = params[m].biga * params[m].epsilon * params[m].powerq *
      pow(params[m].sigma, params[m].powerq);
    params[m].c3 = params[m].biga * params[m].epsilon * params[m].bigb *
      pow(params[m].sigma, params[m].powerp+1.0);
    params[m].c4 = params[m].biga * params[m].epsilon *
      pow(params[m].sigma, params[m].powerq+1.0);
    params[m].c5 = params[m].biga * params[m].epsilon * params[m].bigb *
      pow(params[m].sigma, params[m].powerp);
    params[m].c6 = params[m].biga * params[m].epsilon *
      pow(params[m].sigma, params[m].powerq);
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/*  ----------------------------------------------------------------------  */

void PairSW::twobody(Param *param, double rsq, double &fforce, 
    int eflag, double &eng)
{
  double r, rinvsq, rp, rq, rainv, rainvsq, expsrainv;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r, -param->powerp);
  rq = pow(r, -param->powerq);
  rainv = 1.0 / (r - param->cut);
  rainvsq = rainv * rainv * r;
  expsrainv = exp(param->sigma * rainv);
  fforce = (param->c1 * rp - param->c2 * rq +
      (param->c3 * rp -param->c4 * rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param->c5 * rp - param->c6 * rq) * expsrainv;
}

/*  ----------------------------------------------------------------------  */

void PairSW::threebody(Param *paramij, Param *paramik, Param *paramijk, 
    double rsq1, double rsq2, 
    double *delr1, double *delr2, 
    double *fj, double *fk, int eflag, double &eng)
{
  double r1, rinvsq1, rainv1, gsrainv1, gsrainvsq1, expgsrainv1;
  double r2, rinvsq2, rainv2, gsrainv2, gsrainvsq2, expgsrainv2;
  double rinv12, cs, delcs, delcssq, facexp, facrad, frad1, frad2;
  double facang, facang12, csfacang, csfac1, csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij->cut);
  gsrainv1 = paramij->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1 * rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik->cut);
  gsrainv2 = paramik->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2 * rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1 * r2);
  cs = (delr1[0] * delr2[0] + delr1[1] * delr2[1] + delr1[2] * delr2[2]) * rinv12;
  delcs = cs - paramijk->costheta;
  delcssq = delcs * delcs;

  facexp = expgsrainv1 * expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon * paramik->lambda_epsilon) *
  //          facexp * delcssq;

  facrad = paramijk->lambda_epsilon * facexp * delcssq;
  frad1 = facrad * gsrainvsq1;
  frad2 = facrad * gsrainvsq2;
  facang = paramijk->lambda_epsilon2 * facexp * delcs;
  facang12 = rinv12 * facang;
  csfacang = cs * facang;
  csfac1 = rinvsq1 * csfacang;

  fj[0] = delr1[0] * (frad1+csfac1)-delr2[0] * facang12;
  fj[1] = delr1[1] * (frad1+csfac1)-delr2[1] * facang12;
  fj[2] = delr1[2] * (frad1+csfac1)-delr2[2] * facang12;

  csfac2 = rinvsq2 * csfacang;

  fk[0] = delr2[0] * (frad2+csfac2)-delr1[0] * facang12;
  fk[1] = delr2[1] * (frad2+csfac2)-delr1[1] * facang12;
  fk[2] = delr2[2] * (frad2+csfac2)-delr1[2] * facang12;

  if (eflag) eng = facrad;
}

/*  ----------------------------------------------------------------------
    memory usage of local elem-based arrays
    -------------------------------------------------------------------------  */

double PairSW::memory_usage()
{
  double bytes = nnvamax * 2 * sizeof(int);
  return bytes;
}

