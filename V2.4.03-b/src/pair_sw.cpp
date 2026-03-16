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

/* ---------------------------------------------------------------------- */

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

  maxshort_a = 10;
  maxshort_ia = 10;
  neighshort_a = neighshort_ia = neighshort_ia_index = NULL;

  nniamax = 0;
  shortnia_list = shortnia_flag = NULL;

  nemax = 0;
  intgf = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
   ------------------------------------------------------------------------- */

PairSW::~PairSW()
{
  if (copymode) return;
  memory->destroy(intgf);
  memory->destroy(shortnia_flag); 
  memory->destroy(shortnia_list);

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort_a);
    memory->destroy(neighshort_ia);
    memory->destroy(neighshort_ia_index);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairSW::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inode,jnode,knode,node,inia;
  int inpe;
  int jnum,*jlist,*jindexlist;
  int ictype,jctype,kctype,ietype,jetype,ketype;
  int iintpl,jintpl,kintpl,iintg,jintg,kintg,iintg_local;
  int ijparam,ikparam,ijkparam;
  int numshort_a,numshort_ia;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp;
  double delx,dely,delz,evdwl,fpair;
  double fix,fiy,fiz,fjx,fjy,fjz; // force accumulations for i and j
  double rsq,rsq1,rsq2;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double fj[3],fk[3]; // force component from threebody term

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;
  double **ex = element->x;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  int *nintg = element->nintg;
  int *nintpl = element->nintpl;
  int **i2ia = element->i2ia;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int max_nintg = element->max_nintg;
  double ***shape_array = element->shape_array;
  double ***weighted_shape_array = element->weighted_shape_array;

  int nalocal = element->nlocal;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  tagint *atag = atom->tag;
  tagint *etag = element->tag;

  int newton_pair = force->newton_pair;

  // pointers and parameters from neighbor list

  int ainum = list->ainum;
  int einum = list->einum;
  int niainum = list->niainum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int *e2ilist = list->e2ilist;
  int **ia2nialist = list->ia2nialist;
  int *nia2ialist = list->nia2ialist;
  int *nia2ia_indexlist = list->nia2ia_indexlist;

  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighi2a = list->numneighi2a;
  int *numneighi2ia = list->numneighi2ia;
  int *numneighnia2a = list->numneighnia2a;
  int *numneighnia2ia = list->numneighnia2ia;

  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighi2a = list->firstneighi2a;
  int **firstneighi2ia = list->firstneighi2ia;
  int **firstneighi2ia_index = list->firstneighi2ia_index;
  int **firstneighnia2a = list->firstneighnia2a;
  int **firstneighnia2ia = list->firstneighnia2ia;
  int **firstneighnia2ia_index = list->firstneighnia2ia_index;

  // grow intgf array if necessary
  // need to be at least element->nmax in length

  if (element->nelements && element->nmax > nemax) {
    memory->destroy(intgf);
    nemax = element->nmax;
    memory->create(intgf,nemax,max_nintg,3,"pair:intgf");
  }

  // grow shortnia arrays if necessary
  // need to be at least niainum in length
  // grow to niainum*1.5 to reduce growth frequency

  if (element->nelements && nniamax < niainum) {
    memory->destroy(shortnia_list);
    memory->destroy(shortnia_flag);
    nniamax = static_cast<int> (niainum*1.5);
    memory->create(shortnia_list,nniamax,"pair:shortnia_list");
    memory->create(shortnia_flag,nniamax,"pair:shortnia_flag");
  }

  // zero out force

  size_t nbytes = sizeof(double) * (element->nlocal + element->nghost) * max_nintg * 3;
  if (nbytes) memset(&intgf[0][0][0],0,nbytes);

  // initialztmpe shortnia_flag to 0

  int numshortnia = 0;
  nbytes = sizeof(int) * nniamax;
  if (nbytes) memset(shortnia_flag,0,nbytes);

  //-----------------------------------------
  // Outline:
  // 1. Loop through atom's neighbor to calculate pair component,
  //    extract short neighborlist, mark interpolated atoms in short 
  //    neighbor list that are not integration points as short 
  //    nia (neighboring interpolated atoms, including ghost), then loop through 
  //    short neighbor list of atoms to calculate threebody terms.
  // 2. Loop through integration's neighbor to calculate pair component,
  //    extract short neighborlist, mark interpolated atoms in short 
  //    neighbor list that are not integration points as 
  //    short nia (neighboring interpolated atoms), then loop through
  //    short neighbor list of integration points to calculate threebody terms.
  // 3. Loop through neighbor of short nias from neighbor list,
  //    extract short neighborlist, skip if there is an atom or intg in the short neighborlist,
  //    loop through short neighborlist to calculate threebody terms and 
  //    accumulate forces on its neighboring atoms/integration points
  // ---------------------------------------------

  //--------------------------
  // 1. FORCES ON ATOM
  //--------------------------

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];
    itag = atag[i];

    ictype = map[atype[i]];
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];
    fix = fiy = fiz = 0.0;

    // **************************************
    // two-body interactions
    // **************************************

    // loop over neighboring atoms

    jnum = numneigha2a[i];
    numshort_a = 0;
    if (jnum) {
      jlist = firstneigha2a[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        //j &= NEIGHMASK;

        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        jctype = map[atype[j]];
        ijparam = elem2param[ictype][jctype][jctype];
        if (rsq >= params[ijparam].cutsq) continue;

        if (numshort_a >= maxshort_a) {
          maxshort_a += maxshort_a/2;
          memory->grow(neighshort_a,maxshort_a,"pair:neighshort_a");
        }

        neighshort_a[numshort_a++] = j;
        jtag = atag[j];

        // i and j are distint atoms, only count as one pair
        // odd and even sum of tag condition is to have 
        // the pairs distributed among procs more evenly

        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {

          // if i and j have the same tag, atom j is a ghost image of atom i
          // only consider the pair from upper right side since
          // force contribution from the otherside will be added in reverse_comm
          // see LAMMPS documents

          if (ax[j][2] < ztmp) continue;
          if (ax[j][2] == ztmp && ax[j][1] < ytmp) continue;
          if (ax[j][2] == ztmp && ax[j][1] == ytmp && ax[j][0] < xtmp) continue;
        }

        twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

        fix += delx*fpair;
        fiy += dely*fpair;
        fiz += delz*fpair;
        af[j][0] -= delx*fpair;
        af[j][1] -= dely*fpair;
        af[j][2] -= delz*fpair;

        if (evflag) atom_ev_tally(i,j,nalocal,newton_pair,
            evdwl,0.0,fpair,delx,dely,delz);
      }
    }

    // loop over neighboring interpolated atoms

    jnum = numneigha2ia[i];
    numshort_ia = 0;
    if (jnum) {
      jlist = firstneigha2ia[i];
      jindexlist = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jindexlist[jj];
        jctype = map[ctype[j]];
        jetype = etype[j];
        j &= NEIGHMASK;

        // interpolate positions at site J from nodex

        delx = xtmp; dely = ytmp; delz = ztmp;
        for (node = 0; node < npe[jetype]; node++) {
          delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
          dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
          delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
        }

        rsq = delx*delx + dely*dely + delz*delz;

        ijparam = elem2param[ictype][jctype][jctype];

        if (rsq >= params[ijparam].cutsq) continue;

        if (numshort_ia >= maxshort_ia) {
          maxshort_ia += maxshort_ia/2;
          memory->grow(neighshort_ia,maxshort_ia,"pair:neighshort_ia");
          memory->grow(neighshort_ia_index,maxshort_ia,"pair:neighshort_ia_index");
        }

        neighshort_ia[numshort_ia] = j;
        neighshort_ia_index[numshort_ia++] = jintpl;

        // skip if j is an integration point
        // this will be added when calculating force 
        // on j integration point

        if (ia2i[jetype][jintpl] >= 0) continue;

        // else, mark it as short nia and calculate pair force for atom

        inia = ia2nialist[j][jintpl];  
        if (!shortnia_flag[inia]) {
          shortnia_list[numshortnia++] = inia;
          shortnia_flag[inia] = 1;
        }

        twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

        fix += delx*fpair;
        fiy += dely*fpair;
        fiz += delz*fpair;
        if (evflag) atom_ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }

    // ************************************** 
    // three-body interaction (3 types of triples: a-a-a,a-a-ia,a-ia-ia)
    // ************************************** 

    // loop over short neighboring atoms 
    // include the last atom for a-ia pairs
    // last one will be ignored for a-a pairs

    for (jj = 0; jj < numshort_a; jj++) {
      j = neighshort_a[jj];
      jctype = map[atype[j]];
      ijparam = elem2param[ictype][jctype][jctype];
      delx1 = ax[j][0] - xtmp;
      dely1 = ax[j][1] - ytmp;
      delz1 = ax[j][2] - ztmp;
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      fjx = fjy = fjz = 0.0;

      // loop over short neighboring atoms (a-a-a triples, same as LAMMPS)

      for (kk = jj+1; kk < numshort_a; kk++) {
        k = neighshort_a[kk];
        kctype = map[atype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        delx2 = ax[k][0] - xtmp;
        dely2 = ax[k][1] - ytmp;
        delz2 = ax[k][2] - ztmp;
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        fix -= fj[0] + fk[0];
        fiy -= fj[1] + fk[1];
        fiz -= fj[2] + fk[2];
        fjx += fj[0];
        fjy += fj[1];
        fjz += fj[2];
        af[k][0] += fk[0];
        af[k][1] += fk[1];
        af[k][2] += fk[2];
        if (evflag) atom_ev_tally3(i,j,k,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
      }


      // loop over short neighboring interpolated atoms (a-a-ia triples)

      for (kk = 0; kk < numshort_ia; kk++) {
        k = neighshort_ia[kk];
        kintpl = neighshort_ia_index[kk];
        ketype = etype[k];
        kintg = ia2i[ketype][kintpl];
        kctype = map[ctype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        // del2 = kx - xtmp

        delx2 = -xtmp;
        dely2 = -ytmp;
        delz2 = -ztmp;
        for (node = 0; node < npe[ketype]; node++) {
          delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
          dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
          delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
        }
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        // threebody force contribution from i centered triple

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        fix -= fj[0] + fk[0];
        fiy -= fj[1] + fk[1];
        fiz -= fj[2] + fk[2];
        fjx += fj[0];
        fjy += fj[1];
        fjz += fj[2];

        // if k is an integration point, accumulate force

        if (kintg >= 0) {
          intgf[k][kintg][0] += fk[0];
          intgf[k][kintg][1] += fk[1];
          intgf[k][kintg][2] += fk[2];
        }
        if (evflag) {
          atom_ev_tally3(i,j,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          if (kintg >= 0) {
            knode = i2n[ketype][kintg];
            if (knode >= 0) 
              node_ev_tally3(0,-1,0,-1,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          }
        }

      }
      af[j][0] += fjx;
      af[j][1] += fjy;
      af[j][2] += fjz;
    }

    // loop over short neighboring interpolated atoms (a-ia-ia triple only)

    for (jj = 0; jj < numshort_ia-1; jj++) {
      j = neighshort_ia[jj];
      jintpl = neighshort_ia_index[jj];
      jetype = etype[j];
      jintg = ia2i[jetype][jintpl];
      jctype = map[ctype[j]];
      ijparam = elem2param[ictype][jctype][jctype];
      fjx = fjy = fjz = 0.0;

      delx1 = -xtmp;
      dely1 = -ytmp;
      delz1 = -ztmp;
      for (node = 0; node < npe[jetype]; node++) {
        delx1 += shape_array[jetype][jintpl][node]*nodex[j][node][0];
        dely1 += shape_array[jetype][jintpl][node]*nodex[j][node][1];
        delz1 += shape_array[jetype][jintpl][node]*nodex[j][node][2];
      }
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;

      // loop over short neighboring interpolated atoms (ia-ia pairs)

      for (kk = jj+1; kk < numshort_ia; kk++) {
        k = neighshort_ia[kk];
        kintpl = neighshort_ia_index[kk];
        ketype = etype[k];
        kintg = ia2i[ketype][kintpl];
        kctype = map[ctype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        // del2 = kx - xtmp

        delx2 = -xtmp;
        dely2 = -ytmp;
        delz2 = -ztmp;
        for (node = 0; node < npe[ketype]; node++) {
          delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
          dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
          delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
        }
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        fix -= fj[0] + fk[0];
        fiy -= fj[1] + fk[1];
        fiz -= fj[2] + fk[2];

        // if j or k are integration points, accumulate force

        if (jintg >= 0) {
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];
        } 

        if (kintg >= 0) {
          intgf[k][kintg][0] += fk[0];
          intgf[k][kintg][1] += fk[1];
          intgf[k][kintg][2] += fk[2];
        }

        if (evflag) {
          atom_ev_tally3(i,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          if (jintg >= 0) jnode = i2n[jetype][jintg];
          else jnode = -1;
          if (kintg >= 0) knode = i2n[ketype][kintg];
          else knode = -1;
          if (jnode >= 0 || knode >= 0) node_ev_tally3(0,-1,j,jnode,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);

        }
      }

      // accumulate force if j is an integration point

      if (jintg >= 0) {
        intgf[j][jintg][0] += fjx;
        intgf[j][jintg][1] += fjy;
        intgf[j][jintg][2] += fjz;
      }
    }
    af[i][0] += fix;
    af[i][1] += fiy;
    af[i][2] += fiz;
  }

  //--------------------------
  // 2. FORCES ON INTEGRATION POINTS
  //--------------------------

  // loop over full neighbor list of my integration points

  for (ii = 0; ii < nelocal; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    ictype = map[ctype[i]];
    inpe = npe[ietype];

    for (iintg = 0; iintg < nintg[ietype]; iintg++) {

      fix = fiy = fiz = 0.0;

      // interpolate postions of integration point from node values

      iintpl = i2ia[ietype][iintg];
      inode = i2n[ietype][iintg];
      iintg_local = e2ilist[ii] + iintg;

      xtmp = ytmp = ztmp = 0.0;
      for (node = 0; node < inpe; node++) {
        xtmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
        ytmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
        ztmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
      }

      // **************************************
      // two-body interactions
      // **************************************

      // loop over neighboring atoms

      jnum = numneighi2a[iintg_local];
      numshort_a = 0;
      if (jnum) {
        jlist = firstneighi2a[iintg_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          delx = xtmp - ax[j][0];
          dely = ytmp - ax[j][1];
          delz = ztmp - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          jctype = map[atype[j]];
          ijparam = elem2param[ictype][jctype][jctype];
          if (rsq >= params[ijparam].cutsq) continue;

          if (numshort_a >= maxshort_a) {
            maxshort_a += maxshort_a/2;
            memory->grow(neighshort_a,maxshort_a,"pair:neighshort_a");
          }

          neighshort_a[numshort_a++] = j;

          twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

          fix += delx*fpair;
          fiy += dely*fpair;
          fiz += delz*fpair;
          af[j][0] -= delx*fpair;
          af[j][1] -= dely*fpair;
          af[j][2] -= delz*fpair;

          if (evflag) {
            atom_ev_tally_full(j,evdwl,0.0,fpair,delx,dely,delz);
            node_ev_tally(i,inode,evdwl,0.0,fpair,delx,dely,delz);
          }
        }
      }

      // loop over neighboring interpolated atoms

      jnum = numneighi2ia[iintg_local];
      numshort_ia = 0;
      if (jnum) {
        jlist = firstneighi2ia[iintg_local];
        jindexlist = firstneighi2ia_index[iintg_local];
        for (jj = 0; jj < jnum; jj++) {

          j = jlist[jj];
          jintpl = jindexlist[jj];
          jctype = map[ctype[j]];
          jetype = etype[j];

          // interpolate positions at site J from nodex

          delx = xtmp; dely = ytmp; delz = ztmp;
          for (node = 0; node < npe[jetype]; node++) {
            delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
            dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
            delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
          }

          rsq = delx*delx + dely*dely + delz*delz;

          ijparam = elem2param[ictype][jctype][jctype];
          if (rsq >= params[ijparam].cutsq) continue;

          if (numshort_ia >= maxshort_ia) {
            maxshort_ia += maxshort_ia/2;
            memory->grow(neighshort_ia,maxshort_ia,"pair:neighshort_ia");
            memory->grow(neighshort_ia_index,maxshort_ia,"pair:neighshort_ia_index");
          }
          neighshort_ia[numshort_ia] = j;
          neighshort_ia_index[numshort_ia] = jintpl;
          numshort_ia++;

          // if j is not an integration point, mark it as short nia

          if (ia2i[jetype][jintpl] < 0) {
            inia = ia2nialist[j][jintpl];  
            if (!shortnia_flag[inia]) {
              shortnia_list[numshortnia++] = inia;
              shortnia_flag[inia] = 1;
            }
          }

          twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

          fix += delx*fpair;
          fiy += dely*fpair;
          fiz += delz*fpair;

          if (evflag) node_ev_tally(i,inode,evdwl,0.0,fpair,delx,dely,delz);
        }
      }

      // ************************************** 
      // three-body interaction (3 types of triples: ia-a-a, ia-a-ia, ia-ia-ia)
      // ************************************** 

      // loop over short neighboring atoms 
      // include the last atom for ia-a-ia triples
      // last one will be ignored for ia-a-a triples

      for (jj = 0; jj < numshort_a; jj++) {
        j = neighshort_a[jj];
        jctype = map[atype[j]];
        ijparam = elem2param[ictype][jctype][jctype];
        delx1 = ax[j][0] - xtmp;
        dely1 = ax[j][1] - ytmp;
        delz1 = ax[j][2] - ztmp;
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
        fjx = fjy = fjz = 0.0;

        // loop over short neighboring atoms (ia-a-a triple)

        for (kk = jj+1; kk < numshort_a; kk++) {
          k = neighshort_a[kk];
          kctype = map[atype[k]];
          ikparam = elem2param[ictype][kctype][kctype];
          ijkparam = elem2param[ictype][jctype][kctype];

          delx2 = ax[k][0] - xtmp;
          dely2 = ax[k][1] - ytmp;
          delz2 = ax[k][2] - ztmp;
          rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

          threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
              rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

          fix -= fj[0] + fk[0];
          fiy -= fj[1] + fk[1];
          fiz -= fj[2] + fk[2];
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];
          af[k][0] += fk[0];
          af[k][1] += fk[1];
          af[k][2] += fk[2];
          if (evflag) {
            atom_ev_tally3(j,k,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
            node_ev_tally3(i,inode,0,-1,0,-1,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          }
        }

        // loop over short neighboring interpolated atoms (ia-a-ia triple)

        for (kk = 0; kk < numshort_ia; kk++) {
          k = neighshort_ia[kk];
          kintpl = neighshort_ia_index[kk];
          ketype = etype[k];
          kintg = ia2i[ketype][kintpl];
          kctype = map[ctype[k]];
          ikparam = elem2param[ictype][kctype][kctype];
          ijkparam = elem2param[ictype][jctype][kctype];

          // del2 = kx - xtmp

          delx2 = -xtmp;
          dely2 = -ytmp;
          delz2 = -ztmp;
          for (node = 0; node < npe[ketype]; node++) {
            delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
            dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
            delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
          }
          rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

          threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
              rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

          fix -= fj[0] + fk[0];
          fiy -= fj[1] + fk[1];
          fiz -= fj[2] + fk[2];
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];

          // if k is an integration point, accumulate force

          if (kintg >= 0) {
            intgf[k][kintg][0] += fk[0];
            intgf[k][kintg][1] += fk[1];
            intgf[k][kintg][2] += fk[2];
          }

          if (evflag) {
            atom_ev_tally3(j,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
            if (kintg >= 0) knode = i2n[ketype][kintg];
            else knode = -1;
            node_ev_tally3(i,inode,0,-1,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          }
        }
        af[j][0] += fjx;
        af[j][1] += fjy;
        af[j][2] += fjz;
      }

      // loop over short neighboring interpolated atoms 

      for (jj = 0; jj < numshort_ia-1; jj++) {
        j = neighshort_ia[jj];
        jintpl = neighshort_ia_index[jj];
        jetype = etype[j];
        jintg = ia2i[jetype][jintpl];
        jctype = map[ctype[j]];
        ijparam = elem2param[ictype][jctype][jctype];
        fjx = fjy = fjz = 0.0;

        // del1 = jx - xtmp

        delx1 = -xtmp;
        dely1 = -ytmp;
        delz1 = -ztmp;
        for (node = 0; node < npe[jetype]; node++) {
          delx1 += shape_array[jetype][jintpl][node]*nodex[j][node][0];
          dely1 += shape_array[jetype][jintpl][node]*nodex[j][node][1];
          delz1 += shape_array[jetype][jintpl][node]*nodex[j][node][2];
        }
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;

        // loop over short neighboring interpolated atoms (ia-ia pairs)

        for (kk = jj+1; kk < numshort_ia; kk++) {
          k = neighshort_ia[kk];
          kintpl = neighshort_ia_index[kk];
          ketype = etype[k];
          kintg = ia2i[ketype][kintpl];
          kctype = map[ctype[k]];
          ikparam = elem2param[ictype][kctype][kctype];
          ijkparam = elem2param[ictype][jctype][kctype];

          // del2 = kx - xtmp

          delx2 = -xtmp;
          dely2 = -ytmp;
          delz2 = -ztmp;
          for (node = 0; node < npe[ketype]; node++) {
            delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
            dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
            delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
          }

          rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

          threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
              rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

          fix -= fj[0] + fk[0];
          fiy -= fj[1] + fk[1];
          fiz -= fj[2] + fk[2];

          // if j and/or k are integration points, accumulate force

          if (jintg >= 0) {
            fjx += fj[0];
            fjy += fj[1];
            fjz += fj[2];
          }
          if (kintg >= 0) {
            intgf[k][kintg][0] += fk[0];
            intgf[k][kintg][1] += fk[1];
            intgf[k][kintg][2] += fk[2];
          }

          if (evflag) {
            if (jintg >= 0) jnode = i2n[jetype][jintg];
            else jnode = -1;
            if (kintg >= 0) knode = i2n[ketype][kintg];
            else knode = -1;
            node_ev_tally3(i,inode,j,jnode,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          }
        }
        if (jintg >= 0) {
          intgf[j][jintg][0] += fjx;
          intgf[j][jintg][1] += fjy;
          intgf[j][jintg][2] += fjz;
        }
      }
      intgf[i][iintg][0] += fix;
      intgf[i][iintg][1] += fiy;
      intgf[i][iintg][2] += fiz;
    }
  }

  //--------------------------
  // 3. THREEBODY COMPONENTS FROM SHORT NIA
  //--------------------------
  // no accumulate forces for ghost atoms/elements to avoid double count

  // loop over shortnia
  
  for (ii = 0; ii < numshortnia; ii++) {
    inia = shortnia_list[ii];
    i = nia2ialist[inia];
    iintpl = nia2ia_indexlist[inia];
    ietype = etype[i];
    ictype = map[ctype[i]];

    xtmp = ytmp = ztmp = 0.0;
    for (node = 0; node < npe[ietype]; node++) {
      xtmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
      ytmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
      ztmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
    }

    // loop over its neighbors and extract short neighbor list
    // if I is ghost and has ghost atom/intg as short neighbor, then skip since it will be looped through in other procs

    // int continue_flag = 0;

    // loop over neighboring atoms

    jnum = numneighnia2a[inia];
    numshort_a = 0;

    if (jnum) {
      jlist = firstneighnia2a[inia];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        jctype = map[atype[j]];
        ijparam = elem2param[ictype][jctype][jctype];
        if (rsq >= params[ijparam].cutsq) continue;

        //if (i >= nelocal && j >= nalocal) {
        //  continue_flag = 1;
        //  break;
        //}

        if (numshort_a >= maxshort_a) {
          maxshort_a += maxshort_a/2;
          memory->grow(neighshort_a,maxshort_a,"pair:neighshort_a");
        }
        neighshort_a[numshort_a++] = j;
      }
    }

    //if (continue_flag) continue;

    // loop over neighboring interpolated atoms

    jnum = numneighnia2ia[inia];
    numshort_ia = 0;
    if (jnum) {
      jlist = firstneighnia2ia[inia];
      jindexlist = firstneighnia2ia_index[inia];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jindexlist[jj];
        jctype = map[ctype[j]];
        jetype = etype[j];
        jintg = ia2i[jetype][jintpl];

        // interpolate positions at site J from nodex

        delx = xtmp; dely = ytmp; delz = ztmp;
        for (node = 0; node < npe[jetype]; node++) {
          delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
          dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
          delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
        }

        rsq = delx*delx + dely*dely + delz*delz;

        ijparam = elem2param[ictype][jctype][jctype];
        if (rsq >= params[ijparam].cutsq) continue;

        //if (i >= nelocal && j >= nelocal && jintg >= 0) {
        //  continue_flag = 1;
        //  break;
        //}

        if (numshort_ia >= maxshort_ia) {
          maxshort_ia += maxshort_ia/2;
          memory->grow(neighshort_ia,maxshort_ia,"pair:neighshort_ia");
          memory->grow(neighshort_ia_index,maxshort_ia,"pair:neighshort_ia_index");
        }
        neighshort_ia[numshort_ia] = j;
        neighshort_ia_index[numshort_ia++] = jintpl;
      }
    }

    //if (continue_flag) continue;

    // loop over short neighbor list to calculate threebody term

    for (jj = 0; jj < numshort_a; jj++) {
      j = neighshort_a[jj];
      jctype = map[atype[j]];
      ijparam = elem2param[ictype][jctype][jctype];
      delx1 = ax[j][0] - xtmp;
      dely1 = ax[j][1] - ytmp;
      delz1 = ax[j][2] - ztmp;
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      fjx = fjy = fjz = 0.0;

      // loop over short neighboring atoms (ia-a-a triple)

      for (kk = jj+1; kk < numshort_a; kk++) {
        k = neighshort_a[kk];
        if (j >= nalocal && k >= nalocal) continue;
        kctype = map[atype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];
        delx2 = ax[k][0] - xtmp;
        dely2 = ax[k][1] - ytmp;
        delz2 = ax[k][2] - ztmp;
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        if (j < nalocal) {
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];
        }
        if (k < nalocal) {
          af[k][0] += fk[0];
          af[k][1] += fk[1];
          af[k][2] += fk[2];
        }

        if (evflag) atom_ev_tally3(j,k,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
      }

      // loop over short neighboring interpolated atoms (ia-a-ia triple)

      for (kk = 0; kk < numshort_ia; kk++) {
        k = neighshort_ia[kk];
        kintpl = neighshort_ia_index[kk];
        ketype = etype[k];
        kintg = ia2i[ketype][kintpl];
        if (j >= nalocal && (kintg < 0 || k >= nelocal)) continue;
        kctype = map[ctype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        // del2 = kx - xtmp

        delx2 = -xtmp;
        dely2 = -ytmp;
        delz2 = -ztmp;
        for (node = 0; node < npe[ketype]; node++) {
          delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
          dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
          delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
        }
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        if (j < nalocal) {
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];
        } 

        // if k is an integration point, accumulate force

        if (kintg >= 0 && k < nelocal) {
          intgf[k][kintg][0] += fk[0];
          intgf[k][kintg][1] += fk[1];
          intgf[k][kintg][2] += fk[2];
        }

        if (evflag) {
          atom_ev_tally3(j,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          if (kintg >= 0) {
            knode = i2n[ketype][kintg];
            node_ev_tally3(0,-1,0,-1,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
          }
        }
      }
      af[j][0] += fjx;
      af[j][1] += fjy;
      af[j][2] += fjz;
    }

    // loop over short neighboring interpolated atoms 

    for (jj = 0; jj < numshort_ia-1; jj++) {
      j = neighshort_ia[jj];
      jintpl = neighshort_ia_index[jj];
      jetype = etype[j];
      jintg = ia2i[jetype][jintpl];
      jctype = map[ctype[j]];
      ijparam = elem2param[ictype][jctype][jctype];
      fjx = fjy = fjz = 0.0;

      // del1 = jx - xtmp

      delx1 = -xtmp;
      dely1 = -ytmp;
      delz1 = -ztmp;
      for (node = 0; node < npe[jetype]; node++) {
        delx1 += shape_array[jetype][jintpl][node]*nodex[j][node][0];
        dely1 += shape_array[jetype][jintpl][node]*nodex[j][node][1];
        delz1 += shape_array[jetype][jintpl][node]*nodex[j][node][2];
      }

      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;

      // loop over short neighboring interpolated atoms (ia-ia pairs)

      for (kk = jj+1; kk < numshort_ia; kk++) {
        k = neighshort_ia[kk];
        kintpl = neighshort_ia_index[kk];
        ketype = etype[k];
        kintg = ia2i[ketype][kintpl];

        // skip if both j and k are not integration points

        if ((jintg < 0 || j >= nelocal) && (kintg < 0 || k >= nelocal)) continue;

        kctype = map[ctype[k]];
        ikparam = elem2param[ictype][kctype][kctype];
        ijkparam = elem2param[ictype][jctype][kctype];

        // del2 = kx - xtmp

        delx2 = -xtmp;
        dely2 = -ytmp;
        delz2 = -ztmp;
        for (node = 0; node < npe[ketype]; node++) {
          delx2 += shape_array[ketype][kintpl][node]*nodex[k][node][0];
          dely2 += shape_array[ketype][kintpl][node]*nodex[k][node][1];
          delz2 += shape_array[ketype][kintpl][node]*nodex[k][node][2];
        }

        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
            rsq1,rsq2,delx1,dely1,delz1,delx2,dely2,delz2,fj,fk,eflag,evdwl);

        if (jintg >= 0 && j < nelocal) {
          fjx += fj[0];
          fjy += fj[1];
          fjz += fj[2];
        }
        if (kintg >= 0 && k < nelocal) {
          intgf[k][kintg][0] += fk[0];
          intgf[k][kintg][1] += fk[1];
          intgf[k][kintg][2] += fk[2];
        }
        if (evflag) {
          if (jintg >= 0) jnode = i2n[jetype][jintg];
          else jnode = -1;
          if (kintg >= 0) knode = i2n[ketype][kintg];
          else knode = -1;
          node_ev_tally3(0,-1,j,jnode,k,knode,evdwl,0.0,fj,fk,delx1,dely1,delz1,delx2,dely2,delz2);
        }
      }

      if (jintg >= 0) {
        intgf[j][jintg][0] += fjx;
        intgf[j][jintg][1] += fjy;
        intgf[j][jintg][2] += fjz;
      }
    }
  }

  // distribute forces at integration point to forces at nodes 
  // for all elements (including ghost for reverse comm)

  for (i = 0; i < neall; i++) {
    ietype = etype[i];
    inpe = npe[ietype];
    for (iintg = 0; iintg < nintg[ietype]; iintg++) 
      for (node = 0; node < inpe; node++) {
        nodef[i][node][0] += weighted_shape_array[ietype][iintg][node]*intgf[i][iintg][0];
        nodef[i][node][1] += weighted_shape_array[ietype][iintg][node]*intgf[i][iintg][1];
        nodef[i][node][2] += weighted_shape_array[ietype][iintg][node]*intgf[i][iintg][2];
      }
  }
  if (vflag_fdotr) virial_fdotr_compute();

}

/* ---------------------------------------------------------------------- */

void PairSW::allocate()
{

  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort_a,maxshort_a,"pair:neighshort_a");
  memory->create(neighshort_ia,maxshort_ia,"pair:neighshort_ia");
  memory->create(neighshort_ia_index,maxshort_ia,"pair:neighshort_ia_index");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSW::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSW::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

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
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialztmpe potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairSW::init_style()
{
  if (atom->tag_enable == 0 || element->tag_enable == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires atom/element IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor listand threebody 

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->threebody = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

double PairSW::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

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
      sprintf(str,"Cannot open Stillinger-Weber potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = universe->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = universe->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in Stillinger-Weber potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
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
      error->all(FLERR,"Illegal Stillinger-Weber parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairSW::setup_params()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].sigma*params[m].littlea;

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

    params[m].sigma_gamma = params[m].sigma*params[m].gamma;
    params[m].lambda_epsilon = params[m].lambda*params[m].epsilon;
    params[m].lambda_epsilon2 = 2.0*params[m].lambda*params[m].epsilon;
    params[m].c1 = params[m].biga*params[m].epsilon *
      params[m].powerp*params[m].bigb *
      pow(params[m].sigma,params[m].powerp);
    params[m].c2 = params[m].biga*params[m].epsilon*params[m].powerq *
      pow(params[m].sigma,params[m].powerq);
    params[m].c3 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma,params[m].powerp+1.0);
    params[m].c4 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma,params[m].powerq+1.0);
    params[m].c5 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma,params[m].powerp);
    params[m].c6 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma,params[m].powerq);
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairSW::twobody(Param *param, double rsq, double &fforce,
    int eflag, double &eng)
{
  double r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r,-param->powerp);
  rq = pow(r,-param->powerq);
  rainv = 1.0 / (r - param->cut);
  rainvsq = rainv*rainv*r;
  expsrainv = exp(param->sigma * rainv);
  fforce = (param->c1*rp - param->c2*rq +
      (param->c3*rp -param->c4*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param->c5*rp - param->c6*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

void PairSW::threebody(Param *paramij, Param *paramik, Param *paramijk,
    double rsq1, double rsq2, double delx1, double dely1, 
    double delz1, double delx2, double dely2, double delz2,
    double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
  double facang,facang12,csfacang,csfac1,csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij->cut);
  gsrainv1 = paramij->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik->cut);
  gsrainv2 = paramik->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delx1*delx2 + dely1*dely2 + delz1*delz2) * rinv12;
  delcs = cs - paramijk->costheta;
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
  //          facexp*delcssq;

  facrad = paramijk->lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk->lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delx1*(frad1+csfac1)-delx2*facang12;
  fj[1] = dely1*(frad1+csfac1)-dely2*facang12;
  fj[2] = delz1*(frad1+csfac1)-delz2*facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delx2*(frad2+csfac2)-delx1*facang12;
  fk[1] = dely2*(frad2+csfac2)-dely1*facang12;
  fk[2] = delz2*(frad2+csfac2)-delz1*facang12;

  if (eflag) eng = facrad;
}

/* ----------------------------------------------------------------------
   memory usage of local elem-based arrays
   ------------------------------------------------------------------------- */

double PairSW::memory_usage()
{
  double bytes = nemax * element->max_nintg * 3 * sizeof(double);
  bytes += nniamax * 2 * sizeof(int);
  return bytes;
}
