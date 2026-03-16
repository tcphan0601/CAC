#include "npair_halffull_newton_intg.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "my_page.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairHalffullNewtonIntg::NPairHalffullNewtonIntg(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void NPairHalffullNewtonIntg::build(NeighList *list)
{
  int i,j,ii,jj,jnum,ietype,jetype;
  int *a2aneighptr;
  a2aneighptr = NULL;
  int *i2ianeighptr,*i2ia_indexneighptr;
  i2ianeighptr = i2ia_indexneighptr = NULL;
  int *jlist,*jindexlist;;
  int na2a,ni2ia;
  int iintg,iintg_local,iintpl,jintg,jintpl,node;
  double xtmp,ytmp,ztmp,xjtmp,yjtmp,zjtmp;

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nalocal = atom->nlocal;

  int nelocal = element->nlocal;
  int *npe = element->npe;
  int **i2ia = element->i2ia;
  int **ia2i = element->ia2i;
  int *nintg = element->nintg;
  double ***shape_array = element->shape_array;

  // define pointers from neighbor list

  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int **n2ilist = list->n2ilist;
  int *e2ilist = list->e2ilist;
  int *numneigha2a = list->numneigha2a;
  int *numneighi2ia = list->numneighi2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneighi2ia = list->firstneighi2ia;
  int **firstneighi2ia_index = list->firstneighi2ia_index;

  MyPage<int> *a2apage = list->a2apage;
  MyPage<int> *i2iapage = list->i2iapage;
  MyPage<int> *i2ia_indexpage = list->i2ia_indexpage;

  int *ailist_full = list->listfull->ailist;
  int *eilist_full = list->listfull->eilist;
  int *e2ilist_full = list->listfull->e2ilist;
  int *numneigha2a_full = list->listfull->numneigha2a;
  int *numneighi2ia_full = list->listfull->numneighi2ia;
  int **firstneigha2a_full = list->listfull->firstneigha2a;
  int **firstneighi2ia_full = list->listfull->firstneighi2ia;
  int **firstneighi2ia_index_full = list->listfull->firstneighi2ia_index;
  int ainum_full = list->listfull->ainum;
  int einum_full = list->listfull->einum;

  int ainum = 0;
  int einum = 0;
  int iinum = 0;
  // reset all pages

  a2apage->reset();
  i2iapage->reset();
  i2ia_indexpage->reset();

  // loop over parent full atom list

  for (ii = 0; ii < ainum_full; ii++) {
    na2a = 0;
    a2aneighptr = a2apage->vget();

    i = ailist_full[ii];
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];

    // loop over full atom neighbor list

    jlist = firstneigha2a_full[i];
    jnum = numneigha2a_full[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (j < nalocal) {
        if (i > j) continue;
      } else {
        if (ax[j][2] < ztmp) continue;
        if (ax[j][2] == ztmp) {
          if (ax[j][1] < ytmp) continue;
          if (ax[j][1] == ytmp && ax[j][0] < xtmp) continue;
        }
      }
      a2aneighptr[na2a++] = j;
    }

    firstneigha2a[ainum] = a2aneighptr;
    numneigha2a[ainum] = na2a;
    a2apage->vgot(na2a);
    if (a2apage->status())
      error->one(FLERR,"neighbor list overflow, boost neigh_modify one");

    ailist[ainum++] = i;
  }
  list->ainum = ainum;

  // loop over parent full intg list

  for (ii = 0; ii < MIN(einum_full,nelocal); ii++) {
    i = eilist_full[ii];
    eilist[ii] = i;
    ietype = etype[i];
    e2ilist[i] = iinum;
    xtmp = ex[i][0];
    ytmp = ex[i][1];
    ztmp = ex[i][2];
    for (iintg = 0; iintg < nintg[ietype]; iintg++) {
      iintpl = i2ia[ietype][iintg];

      ni2ia = 0;
      i2ianeighptr = i2iapage->vget();
      i2ia_indexneighptr = i2ia_indexpage->vget();

      // loop over full interpolated atom neighbor list

      iintg_local = e2ilist_full[ii] + iintg;
      jlist = firstneighi2ia_full[iintg_local];
      jindexlist = firstneighi2ia_index_full[iintg_local];
      jnum = numneighi2ia_full[iintg_local];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jindexlist[jj];
        jetype = etype[j];
        jintg = ia2i[jetype][jintpl];

        // only count half if j is intergration point
        // use element center coords to check 
        
        if (jintg >= 0)
          if (j < nelocal) {
            if (i > j) continue;
            else if (i == j && iintg > jintg) continue;
          } else {
            if (ex[j][2] < ztmp) continue;
            if (ex[j][2] == ztmp) {
              if (ex[j][1] < ytmp) continue;
              if (ex[j][1] == ytmp && ex[j][0] < xtmp) continue;
            }
          }

        i2ianeighptr[ni2ia] = j;
        i2ia_indexneighptr[ni2ia++] = jintpl;
      }

      firstneighi2ia[iinum] = i2ianeighptr;
      firstneighi2ia_index[iinum] = i2ia_indexneighptr;
      numneighi2ia[iinum] = ni2ia;
      i2iapage->vgot(ni2ia);
      i2ia_indexpage->vgot(ni2ia);
      if (i2iapage->status() ||
          i2ia_indexpage->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

      iinum++;

    }
    einum++;
  }

  list->einum = einum;
  list->iinum = iinum;

}
