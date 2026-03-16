#include "npair_halffull_newton.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "my_page.h"
#include "error.h"

using namespace CAC_NS;
//#define TEST 2462747
//#define TEST2 2463182

/*  ----------------------------------------------------------------------  */

NPairHalffullNewton::NPairHalffullNewton(CAC *cac) : NPair(cac) {}

/*  ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i, j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
-------------------------------------------------------------------------  */

void NPairHalffullNewton::build(NeighList *list)
{
  int i, ii, j, jj, jnum, k, l, n, ietype, jetype, iindex, jindex, iapc, inpe;
  int *neighptr = NULL;
  int *neighindexptr = NULL;
  int *jlist, *jindexlist;
  int inode, igcell, ibasis, iucell, jucell, jgcell, jbasis, japc, jnode;
  double xtmp, ytmp, ztmp;

  double **ax = atom->x;
  double **ex = element->x;
  double ****nodex = element->nodex;
  int *etype = element->etype;
  int nalocal = atom->nlocal;

  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  int **u2g = element->u2g;
  int **g2n = element->g2n;

  // define pointers from neighbor list

  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  MyPage<int> *ipage = list->ipage;
  MyPage<int> *indexpage = list->indexpage;

  int *ilist_full = list->listfull->ilist;
  int *iindexlist_full = list->listfull->iindexlist;

  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int **firstneighindex_full = list->listfull->firstneighindex;


  int inum_full = list->listfull->inum;
  int ainum, ainum_full;
  int ginum, ginum_full;
  int ninum, ninum_full;
  int vinum, vinum_full;

  ainum = ainum_full = list->listfull->ainum;
  ginum = ginum_full = list->listfull->ginum;
  ninum = ninum_full = list->listfull->ninum;
  vinum = vinum_full = list->listfull->vinum;

  // reset all pages

  ipage->reset();
  indexpage->reset();

  int inum = 0;

  // loop over parent full list

  for (ii = 0; ii < inum_full; ii++) {
    n = 0;
    neighptr = ipage->vget();
    neighindexptr = indexpage->vget();

    i = ilist_full[ii];
    iindex = iindexlist_full[ii];

    inode = igcell = iucell = -1;
    if (ii < ainum_full) {
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
    } else {
      xtmp = ex[i][0];
      ytmp = ex[i][1];
      ztmp = ex[i][2];
      ietype = etype[i];
      inpe = npe[ietype];
      iapc = apc[ietype];
      ibasis = iindex % iapc;
      if (ii < ainum_full + ninum_full) 
        inode = iindex / iapc;
      else if (ii < ainum_full + ninum_full + ginum_full) 
        igcell = iindex / iapc;
      else iucell = iindex / iapc;
    }

    // loop over full neighbor list

    jlist = firstneigh_full[ii];
    jindexlist = firstneighindex_full[ii];
    jnum = numneigh_full[ii];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];

      // j is atom
       
      if (jindex < 0) {

        // skip if i is not atom
        
        if (iindex >= 0) continue;

        // if j is owned atom, only store if j is beyond i (j > i)
        // if j is ghost, only store if j coords are "above and to the right" of i

        if (j < nalocal) {
          if (i > j) continue;
        } else {
          if (ax[j][2] < ztmp) continue;
          if (ax[j][2] == ztmp) {
            if (ax[j][1] < ytmp) continue;
            if (ax[j][1] == ytmp) {
              if (ax[j][0] < xtmp) continue;
              if (ax[j][0] == xtmp && j <= i) continue;
            }
          }
        }
      } 

      // j is virtual atom
      // store all pairs if i is atom
      // otherwise, skip half of the same type

      else { 
        if (iindex >= 0) {
          jetype = etype[j];
          japc = apc[jetype];
          jucell = jindex / japc;
          jbasis = jindex % japc;
          jgcell = u2g[jetype][jucell];

          // skip half if i and j are both of the same type

          if (inode >= 0) {
            if (jgcell < 0) jindex = -1;
            else {
              jnode = g2n[jetype][jgcell];
              if (jnode >= 0) jindex = jnode * japc + jbasis;
            }
          } else if (igcell >= 0) {
            if (jgcell < 0) jindex = -1;
            else jindex = jgcell * japc + jbasis;
          } 

          if (jindex >= 0) {
            if (j < nelocal) {
              if (i > j) continue;
              else if (i == j && iindex >= jindex) continue;
            } else {
              if (ex[j][2] < ztmp) continue;
              if (ex[j][2] == ztmp) {
                if (ex[j][1] < ytmp) continue;
                if (ex[j][1] == ytmp) {
                  if (ex[j][0] < xtmp) continue;
                  if (ex[j][0] == xtmp) {
                    if (i > j) continue;
                    else if (i == j && iindex >= jindex) continue;
                  }
                }
              }
            }
          }
        }
      }

      neighptr[n] = j;
      neighindexptr[n++] = jindexlist[jj];
    }

    firstneigh[inum] = neighptr;
    firstneighindex[inum] = neighindexptr;
    numneigh[inum] = n;
    ipage->vgot(n);
    indexpage->vgot(n);
    if (ipage->status() || indexpage->status())
      error->one(FLERR, "neighbor list overflow, boost neigh_modify one");
    ilist[inum] = i;
    iindexlist[inum++] = iindex;
  }
  list->inum = inum;
  list->ainum = ainum;
  list->ninum = ninum;
  list->ginum = ginum;
  list->vinum = vinum;
}
