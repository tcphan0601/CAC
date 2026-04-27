#include "npair_full_bin_sort_atomonly.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

NPairFullBinSortAtomonly::NPairFullBinSortAtomonly(CAC *cac) : NPair(cac) 
{
  maxbucket = 50;
  memory->create(ibucket, maxbucket, "npair:ibucket");
  memory->create(rbucket, maxbucket, "npair:rbucket");
}

NPairFullBinSortAtomonly::~NPairFullBinSortAtomonly() 
{
  memory->destroy(ibucket);
  memory->destroy(rbucket);
}

/*  ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors for atom-atom pairs only
   every neighbor pair appears in list of both atoms i and j
   neighbors are sorted based on distance
-------------------------------------------------------------------------  */

void NPairFullBinSortAtomonly::build(NeighList *list)
{
  int i, j, k, n, itype, jtype, ibin;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int xbin, ybin, zbin, xbin2, ybin2, zbin2;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;


  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  int sortnum = list->sortnum;
  if (sortnum <= 0) error->all(FLERR, "Sorting all neighbor is not supported yet");

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = atom2bin[i];

    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (i == j) continue;
        jtype = type[j];
        if (exclude && exclusion(itype, jtype, mask[i], mask[j])) continue;
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq <= cutneighsq[itype][jtype]) {
          if (n == maxbucket) {
            maxbucket += maxbucket/2;
            memory->grow(ibucket, maxbucket, "npair:ibucket");
            memory->grow(rbucket, maxbucket, "npair:rbucket");
          }
          ibucket[n] = j;
          rbucket[n] = rsq;
          n++;
        }
      }
    }

    sort_neigh(n, sortnum, neighptr, nullptr, rbucket, ibucket, nullptr);

    firstneigh[inum] = neighptr;
    firstneighindex[inum] = nullptr;
    numneigh[inum] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    ilist[inum] = i;
    iindexlist[inum++] = -1;
  }
  list->inum = inum;
}

