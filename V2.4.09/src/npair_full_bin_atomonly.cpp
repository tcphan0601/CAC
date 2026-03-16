#include "npair_full_bin_atomonly.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinAtomonly::NPairFullBinAtomonly(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors for atom-atom pairs only
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullBinAtomonly::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *ilist = list->ailist;
  int *numneigh = list->numneigha2a;
  int **firstneigh = list->firstneigha2a;
  MyPage<int> *ipage = list->a2apage;

  int inum = 0;
  ipage->reset();

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
        if (exclude && exclusion(itype,jtype,mask[i],mask[j])) continue;
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
      }
    }
    ilist[inum++] = i;
    firstneigh[i] = neighptr;

    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

  }

  list->ainum = inum;
}

