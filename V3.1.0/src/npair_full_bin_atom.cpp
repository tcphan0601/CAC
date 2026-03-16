#include "npair_full_bin_atom.h"
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

/*  ----------------------------------------------------------------------  */

NPairFullBinAtom::NPairFullBinAtom(CAC *cac) : NPair(cac) {}

/*  ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors for both atom-atom and atom-vatom pairs
   every neighbor pair appears in list of both atoms i and j
   -------------------------------------------------------------------------  */

void NPairFullBinAtom::build(NeighList *list)
{
  int i, j, jj, jnum, k, l, n, ictype, jctype, ietype, jetype, ibin, isub, ngcell;
  int *neighptr = nullptr;
  int *neighindexptr = nullptr;
  int *jlist;
  int igcell, ibasis, iucell, jucell, node, jbasis, jnpe, japc;
  double delx, dely, delz, rsq;
  double xtmp, ytmp, ztmp, cutoff, cutoffsq;
  double *iboxlo, *jboxlo, *iboxhi, *jboxhi;

  double **ax = atom->x;
  double **ex = element->x;
  double ****nodex = element->nodex;
  double **element_bound_box = element->element_bound_box;
  double *subelem_size = element->subelem_size;
  int *amask = atom->mask;
  int *atype = atom->type;
  int *emask = element->mask;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int neall = nelocal + neghost;
  int *npe = element->npe;
  int *apc = element->apc;
  int *nsubelem = element->nsubelem;
  int **nucell_subelem = element->nucell_subelem;
  int ***us2u = element->us2u;
  int **g2u = element->g2u;
  int **n2g = element->n2g;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;

  // define pointers from neighbor list

  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  MyPage<int> *ipage = list->ipage;
  MyPage<int> *indexpage = list->indexpage;

  int inum = 0;

  // reset all pages

  ipage->reset();
  indexpage->reset();

  /* ------------------------------------------------------------------- 
  *ATOM NEIGHBORS
  *loop over each atom, storing neighbors for atoms
   *----------------------------------------------------------------- */ 

  for (i = 0; i < nalocal; i++) {
    ictype = atype[i];
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];
    ibin = atom2bin[i];
    n = 0;
    neighptr = ipage->vget();
    neighindexptr = indexpage->vget();

    // loop over all atoms in other atom bins in stencil including self
    // only store pair if i != j

    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (i == j) continue;
        jctype = atype[j];
        if (exclude && exclusion(ictype, jctype, amask[i], amask[j])) continue;
        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq <= cutneighsq[ictype][jctype]) {
          neighptr[n] = j;
          neighindexptr[n++] = -1;
        }
      }
    }

    // loop over all elements in other element bins in stencil including self if there are elements

    if (neall) {
      ibin = nb->coord2elembin(xtmp, ytmp, ztmp);
      for (k = 0; k < nelemstencil; k++) {
        for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

          jetype = etype[j];
          jnpe = npe[jetype];
          japc = apc[jetype];
          cutoff = maxcutneigh;
          jboxlo = element_bound_box[j];
          jboxhi = element_bound_box[j]+3;

          // check if atom I is within element j bounding box extended by cutoff

          if (xtmp < jboxlo[0] - cutoff ||
              xtmp > jboxhi[0] + cutoff) 
            continue;
          if (ytmp < jboxlo[1] - cutoff ||
              ytmp > jboxhi[1] + cutoff) 
            continue;
          if (ztmp < jboxlo[2] - cutoff ||
              ztmp > jboxhi[2] + cutoff) 
            continue;

          // interpolate the centers of the subelements of element j and loop over sub elements of element j

          cutoff += subelem_size[j];
          cutoffsq = cutoff * cutoff;
          for (isub = 0; isub < nsubelem[jetype]; isub++) {
            delx = xtmp; dely = ytmp; delz = ztmp;
            for (node = 0; node < jnpe; node++) 
              for (jbasis = 0; jbasis < japc; jbasis++) {
                delx -= shape_array_center_subelem[jetype][isub][node] * 
                  nodex[j][jbasis][node][0];
                dely -= shape_array_center_subelem[jetype][isub][node] * 
                  nodex[j][jbasis][node][1];
                delz -= shape_array_center_subelem[jetype][isub][node] * 
                  nodex[j][jbasis][node][2];
              }
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutoffsq) {

              // loop over interpolate the atoms inside the subelement isub of element j that is neighbor to atom i

              for (l = 0; l < nucell_subelem[jetype][isub]; l++) {
                jucell = us2u[jetype][isub][l];
                for (jbasis = 0; jbasis < japc; jbasis++) {
                  jctype = ctype[j][jbasis];
                  if (exclude && exclusion(ictype, jctype, amask[i], emask[j])) continue;
                  delx = xtmp; dely = ytmp; delz = ztmp;
                  for (node = 0; node < jnpe; node++) {
                    delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
                    dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
                    delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
                  }
                  rsq = delx * delx + dely * dely + delz * delz;
                  if (rsq <= cutneighsq[ictype][jctype]) {
                    neighptr[n] = j;
                    neighindexptr[n++] = jucell * japc + jbasis;
                  }
                }
              }
            }
          }
        }
      }
    }
    firstneigh[inum] = neighptr;
    firstneighindex[inum] = neighindexptr;
    numneigh[inum] = n;
    ipage->vgot(n);
    indexpage->vgot(n);
    if (ipage->status() ||
        indexpage->status())
      error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    ilist[inum] = i;
    iindexlist[inum] = -1;
    inum++;
  }
  list->inum = inum;
  list->ainum = inum;
}
