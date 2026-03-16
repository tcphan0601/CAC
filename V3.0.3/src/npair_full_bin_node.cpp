#include "npair_full_bin_node.h"
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

NPairFullBinNode::NPairFullBinNode(CAC *cac) : NPair(cac) {}

/*  ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   -------------------------------------------------------------------------  */

void NPairFullBinNode::build(NeighList *list)
{
  int i, j, jj, jnum, k, l, n, ictype, jctype, ietype, jetype, ibin, isub, iapc, inpe, ngcell;
  int *neighptr = NULL;
  int *neighindexptr = NULL;
  int *jlist;
  int inode, ibasis, iucell, jucell, node, jbasis, jnpe, japc;
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
  int **n2u = element->n2u;
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
  list->ainum = inum;

  /* ------------------------------------------------------------------- 
  *NODE NEIGHBORS
  *loop over each node in node cells and store neighbors 
  *for each element, first create a temporary neighboring element list that can be used for all gaussian points
   *----------------------------------------------------------------- */ 

  for (i = 0; i < nelocal; i++) {
    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i]+3;
    ietype = etype[i];
    ibin = elem2bin[i];
    inpe = npe[ietype];
    iapc = apc[ietype];

    int ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++) {
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        // check if bouning box overlap
        // include i itself

        if (i != j) {
          cutoff = maxcutneigh;
          jboxlo = element_bound_box[j];
          jboxhi = element_bound_box[j]+3;

          if (iboxlo[0] > jboxhi[0] + cutoff ||
              iboxhi[0] < jboxlo[0] - cutoff)
            continue;
          if (iboxlo[1] > jboxhi[1] + cutoff ||
              iboxhi[1] < jboxlo[1] - cutoff)
            continue;
          if (iboxlo[2] > jboxhi[2] + cutoff ||
              iboxhi[2] < jboxlo[2] - cutoff)
            continue;
        }

        if (ne2e == maxelemneigh) grow_elemneighlist();
        elemneighlist[ne2e++] = j;
      }
    }

    ngcell = element->ngcell[ietype];

    // loop over node cells in element i

    for (inode = 0; inode < inpe; inode++) {
      iucell = n2u[ietype][inode];
      for (ibasis = 0; ibasis < iapc; ibasis++) {
        ictype = ctype[i][ibasis];
        n = 0;
        neighptr = ipage->vget();
        neighindexptr = indexpage->vget();

        xtmp = nodex[i][ibasis][inode][0];
        ytmp = nodex[i][ibasis][inode][1];
        ztmp = nodex[i][ibasis][inode][2];

        // loop over atoms and store all pairs

        ibin = nb->coord2atombin(xtmp, ytmp, ztmp);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            if (exclude && exclusion(ictype, jctype, emask[i], amask[j])) continue;
            delx = xtmp - ax[j][0];
            dely = ytmp - ax[j][1];
            delz = ztmp - ax[j][2];
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq <= cutneighsq[ictype][jctype]) {
              neighptr[n] = j;
              neighindexptr[n++] = -1;
            }
          }

        // loop over neighbor elements of element i and search through 
        // its sub-elements then search through virtual atoms 
        // within sub-elements

        for (jj = 0; jj < ne2e; jj++) {

          j = elemneighlist[jj];
          jetype = etype[j];
          jnpe = npe[jetype];
          japc = apc[jetype];
          jboxlo = element_bound_box[j];
          jboxhi = element_bound_box[j]+3;
          cutoff = maxcutneigh;

          if (xtmp < jboxlo[0] - cutoff ||
              xtmp > jboxhi[0] + cutoff) 
            continue;
          if (ytmp < jboxlo[1] - cutoff ||
              ytmp > jboxhi[1] + cutoff) 
            continue;
          if (ztmp < jboxlo[2] - cutoff ||
              ztmp > jboxhi[2] + cutoff) 
            continue;

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

              // loop over virtual atoms inside the subelement isub of element j that is neighbor to i

              for (l = 0; l < nucell_subelem[jetype][isub]; l++) {
                jucell = us2u[jetype][isub][l];
                for (jbasis = 0; jbasis < japc; jbasis++) {

                  // skip if i and j are same atoms
                  
                  if (i == j && ibasis == jbasis && iucell == jucell) continue;
                  jctype = ctype[j][jbasis];
                  if (exclude && exclusion(ictype, jctype, emask[i], emask[j])) continue;

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
        firstneigh[inum] = neighptr;
        firstneighindex[inum] = neighindexptr;
        numneigh[inum] = n;
        ipage->vgot(n);
        indexpage->vgot(n);
        if (ipage->status() || indexpage->status())
          error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

        ilist[inum] = i;
        iindexlist[inum] = inode * iapc + ibasis;
        inum++;
      }
    }
  }
  
  list->ninum = inum - list->ainum;
  list->inum = inum;
}
