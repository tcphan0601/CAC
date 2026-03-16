#include "npair_full_bin_elem.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "compute.h"
#include "update.h"
#include "my_page.h"
#include "error.h"
#include "comm.h"
#include "nbin.h"
#include "nstencil.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

NPairFullBinElem::NPairFullBinElem(CAC *cac) : NPair(cac) {}

/*  ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   -------------------------------------------------------------------------  */

void NPairFullBinElem::build(NeighList *list)
{

  int i, j, k, ibin;
  int *neighptr, *neighindexptr;
  neighptr = neighindexptr = NULL;
  int *jlist;
  int n;
  double *iboxlo, *jboxlo, *iboxhi, *jboxhi;
  double cutoff;

  int naall = atom->nlocal + atom->nghost;
  int nelocal = element->nlocal;
  int *atype = atom->type;

  double **ax = atom->x;
  double **element_bound_box = element->element_bound_box;

  // define pointers from neighbor list

  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *emask = element->mask;
  int *amask = atom->mask;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  MyPage<int> *ipage = list->ipage;
  MyPage<int> *indexpage = list->indexpage;


  int inum = 0;

  // reset all pages

  ipage->reset();
  indexpage->reset();


  // invoke compute peratom if filter atom
  // if invoked already, don't need to forward comm since it will be done from the class that called this npair

  double *filter_peratom;
  if (atom_filter) {
    if (atom_filter_compute->invoked_peratom != update->ntimestep) {
      atom_filter_compute->compute_peratom();
      atom_filter_compute->preneighflag = atom_filter_preneighflag;

      comm->forward_comm_compute(atom_filter_compute);
    }
    filter_peratom = atom_filter_compute->vector_atom;
  }

  /* ------------------------------------------------------------------- 
  *ELEMENT NEIGHBORS
  *loop over each element and store neighbors for elements
  *use the same cutoff for element and atom neighbors
   *----------------------------------------------------------------- */ 

  // store neighbor if i!=j

  for (i = 0; i < nelocal; i++) {
    if (group && !(emask[i] & groupbit)) continue;

    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i]+3;
    ibin = elem2bin[i];

    neighptr = ipage->vget();
    neighindexptr = indexpage->vget();
    n = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++) {
      if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->melembins) {
        error->one(FLERR, "TEST"); 
      }

      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        if (group && !(emask[j] & groupbit)) continue;
        // check if bouning box overlap

        if (i == j) continue;


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

        neighptr[n] = j;
        neighindexptr[n++] = 0;
      }
    }
    // loop over all atoms in all atom bins in stencil

    if (naall) {

      // atom stencil is identical to element stencil so use element stencil here
      // filter_atom if needed

      for (k = 0; k < nelemstencil; k++) {
        if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->matombins) {
          error->one(FLERR, "TEST"); 
        }

        for (j = atombinhead[ibin+elemstencil[k]]; j >= 0; j = atombins[j]) {
          cutoff = maxcutneigh;

          if (group && !(amask[j] & groupbit)) continue;

          if (atom_filter && atom_filter_value == filter_peratom[j])
            continue;

          if (iboxlo[0] > ax[j][0] + cutoff ||
              iboxhi[0] < ax[j][0] - cutoff)
            continue;
          if (iboxlo[1] > ax[j][1] + cutoff ||
              iboxhi[1] < ax[j][1] - cutoff)
            continue;
          if (iboxlo[2] > ax[j][2] + cutoff ||
              iboxhi[2] < ax[j][2] - cutoff)
            continue;

          neighptr[n] = j;
          neighindexptr[n++] = -1;
        }
      }
    }
    firstneigh[inum] = neighptr;
    firstneighindex[inum] = neighindexptr;
    numneigh[inum] = n;
    ipage->vgot(n);
    indexpage->vgot(n);
    if (ipage->status() || indexpage->status()) {
      char *errms = new char[100];
      sprintf(errms, "Neighbor list overflow for element ID %d with # = %d, boost neigh_modify one", element->tag[i], n);
      error->one(FLERR, errms);
    }

    ilist[inum++] = i;
    iindexlist[inum++] = 0;
  }
  list->inum = list->einum = inum;
}
