#include "npair_full_bin_elem.h"
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

NPairFullBinElem::NPairFullBinElem(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   ------------------------------------------------------------------------- */

void NPairFullBinElem::build(NeighList *list)
{

  int i,j,k,ibin,ictype,jctype;
  int *e2aneighptr,*e2eneighptr;
  e2aneighptr = e2eneighptr = NULL;
  int *jlist;
  int ne2a,ne2e;
  double *iboxlo,*jboxlo,*iboxhi,*jboxhi;
  double cutoff;

  int naall = atom->nlocal + atom->nghost;
  int nelocal = element->nlocal;
  int *ctype = element->ctype;

  double **element_bound_box = element->element_bound_box;

  // define pointers from neighbor list

  int *eilist = list->eilist;
  int *numneighe2a = list->numneighe2a;
  int *numneighe2e = list->numneighe2e;
  int **firstneighe2a = list->firstneighi2a;
  int **firstneighe2e = list->firstneighe2e;

  MyPage<int> *e2apage = list->e2apage;
  MyPage<int> *e2epage = list->e2epage;


  int einum = 0;

  // reset all pages

  e2apage->reset();
  e2epage->reset();


  /*------------------------------------------------------------------- 
   * ELEMENT NEIGHBORS
   * loop over each element and store neighbors for elements
   * use the same cutoff for element and atom neighbors
   *-----------------------------------------------------------------*/ 

  // store neighbor if i!=j

  for (i = 0; i < nelocal; i++) {
    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i]+3;
    ibin = elem2bin[i];
    ictype = ctype[i];
    e2eneighptr = e2epage->vget();
    ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++)
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        // check if bouning box overlap

        if (i == j) continue;

        jctype = ctype[j];
        cutoff = cutneigh[ictype][jctype];
        jboxlo = element_bound_box[j];
        jboxhi = element_bound_box[j]+3;

        for (int idim = 0; idim < domain->dimension; idim++)
          if (iboxlo[idim] > jboxhi[idim] + cutoff ||
              iboxhi[idim] < jboxlo[idim] - cutoff)
            continue;
        e2eneighptr[ne2e++] = j;
      }
    firstneighe2e[einum] = e2eneighptr;
    numneighe2e[einum] = ne2e;
    e2epage->vgot(ne2e);
    if (e2epage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    eilist[einum] = i;


    // loop over all atoms in all atom bins in stencil
    
    //if (naall) {
    //  ne2a = 0;
    //  e2aneighptr = e2apage->vget();
    //  ibin = nb->coord2atombin(xtmp,ytmp,ztmp);
    //  for (k = 0; k < natomstencil; k++) 
    //    for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
    //      delx = xtmp - ax[j][0];
    //      dely = ytmp - ax[j][1];
    //      delz = ztmp - ax[j][2];
    //      rsq = delx*delx + dely*dely + delz*delz;
    //      if (rsq < cut_allsq) e2aneighptr[ne2a++] = j;
    //    }
    //  firstneighe2a[einum] = e2aneighptr;
    //  numneighe2a[einum] = ne2a;
    //  e2apage->vgot(ne2a);
    //  if (e2apage->status())
    //    error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    //} else {

      // no atoms (local + ghost)

    //  numneighe2a[einum] = 0;
    //  firstneighe2a[einum] = NULL;
    //}
//error->all(FLERR,"TEST"); 
    einum++;
  }

  list->einum = einum;
}
