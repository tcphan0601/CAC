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

/* ---------------------------------------------------------------------- */

NPairFullBinAtom::NPairFullBinAtom(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors for both atom-atom and atom-intpl pairs
   every neighbor pair appears in list of both atoms i and j
   ------------------------------------------------------------------------- */

void NPairFullBinAtom::build(NeighList *list)
{

  int iintg,i,j,jj,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,isub;
  int *a2aneighptr,*a2ianeighptr,*a2ia_indexneighptr;
  a2aneighptr = a2ianeighptr = a2ia_indexneighptr = NULL;
  int *jlist;
  int na2a,na2ia;
  int jintpl,node,jnpe;
  double delx,dely,delz,rsq;
  double xtmp,ytmp,ztmp,cutoff,cutoffsq;
  double *jboxlo,*jboxhi;

  double **ax = atom->x;
  double ***nodex = element->nodex;
  double **element_bound_box = element->element_bound_box;
  double *subelem_size = element->subelem_size;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;

  int neall = element->nlocal + element->nghost;
  int *npe = element->npe;
  int *nsubelem = element->nsubelem;
  int **natom_subelem = element->natom_subelem;
  int ***ias2ia = element->ias2ia;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;

  // define pointers from neighbor list

  int *ailist = list->ailist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;

  MyPage<int> *a2apage = list->a2apage;
  MyPage<int> *a2iapage = list->a2iapage;
  MyPage<int> *a2ia_indexpage = list->a2ia_indexpage;


  int ainum = 0;

  // reset all pages

  a2apage->reset();
  a2iapage->reset();
  a2ia_indexpage->reset();

  /*------------------------------------------------------------------- 
   * ATOM NEIGHBORS
   * loop over each atom, storing neighbors for atoms
   *-----------------------------------------------------------------*/ 

  for (i = 0; i < nalocal; i++) {
    ictype = atype[i];
    xtmp = ax[i][0];
    ytmp = ax[i][1];
    ztmp = ax[i][2];
    ibin = atom2bin[i];

    // loop over all atoms in other atom bins in stencil including self
    // only store pair if i != j

    na2a = 0;
    a2aneighptr = a2apage->vget();
    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (i == j) continue;
        jctype = atype[j];
        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutneighsq[ictype][jctype])
          a2aneighptr[na2a++] = j;
      }
    }
    firstneigha2a[ainum] = a2aneighptr;
    numneigha2a[ainum] = na2a;
    a2apage->vgot(na2a);
    if (a2apage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    // loop over all elements in other element bins in stencil including self if there are elements

    if (neall) {
      na2ia = 0;
      a2ianeighptr = a2iapage->vget();
      a2ia_indexneighptr = a2ia_indexpage->vget();
      ibin = nb->coord2elembin(xtmp,ytmp,ztmp);
      for (k = 0; k < nelemstencil; k++) {
        for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

          jetype = etype[j];
          jctype = ctype[j];
          jnpe = npe[jetype];
          cutoff = cutneigh[ictype][jctype];
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
          cutoffsq = cutoff*cutoff;
          for (isub = 0; isub < nsubelem[jetype]; isub++) {
            delx = xtmp; dely = ytmp; delz = ztmp;
            for (node = 0; node < jnpe; node++) {
              delx -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][0];
              dely -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][1];
              delz -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutoffsq) {

              // loop over interpolate the atoms inside the subelement isub of element j that is neighbor to atom i

              for (l = 0; l < natom_subelem[jetype][isub]; l++) {
                delx = xtmp; dely = ytmp; delz = ztmp;
                jintpl = ias2ia[jetype][isub][l];
                for (node = 0; node < jnpe; node++) {
                  delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
                  dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
                  delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
                }
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutneighsq[ictype][jctype]) {
                  a2ianeighptr[na2ia] = j;
                  a2ia_indexneighptr[na2ia++] = jintpl;
                }
              }
            }
          }
        }
      }
      firstneigha2ia[ainum] = a2ianeighptr;
      firstneigha2ia_index[ainum] = a2ia_indexneighptr;
      numneigha2ia[ainum] = na2ia;
      a2iapage->vgot(na2ia);
      a2ia_indexpage->vgot(na2ia);
      if (a2iapage->status() ||
          a2ia_indexpage->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    } else {

      // no elements (local + ghost)

      numneigha2ia[ainum] = 0;
      firstneigha2ia[ainum] = NULL;
      firstneigha2ia_index[ainum] = NULL;
    }
    ailist[ainum++] = i;
  }
  list->ainum = ainum;

}
