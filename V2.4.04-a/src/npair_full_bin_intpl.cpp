#include "npair_full_bin_intpl.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "nbin.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinIntpl::NPairFullBinIntpl(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   build list for atoms, elements, and interpolated atoms
------------------------------------------------------------------------- */

void NPairFullBinIntpl::build(NeighList *list)
{
  int ii,i,j,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,jsub;
  int *a2aneighptr,*a2ianeighptr,*a2ia_indexneighptr;
  a2aneighptr = a2ianeighptr = a2ia_indexneighptr = NULL;
  int *ia2aneighptr,*e2eneighptr,*ia2ianeighptr,*ia2ia_indexneighptr;
  ia2aneighptr = e2eneighptr = ia2ianeighptr = ia2ia_indexneighptr = NULL;
  int *jlist;
  int na2a,na2ia,nia2a,ne2e,nia2ia;
  int iintpl,jintpl,inode,jnode,jnpe;
  double delx,dely,delz,rsq;
  double xtmp,ytmp,ztmp,cutoff,cutoffsq;
  double *iboxlo,*jboxlo,*iboxhi,*jboxhi;

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  double **element_bound_box = element->element_bound_box;
  double *subelem_size = element->subelem_size;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int neall = nelocal + neghost;
  int *npe = element->npe;
  int *nsubelem = element->nsubelem;
  int *nintpl = element->nintpl;
  int **natom_subelem = element->natom_subelem;
  int ***ias2ia = element->ias2ia;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;

  // define pointers from neighbor list

  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int *e2ilist = list->e2ilist;
  bigint *e2ialist = list->e2ialist;
  int *ia2elist = list->ia2elist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighia2a = list->numneighia2a;
  int *numneighe2e = list->numneighe2e;
  int *numneighia2ia = list->numneighia2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighia2a = list->firstneighia2a;
  int **firstneighe2e = list->firstneighe2e;
  int **firstneighia2ia = list->firstneighia2ia;
  int **firstneighia2ia_index = list->firstneighia2ia_index;

  MyPage<int> *a2apage = list->a2apage;
  MyPage<int> *a2iapage = list->a2iapage;
  MyPage<int> *a2ia_indexpage = list->a2ia_indexpage;

  MyPage<int> *ia2apage = list->ia2apage;
  MyPage<int> *e2epage = list->e2epage;
  MyPage<int> *ia2iapage = list->ia2iapage;
  MyPage<int> *ia2ia_indexpage = list->ia2ia_indexpage;


  int ainum = 0;
  int einum = 0;
  bigint iainum = 0;

  // reset all pages

  a2apage->reset();
  a2iapage->reset();
  a2ia_indexpage->reset();
  ia2apage->reset();
  e2epage->reset();
  ia2iapage->reset();
  ia2ia_indexpage->reset();


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

    // loop over all atoms in other bins in stencil including self
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

          for (jsub = 0; jsub < nsubelem[jetype]; jsub++) {
            delx = xtmp; dely = ytmp; delz = ztmp;
            for (jnode = 0; jnode < jnpe; jnode++) {
              delx -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][0];
              dely -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][1];
              delz -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutoffsq)

              // loop over interpolate the atoms inside the subelement jsub of element j that is neighbor to atom i

              for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
                jintpl = ias2ia[jetype][jsub][l];
                delx = xtmp; dely = ytmp; delz = ztmp;
                for (jnode = 0; jnode < jnpe; jnode++) {
                  delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
                  dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
                  delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
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

  /*------------------------------------------------------------------- 
   * ELEMENT NEIGHBORS
   * loop over each element and store neighbors for elements
   * need for finding neighbors of integration points
   *-----------------------------------------------------------------*/ 

  // store neighbor if i!=j

  for (i = 0; i < nelocal; i++) {
    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i]+3;
    ictype = ctype[i];
    ibin = elem2bin[i];
    e2eneighptr = e2epage->vget();
    ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++)
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        // check if bounding box overlap

        if (i == j) continue;

        jctype = ctype[j];
        cutoff = cutneigh[ictype][jctype];
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

        e2eneighptr[ne2e++] = j;
      }
    firstneighe2e[einum] = e2eneighptr;
    numneighe2e[einum] = ne2e;
    e2epage->vgot(ne2e);
    if (e2epage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    eilist[einum++] = i;

  }
  list->einum = einum;

  /*------------------------------------------------------------------- 
   * INTERPOLATED ATOM NEIGHBORS
   * loop over each interpolated atom and store neighbors 
   * for interpolated atoms
   *-----------------------------------------------------------------*/ 

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    ictype = ctype[i];
    jlist = firstneighe2e[i];
    jnum = numneighe2e[i];

    // store list of mapping between element index and its first interpolated atom index

    e2ialist[i] = iainum; 

    // loop over integration points in element i

    for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {

      // store list of mapping between interpolated atom index and element index

      ia2elist[iainum] = i;

      // interpolating positions for interpolated atom I

      xtmp = ytmp = ztmp = 0.0;
      for (inode = 0; inode < npe[ietype]; inode++) {
        xtmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        ytmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        ztmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
      }

      //loop over all atoms in all atom bins in stencil including self

      if (naall) {

        nia2a = 0;
        ia2aneighptr = ia2apage->vget();
        ibin = nb->coord2atombin(xtmp,ytmp,ztmp);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            delx = xtmp - ax[j][0];
            dely = ytmp - ax[j][1];
            delz = ztmp - ax[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutneighsq[ictype][jctype])
              ia2aneighptr[nia2a++] = j;
          }
        firstneighia2a[iainum] = ia2aneighptr;
        numneighia2a[iainum] = nia2a;
        ia2apage->vgot(nia2a);
        if (ia2apage->status())
          error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
      } else { 

        // no atoms (local + ghost)

        numneighia2a[iainum] = 0;
        firstneighia2a[iainum] = NULL;
      }

      // loop over neighbor elements of element i and search through 
      // its sub-elements then search through interpolated atoms 
      // within sub-element

      nia2ia = 0;
      ia2ianeighptr = ia2iapage->vget();
      ia2ia_indexneighptr = ia2ia_indexpage->vget();
      for (k = 0; k < jnum; k++) {
        j = jlist[k];
        jctype = ctype[j];
        jetype = etype[j];
        jnpe = npe[jetype];
        jboxlo = element_bound_box[j];
        jboxhi = element_bound_box[j]+3;
        cutoff = cutneigh[ictype][jctype];
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
        cutoffsq = cutoff*cutoff;
        for (jsub = 0; jsub < nsubelem[jetype]; jsub++) {
          delx = xtmp; dely = ytmp; delz = ztmp;
          for (jnode = 0; jnode < jnpe; jnode++) {
            delx -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][0];
            dely -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][1];
            delz -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][2];
          }
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutoffsq) 
            for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
              jintpl = ias2ia[jetype][jsub][l];
              delx = xtmp; dely = ytmp; delz = ztmp;
              for (jnode = 0; jnode < jnpe; jnode++) {
                delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
                dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
                delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutneighsq[ictype][jctype]) {
                ia2ianeighptr[nia2ia] = j;
                ia2ia_indexneighptr[nia2ia++] = jintpl;
              }
            }
        }
      }

      // loop over sub-element in element i and search through 
      // interpolated atoms within sub-element

      jctype = ictype;
      jetype = ietype;
      j = i;
      jnpe = npe[jetype];
      cutoff = cutneigh[ictype][jctype] + subelem_size[j];
      cutoffsq = cutoff*cutoff;

      for (jsub = 0; jsub < nsubelem[jetype]; jsub++) {
        delx = xtmp; dely = ytmp; delz = ztmp;
        for (jnode = 0; jnode < jnpe; jnode++) {
          delx -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][0];
          dely -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][1];
          delz -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][2];
        }
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutoffsq)
          for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
            jintpl = ias2ia[jetype][jsub][l];

            // skip if it's the same atom, i.e. iintpl = jintpl 

            if (iintpl == jintpl) continue;

            delx = xtmp; dely = ytmp; delz = ztmp;
            for (jnode = 0; jnode < jnpe; jnode++) {
              delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
              dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
              delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;

            if (rsq < cutneighsq[ictype][ictype]) {
              ia2ianeighptr[nia2ia] = i;
              ia2ia_indexneighptr[nia2ia++] = jintpl;
            }
          }
      }
      firstneighia2ia[iainum] = ia2ianeighptr;
      firstneighia2ia_index[iainum] = ia2ia_indexneighptr;
      numneighia2ia[iainum] = nia2ia;
      ia2iapage->vgot(nia2ia);
      ia2ia_indexpage->vgot(nia2ia);
      if (ia2iapage->status() ||
          ia2ia_indexpage->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
      iainum++;
    }
  }
  list->iainum = iainum;

  // finish e2ialist and ia2elist for ghost elements if needed

  if (neghost) 
    for (i = nelocal; i < neall; i++) {
      ietype = etype[i];
      e2ialist[i] = iainum; 
      for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
        ia2elist[iainum] = i;
        iainum++;
      }
    }
}
