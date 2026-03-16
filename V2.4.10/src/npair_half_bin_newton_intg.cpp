#include "npair_half_bin_newton_intg.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonIntg::NPairHalfBinNewtonIntg(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   ------------------------------------------------------------------------- */

void NPairHalfBinNewtonIntg::build(NeighList *list)
{
  int i,j,jj,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,isub,nintg;
  int *a2aneighptr,*a2ianeighptr,*a2ia_indexneighptr;
  a2aneighptr = a2ianeighptr = a2ia_indexneighptr = NULL;
  int *i2aneighptr,*e2eneighptr,*i2ianeighptr,*i2ia_indexneighptr;
  i2aneighptr = e2eneighptr = i2ianeighptr = i2ia_indexneighptr = NULL;
  int *jlist;
  int na2a,na2ia,ni2a,ne2e,ni2ia;
  int iintg,iintpl,jintpl,node,jnpe;
  double delx,dely,delz,rsq;
  double xtmp,ytmp,ztmp,cutoff,cutoffsq;
  double *iboxlo,*jboxlo,*iboxhi,*jboxhi;

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  double **element_bound_box = element->element_bound_box;
  double *subelem_size = element->subelem_size;
  int *amask = atom->mask;
  int *atype = atom->type;
  int *emask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int neall = nelocal + neghost;
  int *npe = element->npe;
  int *nsubelem = element->nsubelem;
  int **natom_subelem = element->natom_subelem;
  int ***ias2ia = element->ias2ia;
  int **i2ia = element->i2ia;
  int **ia2i = element->ia2i;
  int **n2i = element->n2i;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;

  // define pointers from neighbor list

  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int **n2ilist = list->n2ilist;
  int *e2ilist = list->e2ilist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighi2a = list->numneighi2a;
  int *numneighe2e = list->numneighe2e;
  int *numneighi2ia = list->numneighi2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighi2a = list->firstneighi2a;
  int **firstneighe2e = list->firstneighe2e;
  int **firstneighi2ia = list->firstneighi2ia;
  int **firstneighi2ia_index = list->firstneighi2ia_index;

  MyPage<int> *a2apage = list->a2apage;
  MyPage<int> *a2iapage = list->a2iapage;
  MyPage<int> *a2ia_indexpage = list->a2ia_indexpage;

  MyPage<int> *i2apage = list->i2apage;
  MyPage<int> *e2epage = list->e2epage;
  MyPage<int> *i2iapage = list->i2iapage;
  MyPage<int> *i2ia_indexpage = list->i2ia_indexpage;


  int ainum = 0;
  int einum = 0;
  int iinum = 0;

  // reset all pages

  a2apage->reset();
  a2iapage->reset();
  a2ia_indexpage->reset();
  i2apage->reset();
  e2epage->reset();
  i2iapage->reset();
  i2ia_indexpage->reset();

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
    // only store pair if j coords are "above and to the right" of i

    na2a = 0;
    a2aneighptr = a2apage->vget();
    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (ax[j][2] < ztmp) continue;
        if (ax[j][2] == ztmp) {
          if (ax[j][1] < ytmp) continue;
          if (ax[j][1] == ytmp && ax[j][0] <= xtmp) continue;
        }
        jctype = atype[j];
        if (exclude && exclusion(ictype,jctype,amask[i],amask[j])) continue;
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
          if (exclude && exclusion(ictype,jctype,amask[i],emask[j])) continue;
          jnpe = npe[jetype];
          cutoff = cutneigh[ictype][jctype];
          jboxlo = element_bound_box[j];
          jboxhi = element_bound_box[j]+3;

          // check if atom I is within element J bounding box extended by cutoff

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

    for (k = 0; k < nelemstencil; k++) {
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        // check if bouning box overlap

        if (i == j) continue;

        jctype = ctype[j];
        if (exclude && exclusion(ictype,jctype,emask[i],emask[j])) continue;
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
   * INTEGRATION POINT NEIGHBORS
   * loop over each integration point and store neighbors 
   * for integration points
   *-----------------------------------------------------------------*/ 

  for (i = 0; i < einum; i++) {
    ietype = etype[i];
    ictype = ctype[i];
    nintg = element->nintg[ietype];
    jlist = firstneighe2e[i];
    jnum = numneighe2e[i];

    // store list of mapping between node index and integration point index

    for (node = 0; node < npe[ietype]; node++)
      n2ilist[i][node] = iinum + n2i[ietype][node];
    e2ilist[i] = iinum; 

    // loop over integration points in element i

    for (iintg = 0; iintg < nintg; iintg++) {

      xtmp = ytmp = ztmp = 0.0;
      iintpl = i2ia[ietype][iintg];
      for (node = 0; node < npe[ietype]; node++) {
        xtmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
        ytmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
        ztmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
      }

      // loop over all atoms in all atom bins in stencil including self

      if (naall) {
        ni2a = 0;
        i2aneighptr = i2apage->vget();
        ibin = nb->coord2atombin(xtmp,ytmp,ztmp);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            if (exclude && exclusion(ictype,jctype,emask[i],amask[j])) continue;
            delx = xtmp - ax[j][0];
            dely = ytmp - ax[j][1];
            delz = ztmp - ax[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutneighsq[ictype][jctype])
              i2aneighptr[ni2a++] = j;
          }

        firstneighi2a[iinum] = i2aneighptr;
        numneighi2a[iinum] = ni2a;
        i2apage->vgot(ni2a);
        if (i2apage->status())
          error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
      } else { 

        // no atoms (local + ghost)

        numneighi2a[iinum] = 0;
        firstneighi2a[iinum] = NULL;
      }

      // loop over neighbor elements of element i and search through 
      // its sub-elements then search through interpolated atoms 
      // within sub-elements

      ni2ia = 0;
      i2ianeighptr = i2iapage->vget();
      i2ia_indexneighptr = i2ia_indexpage->vget();
      for (jj = 0; jj < jnum; jj++) {

        j = jlist[jj];
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

        for (isub = 0; isub < nsubelem[jetype]; isub++) {

          delx = xtmp; dely = ytmp; delz = ztmp;
          for (node = 0; node < jnpe; node++) {
            delx -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][0];
            dely -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][1];
            delz -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][2];
          }
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < cutoffsq) {
            for (l = 0; l < natom_subelem[jetype][isub]; l++) {
              jintpl = ias2ia[jetype][isub][l];
              delx = xtmp; dely = ytmp; delz = ztmp;
              for (node = 0; node < jnpe; node++) {
                delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
                dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
                delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutneighsq[ictype][jctype]) {

                i2ianeighptr[ni2ia] = j;
                i2ia_indexneighptr[ni2ia++] = jintpl;
              }
            }
          }
        }
      }

      // loop over sub-element in element i and search through interpolated atoms within sub-element

      jctype = ictype;
      jetype = ietype;
      j = i;
      jnpe = npe[jetype];
      cutoff = cutneigh[ictype][jctype] + subelem_size[j];
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
          for (l = 0; l < natom_subelem[jetype][isub]; l++) {
            jintpl = ias2ia[jetype][isub][l];

            // skip if it's the same atom, i.e. iintpl = jintpl 

            if (iintpl == jintpl) continue;

            delx = xtmp; dely = ytmp; delz = ztmp;
            for (node = 0; node < jnpe; node++) {
              delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
              dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
              delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutneighsq[ictype][jctype]) {
              i2ianeighptr[ni2ia] = j;
              i2ia_indexneighptr[ni2ia++] = jintpl;
            }
          }
        }
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
  }
  list->iinum = iinum;

  // finish n2ilist and e2ilist for ghost elements if needed

  if (neghost)
    for (i = nelocal; i < neall; i++) {
      ietype = etype[i];
      nintg = element->nintg[ietype];
      for (k = 0; k < npe[ietype]; k++)
        n2ilist[i][k] = iinum + n2i[ietype][k];
      e2ilist[i] = iinum; 
      iinum += nintg;
    }
}
