#include "npair_half_bin_newtoff_intg.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtoffIntg::NPairHalfBinNewtoffIntg(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   ------------------------------------------------------------------------- */

void NPairHalfBinNewtoffIntg::build(NeighList *list)
{

  int i,j,jj,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,isub,nintg;
  int *a2aneighptr,*a2ianeighptr,*a2ia_indexneighptr;
  a2aneighptr = a2ianeighptr = a2ia_indexneighptr = NULL;
  int *i2aneighptr,*e2eneighptr,*i2ianeighptr,*i2ia_indexneighptr;
  i2aneighptr = e2eneighptr = i2ianeighptr = i2ia_indexneighptr = NULL;
  int *jlist;
  int na2a,na2ia,ni2a,ne2e,ni2ia;
  int iintg,jintg,iintpl,jintpl,inode;
  double delx,dely,delz,rsq;
  double *ix = new double[3];
  double *jx = new double[3];

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int neall = nelocal + neghost;
  int npe = element->npe;
  int nsubelem = element->nsubelem;
  int **natom_subelem = element->natom_subelem;
  int ***ias2ia = element->ias2ia;
  int **i2ia = element->i2ia;
  int **ia2i = element->ia2i;
  int **n2i = element->n2i;
  double ***shape_array = element->shape_array;
  double **shape_array_center_subelem = element->shape_array_center_subelem;

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
    ix[0] = ax[i][0];
    ix[1] = ax[i][1];
    ix[2] = ax[i][2];
    ibin = atom2bin[i];

    // loop over all atoms in other atom bins in stencil including self
    // only store pair if j > i
    // store own/own pairs only once
    // store own/ghost pairs on both procs

    na2a = 0;
    a2aneighptr = a2apage->vget();
    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {

        if (j <= i) continue;

        jctype = atype[j];
        delx = ix[0] - ax[j][0];
        dely = ix[1] - ax[j][1];
        delz = ix[2] - ax[j][2];
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

    if (neall){
      na2ia = 0;
      a2ianeighptr = a2iapage->vget();
      a2ia_indexneighptr = a2ia_indexpage->vget();
      ibin = nb->coord2elembin(ix);
      for (k = 0; k < nelemstencil; k++) {
        for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

          jetype = etype[j];
          jctype = ctype[j];
          delx = ix[0] - ex[j][0];
          dely = ix[1] - ex[j][1];
          delz = ix[2] - ex[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cut_elemsq) {

            // interpolate the centers of the subelements of element j and loop over sub elements of element j

            for (isub = 0; isub < nsubelem; isub++) {
              jx[0] = jx[1] = jx[2] = 0.0;
              for (inode = 0; inode < npe; inode++) {
                jx[0] += shape_array_center_subelem[isub][inode]*nodex[j][inode][0];
                jx[1] += shape_array_center_subelem[isub][inode]*nodex[j][inode][1];
                jx[2] += shape_array_center_subelem[isub][inode]*nodex[j][inode][2];
              }
              delx = ix[0] - jx[0]; 
              dely = ix[1] - jx[1]; 
              delz = ix[2] - jx[2]; 
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cut_subelemsq) {

                // loop over interpolate the atoms inside the subelement isub of element j that is neighbor to atom i

                for (l = 0; l < natom_subelem[jetype][isub]; l++) {
                  jx[0] = jx[1] = jx[2] = 0.0;
                  jintpl = ias2ia[jetype][isub][l];
                  for (inode = 0; inode < npe; inode++) {
                    jx[0] += shape_array[jetype][jintpl][inode]*nodex[j][inode][0];
                    jx[1] += shape_array[jetype][jintpl][inode]*nodex[j][inode][1];
                    jx[2] += shape_array[jetype][jintpl][inode]*nodex[j][inode][2];
                  }
                  delx = ix[0] - jx[0]; 
                  dely = ix[1] - jx[1]; 
                  delz = ix[2] - jx[2]; 
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

  for (i = 0; i < nelocal; i++){
    ix[0] = ex[i][0];
    ix[1] = ex[i][1];
    ix[2] = ex[i][2];
    ibin = elem2bin[i];
    e2eneighptr = e2epage->vget();
    ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++)
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]){
        if (i == j) continue;
        delx = ix[0] - ex[j][0];
        dely = ix[1] - ex[j][1];
        delz = ix[2] - ex[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_allsq) e2eneighptr[ne2e++]=j;
      }
    firstneighe2e[einum] = e2eneighptr;
    numneighe2e[einum] = ne2e;
    //    fprintf(screen,"i = %d numneigh = %d cutneighsq = %g\n",i,ne2e,cutneighsq[1][1]);
    //    error->all(FLERR,"TEST"); 
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

  for (i = 0; i < einum; i++){
    ietype = etype[i];
    ictype = ctype[i];
    nintg = element->nintg[ietype];

    // store list of mapping between node index and integration point index

    for (inode = 0; inode < npe; inode++)
      n2ilist[i][inode] = iinum + n2i[ietype][inode];
    e2ilist[i] = iinum; 

    // loop over integration points in element i

    for (iintg = 0; iintg < nintg; iintg++){

      ix[0] = ix[1] = ix[2] = 0.0;
      iintpl = i2ia[ietype][iintg];
      for (inode = 0; inode < npe; inode++) {
        ix[0] += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        ix[1] += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        ix[2] += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
      }

      // loop over all atoms in all atom bins in stencil including self

      if (naall) {
        ni2a = 0;
        i2aneighptr = i2apage->vget();
        ibin = nb->coord2atombin(ix);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            delx = ix[0] - ax[j][0];
            dely = ix[1] - ax[j][1];
            delz = ix[2] - ax[j][2];
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
      jlist = firstneighe2e[i];
      jnum = numneighe2e[i];
      for (jj = 0; jj < jnum; jj++){
        j = jlist[jj];
        jctype = ctype[j];
        jetype = etype[j];
        delx = ix[0] - ex[j][0];
        dely = ix[1] - ex[j][1];
        delz = ix[2] - ex[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_elemsq) {

          for (isub = 0; isub < nsubelem; isub++){
            jx[0] = jx[1] = jx[2] = 0.0;
            for (inode = 0; inode < npe; inode++) {
              jx[0] += shape_array_center_subelem[isub][inode]*nodex[j][inode][0];
              jx[1] += shape_array_center_subelem[isub][inode]*nodex[j][inode][1];
              jx[2] += shape_array_center_subelem[isub][inode]*nodex[j][inode][2];
            }
            delx = ix[0] - jx[0]; 
            dely = ix[1] - jx[1]; 
            delz = ix[2] - jx[2]; 
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cut_subelemsq)

              for (l = 0; l < natom_subelem[jetype][isub]; l++){
                jx[0] = jx[1] = jx[2] = 0.0;
                jintpl = ias2ia[jetype][isub][l];
                for (inode = 0; inode < npe; inode++) {
                  jx[0] += shape_array[jetype][jintpl][inode]*nodex[j][inode][0];
                  jx[1] += shape_array[jetype][jintpl][inode]*nodex[j][inode][1];
                  jx[2] += shape_array[jetype][jintpl][inode]*nodex[j][inode][2];
                }
                delx = ix[0] - jx[0]; 
                dely = ix[1] - jx[1]; 
                delz = ix[2] - jx[2]; 
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

      for (isub = 0; isub < nsubelem; isub++){
        jctype = ictype;
        jetype = ietype;
        j = i;
        jx[0] = jx[1] = jx[2] = 0.0;
        for (inode = 0; inode < npe; inode++) {
          jx[0] += shape_array_center_subelem[isub][inode]*nodex[j][inode][0];
          jx[1] += shape_array_center_subelem[isub][inode]*nodex[j][inode][1];
          jx[2] += shape_array_center_subelem[isub][inode]*nodex[j][inode][2];
        }
        delx = ix[0] - jx[0]; 
        dely = ix[1] - jx[1]; 
        delz = ix[2] - jx[2]; 
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_subelemsq)

          for (l = 0; l < natom_subelem[jetype][isub]; l++) {
            jintpl = ias2ia[jetype][isub][l];

            // skip if it's the same atom, i.e. iintpl = jintpl 

            if (iintpl == jintpl) continue;

            jx[0] = jx[1] = jx[2] = 0.0;
            for (inode = 0; inode < npe; inode++) {
              jx[0] += shape_array[jetype][jintpl][inode]*nodex[j][inode][0];
              jx[1] += shape_array[jetype][jintpl][inode]*nodex[j][inode][1];
              jx[2] += shape_array[jetype][jintpl][inode]*nodex[j][inode][2];
            }
            delx = ix[0] - jx[0]; 
            dely = ix[1] - jx[1]; 
            delz = ix[2] - jx[2]; 
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutneighsq[ictype][jctype]) {
              i2ianeighptr[ni2ia] = j;
              i2ia_indexneighptr[ni2ia++] = jintpl;
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
    for (i = nelocal; i < neall; i++){
      ietype = etype[i];
      nintg = element->nintg[ietype];
      for (k = 0; k < npe; k++)
        n2ilist[i][k] = iinum + n2i[ietype][k];
      e2ilist[i] = iinum; 
      iinum +=nintg;
    }

  delete [] ix;
  delete [] jx;
}
