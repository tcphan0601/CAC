#include "npair_full_bin_sort_intpl.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "nbin.h"
#include "memory.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinSortIntpl::NPairFullBinSortIntpl(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   build list for atoms, elements, and interpolated atoms
   neighbor are sorted in nearest order
------------------------------------------------------------------------- */

void NPairFullBinSortIntpl::build(NeighList *list)
{
  int ii,i,j,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,jsub;
  int *aneighptr,*a_indexneighptr;
  aneighptr = a_indexneighptr = NULL;
  int *e2eneighptr,*ianeighptr,*ia_indexneighptr;
  e2eneighptr = ianeighptr = ia_indexneighptr = NULL;
  int *jlist;
  int na,nia,ne2e;
  int iintpl,jintpl,inode,jnode;
  double delx,dely,delz,rsq;
  double ix,iy,iz;
  double jx,jy,jz;

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
  int *nintpl = element->nintpl;
  int **natom_subelem = element->natom_subelem;
  int ***ias2ia = element->ias2ia;
  double ***shape_array = element->shape_array;
  double **shape_array_center_subelem = element->shape_array_center_subelem;
  
  int sortnum = list->sortnum;
  int *ibucket;
  int *indexbucket;
  double *rbucket;
  int maxbucket = 50;
  if (sortnum <= 0) error->all(FLERR,"Sorting all neighbor is not supported yet");

  memory->create(ibucket,maxbucket,"npair:ibucket");
  memory->create(indexbucket,maxbucket,"npair:indexbucket");
  memory->create(rbucket,maxbucket,"npair:rbucket");

  // define pointers from neighbor list

  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int *e2ilist = list->e2ilist;
  bigint *e2ialist = list->e2ialist;
  int *ia2elist = list->ia2elist;
  int *numneigha = list->numneigha;
  int *numneighia = list->numneighia;
  int *numneighe2e = list->numneighe2e;
  int **firstneigha = list->firstneigha;
  int **firstneigha_index = list->firstneigha_index;
  int **firstneighe2e = list->firstneighe2e;
  int **firstneighia = list->firstneighia;
  int **firstneighia_index = list->firstneighia_index;

  MyPage<int> *apage = list->apage;
  MyPage<int> *iapage = list->iapage;
  MyPage<int> *e2epage = list->e2epage;
  MyPage<int> *a_indexpage = list->a_indexpage;
  MyPage<int> *ia_indexpage = list->ia_indexpage;

  int ainum = 0;
  int einum = 0;
  bigint iainum = 0;

  // reset all pages

  apage->reset();
  iapage->reset();
  a_indexpage->reset();
  e2epage->reset();
  ia_indexpage->reset();


  /*------------------------------------------------------------------- 
   * ATOM NEIGHBORS
   * loop over each atom, storing and sorting neighbors for atoms
   *-----------------------------------------------------------------*/ 

  for (i = 0; i < nalocal; i++) {
    ictype = atype[i];
    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];
    ibin = atom2bin[i];
    
    na = 0;
    aneighptr = apage->vget();
    a_indexneighptr = a_indexpage->vget();

    // loop over all atoms in other bins in stencil including self
    // only store pair if i != j
    
    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (i == j) continue;
        jctype = atype[j];
        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutneighsq[ictype][jctype]) {
          if (na == maxbucket) {
            maxbucket += maxbucket/2;
            memory->grow(ibucket,maxbucket,"npair:ibucket");
            memory->grow(indexbucket,maxbucket,"npair:indexbucket");
            memory->grow(rbucket,maxbucket,"npair:rbucket");
          }
          ibucket[na] = j;
          indexbucket[na] = -1;
          rbucket[na] = rsq;
          na++;
        }
      }
    }

    // loop over all elements in other element bins in stencil including self if there are elements

    if (neall) {
      ibin = nb->coord2elembin(ix,iy,iz);
      for (k = 0; k < nelemstencil; k++) {
        for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {
          jetype = etype[j];
          jctype = ctype[j];
          delx = ix - ex[j][0];
          dely = iy - ex[j][1];
          delz = iz - ex[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cut_allsq) {

            // interpolate the centers of the subelements of element j and loop over sub elements of element j

            for (jsub = 0; jsub < nsubelem; jsub++) {
              jx = jy = jz = 0.0;
              for (jnode = 0; jnode < npe; jnode++) {
                jx += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][0];
                jy += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][1];
                jz += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][2];
              }
              delx = ix - jx; 
              dely = iy - jy; 
              delz = iz - jz; 
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cut_subelemsq)

                // loop over interpolate the atoms inside the subelement jsub of element j that is neighbor to atom i

                for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
                  jx = jy = jz = 0.0;
                  jintpl = ias2ia[jetype][jsub][l];
                  for (jnode = 0; jnode < npe; jnode++) {
                    jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
                    jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
                    jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
                  }
                  delx = ix - jx; 
                  dely = iy - jy; 
                  delz = iz - jz; 
                  rsq = delx*delx + dely*dely + delz*delz;
                  if (rsq < cutneighsq[ictype][jctype]) {
                    if (na == maxbucket) {
                      maxbucket += maxbucket/2;
                      memory->grow(ibucket,maxbucket,"npair:ibucket");
                      memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                      memory->grow(rbucket,maxbucket,"npair:rbucket");
                    }
                    ibucket[na] = j;
                    indexbucket[na] = jintpl;
                    rbucket[na] = rsq;
                    na++;
                  }
                }
            }
          }
        }
      }
    }

    int numneighsort;

    // if not enough requested neighbors sort everything

    numneighsort = MIN(na,sortnum);

    // perform bucket sort to find nearest neighbors

    int permute[na];
    for (j = 0; j < na; j++)
      permute[j] = j;
    for (j = 0; j < numneighsort; j++) {

      // after this loop, permute[j] (index in bucket) will be the next nearest neighbor
      // sort from bottom up

      for (k = na-1; k > j; k--) {
        int current = permute[k];
        int next = permute[k-1];

        if (rbucket[current] < rbucket[next]) {
          permute[k] = next;
          permute[k-1] = current;
        }
      }
      aneighptr[j] = ibucket[permute[j]];
      a_indexneighptr[j] = indexbucket[permute[j]];
    }
  
    // fill in the rest of the neighbors without sorting

    for (j = numneighsort; j < na; j++) {
      aneighptr[j] = ibucket[permute[j]];
      a_indexneighptr[j] = indexbucket[permute[j]];
    }

    numneigha[ainum] = na;
    firstneigha[ainum] = aneighptr;
    firstneigha_index[ainum] = a_indexneighptr;
    apage->vgot(na);
    a_indexpage->vgot(na);
    if (apage->status() || a_indexpage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
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
    ix = ex[i][0];
    iy = ex[i][1];
    iz = ex[i][2];
    ibin = elem2bin[i];
    e2eneighptr = e2epage->vget();
    ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++)
      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {
        if (i == j) continue;
        delx = ix - ex[i][0];
        dely = iy - ex[i][1];
        delz = iz - ex[i][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_allsq) e2eneighptr[ne2e++] = j;
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
   * loop over each interpolated atom, 
   * store and sort neighbors for interpolated atoms
   *-----------------------------------------------------------------*/ 

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ibin = elem2bin[i];
    ietype = etype[i];
    ictype = ctype[i];

    // store list of mapping between interpolated atoms and element index

    e2ialist[i] = iainum; 

    // loop over integration points in element i

    for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {

      // store list of mapping between and element index and its first interpolated atom index

      ia2elist[iainum] = i;

      // interpolating positions for interpolated atom I

      ix = iy = iz = 0.0;
      for (inode = 0; inode < npe; inode++) {
        ix += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        iy += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        iz += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
      }

      nia = 0;
      ianeighptr = iapage->vget();
      ia_indexneighptr = ia_indexpage->vget();

      // loop over all atoms in all atom bins in stencil including self

      if (naall) {

        ibin = nb->coord2atombin(ix,iy,iz);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            delx = ix - ax[j][0];
            dely = iy - ax[j][1];
            delz = iz - ax[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutneighsq[ictype][jctype]) {
              if (nia == maxbucket) {
                maxbucket += maxbucket/2;
                memory->grow(ibucket,maxbucket,"npair:ibucket");
                memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                memory->grow(rbucket,maxbucket,"npair:rbucket");
              }
              ibucket[nia] = j;
              indexbucket[nia] = -1;
              rbucket[nia] = rsq;
              nia++;
            }
          }
      }

      // loop over neighbor elements of element i and search through 
      // its sub-elements then search through interpolated atoms 
      // within sub-element

      jlist = firstneighe2e[i];
      jnum = numneighe2e[i];
      for (k = 0; k < jnum; k++) {
        j = jlist[k];
        jctype = ctype[j];
        jetype = etype[j];
        delx = ix - ex[j][0];
        dely = iy - ex[j][1];
        delz = iz - ex[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_elemsq)
          for (jsub = 0; jsub < nsubelem; jsub++) {
            jx = jy = jz = 0.0;
            for (jnode = 0; jnode < npe; jnode++) {
              jx += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][0];
              jy += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][1];
              jz += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][2];
            }
            delx = ix - jx; 
            dely = iy - jy; 
            delz = iz - jz; 
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cut_subelemsq)
              for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
                jx = jy = jz = 0.0;
                jintpl = ias2ia[jetype][jsub][l];
                for (jnode = 0; jnode < npe; jnode++) {
                  jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
                  jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
                  jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
                }
                delx = ix - jx; 
                dely = iy - jy; 
                delz = iz - jz; 
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutneighsq[ictype][jctype]) {
                  if (nia == maxbucket) {
                    maxbucket += maxbucket/2;
                    memory->grow(ibucket,maxbucket,"npair:ibucket");
                    memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                    memory->grow(rbucket,maxbucket,"npair:rbucket");
                  }
                  ibucket[nia] = j;
                  indexbucket[nia] = jintpl;
                  rbucket[nia] = rsq;
                  nia++;
                }
              }
          }
      }

      // loop over sub-element in element i and search through 
      // interpolated atoms within sub-element

      for (jsub = 0; jsub < nsubelem; jsub++) {
        jctype = ictype;
        jetype = ietype;
        j = i;
        jx = jy = jz = 0.0;
        for (jnode = 0; jnode < npe; jnode++) {
          jx += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][0];
          jy += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][1];
          jz += shape_array_center_subelem[jsub][jnode]*nodex[j][jnode][2];
        }
        delx = ix - jx; 
        dely = iy - jy; 
        delz = iz - jz; 
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_subelemsq)
          for (l = 0; l < natom_subelem[jetype][jsub]; l++) {
            jintpl = ias2ia[jetype][jsub][l];

            // skip if it's the same atom, i.e. iintpl = jintpl 

            if (iintpl == jintpl) continue;

            jx = jy = jz = 0.0;
            for (jnode = 0; jnode < npe; jnode++) {
              jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
              jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
              jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
            }
            delx = ix - jx; 
            dely = iy - jy; 
            delz = iz - jz; 
            rsq = delx*delx + dely*dely + delz*delz;

            if (rsq < cutneighsq[ictype][ictype]) {
              if (nia == maxbucket) {
                maxbucket += maxbucket/2;
                memory->grow(ibucket,maxbucket,"npair:ibucket");
                memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                memory->grow(rbucket,maxbucket,"npair:rbucket");
              }
              ibucket[nia] = j;
              indexbucket[nia] = jintpl;
              rbucket[nia] = rsq;
              nia++;
            }
          }
      }

      int numneighsort;

      // if not enough requested neighbors sort all

      numneighsort = MIN(nia,sortnum);

      // perform bucket sort to find nearest neighbors

      int permute[nia];
      for (j = 0; j < nia; j++)
        permute[j] = j;
      for (j = 0; j < numneighsort; j++) {

        // after this loop, permute[j] (index in bucket) will be the next nearest neighbor
        // sort from bottom up

        for (k = nia-1; k > j; k--) {
          int current = permute[k];
          int next = permute[k-1];
          if (rbucket[current] < rbucket[next]) {
            permute[k] = next;
            permute[k-1] = current;
          }
        }
        ianeighptr[j] = ibucket[permute[j]];
        ia_indexneighptr[j] = indexbucket[permute[j]];
      }

      // fill in the rest neighbors from buclet
     
      for (j = numneighsort; j < nia; j++) {
        ianeighptr[j] = ibucket[permute[j]];
        ia_indexneighptr[j] = indexbucket[permute[j]];
      }

      numneighia[iainum] = nia;
      firstneighia[iainum] = ianeighptr;
      firstneighia_index[iainum] = ia_indexneighptr;
      iapage->vgot(nia);
      ia_indexpage->vgot(nia);
      if (iapage->status() || ia_indexpage->status())
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
      for (iintpl = 0; iintpl < nintpl[ietype];iintpl++) {
        ia2elist[iainum] = i;
        iainum++;
      }
    }

  memory->destroy(ibucket);
  memory->destroy(indexbucket);
  memory->destroy(rbucket);
}
