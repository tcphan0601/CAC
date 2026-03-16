#include "npair_full_bin_intg_double.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "my_page.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinIntgDouble::NPairFullBinIntgDouble(CAC *cac) : NPair(cac) 
{
  maxbucket = 50;
  memory->create(ibucket,maxbucket,"npair:ibucket");
  memory->create(indexbucket,maxbucket,"npair:indexbucket");
  memory->create(rbucket,maxbucket,"npair:rbucket");
}

NPairFullBinIntgDouble::~NPairFullBinIntgDouble() 
{
  memory->destroy(ibucket);
  memory->destroy(indexbucket);
  memory->destroy(rbucket);
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
   ------------------------------------------------------------------------- */

void NPairFullBinIntgDouble::build(NeighList *list)
{
  int iintg,i,ii,j,jj,jnum,k,l,n,ictype,jctype,ietype,jetype,ibin,isub,nintg;
  int *a2aneighptr,*a2ianeighptr,*a2ia_indexneighptr;
  a2aneighptr = a2ianeighptr = a2ia_indexneighptr = NULL;
  int *i2aneighptr,*e2eneighptr,*i2ianeighptr,*i2ia_indexneighptr;
  i2aneighptr = e2eneighptr = i2ianeighptr = i2ia_indexneighptr = NULL;
  int *a2aouterneighptr,*a2iaouterneighptr,*a2ia_indexouterneighptr;
  a2aouterneighptr = a2iaouterneighptr = a2ia_indexouterneighptr = NULL;
  int *i2aouterneighptr,*i2iaouterneighptr,*i2ia_indexouterneighptr;
  i2aouterneighptr = i2iaouterneighptr = i2ia_indexouterneighptr = NULL;
  int *nia2aneighptr,*nia2ianeighptr,*nia2ia_indexneighptr;
  nia2aneighptr = nia2ianeighptr = nia2ia_indexneighptr = NULL;

  int *jlist;
  int na2a,na2ia,ni2a,ne2e,ni2ia;
  int na2aouter,na2iaouter,ni2aouter,ne2eouter,ni2iaouter;
  int iintpl,jintpl,node,jnpe;
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
  int **n2i = element->n2i;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;

  // nia
  int nnia2a,nnia2ia;
  int **ia2i = element->ia2i;
  int **ia2nialist = list->ia2nialist;
  int *nia2ialist = list->nia2ialist;
  int *nia2ia_indexlist = list->nia2ia_indexlist;
  for (i = 0; i < neall; i++)
    for (j = 0; j < element->max_nintpl; j++)
      ia2nialist[i][j] = -1;
  MyPage<int> *nia2apage = list->nia2apage;
  MyPage<int> *nia2iapage = list->nia2iapage;
  MyPage<int> *nia2ia_indexpage = list->nia2ia_indexpage;
  int **firstneighnia2a = list->firstneighnia2a;
  int **firstneighnia2ia = list->firstneighnia2ia;
  int **firstneighnia2ia_index = list->firstneighnia2ia_index;
  int *numneighnia2a = list->numneighnia2a;
  int *numneighnia2ia = list->numneighnia2ia;
  nia2apage->reset();
  nia2iapage->reset();
  nia2ia_indexpage->reset();


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

  int *numneigha2a_outer = list->numneigha2a_outer;
  int *numneigha2ia_outer = list->numneigha2ia_outer;
  int *numneighi2ia_outer = list->numneighi2ia_outer;
  int *numneighi2a_outer = list->numneighi2a_outer;
  int **firstneigha2a_outer = list->firstneigha2a_outer;
  int **firstneigha2ia_outer = list->firstneigha2ia_outer;
  int **firstneigha2ia_index_outer = list->firstneigha2ia_index_outer;
  int **firstneighi2a_outer = list->firstneighi2a_outer;
  int **firstneighi2ia_outer = list->firstneighi2ia_outer;
  int **firstneighi2ia_index_outer = list->firstneighi2ia_index_outer;

  MyPage<int> *a2apage = list->a2apage;
  MyPage<int> *a2iapage = list->a2iapage;
  MyPage<int> *a2ia_indexpage = list->a2ia_indexpage;

  MyPage<int> *i2apage = list->i2apage;
  MyPage<int> *e2epage = list->e2epage;
  MyPage<int> *i2iapage = list->i2iapage;
  MyPage<int> *i2ia_indexpage = list->i2ia_indexpage;

  MyPage<int> *a2a_outer_page = list->a2a_outer_page;
  MyPage<int> *a2ia_outer_page = list->a2ia_outer_page;
  MyPage<int> *a2ia_index_outer_page = list->a2ia_index_outer_page;

  MyPage<int> *i2a_outer_page = list->i2a_outer_page;
  MyPage<int> *i2ia_outer_page = list->i2ia_outer_page;
  MyPage<int> *i2ia_index_outer_page = list->i2ia_index_outer_page;

  int ainum = 0;
  int einum = 0;
  int iinum = 0;
  int niainum = 0;

  // reset all pages

  a2apage->reset();
  a2iapage->reset();
  a2ia_indexpage->reset();
  i2apage->reset();
  e2epage->reset();
  i2iapage->reset();
  i2ia_indexpage->reset();

  a2a_outer_page->reset();
  a2ia_outer_page->reset();
  a2ia_index_outer_page->reset();
  i2a_outer_page->reset();
  i2ia_outer_page->reset();
  i2ia_index_outer_page->reset();

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
    na2aouter = 0;
    a2aouterneighptr = a2a_outer_page->vget();
    for (k = 0; k < natomstencil; k++) {
      for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
        if (i == j) continue;
        jctype = atype[j];
        if (exclude && exclusion(ictype,jctype,amask[i],amask[j])) continue;
        delx = xtmp - ax[j][0];
        dely = ytmp - ax[j][1];
        delz = ztmp - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < innercutneighsq[ictype][jctype])
          a2aneighptr[na2a++] = j;
        else if (rsq < cutneighsq[ictype][jctype])
          a2aouterneighptr[na2aouter++] = j;
      }
    }
    firstneigha2a[ainum] = a2aneighptr;
    numneigha2a[ainum] = na2a;
    firstneigha2a_outer[ainum] = a2aouterneighptr;
    numneigha2a_outer[ainum] = na2aouter;
    a2apage->vgot(na2a);
    if (a2apage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    a2a_outer_page->vgot(na2aouter);
    if (a2a_outer_page->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    // loop over all elements in other element bins in stencil including self if there are elements

    if (neall) {
      na2ia = 0;
      a2ianeighptr = a2iapage->vget();
      a2ia_indexneighptr = a2ia_indexpage->vget();
      na2iaouter = 0;
      a2iaouterneighptr = a2ia_outer_page->vget();
      a2ia_indexouterneighptr = a2ia_index_outer_page->vget();

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
                if (rsq < innercutneighsq[ictype][jctype]) {
                  a2ianeighptr[na2ia] = j;
                  a2ia_indexneighptr[na2ia++] = jintpl;
                  if (ia2i[jetype][jintpl] < 0 && ia2nialist[j][jintpl] < 0) {
                    ia2nialist[j][jintpl] = niainum;
                    nia2ialist[niainum] = j;
                    nia2ia_indexlist[niainum] = jintpl;
                    niainum++;
                  }

                } else if (rsq < cutneighsq[ictype][jctype]) {
                  a2iaouterneighptr[na2iaouter] = j;
                  a2ia_indexouterneighptr[na2iaouter++] = jintpl;
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

      firstneigha2ia_outer[ainum] = a2iaouterneighptr;
      firstneigha2ia_index_outer[ainum] = a2ia_indexouterneighptr;
      numneigha2ia_outer[ainum] = na2iaouter;
      a2ia_outer_page->vgot(na2iaouter);
      a2ia_index_outer_page->vgot(na2iaouter);
      if (a2ia_outer_page->status() ||
          a2ia_index_outer_page->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    } else {

      // no elements (local + ghost)

      numneigha2ia[ainum] = 0;
      firstneigha2ia[ainum] = NULL;
      firstneigha2ia_index[ainum] = NULL;
      numneigha2ia_outer[ainum] = 0;
      firstneigha2ia_outer[ainum] = NULL;
      firstneigha2ia_index_outer[ainum] = NULL;
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

  for (i = 0; i < neall; i++) {
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
    firstneighe2e[einum] = e2eneighptr;
    numneighe2e[einum] = ne2e;
    e2epage->vgot(ne2e);
    if (e2epage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    eilist[einum] = i;
    einum++;

  }
  list->einum = nelocal;

  /*------------------------------------------------------------------- 
   * INTEGRATION POINT NEIGHBORS
   * loop over each integration point and store neighbors 
   * for integration points
   *-----------------------------------------------------------------*/ 

  for (i = 0; i < nelocal; i++) {
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
        ni2aouter = 0;
        i2aouterneighptr = i2a_outer_page->vget();

        ibin = nb->coord2atombin(xtmp,ytmp,ztmp);
        for (k = 0; k < natomstencil; k++) 
          for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
            jctype = atype[j];
            if (exclude && exclusion(ictype,jctype,emask[i],amask[j])) continue;
            delx = xtmp - ax[j][0];
            dely = ytmp - ax[j][1];
            delz = ztmp - ax[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < innercutneighsq[ictype][jctype])
              i2aneighptr[ni2a++] = j;
            else if (rsq < cutneighsq[ictype][jctype])
              i2aouterneighptr[ni2aouter++] = j;

          }

        firstneighi2a[iinum] = i2aneighptr;
        numneighi2a[iinum] = ni2a;
        i2apage->vgot(ni2a);
        if (i2apage->status())
          error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

        firstneighi2a_outer[iinum] = i2aouterneighptr;
        numneighi2a_outer[iinum] = ni2aouter;
        i2a_outer_page->vgot(ni2aouter);
        if (i2a_outer_page->status())
          error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
      } else { 

        // no atoms (local + ghost)

        numneighi2a[iinum] = 0;
        firstneighi2a[iinum] = NULL;
        numneighi2a_outer[iinum] = 0;
        firstneighi2a_outer[iinum] = NULL;
      }

      // loop over neighbor elements of element i and search through 
      // its sub-elements then search through interpolated atoms 
      // within sub-elements

      ni2ia = 0;
      i2ianeighptr = i2iapage->vget();
      i2ia_indexneighptr = i2ia_indexpage->vget();
      ni2iaouter = 0;
      i2iaouterneighptr = i2ia_outer_page->vget();
      i2ia_indexouterneighptr = i2ia_index_outer_page->vget();

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
          if (rsq < cutoffsq)

            for (l = 0; l < natom_subelem[jetype][isub]; l++) {
              jintpl = ias2ia[jetype][isub][l];
              delx = xtmp; dely = ytmp; delz = ztmp;
              for (node = 0; node < jnpe; node++) {
                delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
                dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
                delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < innercutneighsq[ictype][jctype]) {
                if (ni2ia == maxbucket) {
                  maxbucket += maxbucket/2;
                  memory->grow(ibucket,maxbucket,"npair:ibucket");
                  memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                  memory->grow(rbucket,maxbucket,"npair:rbucket");
                }
                ibucket[ni2ia] = j;
                indexbucket[ni2ia] = jintpl;
                rbucket[ni2ia] = rsq;
                ni2ia++;

                if (ia2i[jetype][jintpl] < 0 && ia2nialist[j][jintpl] < 0) {
                  ia2nialist[j][jintpl] = niainum;
                  nia2ialist[niainum] = j;
                  nia2ia_indexlist[niainum] = jintpl;
                  niainum++;
                }

              } else if (rsq < cutneighsq[ictype][jctype]) {
                i2iaouterneighptr[ni2iaouter] = j;
                i2ia_indexouterneighptr[ni2iaouter++] = jintpl;
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
            if (rsq < innercutneighsq[ictype][jctype]) {
              if (ni2ia == maxbucket) {
                maxbucket += maxbucket/2;
                memory->grow(ibucket,maxbucket,"npair:ibucket");
                memory->grow(indexbucket,maxbucket,"npair:indexbucket");
                memory->grow(rbucket,maxbucket,"npair:rbucket");
              }
              ibucket[ni2ia] = j;
              indexbucket[ni2ia] = jintpl;
              rbucket[ni2ia] = rsq;
              ni2ia++;
              if (ia2i[jetype][jintpl] < 0 && ia2nialist[j][jintpl] < 0) {
                ia2nialist[j][jintpl] = niainum;
                nia2ialist[niainum] = j;
                nia2ia_indexlist[niainum] = jintpl;
                niainum++;
              }
            } else if (rsq < cutneighsq[ictype][jctype]) {
              i2iaouterneighptr[ni2iaouter] = j;
              i2ia_indexouterneighptr[ni2iaouter++] = jintpl;
            }
          }
        }
      }

      // perform bucket sort to find nearest neighbors

      int numneighsort = MIN(ni2ia,0);
      int permute[ni2ia];
      for (j = 0; j < ni2ia; j++)
        permute[j] = j;
      for (j = 0; j < numneighsort; j++) {

        // after this loop, permute[j] (index in bucket) will be the next nearest neighbor
        // sort from bottom up

        for (k = ni2ia-1; k > j; k--) {
          int current = permute[k];
          int next = permute[k-1];
          if (rbucket[current] < rbucket[next]) {
            permute[k] = next;
            permute[k-1] = current;
          }
        }
        i2ianeighptr[j] = ibucket[permute[j]];
        i2ia_indexneighptr[j] = indexbucket[permute[j]];
      }

      // fill in the rest neighbors from buclet

      for (j = numneighsort; j < ni2ia; j++) {
        i2ianeighptr[j] = ibucket[permute[j]];
        i2ia_indexneighptr[j] = indexbucket[permute[j]];
      }

      firstneighi2ia[iinum] = i2ianeighptr;
      firstneighi2ia_index[iinum] = i2ia_indexneighptr;
      numneighi2ia[iinum] = ni2ia;
      i2iapage->vgot(ni2ia);
      i2ia_indexpage->vgot(ni2ia);
      if (i2iapage->status() ||
          i2ia_indexpage->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

      firstneighi2ia_outer[iinum] = i2iaouterneighptr;
      firstneighi2ia_index_outer[iinum] = i2ia_indexouterneighptr;
      numneighi2ia_outer[iinum] = ni2iaouter;
      i2ia_outer_page->vgot(ni2iaouter);
      i2ia_index_outer_page->vgot(ni2iaouter);
      if (i2ia_outer_page->status() ||
          i2ia_index_outer_page->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

      iinum++;
    }

  }
  list->iinum = iinum;

  list->niainum = niainum; 
  for (ii = 0; ii < niainum; ii++) {
    i = nia2ialist[ii];
    iintpl = nia2ia_indexlist[ii];
    ietype = etype[i];
    ictype = ctype[i];

    xtmp = ytmp = ztmp = 0.0;
    for (node = 0; node < npe[ietype]; node++) {
      xtmp += shape_array[ietype][iintpl][node]*nodex[i][node][0];
      ytmp += shape_array[ietype][iintpl][node]*nodex[i][node][1];
      ztmp += shape_array[ietype][iintpl][node]*nodex[i][node][2];
    } 

    // loop over all atoms in all atom bins in stencil including self

    if (naall) {
      nnia2a = 0;
      nia2aneighptr = nia2apage->vget();
      ibin = nb->coord2atombin(xtmp,ytmp,ztmp);
      for (k = 0; k < natomstencil; k++) 
        for (j = atombinhead[ibin+atomstencil[k]]; j >= 0; j = atombins[j]) {
          jctype = atype[j];
          if (exclude && exclusion(ictype,jctype,emask[i],amask[j])) continue;
          delx = xtmp - ax[j][0];
          dely = ytmp - ax[j][1];
          delz = ztmp - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < innercutneighsq[ictype][jctype])
            nia2aneighptr[nnia2a++] = j;
        }
      firstneighnia2a[ii] = nia2aneighptr;
      numneighnia2a[ii] = nnia2a;
      nia2apage->vgot(nnia2a);
      if (nia2apage->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    } else { 

      // no atoms (local + ghost)

      numneighnia2a[ii] = 0;
      firstneighnia2a[ii] = NULL;
    }

    // loop over neighbor elements of element i and search through 
    // its sub-elements then search through interpolated atoms 
    // within sub-elements

    nnia2ia = 0;

    jnum = numneighe2e[i];
    jlist = firstneighe2e[i];
    nia2ia_indexneighptr = nia2ia_indexpage->vget();
    nia2ianeighptr = nia2iapage->vget();

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jctype = ctype[j];
      jetype = etype[j];
      jnpe = npe[jetype];
      jboxlo = element_bound_box[j];
      jboxhi = element_bound_box[j]+3;
      cutoff = sqrt(innercutneighsq[ictype][jctype]);
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
        if (rsq < cutoffsq)

          for (l = 0; l < natom_subelem[jetype][isub]; l++) {
            jintpl = ias2ia[jetype][isub][l];
            delx = xtmp; dely = ytmp; delz = ztmp;
            for (node = 0; node < jnpe; node++) {
              delx -= shape_array[jetype][jintpl][node]*nodex[j][node][0];
              dely -= shape_array[jetype][jintpl][node]*nodex[j][node][1];
              delz -= shape_array[jetype][jintpl][node]*nodex[j][node][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < innercutneighsq[ictype][jctype]) {

              nia2ianeighptr[nnia2ia] = j;
              nia2ia_indexneighptr[nnia2ia++] = jintpl;
            }
          }

      }

    }

    // loop over sub-element in element i and search through interpolated atoms within sub-element

    jctype = ictype;
    jetype = ietype;
    j = i;
    jnpe = npe[jetype];
    cutoff = sqrt(innercutneighsq[ictype][jctype]) + subelem_size[j];
    cutoffsq = cutoff*cutoff;

    for (isub = 0; isub < nsubelem[jetype]; isub++) {
      delx = xtmp; dely = ytmp; delz = ztmp;
      for (node = 0; node < jnpe; node++) {
        delx -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][0];
        dely -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][1];
        delz -= shape_array_center_subelem[jetype][isub][node]*nodex[j][node][2];
      }
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutoffsq)

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
          if (rsq < innercutneighsq[ictype][jctype]) {
            nia2ianeighptr[nnia2ia] = j;
            nia2ia_indexneighptr[nnia2ia++] = jintpl;
          }
        }
    }

    firstneighnia2ia[ii] = nia2ianeighptr;
    firstneighnia2ia_index[ii] = nia2ia_indexneighptr;
    numneighnia2ia[ii] = nnia2ia;
    nia2iapage->vgot(nnia2ia);
    if (nia2iapage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    nia2ia_indexpage->vgot(nnia2ia);
    if (nia2ia_indexpage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

  }

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
