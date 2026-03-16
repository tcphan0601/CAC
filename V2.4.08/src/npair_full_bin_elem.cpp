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
  int *atype = atom->type;

  double **ax = atom->x;
  double **element_bound_box = element->element_bound_box;

  // define pointers from neighbor list

  int *eilist = list->eilist;
  int *emask = element->mask;
  int *amask = atom->mask;
  int *numneighe2a = list->numneighe2a;
  int *numneighe2e = list->numneighe2e;
  int **firstneighe2a = list->firstneighe2a;
  int **firstneighe2e = list->firstneighe2e;

  MyPage<int> *e2apage = list->e2apage;
  MyPage<int> *e2epage = list->e2epage;


  int einum = 0;

  // reset all pages

  e2apage->reset();
  e2epage->reset();


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

  /*------------------------------------------------------------------- 
   * ELEMENT NEIGHBORS
   * loop over each element and store neighbors for elements
   * use the same cutoff for element and atom neighbors
   *-----------------------------------------------------------------*/ 

  // store neighbor if i!=j

  for (i = 0; i < nelocal; i++) {
    if (group && !(emask[i] & groupbit)) continue;

    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i]+3;
    ibin = elem2bin[i];
    ictype = ctype[i];

    e2eneighptr = e2epage->vget();
    ne2e = 0;

    // loop over all elements in other element bins in stencil including self

    for (k = 0; k < nelemstencil; k++) {
      if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->melembins) {
        printf("Step = %d me = %d ibin = %d elemstencil[%d] = %d melembins = %d\n",
            update->ntimestep,comm->me,ibin,k,elemstencil[k],nb->melembins);
        printf("me = %04d i = %08d bound_box = %g/%g %g/%g %g/%g\n     subbox = %g/%g %g/%g %g/%g\n"
            ,comm->me,i,iboxlo[0],iboxhi[0],iboxlo[1],iboxhi[1],iboxlo[2],iboxhi[2],
            domain->sublo[0],domain->subhi[0],
            domain->sublo[1],domain->subhi[1],
            domain->sublo[2],domain->subhi[2]
            );
        printf("x = %g %g %g bbox = %g/%g %g/%g %g/%g\n     elembininv = %g %g %g nelembin = %d %d %d melembin = %d %d %d melembinlo = %d %d %d\n ibinstyle = %d\n  from stencil melembin = %d %d %d ibinstyle = %d nelembin = %d %d %d melembinlo = %d %d %d\n"
            ,element->x[i][0]
            ,element->x[i][1]
            ,element->x[i][2]
            ,neighbor->bboxlo[0],neighbor->bboxhi[0]
            ,neighbor->bboxlo[1],neighbor->bboxhi[1]
            ,neighbor->bboxlo[2],neighbor->bboxhi[2]
            ,elembininvx,elembininvy,elembininvz
            ,nelembinx,nelembiny,nelembinz
            ,melembinx,melembiny,melembinz
            ,melembinxlo,melembinylo,melembinzlo
            ,nb->istyle
            ,ns->nb->melembinx
            ,ns->nb->melembiny
            ,ns->nb->melembinz
            ,ns->nb->istyle
            ,ns->nb->nelembinx
            ,ns->nb->nelembiny
            ,ns->nb->nelembinz
            ,ns->nb->melembinxlo
            ,ns->nb->melembinylo
            ,ns->nb->melembinzlo
            );
        error->one(FLERR,"TEST"); 
      }

      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        if (j > element->nlocal+element->nghost) {
          printf("Step = %d me = %d element j = %d nlocal = %d nghost = %d\n",
              update->ntimestep,comm->me,j,element->nlocal,element->nghost);
          error->one(FLERR,"TEST"); 
        }

        if (group && !(emask[j] & groupbit)) continue;
        // check if bouning box overlap

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
    }
    firstneighe2e[einum] = e2eneighptr;
    numneighe2e[einum] = ne2e;
    e2epage->vgot(ne2e);
    if (e2epage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    // loop over all atoms in all atom bins in stencil

    if (naall) {
      ne2a = 0;
      e2aneighptr = e2apage->vget();

      // atom stencil is identical to element stencil so use element stencil here
      // filter_atom if needed

      //if ((element->tag[i] == 382 || element->tag[i] == 381) && update->ntimestep > 11700) printf("atoms in stencil for element %d timestep  = %d me = %d cutoff = %g  box = \n%g %g\n%g %g\n%g %g\n",element->tag[i],update->ntimestep,comm->me,cutoff,
      //    iboxlo[0]-cutoff,iboxhi[0]+cutoff,
      //    iboxlo[1]-cutoff,iboxhi[1]+cutoff,
      //    iboxlo[2]-cutoff,iboxhi[2]+cutoff);
      for (k = 0; k < nelemstencil; k++) {
        if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->matombins) {
          printf("Step = %d me = %d ibin = %d elemstencil[%d] = %d matombins = %d\n",
              update->ntimestep,comm->me,ibin,k,elemstencil[k],nb->matombins);
          error->one(FLERR,"TEST"); 
        }

        for (j = atombinhead[ibin+elemstencil[k]]; j >= 0; j = atombins[j]) {
          if (j > atom->nlocal+atom->nghost) {
            printf("Step = %d me = %d atom j = %d nlocal = %d nghost = %d\n",
                update->ntimestep,comm->me,j,atom->nlocal,atom->nghost);
            error->one(FLERR,"TEST"); 

          }

          jctype = atype[j];
          cutoff = cutneigh[ictype][jctype];
          //if ((element->tag[i] == 382 || element->tag[i] == 381) && update->ntimestep > 11700) 
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("%d %d %g %g %g structure = %g group check = %d eid = %d cutoff\n",
          //      atom->tag[j],atom->type[j],atom->x[j][0],atom->x[j][1],atom->x[j][2],filter_peratom[j],!(amask[j] & groupbit),element->tag[i],cutoff);
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("%d %d %g %g %g structure = %g group check = %d eid = %d cutoff\n",
          //      atom->tag[j],atom->type[j],atom->x[j][0],atom->x[j][1],atom->x[j][2],filter_peratom[j],!(amask[j] & groupbit),element->tag[i],cutoff);

          if (group && !(amask[j] & groupbit)) continue;
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("Test 1 381\n");
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("Test 1 382\n");

          if (atom_filter && atom_filter_value == filter_peratom[j])
            continue;
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("Test 2 381\n");
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("Test 2 382\n");

          if (iboxlo[0] > ax[j][0] + cutoff ||
              iboxhi[0] < ax[j][0] - cutoff)
            continue;
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("Test 3 381\n");
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("Test 3 382\n");

          if (iboxlo[1] > ax[j][1] + cutoff ||
              iboxhi[1] < ax[j][1] - cutoff)
            continue;
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("Test 4 381\n");
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("Test 4 382\n");

          if (iboxlo[2] > ax[j][2] + cutoff ||
              iboxhi[2] < ax[j][2] - cutoff)
            continue;
          //if (atom->tag[j] == 51666 && element->tag[i] == 381 && update->ntimestep > 11700) 
          //  printf("Test 5 381\n");
          //if (atom->tag[j] == 54606 && element->tag[i] == 382 && update->ntimestep > 11700) 
          //  printf("Test 5 382\n");

          e2aneighptr[ne2a++] = j;
          //printf("j = %d structure = %g\n",atom->tag[j],filter_peratom[j]);
        }
      }
      firstneighe2a[einum] = e2aneighptr;
      //if (element->tag[i] == 381 || element->tag[i] == 382) 
      //  printf("Number of neighbor for element %d is %d\n",element->tag[i],ne2a);
      numneighe2a[einum] = ne2a;
      e2apage->vgot(ne2a);
      if (e2apage->status()) {
        char *errms = new char[100];
        sprintf(errms,"Neighbor list overflow for element ID %d with # = %d, boost neigh_modify one",element->tag[i],ne2a);
        error->one(FLERR,errms);
      }
    } else {

      // no atoms (local + ghost)

      numneighe2a[einum] = 0;
      firstneighe2a[einum] = NULL;
    }

    eilist[einum] = i;
    einum++;
  }

  list->einum = einum;
}
