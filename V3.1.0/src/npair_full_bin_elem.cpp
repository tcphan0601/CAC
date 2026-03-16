#include "npair_full_bin_elem.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
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
  neighptr = neighindexptr = nullptr;
  int *jlist;
  int n;
  double *iboxlo, *jboxlo, *iboxhi, *jboxhi;
  double xi[3];
  double cutoff, *length_scale;

  int naall = atom->nlocal + atom->nghost;
  int nelocal = element->nlocal;
  int *atype = atom->type;

  double **ax = atom->x;
  double **element_bound_box = element->element_bound_box;
  double **element_box2xi_scale = element->element_box2xi_scale;

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
    //if (group && !(emask[i] & groupbit)) continue;

    iboxlo = element_bound_box[i];
    iboxhi = element_bound_box[i] + 3;
    ibin = elem2bin[i];

    neighptr = ipage->vget();
    neighindexptr = indexpage->vget();
    n = 0;

    // loop over all atoms in all atom bins in stencil
   /*
      FILE *fp = fopen("dump_test_508.atom","w");
      fprintf(fp, "ITEM: TIMESTEP\n");

      fprintf(fp, "%d\n", comm->me);

      fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
      fprintf(fp, "%d\n", 24258);
      char boundstr[9];
      domain->boundary_string(boundstr);
      if (domain->triclinic == 0) {
        fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
        fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
        fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
        fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
      } else {
        fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
        fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0], domain->boxhi_bound[0], domain->xy);
        fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1], domain->boxhi_bound[1], domain->xz);
        fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2], domain->boxhi_bound[2], domain->yz);
      }
      fprintf(fp, "ITEM: ATOMS id type x y z structure_type tag jindex j\n");
      double coord[3];
      int count = 0;
      for (int ibasis = 0; ibasis < element->apc[element->etype[i]]; ibasis++) {
        for (int iucell = 0; iucell < element->nucell[element->etype[i]]; iucell++) {
          element->evec->interpolate(coord,element->nodex,i,ibasis,iucell,3);
          fprintf(fp,"%d %d %g %g %g %d 0 -1 -1\n",++count,1,coord[0],coord[1],coord[2],0);
        }
      }
   */ 

    if (naall) {

      // atom stencil is identical to element stencil so use element stencil here
      // filter_atom if needed
      
      cutoff = maxcutneigh;
      length_scale = element_box2xi_scale[i]; 
      for (k = 0; k < nelemstencil; k++) {
        if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->matombins) {
          error->one(FLERR, "TEST"); 
        }

        for (j = atombinhead[ibin+elemstencil[k]]; j >= 0; j = atombins[j]) {

          if (group && !(amask[j] & groupbit)) continue;
          //if (element->tag[i]==893)
          //  printf("i = %d atom_filter = %d atom_filter_value = %d filter_peratom[%d] = %g\n",element->tag[i],atom_filter,atom_filter_value,atom->tag[j],filter_peratom[j]);
          if (atom_filter && atom_filter_value == (int) filter_peratom[j])
            continue;

          // convert atom coords to element I xi coordinates 

          element->box2natural(ax[j], xi, i);
          int jindex;
          int tmp[3];
          if (xi[0] + cutoff * length_scale[0] < -1) continue;
          if (xi[0] - cutoff * length_scale[0] > 1) continue;
          if (xi[1] + cutoff * length_scale[1] < -1) continue;
          if (xi[1] - cutoff * length_scale[1] > 1) continue;
          if (xi[2] + cutoff * length_scale[2] < -1) continue;
          if (xi[2] - cutoff * length_scale[2] > 1) continue;

          for (int dim = 0; dim < 3; dim++) {
            if (xi[dim] < -1) tmp[dim] = 0;
            else if (xi[dim] <= 1) tmp[dim] = 1;
            else tmp[dim] = 2;
          }
          jindex = tmp[2] * 9 + tmp[1] * 3 + tmp[0];
//          fprintf(fp,"%d %d %g %g %g %g %d %d %d\n",++count,2,ax[j][0],ax[j][1],ax[j][2],filter_peratom[j],atom->tag[j],jindex,j);
//          printf(p,"%d %d %g %g %g %g %d %d %d\n",++count,2,ax[j][0],ax[j][1],ax[j][2],filter_peratom[j],atom->tag[j],jindex,j);
 //         if (jindex == 0 || jindex == 2 || jindex == 6 ||
 //             jindex == 8 || jindex == 13 || jindex == 18 || 
 //             jindex == 20 || jindex == 24 || jindex == 26) 
 //           continue;
          neighptr[n] = j;
          neighindexptr[n++] = jindex;
        }
      }
    }

    //fclose(fp);

    // loop over all elements in other element bins in stencil including self

    cutoff = maxcutneigh;
    for (k = 0; k < nelemstencil; k++) {
      if (ibin+elemstencil[k] < 0 || ibin+elemstencil[k] >= nb->melembins) {
        error->one(FLERR, "TEST"); 
      }

      for (j = elembinhead[ibin+elemstencil[k]]; j >= 0; j = elembins[j]) {

        if (group && !(emask[j] & groupbit)) continue;
        // check if bouning box overlap

        if (i == j) continue;

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
        neighindexptr[n++] = -1;
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

    ilist[inum] = i;
    iindexlist[inum++] = 0;
  }
  list->inum = list->einum = inum;
}
