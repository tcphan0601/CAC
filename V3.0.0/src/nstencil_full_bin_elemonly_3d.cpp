#include "nstencil_full_bin_elemonly_3d.h"
#include "atom.h"
#include "element.h"
#include "memory.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "nbin.h"

using namespace CAC_NS;
#define BIG           1e30

/*  ----------------------------------------------------------------------  */

NStencilFullBinElemonly3d::NStencilFullBinElemonly3d(CAC *cac) : NStencil(cac) {}

/*  ----------------------------------------------------------------------
   create element stencil based on bin geometry and cutoff
-------------------------------------------------------------------------  */

void NStencilFullBinElemonly3d::create()
{
  int i, j, k;

  // create element stencil

  if (element->nelements) {
    nelemstencil = 0;
    for (k = -elemsz; k <= elemsz; k++)
      for (j = -elemsy; j <= elemsy; j++)
        for (i = -elemsx; i <= elemsx; i++)
          if (elembin_distance(i, 0) < cutelemx &&
              elembin_distance(j, 1) < cutelemy &&
              elembin_distance(k, 2) < cutelemz) 
            elemstencil[nelemstencil++] = k * melembiny * melembinx + j * melembinx + i;
          
  }
}

/*  ----------------------------------------------------------------------
   specific setup for elemonly 
   create element stencil only
   insure NBin data is current
   insure stencils are allocated large enough
   -------------------------------------------------------------------------  */

void NStencilFullBinElemonly3d::create_setup()
{
  if (nb) copy_bin_info();
  last_stencil = update->ntimestep;

  // if triclinic box, recalculate max_box in box coords since element->max_element_bound_box_size is in lamda coords
  
  if (domain->triclinic) {
    domain->lamda2nodex(element->nlocal, element->nodex);
    double maxcoord, mincoord, max[3], max_box[3];
    max[2] = 0.0;
    for (int i = 0; i < element->nlocal; i++) {
      double ***inodex = element->nodex[i];
      int inpe = element->npe[element->etype[i]];
      int iapc = element->apc[element->etype[i]];
      for (int dim = 0; dim < dimension; dim++) {
        maxcoord = -BIG;
        mincoord = BIG;
        for (int j = 0; j < iapc; j++) 
          for (int k = 0; k < inpe; k++) {
            maxcoord = MAX(maxcoord, inodex[j][k][dim]);
            mincoord = MIN(mincoord, inodex[j][k][dim]);
          }
        max[dim] = MAX(maxcoord-mincoord, max[dim]);
      }
    }
    MPI_Allreduce(max, max_box, 3, MPI_DOUBLE, MPI_MAX, world);
    domain->nodex2lamda(element->nlocal, element->nodex);
    cutelemx = cutneighmaxatom + max_box[0];
    cutelemy = cutneighmaxatom + max_box[1];
    cutelemz = cutneighmaxatom + max_box[2];
  } else {
    cutelemx = cutneighmaxatom + element->max_element_bound_box_size[0];
    cutelemy = cutneighmaxatom + element->max_element_bound_box_size[1];
    cutelemz = cutneighmaxatom + element->max_element_bound_box_size[2];
  }

  // sx, sy, sz = max range of atom stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil will be empty if cutneighmax = 0.0

  // element stencil

  if (element->nelements) {
    elemsx = static_cast<int> (cutelemx * elembininvx);
    if (elemsx * elembinsizex < cutelemx) elemsx++;
    elemsy = static_cast<int> (cutelemy * elembininvy);
    if (elemsy * elembinsizey < cutelemy) elemsy++;

    if (dimension == 3) {
      elemsz = static_cast<int> (cutelemz * elembininvz);
      if (elemsz * elembinsizez < cutelemz) elemsz++;
    } else elemsz = 0;

    int elemsmax = (2 * elemsx+1) * (2 * elemsy+1) * (2 * elemsz+1);

    // reallocate stencil structs if necessary
    // for BIN and MULTI styles

    if (neighstyle == Neighbor::BIN) {
      if (elemsmax > maxelemstencil) {
        maxelemstencil = elemsmax;
        memory->destroy(elemstencil);
        memory->create(elemstencil, maxelemstencil, "neighstencil:elemstencil");
        if (xyzflag) {
          memory->destroy(elemstencilxyz);
          memory->create(elemstencilxyz, maxelemstencil, 3, "neighstencil:elemstencilxyz");
        }
      }
    } else {
      int i;
      int n = atom->ntypes;
      if (maxelemstencil_multi == 0) {
        nelemstencil_multi = new int[n+1];
        elemstencil_multi = new int * [n+1];
        elemdistsq_multi = new double * [n+1];
        for (i = 1; i <= n; i++) {
          nelemstencil_multi[i] = 0;
          elemstencil_multi[i] = NULL;
          elemdistsq_multi[i] = NULL;
        }
      }
      if (elemsmax > maxelemstencil_multi) {
        maxelemstencil_multi = elemsmax;
        for (i = 1; i <= n; i++) {
          memory->destroy(elemstencil_multi[i]);
          memory->destroy(elemdistsq_multi[i]);
          memory->create(elemstencil_multi[i], maxelemstencil_multi, 
              "neighstencil:elemstencil_multi");
          memory->create(elemdistsq_multi[i], maxelemstencil_multi, 
              "neighstencil:elemdistsq_multi");
        }
      }
    }
  }
}

