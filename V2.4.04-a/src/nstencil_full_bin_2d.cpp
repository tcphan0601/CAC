#include "nstencil_full_bin_2d.h"
#include "neighbor.h"
#include "atom.h"
#include "element.h"
#include "neigh_list.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NStencilFullBin2d::NStencilFullBin2d(CAC *cac) : NStencil(cac) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullBin2d::create()
{
  int i,j;

  // create atom stencil
  
  if (atom->natoms) {
  natomstencil = 0;

  for (j = -atomsy; j <= atomsy; j++)
    for (i = -atomsx; i <= atomsx; i++)
      if (atombin_distance(i,j,0) < cutneighmaxatomsq)
        atomstencil[natomstencil++] = j*matombinx + i;
  }

  // create element stencil
  
  if (element->nelements) { 
  nelemstencil = 0;

  for (j = -elemsy; j <= elemsy; j++)
    for (i = -elemsx; i <= elemsx; i++)
      if (elembin_distance(i,0) < cutelemx && 
          elembin_distance(j,1) < cutelemy) 
        elemstencil[nelemstencil++] = j*melembinx + i;
  }
}
