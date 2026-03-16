#include "nstencil_half_bin_2d.h"
#include "atom.h"
#include "element.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

NStencilHalfBin2d::NStencilHalfBin2d(CAC *cac) : NStencil(cac) {}

/*  ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
-------------------------------------------------------------------------  */

void NStencilHalfBin2d::create()
{
  int i, j;

  // create atom stencil

  if (atom->natoms) {
    natomstencil = 0;
    for (j = 0; j <= atomsy; j++)
      for (i = -atomsx; i <= atomsx; i++)
        if (j > 0 || (j == 0 && i > 0))
          if (atombin_distance(i, j, 0) < cutneighmaxatomsq)
            atomstencil[natomstencil++] = j * matombinx + i;
  }

  // create element stencil

  if (!atomonly && element->nelements) { 
    nelemstencil = 0;

    for (j = -elemsy; j <= elemsy; j++)
      for (i = -elemsx; i <= elemsx; i++)
        if (elembin_distance(i, 0) < cutelemx && 
            elembin_distance(j, 1) < cutelemy) 
          elemstencil[nelemstencil++] = j * melembinx + i;
  }
}
