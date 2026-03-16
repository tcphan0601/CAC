#include "nstencil_full_bin_3d.h"
#include "atom.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NStencilFullBin3d::NStencilFullBin3d(CAC *cac) : NStencil(cac) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullBin3d::create()
{
  int i,j,k;

  // create atom stencil

  if (atom->natoms) {
    natomstencil = 0;

    for (k = -atomsz; k <= atomsz; k++)
      for (j = -atomsy; j <= atomsy; j++)
        for (i = -atomsx; i <= atomsx; i++)
          if (atombin_distance(i,j,k) < cutneighmaxatomsq)
            atomstencil[natomstencil++] = k*matombiny*matombinx + j*matombinx + i;
  }

  // create element stencil

  if (!atomonly && element->nelements) {
    nelemstencil = 0;

    for (k = -elemsz; k <= elemsz; k++)
      for (j = -elemsy; j <= elemsy; j++)
        for (i = -elemsx; i <= elemsx; i++)
          if (elembin_distance(i,0) < cutelemx &&
              elembin_distance(j,1) < cutelemy &&
              elembin_distance(k,2) < cutelemz)
            elemstencil[nelemstencil++] = k*melembiny*melembinx + j*melembinx + i;
  }
}
