#include "nstencil_full_ghost_bin_3d.h"
#include "atom.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NStencilFullGhostBin3d::NStencilFullGhostBin3d(CAC *cac) : NStencil(cac)
{
  xyzflag = 1;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullGhostBin3d::create()
{
  int i,j,k;

  // create atom stencil

  if (atom->natoms) {
    natomstencil = 0;

    for (k = -atomsz; k <= atomsz; k++)
      for (j = -atomsy; j <= atomsy; j++)
        for (i = -atomsx; i <= atomsx; i++)
          if (atombin_distance(i,j,k) < cutneighmaxatomsq) {
            atomstencilxyz[natomstencil][0] = i;
            atomstencilxyz[natomstencil][1] = j;
            atomstencilxyz[natomstencil][2] = k;
            atomstencil[natomstencil++] = k*matombiny*matombinx + j*matombinx + i;
          }
  }

  // create element stencil

  if (element->nelements) {
    nelemstencil = 0;

    for (k = -elemsz; k <= elemsz; k++)
      for (j = -elemsy; j <= elemsy; j++)
        for (i = -elemsx; i <= elemsx; i++)
          if (elembin_distance(i,0) < cutelemx &&
              elembin_distance(j,1) < cutelemy &&
              elembin_distance(k,2) < cutelemz) {
            elemstencilxyz[nelemstencil][0] = i;
            elemstencilxyz[nelemstencil][1] = j;
            elemstencilxyz[nelemstencil][2] = k;
            elemstencil[nelemstencil++] = k*melembiny*melembinx + j*melembinx + i;
          }
  }

}
