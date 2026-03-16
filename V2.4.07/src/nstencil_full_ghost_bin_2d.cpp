#include "nstencil_full_ghost_bin_2d.h"
#include "atom.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NStencilFullGhostBin2d::NStencilFullGhostBin2d(CAC *cac) : NStencil(cac)
{
  xyzflag = 1;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullGhostBin2d::create()
{
  int i,j;

  // create atom stencil

  if (atom->natoms) {
    natomstencil = 0;

    for (j = -atomsy; j <= atomsy; j++)
      for (i = -atomsx; i <= atomsx; i++)
        if (atombin_distance(i,j,0) < cutneighmaxatomsq) {
          atomstencilxyz[natomstencil][0] = i;
          atomstencilxyz[natomstencil][1] = j;
          atomstencilxyz[natomstencil][2] = 0;
          atomstencil[natomstencil++] = j*matombinx + i;
        }
  }

  // create element stencil

  if (element->nelements) { 
    nelemstencil = 0;

    for (j = -elemsy; j <= elemsy; j++)
      for (i = -elemsx; i <= elemsx; i++)
        if (elembin_distance(i,0) < cutelemx && 
            elembin_distance(j,1) < cutelemy) {
          elemstencilxyz[nelemstencil][0] = i;
          elemstencilxyz[nelemstencil][1] = j;
          elemstencilxyz[nelemstencil][2] = 0;
          elemstencil[nelemstencil++] = j*melembinx + i;
        }
  }

}
