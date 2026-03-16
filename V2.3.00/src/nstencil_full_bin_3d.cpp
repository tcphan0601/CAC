#include "nstencil_full_bin_3d.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NStencilFullBin3d::NStencilFullBin3d(CAC *cac) : NStencil(cac) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullBin3d::create()
{
  int i,j,k;
  nstencil = 0;
  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}
