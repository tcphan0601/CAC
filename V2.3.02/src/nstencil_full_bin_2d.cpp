#include "nstencil_full_bin_2d.h"
#include "neighbor.h"
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

  nstencil = 0;

  for (j = -sy; j <= sy; j++)
    for (i = -sx; i <= sx; i++)
      if (bin_distance(i,j,0) < cutneighmaxsq)
        stencil[nstencil++] = j*mbinx + i;
}
