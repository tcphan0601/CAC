#include <math.h>
#include "npair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "nstencil.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "my_page.h"

using namespace CAC_NS;
#define EPSILON 1.0e-4
/* ---------------------------------------------------------------------- */

NPair::NPair(CAC *cac)
  : Pointers(cac), nb(NULL), ns(NULL), atombins(NULL), elembins(NULL), stencil(NULL)
{
  last_build = -1;
  mycutneighsq = NULL;
}

/* ---------------------------------------------------------------------- */

NPair::~NPair()
{
  memory->destroy(mycutneighsq);
}

/* ---------------------------------------------------------------------- */

void NPair::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   done once per run
------------------------------------------------------------------------- */

void NPair::copy_neighbor_info()
{

  // general params

  skin = neighbor->skin;
  cutneighsq = neighbor->cutneighsq;
  cutneighghostsq = neighbor->cutneighghostsq;

  cut_allsq = neighbor->cut_allsq;
  cut_elemsq = neighbor->cut_elemsq;
  cut_subelemsq = neighbor->cut_subelemsq;

  zeroes = neighbor->zeroes;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // special info

  special_flag = neighbor->special_flag;

  // overwrite per-type Neighbor cutoffs with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {
    memory->destroy(mycutneighsq);
    int n = atom->ntypes;
    memory->create(mycutneighsq,n+1,n+1,"npair:cutneighsq");
    int i,j;
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        mycutneighsq[i][j] = cutoff_custom * cutoff_custom;
    cutneighsq = mycutneighsq;
  }
}

/* ----------------------------------------------------------------------
   copy info from NBin class to this build class
------------------------------------------------------------------------- */

void NPair::copy_bin_info()
{
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;

  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;

  atom2bin = nb->atom2bin;
  atombins = nb->atombins;
  atombinhead = nb->atombinhead;
  elem2bin = nb->elem2bin;
  elembins = nb->elembins;
  elembinhead = nb->elembinhead;

}

/* ----------------------------------------------------------------------
   copy info from NStencil class to this build class
------------------------------------------------------------------------- */

void NPair::copy_stencil_info()
{
  nstencil = ns->nstencil;
  stencil = ns->stencil;
  stencilxyz = ns->stencilxyz;

  nstencil_multi = ns->nstencil_multi;
  stencil_multi = ns->stencil_multi;
  distsq_multi = ns->distsq_multi;
  
}

/* ----------------------------------------------------------------------
   copy info from NBin and NStencil classes to this build class
------------------------------------------------------------------------- */

void NPair::build_setup()
{
  if (nb) copy_bin_info();
  if (ns) copy_stencil_info();

  // set here, since build_setup() always called before build()

  last_build = update->ntimestep;
}

/* ----------------------------------------------------------------------
   same as coord2bin in Nbin, but also return ix,iy,iz offsets in each dim
   used by some of the ghost neighbor lists
------------------------------------------------------------------------- */

int NPair::coord2bin(double *x)
{
  int ix,iy,iz;
  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  return iz*mbiny*mbinx + iy*mbinx + ix;
}

