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
  : Pointers(cac), nb(NULL), ns(NULL), atombins(NULL), elembins(NULL), atomstencil(NULL), elemstencil(NULL)
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
  natomstencil = ns->natomstencil;
  atomstencil = ns->atomstencil;
  atomstencilxyz = ns->atomstencilxyz;
  natomstencil_multi = ns->natomstencil_multi;
  atomstencil_multi = ns->atomstencil_multi;
  atomdistsq_multi = ns->atomdistsq_multi;

  nelemstencil = ns->nelemstencil;
  elemstencil = ns->elemstencil;
  elemstencilxyz = ns->elemstencilxyz;
  nelemstencil_multi = ns->nelemstencil_multi;
  elemstencil_multi = ns->elemstencil_multi;
  elemdistsq_multi = ns->elemdistsq_multi;


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


