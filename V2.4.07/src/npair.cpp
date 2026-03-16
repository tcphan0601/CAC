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
  groupbit = 0;
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
  if (nrq->group) groupbit = nrq->groupbit;
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
  cutneigh = neighbor->cutneigh;

  //cutneighghostsq = neighbor->cutneighghostsq;

  //cut_allsq = neighbor->cut_allsq;
  //cut_elemsq = neighbor->cut_elemsq;
  //cut_subelemsq = neighbor->cut_subelemsq;
  
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

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
    cutneighghostsq = mycutneighsq;
  }
}

/* ----------------------------------------------------------------------
   copy info from NBin class to this build class
   ------------------------------------------------------------------------- */

void NPair::copy_bin_info()
{
  natombinx = nb->natombinx;
  natombiny = nb->natombiny;
  natombinz = nb->natombinz;
  matombins = nb->matombins;
  matombinx = nb->matombinx;
  matombiny = nb->matombiny;
  matombinz = nb->matombinz;
  matombinxlo = nb->matombinxlo;
  matombinylo = nb->matombinylo;
  matombinzlo = nb->matombinzlo;

  atombininvx = nb->atombininvx;
  atombininvy = nb->atombininvy;
  atombininvz = nb->atombininvz;

  atom2bin = nb->atom2bin;
  atombins = nb->atombins;
  atombinhead = nb->atombinhead;

  nelembinx = nb->nelembinx;
  nelembiny = nb->nelembiny;
  nelembinz = nb->nelembinz;
  melembins = nb->melembins;
  melembinx = nb->melembinx;
  melembiny = nb->melembiny;
  melembinz = nb->melembinz;
  melembinxlo = nb->melembinxlo;
  melembinylo = nb->melembinylo;
  melembinzlo = nb->melembinzlo;

  elembininvx = nb->elembininvx;
  elembininvy = nb->elembininvy;
  elembininvz = nb->elembininvz;

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

/* ----------------------------------------------------------------------
   same as coord2atombin in Nbin, but also return ix,iy,iz offsets in each dim
   used by some of the ghost neighbor lists
   ------------------------------------------------------------------------- */

int NPair::coord2atombin(double *x, int &ix, int &iy, int &iz)
{
  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*atombininvx) + natombinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*atombininvx);
    ix = MIN(ix,natombinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*atombininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*atombininvy) + natombiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*atombininvy);
    iy = MIN(iy,natombiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*atombininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*atombininvz) + natombinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*atombininvz);
    iz = MIN(iz,natombinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*atombininvz) - 1;

  ix -= matombinxlo;
  iy -= matombinylo;
  iz -= matombinzlo;
  return iz*matombiny*matombinx + iy*matombinx + ix;
}

/* ----------------------------------------------------------------------
   same as coord2elembin in Nbin, but also return ix,iy,iz offsets in each dim
   used by some of the ghost neighbor lists
   ------------------------------------------------------------------------- */

int NPair::coord2elembin(double *x, int &ix, int &iy, int &iz)
{
  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*elembininvx) + nelembinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*elembininvx);
    ix = MIN(ix,nelembinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*elembininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*elembininvy) + nelembiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*elembininvy);
    iy = MIN(iy,nelembiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*elembininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*elembininvz) + nelembinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*elembininvz);
    iz = MIN(iz,nelembinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*elembininvz) - 1;

  ix -= melembinxlo;
  iy -= melembinylo;
  iz -= melembinzlo;
  return iz*melembiny*melembinx + iy*melembinx + ix;
}

