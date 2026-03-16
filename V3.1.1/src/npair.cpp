#include <math.h>
#include "npair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "nstencil.h"
#include "atom.h"
#include "element.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "my_page.h"

using namespace CAC_NS;
#define EPSILON 1.0e-4
#define BIG 1.0e30
/*  ----------------------------------------------------------------------  */

NPair::NPair(CAC *cac)
: Pointers(cac), nb(nullptr), ns(nullptr), atombins(nullptr), elembins(nullptr), atomstencil(nullptr), elemstencil(nullptr)
{
  last_build = -1;
  mycutneigh = nullptr;
  mycutneighsq = nullptr;
  myinnercutneighsq = nullptr;
  groupbit = 0;
  group = 0;
  noinner = 0;
  force_rebuild = 0;
  custom_domain = 0;

  atom_filter = 0;
  atom_filter_value = -1;
  atom_filter_preneighflag = 0;
  atom_filter_compute = nullptr;

  maxelemneigh = 20;
  memory->create(elemneighlist, maxelemneigh, "npair:elemneighlist");
}

/*  ----------------------------------------------------------------------  */

NPair::~NPair()
{
  memory->destroy(mycutneigh);
  memory->destroy(mycutneighsq);
  memory->destroy(myinnercutneighsq);
  memory->destroy(elemneighlist);
}

/*  ----------------------------------------------------------------------  */

void NPair::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) {
    cutoff_custom = nrq->cutoff;
  }
  if (nrq->group) {
    group = nrq->group;
    groupbit = nrq->groupbit;
  }
  noinner = nrq->noinner; 
  force_rebuild = nrq->force_rebuild;
  custom_domain = nrq->custom_domain;
  if (nrq->atom_filter) {
    atom_filter = 1;
    atom_filter_value = nrq->atom_filter_value;
    atom_filter_compute = modify->compute[nrq->atom_filter_icompute];
    atom_filter_preneighflag = nrq->atom_filter_preneighflag;
  }
}

/*  ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   done once per run
   -------------------------------------------------------------------------  */

void NPair::copy_neighbor_info()
{
  int n = atom->ntypes;
  int i, j;


  // general params

  exclude = neighbor->exclude;
  skin = neighbor->skin;
  cutneigh = neighbor->cutneigh;
  cutneighsq = neighbor->cutneighsq;
  cutneighghostsq = neighbor->cutneighghostsq;


  // exclusion info

  nex_type = neighbor->nex_type;
  ex1_type = neighbor->ex1_type;
  ex2_type = neighbor->ex2_type;
  ex_type = neighbor->ex_type;

  nex_group = neighbor->nex_group;
  ex1_group = neighbor->ex1_group;
  ex2_group = neighbor->ex2_group;
  ex1_bit = neighbor->ex1_bit;
  ex2_bit = neighbor->ex2_bit;

  //cutneighghostsq = neighbor->cutneighghostsq;

  //cut_allsq = neighbor->cut_allsq;
  //cut_elemsq = neighbor->cut_elemsq;
  //cut_subelemsq = neighbor->cut_subelemsq;

  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // overwrite per-type Neighbor cutoffs with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {

    memory->destroy(mycutneigh);
    memory->destroy(mycutneighsq);
    memory->create(mycutneigh, n+1, n+1, "npair:cutneigh");
    memory->create(mycutneighsq, n+1, n+1, "npair:cutneighsq");
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++) {
        mycutneigh[i][j] = cutoff_custom;
        mycutneighsq[i][j] = cutoff_custom * cutoff_custom;
      }
    cutneigh = mycutneigh;
    cutneighsq = mycutneighsq;
    //cutneighghostsq = mycutneighsq;
  }

  maxcutneigh = -BIG;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      maxcutneigh = MAX(maxcutneigh, cutneigh[i][j]);

}

/*  ----------------------------------------------------------------------
   copy info from NBin class to this build class
   -------------------------------------------------------------------------  */

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

/*  ----------------------------------------------------------------------
   copy info from NStencil class to this build class
   -------------------------------------------------------------------------  */

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

/*  ----------------------------------------------------------------------
   copy info from NBin and NStencil classes to this build class
   -------------------------------------------------------------------------  */

void NPair::build_setup()
{
  if (nb) copy_bin_info();
  if (ns) copy_stencil_info();

  // set here, since build_setup() always called before build()

  last_build = update->ntimestep;
}

/*  ----------------------------------------------------------------------
   test if atom pair i, j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
   -------------------------------------------------------------------------  */

int NPair::exclusion(int itype, int jtype, int imask, int jmask) const 
{
  if (nex_type && ex_type[itype][jtype]) return 1;

  if (nex_group) {
    for (int m = 0; m < nex_group; m++) {
      if (imask & ex1_bit[m] && jmask & ex2_bit[m]) return 1;
      if (imask & ex2_bit[m] && jmask & ex1_bit[m]) return 1;
    }
  }

  return 0;
}

/*  ----------------------------------------------------------------------
   same as coord2atombin in Nbin, but also return ix, iy, iz offsets in each dim
   used by some of the ghost neighbor lists
   -------------------------------------------------------------------------  */

int NPair::coord2atombin(double *x, int &ix, int &iy, int &iz)
{
  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
    error->one(FLERR, "Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0]) * atombininvx) + natombinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0]) * atombininvx);
    ix = MIN(ix, natombinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0]) * atombininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1]) * atombininvy) + natombiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1]) * atombininvy);
    iy = MIN(iy, natombiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1]) * atombininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2]) * atombininvz) + natombinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2]) * atombininvz);
    iz = MIN(iz, natombinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2]) * atombininvz) - 1;

  ix -= matombinxlo;
  iy -= matombinylo;
  iz -= matombinzlo;
  return iz * matombiny * matombinx + iy * matombinx + ix;
}

/*  ----------------------------------------------------------------------
   same as coord2elembin in Nbin, but also return ix, iy, iz offsets in each dim
   used by some of the ghost neighbor lists
   -------------------------------------------------------------------------  */

int NPair::coord2elembin(double *x, int &ix, int &iy, int &iz)
{
  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
    error->one(FLERR, "Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0]) * elembininvx) + nelembinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0]) * elembininvx);
    ix = MIN(ix, nelembinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0]) * elembininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1]) * elembininvy) + nelembiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1]) * elembininvy);
    iy = MIN(iy, nelembiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1]) * elembininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2]) * elembininvz) + nelembinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2]) * elembininvz);
    iz = MIN(iz, nelembinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2]) * elembininvz) - 1;

  ix -= melembinxlo;
  iy -= melembinylo;
  iz -= melembinzlo;
  return iz * melembiny * melembinx + iy * melembinx + ix;
}

/*  ----------------------------------------------------------------------
   grow temporary element neighbor list 
   -------------------------------------------------------------------------  */

void NPair::grow_elemneighlist()
{
  maxelemneigh += maxelemneigh / 2;
  memory->grow(elemneighlist, maxelemneigh, "npair:elemneighlist");
}

/*  ----------------------------------------------------------------------
   perform bucket sort on neighbor list up to sortnum
   -------------------------------------------------------------------------  */

void NPair::sort_neigh(int n, int sortnum, int *neighptr, int *neighindexptr, 
    double *rbucket, int *ibucket, int *indexbucket)
{
  int i, j;

  // if not enough requested neighbors sort everything

  int numneighsort = MIN(n, sortnum);

  // perform bucket sort to find nearest neighbors

  int permute[n];
  for (i = 0; i < n; i++)
    permute[i] = i;
  for (i = 0; i < numneighsort; i++) {

    // after this loop, permute[j] (index in bucket) will be the next nearest neighbor
    // sort from bottom up

    for (j = n-1; j > i; j--) {
      int current = permute[j];
      int next = permute[j-1];

      if (rbucket[current] < rbucket[next]) {
        permute[j] = next;
        permute[j-1] = current;
      }
    }
    neighptr[i] = ibucket[permute[i]];
    if (neighindexptr) neighindexptr[i] = indexbucket[permute[i]];
  }

  // fill in the rest of the neighbors without sorting

  for (i = numneighsort; i < n; i++) {
    neighptr[i] = ibucket[permute[i]];
    if (neighindexptr) neighindexptr[i] = indexbucket[permute[i]];
  }
}
