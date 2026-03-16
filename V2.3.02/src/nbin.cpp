#include "nbin.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "domain.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NBin::NBin(CAC *cac) : Pointers(cac)
{
  last_bin = -1;
  maxbin = maxatom = maxelem = 0;
  atombinhead = elembinhead = NULL;
  atombins = elembins = NULL;
  atom2bin = elem2bin = NULL;

  // geometry settings

  dimension = domain->dimension;
  triclinic = domain->triclinic;
}

/* ---------------------------------------------------------------------- */

NBin::~NBin()
{
  memory->destroy(atombinhead);
  memory->destroy(elembinhead);
  memory->destroy(atombins);
  memory->destroy(elembins);
  memory->destroy(atom2bin);
  memory->destroy(elem2bin);
}

/* ---------------------------------------------------------------------- */

void NBin::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class
------------------------------------------------------------------------- */

void NBin::copy_neighbor_info()
{
  cutneighmin = neighbor->cutneighmin;
  cutneighmax = neighbor->cutneighmax;
  binsizeflag = neighbor->binsizeflag;
  binsize_user = neighbor->binsize_user;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // overwrite Neighbor cutoff with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) cutneighmax = cutoff_custom;
}

/* ----------------------------------------------------------------------
   setup for bin_all()
------------------------------------------------------------------------- */

void NBin::bin_setup(int naall,int neall)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for USER-INTEL package

  if (mbins > maxbin) {
    maxbin = mbins;
    memory->destroy(atombinhead);
    memory->destroy(elembinhead);
    memory->create(atombinhead,maxbin,"neigh:atombinhead");
    memory->create(elembinhead,maxbin,"neigh:elebinhead");
  }

  // atombins and atom2bin = per-atom vectors
  // for both local and ghost atoms

  if (naall > maxatom) {
    maxatom = naall;
    memory->destroy(atombins);
    memory->create(atombins,maxatom,"neigh:atombins");
    memory->destroy(atom2bin);
    memory->create(atom2bin,maxatom,"neigh:atom2bin");
  }

  // elembins and elem2bin = per-elem vectors
  // for both local and ghost elements

  if (neall > maxelem) {
    maxelem = neall;
    memory->destroy(elembins);
    memory->create(elembins,maxelem,"neigh:elembins");
    memory->destroy(elem2bin);
    memory->create(elem2bin,maxelem,"neigh:elem2bin");
  }

}

/* ----------------------------------------------------------------------
   convert obj coords into local bin #
   for orthogonal, only ghost objs will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost objs = nbin,nbin+1,etc
     owned objs = 0 to nbin-1
     lo ghost objs = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of objs straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NBin::coord2bin(double *x)
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

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/* ---------------------------------------------------------------------- */

bigint NBin::memory_usage()
{
  bigint bytes = 0;
  bytes += 2*maxbin*sizeof(int);
  bytes += 2*maxatom*sizeof(int);
  bytes += 2*maxelem*sizeof(int);
  return bytes;
}
