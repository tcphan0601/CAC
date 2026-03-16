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
  maxatombin = maxelembin =  maxatom = maxelem = 0;
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
  atombinsizeflag = neighbor->atombinsizeflag;
  atombinsize_user = neighbor->atombinsize_user;
  elembinsizeflag = neighbor->elembinsizeflag;
  elembinsize_user = neighbor->elembinsize_user;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // overwrite Neighbor cutoff with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) cutneighmax = cutoff_custom;
}

/* ----------------------------------------------------------------------
   setup for bin_all()
------------------------------------------------------------------------- */

void NBin::bin_setup(int naall, int neall)
{
  // atombinhead = per-atombin vector, matombins in length

  if (matombins > maxatombin) {
    maxatombin = matombins;
    memory->destroy(atombinhead);
    memory->create(atombinhead,maxatombin,"neigh:atombinhead");
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

  // elembinhead = per-elembin vector, melembins in length
  
  if (melembins > maxelembin) {
    maxelembin = melembins;
    memory->destroy(elembinhead);
    memory->create(elembinhead,maxelembin,"neigh:elebinhead");
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
   convert obj coords into local atom bin #
   for orthogonal, only ghost objs will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost objs = nbin,nbin+1,etc
     owned objs = 0 to nbin-1
     lo ghost objs = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of objs straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NBin::coord2atombin(double x, double y, double z)
{
  int ix,iy,iz;

  if (!ISFINITE(x) || !ISFINITE(y) || !ISFINITE(z))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x >= bboxhi[0])
    ix = static_cast<int> ((x-bboxhi[0])*atombininvx) + natombinx;
  else if (x >= bboxlo[0]) {
    ix = static_cast<int> ((x-bboxlo[0])*atombininvx);
    ix = MIN(ix,natombinx-1);
  } else
    ix = static_cast<int> ((x-bboxlo[0])*atombininvx) - 1;

  if (y >= bboxhi[1])
    iy = static_cast<int> ((y-bboxhi[1])*atombininvy) + natombiny;
  else if (y >= bboxlo[1]) {
    iy = static_cast<int> ((y-bboxlo[1])*atombininvy);
    iy = MIN(iy,natombiny-1);
  } else
    iy = static_cast<int> ((y-bboxlo[1])*atombininvy) - 1;

  if (z >= bboxhi[2])
    iz = static_cast<int> ((z-bboxhi[2])*atombininvz) + natombinz;
  else if (z >= bboxlo[2]) {
    iz = static_cast<int> ((z-bboxlo[2])*atombininvz);
    iz = MIN(iz,natombinz-1);
  } else
    iz = static_cast<int> ((z-bboxlo[2])*atombininvz) - 1;

  return (iz-matombinzlo)*matombiny*matombinx + (iy-matombinylo)*matombinx + (ix-matombinxlo);
}

/* ----------------------------------------------------------------------
   convert obj coords into local element bin #
   for orthogonal, only ghost objs will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost objs = nbin,nbin+1,etc
     owned objs = 0 to nbin-1
     lo ghost objs = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of objs straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NBin::coord2elembin(double x, double y, double z)
{
  int ix,iy,iz;

  if (!ISFINITE(x) || !ISFINITE(y) || !ISFINITE(z))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x >= bboxhi[0])
    ix = static_cast<int> ((x-bboxhi[0])*elembininvx) + nelembinx;
  else if (x >= bboxlo[0]) {
    ix = static_cast<int> ((x-bboxlo[0])*elembininvx);
    ix = MIN(ix,nelembinx-1);
  } else
    ix = static_cast<int> ((x-bboxlo[0])*elembininvx) - 1;

  if (y >= bboxhi[1])
    iy = static_cast<int> ((y-bboxhi[1])*elembininvy) + nelembiny;
  else if (y >= bboxlo[1]) {
    iy = static_cast<int> ((y-bboxlo[1])*elembininvy);
    iy = MIN(iy,nelembiny-1);
  } else
    iy = static_cast<int> ((y-bboxlo[1])*elembininvy) - 1;

  if (z >= bboxhi[2])
    iz = static_cast<int> ((z-bboxhi[2])*elembininvz) + nelembinz;
  else if (z >= bboxlo[2]) {
    iz = static_cast<int> ((z-bboxlo[2])*elembininvz);
    iz = MIN(iz,nelembinz-1);
  } else
    iz = static_cast<int> ((z-bboxlo[2])*elembininvz) - 1;
  
  return (iz-melembinzlo)*melembiny*melembinx + (iy-melembinylo)*melembinx + (ix-melembinxlo);
}


/* ---------------------------------------------------------------------- */

bigint NBin::memory_usage()
{
  bigint bytes = 0;
  bytes += maxatombin*sizeof(int);
  bytes += 2*maxatom*sizeof(int);
  bytes += maxelembin*sizeof(int);
  bytes += 2*maxelem*sizeof(int);
  return bytes;
}
