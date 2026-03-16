#include "nstencil.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace CAC_NS;

enum{NSQ,BIN,MULTI};     // also in Neighbor
enum{ATOM,ELEMENT};

/* ----------------------------------------------------------------------
   NStencil classes
   each has method to create a stencil = list of bin offsets
     invoked each time simulation box size/shape changes
     since induces change in bins
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
     calculated below in create_setup()
   3d creates xyz stencil, 2d creates xy stencil
   for half list with newton off:
     stencil is all surrounding bins including self
     regardless of triclinic
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
     no versions that allow ghost on (no callers need it?)
   for half list with newton on and triclinic:
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
     no versions that allow ghost on (no callers need it?)
   for full list:
     stencil is all surrounding bins including self
     regardless of newton on/off or triclinic
   for multi:
     create one stencil for each atom type
     stencil follows same rules for half/full, newton on/off, triclinic
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
     no versions that allow ghost on (any need for it?)
------------------------------------------------------------------------- */

NStencil::NStencil(CAC *cac) : Pointers(cac)
{
  last_stencil = -1;

  xyzflag = 0;
  maxatomstencil = maxatomstencil_multi = 0;
  atomstencil = NULL;
  atomstencilxyz = NULL;
  natomstencil_multi = NULL;
  atomstencil_multi = NULL;
  atomdistsq_multi = NULL;

  maxelemstencil = maxelemstencil_multi = 0;
  elemstencil = NULL;
  elemstencilxyz = NULL;
  nelemstencil_multi = NULL;
  elemstencil_multi = NULL;
  elemdistsq_multi = NULL;

  dimension = domain->dimension;
}

/* ---------------------------------------------------------------------- */

NStencil::~NStencil()
{
  memory->destroy(atomstencil);
  memory->destroy(atomstencilxyz);
  memory->destroy(elemstencil);
  memory->destroy(elemstencilxyz);

  int n = atom->ntypes;

  if (atomstencil_multi) {
    for (int i = 1; i <= n; i++) {
      memory->destroy(atomstencil_multi[i]);
      memory->destroy(atomdistsq_multi[i]);
    }
    delete [] natomstencil_multi;
    delete [] atomstencil_multi;
    delete [] atomdistsq_multi;
  }

  if (elemstencil_multi) {
    for (int i = 1; i <= n; i++) {
      memory->destroy(elemstencil_multi[i]);
      memory->destroy(elemdistsq_multi[i]);
    }
    delete [] nelemstencil_multi;
    delete [] elemstencil_multi;
    delete [] elemdistsq_multi;
  }
}

/* ---------------------------------------------------------------------- */

void NStencil::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this stencil class
   ------------------------------------------------------------------------- */

void NStencil::copy_neighbor_info()
{
  neighstyle = neighbor->style;
  cutneighmaxatom = neighbor->cutneighmax;
  cutneighmaxatomsq = neighbor->cutneighmaxsq;
  cuttypesq = neighbor->cuttypesq;

  // overwrite Neighbor cutoff with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {
    cutneighmaxatom = cutoff_custom;
    cutneighmaxatomsq = cutneighmaxatom * cutneighmaxatom;
  }

  cutelemx = cutneighmaxatom + element->max_bound_box_size[0];
  cutelemy = cutneighmaxatom + element->max_bound_box_size[1];
  cutelemz = cutneighmaxatom + element->max_bound_box_size[2];
}

/* ----------------------------------------------------------------------
   copy needed info from NBin class to this stencil class
   ------------------------------------------------------------------------- */

void NStencil::copy_bin_info()
{
  matombinx = nb->matombinx;
  matombiny = nb->matombiny;
  matombinz = nb->matombinz;
  atombinsizex = nb->atombinsizex;
  atombinsizey = nb->atombinsizey;
  atombinsizez = nb->atombinsizez;
  atombininvx = nb->atombininvx;
  atombininvy = nb->atombininvy;
  atombininvz = nb->atombininvz;

  melembinx = nb->melembinx;
  melembiny = nb->melembiny;
  melembinz = nb->melembinz;
  elembinsizex = nb->elembinsizex;
  elembinsizey = nb->elembinsizey;
  elembinsizez = nb->elembinsizez;
  elembininvx = nb->elembininvx;
  elembininvy = nb->elembininvy;
  elembininvz = nb->elembininvz;
}

/* ----------------------------------------------------------------------
   insure NBin data is current
   insure stencils are allocated large enough
   ------------------------------------------------------------------------- */

void NStencil::create_setup()
{
  if (nb) copy_bin_info();
  last_stencil = update->ntimestep;

  // sx,sy,sz = max range of atom stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil will be empty if cutneighmax = 0.0

  // atom stencil

  if (atom->natoms) {
    atomsx = static_cast<int> (cutneighmaxatom*atombininvx);
    if (atomsx*atombinsizex < cutneighmaxatom) atomsx++;
    atomsy = static_cast<int> (cutneighmaxatom*atombininvy);
    if (atomsy*atombinsizey < cutneighmaxatom) atomsy++;

    if (dimension == 3) {
      atomsz = static_cast<int> (cutneighmaxatom*atombininvz);
      if (atomsz*atombinsizez < cutneighmaxatom) atomsz++;
    } else atomsz = 0;

    int atomsmax = (2*atomsx+1) * (2*atomsy+1) * (2*atomsz+1);

    // reallocate stencil structs if necessary
    // for BIN and MULTI styles

    if (neighstyle == BIN) {
      if (atomsmax > maxatomstencil) {
        maxatomstencil = atomsmax;
        memory->destroy(atomstencil);
        memory->create(atomstencil,maxatomstencil,"neighstencil:atomstencil");
        if (xyzflag) {
          memory->destroy(atomstencilxyz);
          memory->create(atomstencilxyz,maxatomstencil,3,"neighstencil:atomstencilxyz");
        }
      }
    } else {
      int i;
      int n = atom->ntypes;
      if (maxatomstencil_multi == 0) {
        natomstencil_multi = new int[n+1];
        atomstencil_multi = new int*[n+1];
        atomdistsq_multi = new double*[n+1];
        for (i = 1; i <= n; i++) {
          natomstencil_multi[i] = 0;
          atomstencil_multi[i] = NULL;
          atomdistsq_multi[i] = NULL;
        }
      }
      if (atomsmax > maxatomstencil_multi) {
        maxatomstencil_multi = atomsmax;
        for (i = 1; i <= n; i++) {
          memory->destroy(atomstencil_multi[i]);
          memory->destroy(atomdistsq_multi[i]);
          memory->create(atomstencil_multi[i],maxatomstencil_multi,
              "neighstencil:atomstencil_multi");
          memory->create(atomdistsq_multi[i],maxatomstencil_multi,
              "neighstencil:atomdistsq_multi");
        }
      }
    }
  }

  // element stencil

  if (element->nelements) {
    elemsx = static_cast<int> (cutelemx*elembininvx);
    if (elemsx*elembinsizex < cutelemx) elemsx++;
    elemsy = static_cast<int> (cutelemy*elembininvy);
    if (elemsy*elembinsizey < cutelemy) elemsy++;

    if (dimension == 3) {
      elemsz = static_cast<int> (cutelemz*elembininvz);
      if (elemsz*elembinsizez < cutelemz) elemsz++;
    } else elemsz = 0;

    int elemsmax = (2*elemsx+1) * (2*elemsy+1) * (2*elemsz+1);

    // reallocate stencil structs if necessary
    // for BIN and MULTI styles

    if (neighstyle == BIN) {
      if (elemsmax > maxelemstencil) {
        maxelemstencil = elemsmax;
        memory->destroy(elemstencil);
        memory->create(elemstencil,maxelemstencil,"neighstencil:elemstencil");
        if (xyzflag) {
          memory->destroy(elemstencilxyz);
          memory->create(elemstencilxyz,maxelemstencil,3,"neighstencil:elemstencilxyz");
        }
      }
    } else {
      int i;
      int n = atom->ntypes;
      if (maxelemstencil_multi == 0) {
        nelemstencil_multi = new int[n+1];
        elemstencil_multi = new int*[n+1];
        elemdistsq_multi = new double*[n+1];
        for (i = 1; i <= n; i++) {
          nelemstencil_multi[i] = 0;
          elemstencil_multi[i] = NULL;
          elemdistsq_multi[i] = NULL;
        }
      }
      if (elemsmax > maxelemstencil_multi) {
        maxelemstencil_multi = elemsmax;
        for (i = 1; i <= n; i++) {
          memory->destroy(elemstencil_multi[i]);
          memory->destroy(elemdistsq_multi[i]);
          memory->create(elemstencil_multi[i],maxelemstencil_multi,
              "neighstencil:elemstencil_multi");
          memory->create(elemdistsq_multi[i],maxelemstencil_multi,
              "neighstencil:elemdistsq_multi");
        }
      }
    }
  }
}
/* ----------------------------------------------------------------------
   compute closest distance between central atom bin (0,0,0) and atom bin (i,j,k)
   ------------------------------------------------------------------------- */

double NStencil::atombin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*atombinsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*atombinsizex;

  if (j > 0) dely = (j-1)*atombinsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*atombinsizey;

  if (k > 0) delz = (k-1)*atombinsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*atombinsizez;

  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   compute closest distance between central element bin (0,0,0) and element bin (i,j,k)
   ------------------------------------------------------------------------- */

double NStencil::elembin_distance(int i, int dim)
{
  double distance;
  double binsize;
  if (dim == 0) binsize = elembinsizex;
  if (dim == 1) binsize = elembinsizey;
  if (dim == 2) binsize = elembinsizez;


  if (i > 0) distance = (i-1)*binsize;
  else if (i == 0) distance = 0.0;
  else distance = (i+1)*binsize;

  return distance;
}

/* ---------------------------------------------------------------------- */

bigint NStencil::memory_usage()
{
  bigint bytes = 0;
  if (neighstyle == BIN) {
    bytes += memory->usage(atomstencil,maxatomstencil);
    bytes += memory->usage(elemstencil,maxatomstencil);
    if (xyzflag) {
      bytes += memory->usage(atomstencilxyz,maxatomstencil,3);
      bytes += memory->usage(elemstencilxyz,maxatomstencil,3);
    }
  } else if (neighstyle == MULTI) {
    bytes += atom->ntypes*maxatomstencil_multi * sizeof(int);
    bytes += atom->ntypes*maxatomstencil_multi * sizeof(double);
    bytes += atom->ntypes*maxelemstencil_multi * sizeof(int);
    bytes += atom->ntypes*maxelemstencil_multi * sizeof(double);
  }
  return bytes;
}
