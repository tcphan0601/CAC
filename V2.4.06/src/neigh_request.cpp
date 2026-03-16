#include "neigh_request.h"
#include "atom.h"
#include "memory.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(CAC *cac) : Pointers(cac)
{
  // default ID = 0

  id = 0;

  // class user of list: default is pair request
  // only one is set to 1

  pair = 1;
  fix = compute = command = neigh = 0;

  // kind of list: default is half neighbor list
  // only one is set to 1

  half = 1;
  full = 0;

  // which list to include: defaults are atom, element, and integration point neighbor lists
  // NOTE: atomlist and atomonlylist can't be 1 at the same time
  
  atomlist = 1;
  atomonlylist = 0;
  elemlist = 1;
  intglist = 1;
  intpllist = 0;
  
  // attribute flags, mutiple can be set to 1
  // default is every reneighboring, not occasional
  // default is use newton_pair setting in force
  // default is no neighbors of ghosts
  // default is use cutoffs, not size of particles
  // default is no list-specific cutoff
  // default is no storage of auxiliary floating point values
  // default is no threebody
  // default is no sorting

  occasional = 0;
  newton = 0;
  ghost = 0;
  size = 0;
  cut = 0;
  cutoff = 0.0;
  threebody = 0.0;
  sort = 0;
  sortnum = 0;

  dnum = 0;

  command_style = NULL;

  // info set by Neighbor class when morphing original requests

  off2on = 0;
  copy = 0;
  copylist = -1;
  unique = 0;

  // internal settings

  index_bin = index_stencil = index_pair = -1;
}

/* ---------------------------------------------------------------------- */

NeighRequest::~NeighRequest()
{
//  delete [] iskip;
//  memory->destroy(ijskip);
}

/* ----------------------------------------------------------------------
   compare this request to other request
   identical means requestor identity and all params it sets are the same
   do not check other params that Neighbor can change after requests are made
   return 1 if identical, 0 if not
------------------------------------------------------------------------- */

int NeighRequest::identical(NeighRequest *other)
{
  int same = 1;

  // check for match of requestor_instance and instance counter
  // prevents an old fix from being unfix/refix in same memory location
  //   stored in requestor, and thus appearing old, when really new
  // only needed for classes with persistent neigh lists: Pair, Fix, Compute

  if (requestor != other->requestor) same = 0;
  if (requestor_instance != other->requestor_instance) same = 0;
  if (id != other->id) same = 0;

  // only compare settings made by requestors
  // not settings made later by Neighbor class
  
  if (pair != other->pair) same = 0;
  if (fix != other->fix) same = 0;
  if (compute != other->compute) same = 0;
  if (command != other->command) same = 0;
  if (neigh != other->neigh) same = 0;

  if (half != other->half) same = 0;
  if (full != other->full) same = 0;

  if (atomlist != other->atomlist) same = 0;
  if (elemlist != other->elemlist) same = 0;
  if (intglist != other->intglist) same = 0;
  if (intpllist != other->intpllist) same = 0;

  if (occasional != other->occasional) same = 0;
  if (newton != other->newton) same = 0;
  if (ghost != other->ghost) same = 0;
  if (size != other->size) same = 0;
//  if (history != other->history) same = 0;
//  if (granonesided != other->granonesided) same = 0;
  if (copy != other->copy) same = 0;
  if (cut != other->cut) same = 0;
  else if (cutoff != other->cutoff) same = 0;

  // if requested sortnum is less than other request, can still reuse it
  
  if (sort != other->sort) same = 0;
  else if (sortnum > other->sortnum) same = 0;

  if (dnum != other->dnum) same = 0;

//  if (skip != other->skip) same = 0;
//  if (skip) same = same_skip(other);

  return same;
}

/* ----------------------------------------------------------------------
   compare same info for two requests that have skip = 1
   return 1 if identical, 0 if not
------------------------------------------------------------------------- */
/*
int NeighRequest::same_skip(NeighRequest *other)
{
  int i,j;

  int ntypes = atom->ntypes;
  int same = 1;

  for (i = 1; i <= ntypes; i++)
    if (iskip[i] != other->iskip[i]) same = 0;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      if (ijskip[i][j] != other->ijskip[i][j]) same = 0;

  return same;
}
*/
/* ----------------------------------------------------------------------
   set params in this request to those of other request
   copy same fields that are checked in identical()
   purpose is to allow comparison of new requests to old requests
   skipflag = 1 to copy skip vector/array
------------------------------------------------------------------------- */

void NeighRequest::copy_request(NeighRequest *other)
{
  requestor = other->requestor;
  requestor_instance = other->requestor_instance;
  id = other->id;

  pair = other->pair;
  fix = other->fix;
  compute = other->compute;
  command = other->command;

  half = other->half;
  full = other->full;

  atomlist = other->atomlist;
  elemlist = other->elemlist;
  intglist = other->intglist;
  intpllist = other->intpllist;


  occasional = other->occasional;
  newton = other->newton;
  ghost = other->ghost;
  size = other->size;
//  history = other->history;
//  granonesided = other->granonesided;
  cut = other->cut;
  cutoff = other->cutoff;
  sort = other->sort;
  sortnum = other->sortnum;

  dnum = other->dnum;
/*
  iskip = NULL;
  ijskip = NULL;

  if (!skipflag) return;

  int i,j;
  int ntypes = atom->ntypes;

  if (other->iskip) {
    iskip = new int[ntypes+1];
    for (i = 1; i <= ntypes; i++)
      iskip[i] = other->iskip[i];
  }

  if (other->ijskip) {
    memory->create(ijskip,ntypes+1,ntypes+1,"neigh_request:ijskip");
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = other->ijskip[i][j];
  }*/
}

