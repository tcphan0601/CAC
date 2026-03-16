#include "neigh_list.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "update.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define PGDELTA 1

/*  ----------------------------------------------------------------------  */

NeighList::NeighList(CAC *cac) : Pointers(cac)
{
  // initializations

  maxall = maxnva = 0; 

  inum = ainum = einum = ninum = ginum = vinum = nvainum = gnum = 0;
  ilist = iindexlist = NULL;
  nvailist = nvaindexlist = NULL;
  va2nvalist = NULL;

  numneigh = NULL;
  firstneigh = firstneighindex = NULL;
  ipage = indexpage = NULL;
  nvanumneigh = NULL;
  nvafirstneigh = nvafirstneighindex = NULL;
  nvaipage = nvaindexpage = NULL;

  // defaults, but may be reset by post_constructor()

  atomlist = 1;
  atomonlylist = 0;
  elemlist = 0;
  nodelist = 0;
  gausslist = 1;
  vatomlist = 0;

  occasional = 0;
  ghost = 0;
  copy = 0;
  sort = 0;
  sortnum = 0;
  dnum = 0;
  threebody = 0;

  group = 0;

  maxelem = 0;
  maxnva_elem = 0;
  neighneigh = 0;
  element_firstindex = NULL;
  // ptrs

  //iskip = NULL;
  //ijskip = NULL;

  listcopy = NULL;
  //listskip = NULL;
  listfull = NULL;

  //listhistory = NULL;
  //fix_history = NULL;

  //respamiddle = 0;
  //listinner = NULL;
  //listmiddle = NULL;

  //fix_bond = NULL;


  dpage = NULL;
}

/*  ----------------------------------------------------------------------  */

NeighList::~NeighList()
{

  //if (copymode) return;
  if (copy) return;
  memory->destroy(ilist);
  memory->destroy(iindexlist);

  memory->destroy(numneigh);
  memory->sfree(firstneigh);
  memory->sfree(firstneighindex);
  delete [] ipage;
  delete [] indexpage;

  if (threebody) {
    memory->destroy(nvailist);
    memory->destroy(nvaindexlist);
    memory->destroy(va2nvalist);
    memory->destroy(nvanumneigh);
    memory->sfree(nvafirstneigh);
    memory->sfree(nvafirstneighindex);
    delete [] nvaipage;
    delete [] nvaindexpage;
  }
}

/*  ----------------------------------------------------------------------
   adjust settings to match corresponding NeighRequest
   cannot do this in constructor b/c not all NeighLists are allocated yet
   copy -> set listcopy for list to copy from
   skip -> set listskip for list to skip from, create copy of itype, ijtype
   halffull -> set listfull for full list to derive from
   history -> set LH and FH ptrs in partner list that uses the history info
   respaouter -> set listinner/listmiddle for other rRESPA lists
   bond -> set fix_bond to Fix that made the request
   -------------------------------------------------------------------------  */

void NeighList::post_constructor(NeighRequest *rq)
{
  // copy request settings used by list itself

  occasional = rq->occasional;
  ghost = rq->ghost;
  atomlist = rq->atomlist;
  atomonlylist = rq->atomonlylist;
  elemlist = rq->elemlist;
  nodelist = rq->nodelist;
  gausslist = rq->gausslist;
  vatomlist = rq->vatomlist;

  threebody = rq->threebody;
  copy = rq->copy;
  sort = rq->sort;
  sortnum = rq->sortnum;
  group = rq->group;
  dnum = rq->dnum;

  if (rq->copy)
    listcopy = neighbor->lists[rq->copylist];

  if (rq->halffull)
    listfull = neighbor->lists[rq->halffulllist];

  //if (rq->skip) {
  //  listskip = neighbor->lists[rq->skiplist];
  //  int ntypes = atom->ntypes;
  //  iskip = new int[ntypes+1];
  //  memory->create(ijskip, ntypes+1, ntypes+1, "neigh_list:ijskip");
  //  int i, j;
  //  for (i = 1; i <= ntypes; i++) iskip[i] = rq->iskip[i];
  //  for (i = 1; i <= ntypes; i++)
  //    for (j = 1; j <= ntypes; j++)
  //      ijskip[i][j] = rq->ijskip[i][j];
  //}

}

/*  ----------------------------------------------------------------------  */

void NeighList::setup_pages(int pgsize_caller, int oneatom_caller)
{
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;

  int nmypage = comm->nthreads;
  ipage = new MyPage<int>[nmypage];
  indexpage = new MyPage<int>[nmypage];

  for (int i = 0; i < nmypage; i++) {
    ipage[i].init(oneatom, pgsize, PGDELTA);
    indexpage[i].init(oneatom, pgsize, PGDELTA);
  }

  if (threebody) {
    nvaipage = new MyPage<int>[nmypage];
    nvaindexpage = new MyPage<int>[nmypage];

    for (int i = 0; i < nmypage; i++) {
      nvaipage[i].init(oneatom, pgsize, PGDELTA);
      nvaindexpage[i].init(oneatom, pgsize, PGDELTA);
    }
  }

  if (dnum) {
    dpage = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++)
      dpage[i].init(dnum * oneatom, dnum * pgsize, PGDELTA);
  } else dpage = NULL;
}

/*  ----------------------------------------------------------------------
   grow per-atom data to allow for nlocal/nall atoms/gcells
   triggered by neighbor list build
   not called if a copy list
   -------------------------------------------------------------------------  */

void NeighList::grow(int nlocal, int nall)
{
  // grow element_firstindex if needed
  // check this before checking maxall
  
  if (neighneigh) 
    if (maxelem < element->nlocal) {
      maxelem = element->nmax;
      memory->destroy(element_firstindex);
      memory->create(element_firstindex, maxelem, "neighlist:element_firstindex");
    }

  if (ghost) {
    if (nall <= maxall) return;
  } else {
    if (nlocal <= maxall) return;
  } 
  maxall = 0;
  if (atomlist || atomonlylist) maxall += atom->nmax;
  if (elemlist) maxall += element->nmax;
  if (nodelist) maxall += element->nmax * element->maxnpe * element->maxapc;
  if (gausslist) maxall += element->nmax * element->maxgcell * element->maxapc;
  if (vatomlist) maxall += element->nmax * element->maxucell * element->maxapc;

  memory->destroy(ilist);
  memory->destroy(iindexlist);
  memory->create(ilist, maxall, "neighlist:ilist");
  memory->create(iindexlist, maxall, "neighlist:iindexlist");

  memory->destroy(numneigh);
  memory->create(numneigh, maxall, "neighlist:numneigh");
  memory->sfree(firstneigh);
  memory->sfree(firstneighindex);
  firstneigh = (int **) memory->smalloc(maxall * sizeof(int *), "neighlist:firstneigh");
  firstneighindex = (int **) memory->smalloc(maxall * sizeof(int *), "neighlist:firstneighindex");

}


/*  ----------------------------------------------------------------------
   grow per-nva data to allow for n nva
   triggered by neighbor list build with threebody 
   not called if a copy list
   -------------------------------------------------------------------------  */

void NeighList::grow_nva(int nall, int neall)
{
  if (nall > maxnva) {
    maxnva = element->nmax * (element->maxucell - element->maxgcell) * element->maxapc;
    memory->destroy(nvanumneigh);
    memory->destroy(nvailist);
    memory->destroy(nvaindexlist);
    memory->sfree(nvafirstneigh);
    memory->sfree(nvafirstneighindex);

    memory->create(nvanumneigh, maxnva, "neighlist:nvanumneigh");
    memory->create(nvailist, maxnva, "neighlist:nvailist");
    memory->create(nvaindexlist, maxnva, "neighlist:nvaindexlist");
    nvafirstneigh = (int **) memory->smalloc(maxnva * sizeof(int *), "neighlist:nvafirstneigh");
    nvafirstneighindex = (int **) memory->smalloc(maxnva * sizeof(int *), "neighlist:nvafirstneighindex");
  }
  if (neall > maxnva_elem) {
    maxnva_elem = element->nmax;
    memory->destroy(va2nvalist);
    memory->create(va2nvalist, maxnva_elem, element->maxapc, element->maxucell, "neighlist:va2vnalist");
  }

}

/*  ----------------------------------------------------------------------
   print attributes of this list and associated request
   -------------------------------------------------------------------------  */

void NeighList::print_attributes()
{
  if (comm->me != 0) return;

  NeighRequest *rq = neighbor->requests[index];

  printf("Neighbor list/request %d:\n", index);
  printf("  %p = requestor ptr (instance %d id %d)\n", 
      rq->requestor, rq->requestor_instance, rq->id);
  printf("  %d = pair\n", rq->pair);
  printf("  %d = fix\n", rq->fix);
  printf("  %d = compute\n", rq->compute);
  printf("  %d = command\n", rq->command);
  printf("  %d = neigh\n", rq->neigh);
  printf("\n");
  printf("  %d = half\n", rq->half);
  printf("  %d = full\n", rq->full);
  printf("\n");
  printf("  %d = occasional\n", occasional);
  printf("  %d = newton\n", rq->newton);
  printf("  %d = ghost flag\n", ghost);
  printf("  %d = size\n", rq->size);
  printf("  %d = dnum\n", dnum);
  printf("\n");
  //printf("  %d = skip flag\n", rq->skip);
  //printf("  %d = off2on\n", rq->off2on);
  printf("  %d = copy flag\n", rq->copy);
  printf("\n");
}

/*  ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxatom & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
   -------------------------------------------------------------------------  */

bigint NeighList::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(ilist, maxall);
  bytes += memory->usage(iindexlist, maxall);
  bytes += memory->usage(numneigh, maxall);
  bytes += 2 * maxall * sizeof(int * );
  if (threebody) {
    bytes += memory->usage(nvailist, maxnva);
    bytes += memory->usage(nvaindexlist, maxnva);
    bytes += memory->usage(va2nvalist, maxnva_elem, element->maxapc, element->max_nucell);
    bytes += memory->usage(nvanumneigh, maxnva);
    bytes += 2 * maxnva * sizeof(int * );
  }
  int nmypage = comm->nthreads;

  if (ipage) for (int i = 0; i < nmypage; i++) bytes += ipage[i].size();
  if (indexpage) for (int i = 0; i < nmypage; i++) bytes += indexpage[i].size();
  if (nvaipage) for (int i = 0; i < nmypage; i++) bytes += nvaipage[i].size();
  if (nvaindexpage) for (int i = 0; i < nmypage; i++) bytes += nvaindexpage[i].size();

  //if (dnum && dpage) {
  //  for (int i = 0; i < nmypage; i++) {
  //    bytes += maxatom * sizeof(double *);
  //    bytes += dpage[i].size();
  //  }
  //}

  return bytes;
}
