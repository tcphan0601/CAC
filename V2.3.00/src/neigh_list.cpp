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

enum{NSQ,BIN,MULTI};     // also in Neighbor

/* ---------------------------------------------------------------------- */

NeighList::NeighList(CAC *cac) : Pointers(cac)
{
  // initializations

  maxatom = 0;
  maxelem = 0;
  maxelem_conv = 0;
  maxintg = 0;
  maxintpl = 0;
  maxintpl_conv = 0;

  ainum = einum = iinum = iainum = gnum = 0;
  ailist = eilist = e2ilist = ia2elist = NULL;
  e2ialist = NULL;
  n2ilist = NULL;
  numneigha2a = NULL;
  numneigha2e = NULL;
  numneigha2ia = NULL;
  numneighi2a = NULL;
  numneighe2e = NULL;
  numneighi2ia = NULL;
  numneighia2ia = NULL;
  numneighia2a = NULL;

  firstneigha2a = NULL;
  firstneigha2e = NULL;
  firstneigha2ia = NULL;
  firstneigha2ia_index = NULL;
  firstneighi2a = NULL;
  firstneighia2a = NULL;
  firstneighe2e = NULL;
  firstneighi2ia = NULL;
  firstneighia2ia = NULL;
  firstneighi2ia_index = NULL;
  firstneighia2ia_index = NULL;

  // defaults, but may be reset by post_constructor()
  atomlist = 1;
  elemlist = 1;
  intpllist = 0;
  intglist = 1;
  occasional = 0;
  ghost = 0;
  copy = 0;
  dnum = 0;

  // ptrs

//  iskip = NULL;
//  ijskip = NULL;

  listcopy = NULL;
//  listskip = NULL;
  listfull = NULL;
/*
  listhistory = NULL;
  fix_history = NULL;

  respamiddle = 0;
  listinner = NULL;
  listmiddle = NULL;

  fix_bond = NULL;
*/
  a2apage = NULL;
  a2epage = NULL;
  a2iapage = NULL;
  a2ia_indexpage = NULL;
  i2apage = NULL;
  ia2apage = NULL;
  e2epage = NULL;
  i2iapage = NULL;
  ia2iapage = NULL;
  i2ia_indexpage = NULL;
  ia2ia_indexpage = NULL;
  dpage = NULL;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
//  if (copymode) return;
  if (!copy) {
    memory->destroy(ailist);
    memory->destroy(eilist);
    memory->destroy(n2ilist);
    memory->destroy(e2ilist);
    memory->destroy(e2ialist);
    memory->destroy(ia2elist);
    memory->destroy(numneigha2a);
    memory->destroy(numneigha2e);
    memory->destroy(numneigha2ia);
    memory->destroy(numneighi2a);
    memory->destroy(numneighia2a);
    memory->destroy(numneighe2e);
    memory->destroy(numneighi2ia);
    memory->destroy(numneighia2ia);
    memory->sfree(firstneigha2a);
    memory->sfree(firstneigha2e);
    memory->sfree(firstneigha2ia);
    memory->sfree(firstneigha2ia_index);
    memory->sfree(firstneighi2a);
    memory->sfree(firstneighia2a);
    memory->sfree(firstneighe2e);
    memory->sfree(firstneighi2ia);
    memory->sfree(firstneighia2ia);
    memory->sfree(firstneighi2ia_index);
    memory->sfree(firstneighia2ia_index);

    delete [] a2apage;
    delete [] a2epage;
    delete [] a2iapage;
    delete [] a2ia_indexpage;
    delete [] i2apage;
    delete [] ia2apage;
    delete [] e2epage;
    delete [] i2iapage;
    delete [] ia2iapage;
    delete [] i2ia_indexpage;
    delete [] ia2ia_indexpage;
  }
}

/* ----------------------------------------------------------------------
   adjust settings to match corresponding NeighRequest
   cannot do this in constructor b/c not all NeighLists are allocated yet
   copy -> set listcopy for list to copy from
   skip -> set listskip for list to skip from, create copy of itype,ijtype
   halffull -> set listfull for full list to derive from
   history -> set LH and FH ptrs in partner list that uses the history info
   respaouter -> set listinner/listmiddle for other rRESPA lists
   bond -> set fix_bond to Fix that made the request
------------------------------------------------------------------------- */

void NeighList::post_constructor(NeighRequest *rq)
{
  // copy request settings used by list itself
  
  occasional = rq->occasional;
  ghost = rq->ghost;
  atomlist = rq->atomlist;
  elemlist = rq->elemlist;
  intglist = rq->intglist;
  intpllist = rq->intpllist;
  copy = rq->copy;
  dnum = rq->dnum;

  if (rq->copy)
    listcopy = neighbor->lists[rq->copylist];
/*
  if (rq->skip) {
    listskip = neighbor->lists[rq->skiplist];
    int ntypes = atom->ntypes;
    iskip = new int[ntypes+1];
    memory->create(ijskip,ntypes+1,ntypes+1,"neigh_list:ijskip");
    int i,j;
    for (i = 1; i <= ntypes; i++) iskip[i] = rq->iskip[i];
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = rq->ijskip[i][j];
  }
*/
}

/* ---------------------------------------------------------------------- */

void NeighList::setup_pages(int pgsize_caller, int oneatom_caller)
{
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;

  int nmypage = comm->nthreads;
  a2apage = new MyPage<int>[nmypage];
  a2epage = new MyPage<int>[nmypage];
  a2iapage = new MyPage<int>[nmypage];
  a2ia_indexpage = new MyPage<int>[nmypage];
  i2apage = new MyPage<int>[nmypage];
  ia2apage = new MyPage<int>[nmypage];
  e2epage = new MyPage<int>[nmypage];
  i2iapage = new MyPage<int>[nmypage];
  ia2iapage = new MyPage<int>[nmypage];
  i2ia_indexpage = new MyPage<int>[nmypage];
  ia2ia_indexpage = new MyPage<int>[nmypage];
  for (int i = 0; i < nmypage; i++){
    a2apage[i].init(oneatom,pgsize,PGDELTA);
    a2epage[i].init(oneatom,pgsize,PGDELTA);
    a2iapage[i].init(oneatom,pgsize,PGDELTA);
    a2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    i2apage[i].init(oneatom,pgsize,PGDELTA);
    i2iapage[i].init(oneatom,pgsize,PGDELTA);
    i2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    e2epage[i].init(oneatom,pgsize,PGDELTA);
    ia2apage[i].init(oneatom,pgsize,PGDELTA);
    ia2iapage[i].init(oneatom,pgsize,PGDELTA);
    ia2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
  }

  if (dnum) {
    dpage = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++)
      dpage[i].init(dnum*oneatom,dnum*pgsize,PGDELTA);
  } else dpage = NULL;
}

/* ----------------------------------------------------------------------
   grow per-atom data to allow for nlocal/nall atoms
   triggered by neighbor list build
   not called if a copy list
------------------------------------------------------------------------- */

void NeighList::grow_atom(int nlocal, int nall)
{
  // skip if data structs are already big enough

  if (ghost) {
    if (nall <= maxatom) return;
  } else {
    if (nlocal <= maxatom) return;
  }

  maxatom = atom->nmax;

  memory->destroy(ailist);
  memory->destroy(numneigha2a);
  memory->destroy(numneigha2e);
  memory->destroy(numneigha2ia);

  memory->sfree(firstneigha2a);
  memory->sfree(firstneigha2e);
  memory->sfree(firstneigha2ia);
  memory->sfree(firstneigha2ia_index);

  memory->create(ailist,maxatom,"neighlist:ailist");
  memory->create(numneigha2a,maxatom,"neighlist:numneigha2a");
  memory->create(numneigha2e,maxatom,"neighlist:numneigha2e");
  memory->create(numneigha2ia,maxatom,"neighlist:numneigha2ia");

  firstneigha2a = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2a");
  firstneigha2e = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2e");
  firstneigha2ia = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia");
  firstneigha2ia_index = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia_index");
/*
   if (dnum) {
    memory->sfree(firstdouble);
    firstdouble = (double **) memory->smalloc(maxatom*sizeof(double *),
                                              "neighlist:firstdouble");
  }
*/
}

/* ----------------------------------------------------------------------
   grow per-element data to allow for nlocal/nall elements
   triggered by neighbor list build
   not called if a copy list
------------------------------------------------------------------------- */

void NeighList::grow_elem(int nlocal, int nall)
{
  // ghost elements conversion lists might be required 
  // even if pairs for ghost elements are not stored
  int npe = element->npe; 
  if (maxelem_conv < nall) {
    maxelem_conv = nall;
    memory->destroy(n2ilist);
    memory->destroy(e2ilist);
    memory->destroy(e2ialist);
    if (intpllist)
      memory->create(e2ialist,maxelem_conv,"neighlist:e2ialist");
    if (intglist) {
      memory->create(n2ilist,maxelem_conv,npe,"neighlist:n2ilist");
      memory->create(e2ilist,maxelem_conv,"neighlist:e2ilist");
    }
  }

  // skip if data structs are already big enough
  if (ghost) {
    if (nall <= maxelem) return;
  } else {
    if (nlocal <= maxelem) return;
  }

  maxelem = element->nmax;
  memory->destroy(eilist);
  memory->destroy(numneighe2e);
  memory->sfree(firstneighe2e);

  memory->create(eilist,maxelem,"neighlist:eilist");
  memory->create(numneighe2e,maxelem,"neighlist:numneighe2e");
  firstneighe2e = (int **) memory->smalloc(maxelem*sizeof(int *),"neighlist:firstneighe2e"); 

}

/* ----------------------------------------------------------------------
   grow per-intg data to allow for nlocal/nall integration points
   triggered by neighbor list build
   not called if a copy list
   ------------------------------------------------------------------------- */

void NeighList::grow_intg(int nlocal, int nall)
{
  // skip if data structs are already big enough
  
  if (ghost) {
    if (nall <= maxintg) return;
  } else {
    if (nlocal <= maxintg) return;
  }

  maxintg = element->nmaxintg;

  memory->destroy(numneighi2a);
  memory->destroy(numneighi2ia);

  memory->sfree(firstneighi2a);
  memory->sfree(firstneighi2ia);
  memory->sfree(firstneighi2ia_index);

  memory->create(numneighi2a,maxintg,"neighlist:numneighi2a");
  memory->create(numneighi2ia,maxintg,"neighlist:numneighi2ia");

  firstneighi2a = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2a");
  firstneighi2ia = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2ia");
  firstneighi2ia_index = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2ia_index");
}

/* ----------------------------------------------------------------------
   grow per-intpl data to allow for nlocal/nall interpolated atoms
   triggered by neighbor list build
   not called if a copy list
   ------------------------------------------------------------------------- */

void NeighList::grow_intpl(bigint nlocal, bigint nall)
{
  // ghost elements conversion lists might be required 
  // even if pairs for ghost elements are not stored
  if (maxintpl_conv < nall) {
    maxintpl_conv = element->nmaxintpl;
    memory->destroy(ia2elist); 
    memory->create(ia2elist,maxintpl_conv,"neighlist:ia2elist");
  } 

  // skip if data structs are already big enough
  if (ghost) {
    if (nall <= maxintpl) return;
  } else {
    if (nlocal <= maxintpl) return;
  }

  maxintpl = element->nmaxintpl;
  memory->destroy(numneighia2a);
  memory->destroy(numneighia2ia);
  memory->sfree(firstneighia2a);
  memory->sfree(firstneighia2ia);
  memory->sfree(firstneighia2ia_index);
  memory->create(numneighia2a,maxintpl,"neighlist:numneighia2a");
  memory->create(numneighia2ia,maxintpl,"neighlist:numneighia2ia");

  firstneighia2a = (int **) memory->smalloc(maxintpl*sizeof(int *),"neighlist:firstneighia2a");
  firstneighia2ia = (int **) memory->smalloc(maxintpl*sizeof(int *),"neighlist:firstneighia2ia");
  firstneighia2ia_index = (int **) memory->smalloc(maxintpl*sizeof(int *),"neighlist:firstneighia2ia_index");
}


/* ----------------------------------------------------------------------
   print attributes of this list and associated request
   ------------------------------------------------------------------------- */

void NeighList::print_attributes()
{
  if (comm->me != 0) return;

  NeighRequest *rq = neighbor->requests[index];

  printf("Neighbor list/request %d:\n",index);
  printf("  %p = requestor ptr (instance %d id %d)\n",
      rq->requestor,rq->requestor_instance,rq->id);
  printf("  %d = pair\n",rq->pair);
  printf("  %d = fix\n",rq->fix);
  printf("  %d = compute\n",rq->compute);
  printf("  %d = command\n",rq->command);
  printf("  %d = neigh\n",rq->neigh);
  printf("\n");
  printf("  %d = half\n",rq->half);
  printf("  %d = full\n",rq->full);
  printf("\n");
  printf("  %d = occasional\n",occasional);
  printf("  %d = newton\n",rq->newton);
  printf("  %d = ghost flag\n",ghost);
  printf("  %d = size\n",rq->size);
  printf("  %d = dnum\n",dnum);
  printf("\n");
//  printf("  %d = skip flag\n",rq->skip);
//  printf("  %d = off2on\n",rq->off2on);
  printf("  %d = copy flag\n",rq->copy);
  printf("\n");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxatom & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
   ------------------------------------------------------------------------- */

bigint NeighList::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(ailist,maxatom);
  bytes += memory->usage(numneigha2a,maxatom);
  bytes += memory->usage(numneigha2e,maxatom);
  bytes += memory->usage(numneigha2ia,maxatom);
  bytes += 4*maxatom*sizeof(int*);
  int npe = element->npe;
  bytes += memory->usage(eilist,maxelem);
  if (intglist) {
    bytes += memory->usage(n2ilist,maxelem_conv,npe);
    bytes += memory->usage(e2ilist,maxelem_conv);
  }
  if (intpllist) { 
    bytes += memory->usage(e2ialist,maxelem_conv);
  }
  bytes += memory->usage(numneighe2e,maxelem);
  bytes += maxelem*sizeof(int*);

  bytes += memory->usage(numneighi2a,maxintg);
  bytes += memory->usage(numneighi2ia,maxintg);
  bytes += 3*maxintg*sizeof(int*);

  bytes += memory->usage(ia2elist,maxintpl_conv);
  bytes += memory->usage(numneighia2a,maxintpl);
  bytes += memory->usage(numneighia2ia,maxintpl);
  bytes += 3*maxintpl*sizeof(int*);

  int nmypage = comm->nthreads;

  if (a2apage) for (int i = 0; i < nmypage; i++) bytes += a2apage[i].size();
  if (a2epage) for (int i = 0; i < nmypage; i++) bytes += a2epage[i].size();
  if (a2iapage) for (int i = 0; i < nmypage; i++) bytes += a2iapage[i].size();
  if (a2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += a2ia_indexpage[i].size();
  if (i2apage) for (int i = 0; i < nmypage; i++) bytes += i2apage[i].size();
  if (ia2apage) for (int i = 0; i < nmypage; i++) bytes += ia2apage[i].size();
  if (e2epage) for (int i = 0; i < nmypage; i++) bytes += e2epage[i].size();
  if (i2iapage) for (int i = 0; i < nmypage; i++) bytes += i2iapage[i].size();
  if (ia2iapage) for (int i = 0; i < nmypage; i++) bytes += ia2iapage[i].size();
  if (i2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += i2ia_indexpage[i].size();
  if (ia2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += ia2ia_indexpage[i].size();
  /*
     if (dnum && dpage) {
     for (int i = 0; i < nmypage; i++) {
     bytes += maxatom * sizeof(double *);
     bytes += dpage[i].size();
     }
     }
     */
  return bytes;
}
