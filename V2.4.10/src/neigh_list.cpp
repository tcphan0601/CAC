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

/* ---------------------------------------------------------------------- */

NeighList::NeighList(CAC *cac) : Pointers(cac)
{
  // initializations

  maxatom = 0;
  maxelem = 0;
  maxelem_conv = 0;
  maxintg = 0;
  maxnode = 0;
  maxintpl = 0;
  maxintpl_conv = 0;
  maxnia = 0; 

  ainum = einum = iinum = iainum = niainum = agnum = 0;
  ailist = eilist = e2ilist = ia2elist = NULL;
  e2ialist = NULL;
  n2ilist = ia2nialist = NULL;
  nia2ialist = nia2ia_indexlist = NULL;
  numneigha = NULL;
  numneighia = NULL;
  numneigha2a = NULL;
  numneigha2ia = NULL;
  numneighi2a = NULL;
  numneighn2a = NULL;
  numneighe2a = NULL;
  numneighe2e = NULL;
  numneighi2ia = NULL;
  numneighn2ia = NULL;
  numneighia2ia = NULL;
  numneighia2a = NULL;
  numneighnia2a = NULL;
  numneighnia2ia = NULL;

  firstneigha = NULL;
  firstneigha_index = NULL;
  firstneighia = NULL;
  firstneighia_index = NULL;
  firstneigha2a = NULL;
  firstneigha2ia = NULL;
  firstneigha2ia_index = NULL;
  firstneighi2a = NULL;
  firstneighn2a = NULL;
  firstneighia2a = NULL;
  firstneighe2a = NULL;
  firstneighe2e = NULL;
  firstneighi2ia = NULL;
  firstneighn2ia = NULL;
  firstneighia2ia = NULL;
  firstneighi2ia_index = NULL;
  firstneighn2ia_index = NULL;
  firstneighia2ia_index = NULL;
  firstneighnia2a = NULL;
  firstneighnia2ia = NULL;
  firstneighnia2ia_index = NULL;

  numneigha2a_outer = NULL;
  numneigha2ia_outer = NULL;
  numneighi2a_outer = NULL;
  numneighi2ia_outer = NULL;
  firstneigha2a_outer = NULL;
  firstneigha2ia_outer = NULL;
  firstneigha2ia_index_outer = NULL;
  firstneighi2a_outer = NULL;
  firstneighi2ia_outer = NULL;
  firstneighi2ia_index_outer = NULL;

  // defaults, but may be reset by post_constructor()

  atomlist = 1;
  atomonlylist = 0;
  elemlist = 1;
  nodelist = 0;
  intglist = 1;
  intpllist = 0;
  doublelist = 0;

  occasional = 0;
  ghost = 0;
  copy = 0;
  sort = 0;
  sortnum = 0;
  dnum = 0;
  threebody = 0;

  group = 0;

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

  apage = NULL;
  a_indexpage = NULL;
  iapage = NULL;
  ia_indexpage = NULL;
  a2apage = NULL;
  a2iapage = NULL;
  a2ia_indexpage = NULL;
  i2apage = NULL;
  n2apage = NULL;
  ia2apage = NULL;
  e2apage = NULL;
  e2epage = NULL;
  i2iapage = NULL;
  n2iapage = NULL;
  ia2iapage = NULL;
  i2ia_indexpage = NULL;
  n2ia_indexpage = NULL;
  ia2ia_indexpage = NULL;
  nia2apage = NULL;
  nia2iapage = NULL;
  nia2ia_indexpage = NULL;

  a2a_outer_page = NULL;
  a2ia_outer_page = NULL;
  a2ia_index_outer_page = NULL;
  i2a_outer_page = NULL;
  i2ia_outer_page = NULL;
  i2ia_index_outer_page = NULL;

  dpage = NULL;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{

  //if (copymode) return;
  if (copy) return;
  memory->destroy(ailist);
  memory->destroy(eilist);
  memory->destroy(n2ilist);
  memory->destroy(e2ilist);
  memory->destroy(e2ialist);
  memory->destroy(ia2elist);
  memory->destroy(ia2nialist);
  memory->destroy(nia2ialist);
  memory->destroy(nia2ia_indexlist);

  memory->destroy(numneigha2a);
  memory->destroy(numneigha2ia);
  memory->sfree(firstneigha2a);
  memory->sfree(firstneigha2ia);
  memory->sfree(firstneigha2ia_index);

  memory->destroy(numneigha2a_outer);
  memory->destroy(numneigha2ia_outer);
  memory->sfree(firstneigha2a_outer);
  memory->sfree(firstneigha2ia_outer);
  memory->sfree(firstneigha2ia_index_outer);

  memory->destroy(numneighia2a);
  memory->destroy(numneighia2ia);
  memory->sfree(firstneighia2a);
  memory->sfree(firstneighia2ia);
  memory->sfree(firstneighia2ia_index);

  memory->destroy(numneighe2a);
  memory->destroy(numneighe2e);
  memory->sfree(firstneighe2a);
  memory->sfree(firstneighe2e);

  memory->destroy(numneighi2a);
  memory->destroy(numneighi2ia);
  memory->sfree(firstneighi2a);
  memory->sfree(firstneighi2ia);
  memory->sfree(firstneighi2ia_index);

  memory->destroy(numneighi2a_outer);
  memory->destroy(numneighi2ia_outer);
  memory->sfree(firstneighi2a_outer);
  memory->sfree(firstneighi2ia_outer);
  memory->sfree(firstneighi2ia_index_outer);

  memory->destroy(numneighn2a);
  memory->destroy(numneighn2ia);
  memory->sfree(firstneighn2a);
  memory->sfree(firstneighn2ia);
  memory->sfree(firstneighn2ia_index);

  memory->destroy(numneighnia2a);
  memory->destroy(numneighnia2ia);
  memory->sfree(firstneighnia2a);
  memory->sfree(firstneighnia2ia);
  memory->sfree(firstneighnia2ia_index);

  memory->destroy(numneigha);
  memory->sfree(firstneigha);
  memory->sfree(firstneigha_index);

  memory->destroy(numneighia);
  memory->sfree(firstneighia);
  memory->sfree(firstneighia_index);

  delete [] apage;
  delete [] a_indexpage;
  delete [] iapage;
  delete [] ia_indexpage;
  delete [] a2apage;
  delete [] a2iapage;
  delete [] a2ia_indexpage;
  delete [] ia2apage;
  delete [] ia2iapage;
  delete [] ia2ia_indexpage;
  delete [] e2apage;
  delete [] e2epage;
  delete [] i2apage;
  delete [] i2iapage;
  delete [] i2ia_indexpage;
  delete [] n2apage;
  delete [] n2iapage;
  delete [] n2ia_indexpage;
  delete [] nia2apage;
  delete [] nia2iapage;
  delete [] nia2ia_indexpage;

  delete [] a2a_outer_page;
  delete [] a2ia_outer_page;
  delete [] a2ia_index_outer_page;
  delete [] i2a_outer_page;
  delete [] i2ia_outer_page;
  delete [] i2ia_index_outer_page;

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
  atomonlylist = rq->atomonlylist;
  elemlist = rq->elemlist;
  nodelist = rq->nodelist;
  intglist = rq->intglist;
  intpllist = rq->intpllist;
  doublelist = rq->doublelist;

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
  //  memory->create(ijskip,ntypes+1,ntypes+1,"neigh_list:ijskip");
  //  int i,j;
  //  for (i = 1; i <= ntypes; i++) iskip[i] = rq->iskip[i];
  //  for (i = 1; i <= ntypes; i++)
  //    for (j = 1; j <= ntypes; j++)
  //      ijskip[i][j] = rq->ijskip[i][j];
  //}

}

/* ---------------------------------------------------------------------- */

void NeighList::setup_pages(int pgsize_caller, int oneatom_caller)
{
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;

  int nmypage = comm->nthreads;
  apage = new MyPage<int>[nmypage];
  a_indexpage = new MyPage<int>[nmypage];
  iapage = new MyPage<int>[nmypage];
  ia_indexpage = new MyPage<int>[nmypage];
  a2apage = new MyPage<int>[nmypage];
  a2iapage = new MyPage<int>[nmypage];
  a2ia_indexpage = new MyPage<int>[nmypage];
  ia2apage = new MyPage<int>[nmypage];
  ia2iapage = new MyPage<int>[nmypage];
  ia2ia_indexpage = new MyPage<int>[nmypage];
  e2apage = new MyPage<int>[nmypage];
  e2epage = new MyPage<int>[nmypage];
  i2apage = new MyPage<int>[nmypage];
  i2iapage = new MyPage<int>[nmypage];
  i2ia_indexpage = new MyPage<int>[nmypage];
  n2apage = new MyPage<int>[nmypage];
  n2iapage = new MyPage<int>[nmypage];
  n2ia_indexpage = new MyPage<int>[nmypage];
  nia2apage = new MyPage<int>[nmypage];
  nia2iapage = new MyPage<int>[nmypage];
  nia2ia_indexpage = new MyPage<int>[nmypage];

  a2a_outer_page = new MyPage<int>[nmypage];
  a2ia_outer_page = new MyPage<int>[nmypage];
  a2ia_index_outer_page = new MyPage<int>[nmypage];
  i2a_outer_page = new MyPage<int>[nmypage];
  i2ia_outer_page = new MyPage<int>[nmypage];
  i2ia_index_outer_page = new MyPage<int>[nmypage];

  for (int i = 0; i < nmypage; i++) {
    apage[i].init(oneatom,pgsize,PGDELTA);
    a_indexpage[i].init(oneatom,pgsize,PGDELTA);
    iapage[i].init(oneatom,pgsize,PGDELTA);
    ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    a2apage[i].init(oneatom,pgsize,PGDELTA);
    a2iapage[i].init(oneatom,pgsize,PGDELTA);
    a2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    ia2apage[i].init(oneatom,pgsize,PGDELTA);
    ia2iapage[i].init(oneatom,pgsize,PGDELTA);
    ia2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    i2apage[i].init(oneatom,pgsize,PGDELTA);
    i2iapage[i].init(oneatom,pgsize,PGDELTA);
    i2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    n2apage[i].init(oneatom,pgsize,PGDELTA);
    n2iapage[i].init(oneatom,pgsize,PGDELTA);
    n2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);
    e2apage[i].init(oneatom,pgsize,PGDELTA);
    e2epage[i].init(oneatom,pgsize,PGDELTA);
    nia2apage[i].init(oneatom,pgsize,PGDELTA);
    nia2iapage[i].init(oneatom,pgsize,PGDELTA);
    nia2ia_indexpage[i].init(oneatom,pgsize,PGDELTA);

    a2a_outer_page[i].init(oneatom,pgsize,PGDELTA);
    a2ia_outer_page[i].init(oneatom,pgsize,PGDELTA);
    a2ia_index_outer_page[i].init(oneatom,pgsize,PGDELTA);
    i2a_outer_page[i].init(oneatom,pgsize,PGDELTA);
    i2ia_outer_page[i].init(oneatom,pgsize,PGDELTA);
    i2ia_index_outer_page[i].init(oneatom,pgsize,PGDELTA);
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
  memory->create(ailist,maxatom,"neighlist:ailist");

  if (sort) {
    memory->destroy(numneigha);
    memory->create(numneigha,maxatom,"neighlist:numneigha2a");
    memory->sfree(firstneigha);
    memory->sfree(firstneigha_index);
    firstneigha = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha");
    firstneigha_index = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha_index");
  } else {
    memory->destroy(numneigha2a);
    memory->destroy(numneigha2ia);
    memory->sfree(firstneigha2a);
    memory->sfree(firstneigha2ia);
    memory->sfree(firstneigha2ia_index);

    memory->create(numneigha2a,maxatom,"neighlist:numneigha2a");
    memory->create(numneigha2ia,maxatom,"neighlist:numneigha2ia");
    firstneigha2a = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2a");
    firstneigha2ia = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia");
    firstneigha2ia_index = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia_index");

    if (doublelist) {
      memory->destroy(numneigha2a_outer);
      memory->destroy(numneigha2ia_outer);
      memory->sfree(firstneigha2a_outer);
      memory->sfree(firstneigha2ia_outer);
      memory->sfree(firstneigha2ia_index_outer);

      memory->create(numneigha2a_outer,maxatom,"neighlist:numneigha2a_outer");
      memory->create(numneigha2ia_outer,maxatom,"neighlist:numneigha2ia_outer");
      firstneigha2a_outer = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2a_outer");
      firstneigha2ia_outer = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia_outer");
      firstneigha2ia_index_outer = (int **) memory->smalloc(maxatom*sizeof(int *),"neighlist:firstneigha2ia_index_outer");
    }
  }
  //if (dnum) {
  //  memory->sfree(firstdouble);
  //  firstdouble = (double **) memory->smalloc(maxatom*sizeof(double *),
  //      "neighlist:firstdouble");
  //}

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

  if (maxelem_conv < nall) {
    maxelem_conv = nall;
    if (intpllist) {
      memory->destroy(e2ialist);
      memory->create(e2ialist,maxelem_conv,"neighlist:e2ialist");
    }
    if (intglist) {
      memory->destroy(n2ilist);
      memory->destroy(e2ilist);
      memory->create(n2ilist,maxelem_conv,element->max_npe,"neighlist:n2ilist");
      memory->create(e2ilist,maxelem_conv,"neighlist:e2ilist");
    }
    if (nodelist) 
      memory->create(e2nlist,maxelem_conv,"neighlist:e2nlist");
  }

  // skip if data structs are already big enough

  if (ghost || threebody) {
    if (nall <= maxelem) return;
  } else {
    if (nlocal <= maxelem) return;
  }
  maxelem = element->nmax;

  // only owned interpolated atom can be nia

  if (threebody) {
    memory->destroy(ia2nialist);
    memory->create(ia2nialist,maxelem,element->max_nintpl,"neighlist:ia2nialist");
  }

  memory->destroy(eilist);
  memory->destroy(numneighe2a);
  memory->destroy(numneighe2e);
  memory->sfree(firstneighe2a);
  memory->sfree(firstneighe2e);

  memory->create(eilist,maxelem,"neighlist:eilist");
  memory->create(numneighe2a,maxelem,"neighlist:numneighe2a");
  memory->create(numneighe2e,maxelem,"neighlist:numneighe2e");
  firstneighe2a = (int **) memory->smalloc(maxelem*sizeof(int *),"neighlist:firstneighe2a"); 
  firstneighe2e = (int **) memory->smalloc(maxelem*sizeof(int *),"neighlist:firstneighe2e"); 

}

/* ----------------------------------------------------------------------
   grow per-intg data to allow for nlocal/nall nodes
   triggered by neighbor list build
   not called if a copy list
   ------------------------------------------------------------------------- */

void NeighList::grow_node(int nlocal, int nall)
{
  // skip if data structs are already big enough

  if (ghost) {
    if (nall <= maxnode) return;
  } else {
    if (nlocal <= maxnode) return;
  }

  maxnode = element->nmaxnode;

  memory->destroy(numneighn2a);
  memory->destroy(numneighn2ia);

  memory->sfree(firstneighn2a);
  memory->sfree(firstneighn2ia);
  memory->sfree(firstneighn2ia_index);

  memory->create(numneighn2a,maxnode,"neighlist:numneighn2a");
  memory->create(numneighn2ia,maxnode,"neighlist:numneighn2ia");

  firstneighn2a = (int **) memory->smalloc(maxnode*sizeof(int *),"neighlist:firstneighn2a");
  firstneighn2ia = (int **) memory->smalloc(maxnode*sizeof(int *),"neighlist:firstneighn2ia");
  firstneighn2ia_index = (int **) memory->smalloc(maxnode*sizeof(int *),"neighlist:firstneighn2ia_index");
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

  if (doublelist) {
    memory->destroy(numneighi2a_outer);
    memory->destroy(numneighi2ia_outer);
    memory->sfree(firstneighi2a_outer);
    memory->sfree(firstneighi2ia_outer);
    memory->sfree(firstneighi2ia_index_outer);

    memory->create(numneighi2a_outer,maxintg,"neighlist:numneighi2a_outer");
    memory->create(numneighi2ia_outer,maxintg,"neighlist:numneighi2ia_outer");
    firstneighi2a_outer = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2a_outer");
    firstneighi2ia_outer = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2ia_outer");
    firstneighi2ia_index_outer = (int **) memory->smalloc(maxintg*sizeof(int *),"neighlist:firstneighi2ia_index_outer");
  }
}

/* ----------------------------------------------------------------------
   grow per-nia data to allow for n nia
   triggered by neighbor list build
   not called if a copy list
   ------------------------------------------------------------------------- */

void NeighList::grow_nia(bigint n)
{
  if (n <= maxnia) return;

  maxnia = element->nmaxintpl;
  memory->destroy(numneighnia2a);
  memory->destroy(numneighnia2ia);
  memory->destroy(nia2ialist);
  memory->destroy(nia2ia_indexlist);
  memory->sfree(firstneighnia2a);
  memory->sfree(firstneighnia2ia);
  memory->sfree(firstneighnia2ia_index);

  memory->create(numneighnia2a,maxnia,"neighlist:numneighnia2a");
  memory->create(numneighnia2ia,maxnia,"neighlist:numneighnia2ia");
  memory->create(nia2ialist,maxnia,"neighlist:nia2ialist");
  memory->create(nia2ia_indexlist,maxnia,"neighlist:nia2ia_indexlist");

  firstneighnia2a = (int **) memory->smalloc(maxnia*sizeof(int *),"neighlist:firstneighnia2a");
  firstneighnia2ia = (int **) memory->smalloc(maxnia*sizeof(int *),"neighlist:firstneighnia2ia");
  firstneighnia2ia_index = (int **) memory->smalloc(maxnia*sizeof(int *),"neighlist:firstneighnia2ia_index");

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

  if (sort) {
    memory->destroy(numneighia);
    memory->create(numneighia,maxintpl,"neighlist:numneighia");
    memory->sfree(firstneighia);
    memory->sfree(firstneighia_index);
    firstneighia = (int **) memory->smalloc(maxintpl*sizeof(int *),"neighlist:firstneighia");
    firstneighia_index = (int **) memory->smalloc(maxintpl*sizeof(int *),"neighlist:firstneighia_index");
  } else {

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
  //printf("  %d = skip flag\n",rq->skip);
  //printf("  %d = off2on\n",rq->off2on);
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
  if (sort) {
    bytes += memory->usage(numneigha,maxatom);
    bytes += memory->usage(numneighia,maxintpl);
    bytes += 2*maxatom*sizeof(int*);
    bytes += 2*maxintpl*sizeof(int*);
  } else {
    bytes += memory->usage(numneigha2a,maxatom);
    bytes += memory->usage(numneigha2ia,maxatom);
    bytes += 3*maxatom*sizeof(int*);
    bytes += memory->usage(numneighia2a,maxintpl);
    bytes += memory->usage(numneighia2ia,maxintpl);
    bytes += 3*maxintpl*sizeof(int*);
    if (doublelist) {
      bytes += memory->usage(numneigha2a_outer,maxatom);
      bytes += memory->usage(numneigha2ia_outer,maxatom);
      bytes += 3*maxatom*sizeof(int*);
    }
  }
  bytes += memory->usage(eilist,maxelem);
  if (intglist) {
    bytes += memory->usage(n2ilist,maxelem_conv,element->max_npe);
    bytes += memory->usage(e2ilist,maxelem_conv);
    bytes += memory->usage(numneighi2a,maxintg);
    bytes += memory->usage(numneighi2ia,maxintg);
    bytes += 3*maxintg*sizeof(int*);
    if (doublelist) {
      bytes += memory->usage(numneighi2a_outer,maxatom);
      bytes += memory->usage(numneighi2ia_outer,maxatom);
      bytes += 3*maxatom*sizeof(int*);
    }
  }
  if (nodelist) {
    bytes += memory->usage(e2nlist,maxelem_conv);
    bytes += memory->usage(numneighn2a,maxnode);
    bytes += memory->usage(numneighn2ia,maxnode);
    bytes += 3*maxintg*sizeof(int*);
  }

  if (intpllist) 
    bytes += memory->usage(e2ialist,maxelem_conv);

  if (threebody) {
    bytes += memory->usage(ia2nialist,maxelem,element->max_nintpl);
    bytes += memory->usage(numneighe2e,maxelem);
    bytes += memory->usage(numneighe2a,maxelem);
    bytes += 2*maxelem*sizeof(int*);
    bytes += memory->usage(ia2elist,maxintpl_conv);
    bytes += memory->usage(numneighnia2a,maxnia);
    bytes += memory->usage(numneighnia2ia,maxnia);
    bytes += 3*maxnia*sizeof(int*);
  }
  int nmypage = comm->nthreads;

  if (apage) for (int i = 0; i < nmypage; i++) bytes += apage[i].size();
  if (a_indexpage) for (int i = 0; i < nmypage; i++) bytes += a_indexpage[i].size();
  if (iapage) for (int i = 0; i < nmypage; i++) bytes += iapage[i].size();
  if (ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += ia_indexpage[i].size();
  if (a2apage) for (int i = 0; i < nmypage; i++) bytes += a2apage[i].size();
  if (a2iapage) for (int i = 0; i < nmypage; i++) bytes += a2iapage[i].size();
  if (a2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += a2ia_indexpage[i].size();
  if (i2apage) for (int i = 0; i < nmypage; i++) bytes += i2apage[i].size();
  if (n2apage) for (int i = 0; i < nmypage; i++) bytes += n2apage[i].size();
  if (ia2apage) for (int i = 0; i < nmypage; i++) bytes += ia2apage[i].size();
  if (e2apage) for (int i = 0; i < nmypage; i++) bytes += e2apage[i].size();
  if (e2epage) for (int i = 0; i < nmypage; i++) bytes += e2epage[i].size();
  if (i2iapage) for (int i = 0; i < nmypage; i++) bytes += i2iapage[i].size();
  if (n2iapage) for (int i = 0; i < nmypage; i++) bytes += n2iapage[i].size();
  if (ia2iapage) for (int i = 0; i < nmypage; i++) bytes += ia2iapage[i].size();
  if (i2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += i2ia_indexpage[i].size();
  if (n2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += n2ia_indexpage[i].size();
  if (ia2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += ia2ia_indexpage[i].size();
  if (nia2apage) for (int i = 0; i < nmypage; i++) bytes += nia2apage[i].size();
  if (nia2iapage) for (int i = 0; i < nmypage; i++) bytes += nia2iapage[i].size();
  if (nia2ia_indexpage) for (int i = 0; i < nmypage; i++) bytes += nia2ia_indexpage[i].size();

  if (a2a_outer_page) for (int i = 0; i < nmypage; i++) bytes += a2a_outer_page[i].size();
  if (a2ia_outer_page) for (int i = 0; i < nmypage; i++) bytes += a2ia_outer_page[i].size();
  if (a2ia_index_outer_page) for (int i = 0; i < nmypage; i++) bytes += a2ia_index_outer_page[i].size();
  if (i2a_outer_page) for (int i = 0; i < nmypage; i++) bytes += i2a_outer_page[i].size();
  if (i2ia_outer_page) for (int i = 0; i < nmypage; i++) bytes += i2ia_outer_page[i].size();
  if (i2ia_index_outer_page) for (int i = 0; i < nmypage; i++) bytes += i2ia_index_outer_page[i].size();

  //if (dnum && dpage) {
  //  for (int i = 0; i < nmypage; i++) {
  //    bytes += maxatom * sizeof(double *);
  //    bytes += dpage[i].size();
  //  }
  //}

  return bytes;
}
