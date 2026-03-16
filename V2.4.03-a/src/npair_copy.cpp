#include "npair_copy.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "my_page.h"


using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

NPairCopy::NPairCopy(CAC *cac) : NPair(cac) {}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
   update this list if adding new neighbor lists
------------------------------------------------------------------------- */

void NPairCopy::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->ainum = listcopy->ainum;
  list->einum = listcopy->einum;
  list->iinum = listcopy->iinum;
  list->iainum = listcopy->iainum;
  list->niainum = listcopy->niainum;
  list->gnum = listcopy->gnum;
  list->ailist = listcopy->ailist;
  list->eilist = listcopy->eilist;
  list->n2ilist = listcopy->n2ilist;
  list->e2ilist = listcopy->e2ilist;
  list->e2ialist = listcopy->e2ialist;
  list->ia2elist = listcopy->ia2elist;
  list->ia2nialist = listcopy->ia2nialist;
  list->nia2ialist = listcopy->nia2ialist;
  list->nia2ia_indexlist = listcopy->nia2ia_indexlist;
  
  list->numneigha2a = listcopy->numneigha2a;
  list->numneigha2ia = listcopy->numneigha2ia;
  list->numneighi2a = listcopy->numneighi2a;
  list->numneighia2a = listcopy->numneighia2a;
  list->numneighe2e = listcopy->numneighe2e;
  list->numneighi2ia = listcopy->numneighi2ia;
  list->numneighia2ia = listcopy->numneighia2ia;
  list->numneighnia2a = listcopy->numneighnia2a;
  list->numneighnia2ia = listcopy->numneighnia2ia;

  list->firstneigha2a = listcopy->firstneigha2a;
  list->firstneigha2ia = listcopy->firstneigha2ia;
  list->firstneigha2ia_index = listcopy->firstneigha2ia_index;
  list->firstneighi2a = listcopy->firstneighi2a;
  list->firstneighia2a = listcopy->firstneighia2a;
  list->firstneighe2e = listcopy->firstneighe2e;
  list->firstneighi2ia = listcopy->firstneighi2ia;
  list->firstneighi2ia_index = listcopy->firstneighi2ia_index;
  list->firstneighia2ia = listcopy->firstneighia2ia;
  list->firstneighia2ia_index = listcopy->firstneighia2ia_index;
  list->firstneighnia2a = listcopy->firstneighnia2a;
  list->firstneighnia2ia = listcopy->firstneighnia2ia;
  list->firstneighnia2ia_index = listcopy->firstneighnia2ia_index;

  list->a2apage = listcopy->a2apage;
  list->a2iapage = listcopy->a2iapage;
  list->a2ia_indexpage = listcopy->a2ia_indexpage;
  list->i2apage = listcopy->i2apage;
  list->ia2apage = listcopy->ia2apage;
  list->e2epage = listcopy->e2epage;
  list->i2iapage = listcopy->i2iapage;
  list->i2ia_indexpage = listcopy->i2ia_indexpage;
  list->ia2iapage = listcopy->ia2iapage;
  list->ia2ia_indexpage = listcopy->ia2ia_indexpage;
  list->nia2apage = listcopy->nia2apage;
  list->nia2iapage = listcopy->nia2iapage;
  list->nia2ia_indexpage = listcopy->nia2ia_indexpage;

}
