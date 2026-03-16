#include "npair_copy.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "my_page.h"


using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

NPairCopy::NPairCopy(CAC *cac) : NPair(cac) {}

/*  ----------------------------------------------------------------------
   create list which is simply a copy of parent list
   update this list if adding new neighbor lists
-------------------------------------------------------------------------  */

void NPairCopy::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->ainum = listcopy->ainum;
  list->einum = listcopy->einum;
  list->ginum = listcopy->ginum;
  list->vinum = listcopy->vinum;
  list->nvainum = listcopy->nvainum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->iindexlist = listcopy->iindexlist;
  list->va2nvalist = listcopy->va2nvalist;
  list->nvailist = listcopy->nvailist;
  list->nvaindexlist = listcopy->nvaindexlist;
  
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->firstneighindex = listcopy->firstneighindex;
  list->ipage = listcopy->ipage;
  list->indexpage = listcopy->indexpage;

  list->nvanumneigh = listcopy->nvanumneigh;
  list->nvafirstneigh = listcopy->nvafirstneigh;
  list->nvafirstneighindex = listcopy->nvafirstneighindex;
  list->nvaipage = listcopy->nvaipage;
  list->nvaindexpage = listcopy->nvaindexpage;

}
