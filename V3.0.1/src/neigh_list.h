#ifndef CAC_NEIGH_LIST_H
#define CAC_NEIGH_LIST_H

#include "pointers.h"
#include "my_page.h"

namespace CAC_NS {

class NeighList : protected Pointers {
 public:
  int index;                   // index of which neigh list this is
                               // also indexes the request it came from
                               // and the npair list of NPair classes

  int bin_method;        // 0 if no binning, else 1-N index into binnames
  int stencil_method;    // 0 if no stencil, else 1-N index into stencilnames
  int pair_method;       // 0 if no pair, else 1-N index into pairnames
  int threebody;         // 1 if for threebody pair, 0 otherwise

  // settings from NeighRequest
  int atomlist;              // 1 if include atom neighbor list
  int atomonlylist;          // 1 if include atom-atom pair only neighbor list
  int elemlist;              // 1 if include element neighbor list
  int nodelist;              // 1 if include node neighbor list
  int vatomlist;             // 1 if include virtual atom neighbor list
  int gausslist;              // 1 if include gaussian cell neighbor list

  int occasional;                  // 0 if build every reneighbor, 1 if not
  int ghost;                       // 1 if list stores neighbors of ghosts
//  int ssa;                         // 1 if list stores Shardlow data
  int copy;                        // 1 if this list is (host) copied from another list
//  int copymode;                    // 1 if this is a Kokkos on-device copy
  int sort;                        // 1 if perform neighbor sorting
  int sortnum;                     // number of nearest neighbor sorted and stored 

  int neighneigh;            // 1 if need to access other vatoms/nodes/gauss points neighbor list during inum loop
  int *element_firstindex;   // list of indices of the element first vatom/node/gauss point in ilist
  int maxelem;
                                   // will throw error if not enough neighbor founded within cutoff
                                   // = 0 if just sort all neighbors within cutoff
                                   
  int group;                 // 1 if only find neighbor in group
  int groubit;

                                   
  int dnum;                        // # of doubles per neighbor, 0 if none

  // data structs for storing neighbor pairs I, J and associated values
  // a = atom
  // g = gauss point
  // n = node
  // va = virtual atom
  // nva = neighboring virtual atom

  int ainum;                  // # of I atoms that neighbors are stored for
  int einum;                  // # of I elements that neighbors are stored for
  int ginum;                  // # of I gaussian cells that neighbors are stored for
  int ninum;                  // # of I nodes that neighbors are stored for
  int inum;                   // # of I in the entire list
  bigint vinum;               // # of I virtual atoms that neighbors are stored for
  int nvainum;                // # of I nia that neighbors are stored for

  int gnum;                   // # of ghost atoms neighbors are stored for

  int *ilist;                 // local indices of I atoms
  int *iindexlist;            // local indices of I atoms


  // a combined neighbor list for everything
  // the order is as follow: 
  // 1. atom
  // 2. node
  // 3. gauss
  // 4. virual atom
  // For neighboring atoms, index = -1
  // For neighboring virtual atoms, index = iucell * apc + ibasis
  //                                iucell = index / apc 
  //                                ibasis = index % apc
  
  int *numneigh;
  int **firstneigh;
  int **firstneighindex;
  MyPage<int> *ipage;     
  MyPage<int> *indexpage; 

  // neighbor list for neighboring virtual atoms (nva)
  
  int *nvailist;
  int *nvaindexlist;
  int ***va2nvalist;          // mapping from virtual atoms (va) to index in neighbor list
                              // = -1 if not a nva
  int *nvanumneigh;
  int **nvafirstneigh;
  int **nvafirstneighindex;
  MyPage<int> *nvaipage;     
  MyPage<int> *nvaindexpage; 


  double **firstdouble;           // ptr to 1st J double value of each I atom
  bigint maxall;                     // size of allocated combined arrays
  bigint maxnva;                  // size of allocated per-nia arrays
  int maxnva_elem;

  int pgsize;                     // size of each page
  int oneatom;                    // max size for one atom
                              
  MyPage<double> *dpage;           // pages of neighbor doubles, if dnum > 0

  // atom types to skip when building list
  // copied info from corresponding request into realloced vec/array

//  int *iskip;         // iskip[i] = 1 if atoms of type I are not in list
//  int **ijskip;       // ijskip[i][j] = 1 if pairs of type I, J are not in list

  // settings and pointers for related neighbor lists and fixes

  NeighList *listcopy;          // me = copy list, point to list I copy from
//  NeighList *listskip;          // me = skip list, point to list I skip from
  NeighList *listfull;          // me = half list, point to full I derive from
/* 
  NeighList *listhistory;       // list storing neigh history
  class Fix *fix_history;       // fix that stores history info

  int respamiddle;              // 1 if this respaouter has middle list
  NeighList *listinner;         // me = respaouter, point to respainner
  NeighList *listmiddle;        // me = respaouter, point to respamiddle

  class Fix *fix_bond;          // fix that stores bond info

  // Kokkos package

  int kokkos;                   // 1 if list stores Kokkos data
  ExecutionSpace execution_space;

  // USER-DPD package and Shardlow Splitting Algorithm (SSA) support

  uint16_t ( * ndxAIR_ssa)[8]; // for each atom, last neighbor index of each AIR
 */
  // methods

  NeighList(class CAC *);
  virtual ~NeighList();
  void post_constructor(class NeighRequest *);
  void setup_pages(int, int);           // setup page data structures
  void grow(int, int);
  void grow_nva(int, int);
  void print_attributes();              // debug routine
  bigint memory_usage();
};

}

#endif
