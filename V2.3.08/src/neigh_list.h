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
  int elemlist;           // 1 if include element neighbor list
  int intpllist;             // 1 if include interpolated atom neighbor list
  int intglist;              // 1 if include integration point neighbor list

  int occasional;                  // 0 if build every reneighbor, 1 if not
  int ghost;                       // 1 if list stores neighbors of ghosts
//  int ssa;                         // 1 if list stores Shardlow data
  int copy;                        // 1 if this list is (host) copied from another list
//  int copymode;                    // 1 if this is a Kokkos on-device copy
  int sort;                        // 1 if perform neighbor sorting
  int sortnum;                     // number of nearest neighbor sorted and stored 
                                   // will throw error if not enough neighbor founded within cutoff
                                   // = 0 if just sort all neighbors within cutoff
  int dnum;                        // # of doubles per neighbor, 0 if none

  // data structs to store neighbor pairs I,J and associated values
  // a = atom
  // i = integration point
  // ia = interpolated atom
  // nia = neighboring interpolated atom

  int ainum;                  // # of I atoms that neighbors are stored for
  int einum;                  // # of I elements that neighbors are stored for
  int iinum;                  // # of I integration points that neighbors are stored for
  bigint iainum;              // # of I interpolated atoms that neighbors are stored for
  int niainum;                // # of I nia that neighbors are stored for
  int gnum;                        // # of ghost atoms neighbors are stored for

  int *ailist;                // local indices of I atoms
  int *eilist;                // local indices of I elements
  int **n2ilist;              // local interpolated indices of element node K in element I
  int *e2ilist;               // local indices of first integration point of element I
  bigint *e2ialist;           // local indices of first interpolated point of element I
  int *ia2elist;              // local indices of parent element of interpolated atoms I
  int **ia2nialist;           // mapping from interpolated atoms to index in nia neighbor list
                              // = -1 if not a nia

  int *nia2ialist;
  int *nia2ia_indexlist;

  int *numneigha2a;           // # of J atom neighbors for each 
                              // I atom (size of array should
                              // be # of atoms: ainum)
                              
  int *numneigha2e;           // # of J element neighbors for each 
                              // I atom (size of array should 
                              // be # of atoms: ainum)
                              
  int *numneigha2ia;          // # of J interpolated atom neighbors
                              // for each I atom (size of array
                              // should be # of atoms: ainum)
                              
  int *numneighe2e;           // # of J element neighbors for each
                              // I element (size of array should
                              // be # of elements: einum)
                              
  int *numneighi2ia;          // # of J interpolated atom neighbors
                              // for each I integration point 
                              // (size of array should be
                              // # of integration points: iinum)
                              
  int *numneighi2a;           // # of J "real" atom neighbors 
                              // for each I integration point 
                              // (size of array should be 
                              // # of integration points: iinum)
                              //
  int *numneighia2ia;         // # of J interpolated atom neighbors
                              // for each I interpolated atoms
                              // (size of array should be
                              // # of integration points: iinum)
                              
  int *numneighia2a;          // # of J "real" atom neighbors 
                              // for each I interpolated atoms
                              // (size of array should be 
 
  int **firstneigha2a;        // ptr to array (size numneigha2a) 
                              // of J atom tag of each I atom 
                              // (ainum pointers)
                              
  int **firstneigha2e;        // ptr to array (size numneigha2e) 
                              // of J element tag of each I atom h
                              // (ainum pointers)
                                      
  int **firstneigha2ia;       // ptr to array (size numneigha2ia)
                              // of J interpolated point's parent
                              // elelemt tag of each I atom 
                              // (ainum pointers)
  
  int **firstneigha2ia_index; // ptr to array (size numneigha2ia) 
                              // of J interpolated point's index 
                              // within its parent elelemt of 
                              // each I atom (ainum pointers)
  
  int **firstneighe2e;        // ptr to array (size numneighe2e) 
                              // of J element tag of each I element 
                              // (einum pointers)

  int **firstneighi2a;        // ptr to array (size numneighi2a) 
                              // of J atom tag of each I 
                              // integration point (iinum pointers)

  int **firstneighi2ia;       // ptr to array (size numneighi2ia)
                              // of J interpolated atom's parent 
                              // element tag of each I integration
                              // point (iinum pointers)

  int **firstneighi2ia_index; // ptr to array (size numneighi2ia) 
                              // of J interpolated atom's index 
                              // within its parent element of each 
                              // I integration point (iinum pointers)
                              
  int **firstneighia2a;       // ptr to array (size numneighi2a) 
                              // of J atom tag of each I 
                              // interpolated atoms (iainum pointers)

  int **firstneighia2ia;      // ptr to array (size numneighi2ia)
                              // of J interpolated atom's parent 
                              // element tag of each I interpolated
                              // atoms (iainum pointers)

  int **firstneighia2ia_index;// ptr to array (size numneighi2ia) 
                              // of J interpolated atom's index 
                              // within its parent element of each 
                              //
  // for neighbor of interpolated atom neighbors
  // used for threebody potential

  int *numneighnia2a;
  int *numneighnia2ia;
  int **firstneighnia2a;
  int **firstneighnia2ia;
  int **firstneighnia2ia_index;
  MyPage<int> *nia2apage;         // pages of nia to atom neighbor indices
  MyPage<int> *nia2iapage;        // pages of nia to interpolated atom
  MyPage<int> *nia2ia_indexpage;  // pages of nia to interpolated atom 
 
  // for sorted neighbor list
  
  int *numneigha;
  int *numneighia;
  int **firstneigha;
  int **firstneigha_index;
  int **firstneighia;
  int **firstneighia_index;
  MyPage<int> *apage;         // pages of atom neighbor indices
  MyPage<int> *a_indexpage;   // pages of atom neighbor intpl indices (-1 for atom)
  MyPage<int> *iapage;        // pages of interpolated atom neighbor indices
  MyPage<int> *ia_indexpage;  // pages of interpolated atom neighbor intpl indices (-1 for atom)
 
  double **firstdouble;           // ptr to 1st J double value of each I atom
  int maxatom;                    // size of allocated per-atom arrays
  int maxelem;                    // size of allocated per-element arrays
  int maxelem_conv;               // size of allocated per-element convert to intg or intpl arrays
  int maxintg;                    // size of allocated per-intg arrays
  bigint maxnia;                  // size of allocated per-nia arrays
  bigint maxintpl;                // size of allocated per-intpl arrays
  bigint maxintpl_conv;           // size of allocated per-intpl convert to elem arrays

  int pgsize;                     // size of each page
  int oneatom;                    // max size for one atom
  MyPage<int> *a2apage;           // pages of atom to atom neighbor indices
  MyPage<int> *a2iapage;          // pages of atom to interpolated atom
                                  // neighbor's parent element indices
  MyPage<int> *a2ia_indexpage;    // pages of atom to interpolated atom 
                                  // neighbor's indices within its parent 
                                  // element
  MyPage<int> *i2apage;           // pages of integration point to atom 
                                  // neighbor indices
  MyPage<int> *e2epage;           // pages of element to element neighbor
                                  // indices
  MyPage<int> *i2iapage;          // pages of integration point to interpolated
                                  // atom neighbor's parent element indices
  MyPage<int> *i2ia_indexpage;    // pages of integration point to 
                                  // interpolated atom neighbor's indices 
                                  // within its parent element
  MyPage<int> *ia2apage;          // pages of interpolated atom to atom 
                                  // neighbor indices
                              
  MyPage<int> *ia2iapage;         // pages of interpolated atom to interpolated
                                  // atom neighbor's parent element indices
  MyPage<int> *ia2ia_indexpage;   // pages of interpolated atom to 
                                  // interpolated atom neighbor's indices 
                                  // within its parent element
                            
  MyPage<double> *dpage;           // pages of neighbor doubles, if dnum > 0

  // atom types to skip when building list
  // copied info from corresponding request into realloced vec/array

//  int *iskip;         // iskip[i] = 1 if atoms of type I are not in list
//  int **ijskip;       // ijskip[i][j] = 1 if pairs of type I,J are not in list

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

  uint16_t (*ndxAIR_ssa)[8]; // for each atom, last neighbor index of each AIR
*/
  // methods

  NeighList(class CAC *);
  virtual ~NeighList();
  void post_constructor(class NeighRequest *);
  void setup_pages(int, int);           // setup page data structures
  void grow_atom(int, int);                       // grow size for per atoms array
  void grow_elem(int, int);                       // grow size for per elements array
  void grow_intg(int, int);                       // grow size for per integration points array
  void grow_intpl(bigint, bigint);                       // grow size for per interpolated atoms array
  void grow_nia(bigint);
  void print_attributes();              // debug routine
  bigint memory_usage();
};

}

#endif
