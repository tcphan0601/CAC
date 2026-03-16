#ifndef CAC_NEIGH_REQUEST_H
#define CAC_NEIGH_REQUEST_H

#include "pointers.h"

namespace CAC_NS {

class NeighRequest : protected Pointers {
 public:
  int index;                // index of which neigh request this is
  void *requestor;          // class that made request
  int requestor_instance;   // instance of that class (only Fix, Compute, Pair)
  int id;                   // ID of request as stored by requestor
                            // used to track multiple requests from one class

  // -----------------------------  
  // flags set by requesting class for attributes of neighor list they need
  // all must be set appropriately, all have defaults
  // -----------------------------  

  // which class style requests the list
  // one flag is 1, others are 0

  int pair;              // pair is set by default
  int fix;
  int compute;
  int command;
  int neigh;

  // half/full setting, determines which neighbors appear in list
  // one flag is 1, other is 0

  int half;              // half neigh list (set by default)
  int full;              // full neigh list
  int atomonlylist;      // 1 if include atom neighbor list (for only atom pairs)
  int elemonlylist;      // 1 if include element neighbor list (atom/elemenet neighbors of elements, not intg)
  int atomlist;          // 1 if include atom neighbor list (for both atom-atom pairs and atom-intpl pairs)
  int elemlist;          // 1 if include element neighbor list
  int nodelist;          // 1 if include node neighbor list
  int intglist;          // 1 if include integration point neighbor list
  int intpllist;         // 1 if include interpolated atom neighbor list
  int doublelist;             // 1 if include an additional outer skin neighbor list

  // attribute flags, all are 0 by default

  int occasional;        // how often list is built
                         // 0 if needed every reneighboring during run
                         // 1 if only occasionally needed by a fix, compute, etc
                         // 2 if only needed once for a command

  int newton;            // which owned/ghost pairs appear in list
                         // 0 if use force::newton_pair setting
                         // 1 if override with pair newton on
                         // 2 if override with pair newton off
  
  int ghost;             // 1 if includes ghost atom neighbors
  int size;              // 1 if pair cutoff set by particle radius
//  int history;         // 1 if stores neighbor history info
//  int granonesided;    // 1 if one-sided granular list for 
                         //   sphere/surf interactions
  int cut;               // 1 if use a non-standard cutoff length
  double cutoff;         // special cutoff distance for this list
  int sort;              // 1 if require sorting neighbors
  int sortnum;           // number of nearest neighbor needed (throw error if 
                         // number of neighbors is less than this number)
                         // = 0 if sort all neighbors within cutoff
  int group;             // 1 if only find neighbor for specific group
  int groupbit;          // group bit if only find neighbor for specific group

  int force_rebuild;     // 1 if need to build regardless of last re-neighboring
  int custom_domain;     // > 0 if using custom domain size (not from Domain class

  int noinner;           // 1 if exclude neighbor of inner interpolated atoms

  int dnum;              // # of extra floating point values stored in list

  int threebody;         // 1 if include neighbors of neighboring interpolated atoms
                         // for threebody potentials

  int atom_filter;       // 1 if neighbor atoms are filtered through compute values
  int atom_filter_icompute;
  int atom_filter_value;  
  int atom_filter_preneighflag;

  // flags set by pair hybrid

//  int skip;              // 1 if this list skips atom types from another list
//  int *iskip;            // iskip[i] if atoms of type I are not in list
//  int **ijskip;          // ijskip[i][j] if pairs of type I,J are not in list

  // command_style only set if command = 1
  // allows print_pair_info() to access command name

  const char *command_style;

  // -----------------------------  
  // flags set by Neighbor class to morph original requests
  // -----------------------------  

//  int skiplist;          // index of list to skip from
  int off2on;            // 1 if this is newton on list, but skips from off list

  int copy;              // 1 if this list copied from another list
  int copylist;          // index of list to copy from

  int halffull;          // 1 if half list computed from another full list
  int halffulllist;      // index of full list to derive half from

//  int history_partner;   // 1 if this list partners with a history list
//  int historylist;       // index of the associated history list
                         // for history = 1, index of the non-history partner
                         
  int unique;            // 1 if this list requires its own
                         // NStencil, Nbin class - because of requestor cutoff

  // pointer to FSH class, set by requestor class (not by Neighbor)

//  class FixShearHistory *fix_history;  // fix that stores per-atom history info

  // -----------------------------  
  // internal settings made by Neighbor class
  // -----------------------------  

  int index_bin;         // index of NBin class assigned to this request
  int index_stencil;     // index of NStencil class assigned to this request
  int index_pair;        // index of NPair class assigned to this request

  // methods

  NeighRequest(class CAC *);
  ~NeighRequest();
  void archive();
  int identical(NeighRequest *);
//  int same_skip(NeighRequest *);
  void copy_request(NeighRequest *);
};

}

#endif
