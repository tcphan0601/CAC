#ifdef NPAIR_CLASS

NPairStyle(full/bin/sort/atomonly/ghost, 
           NPairFullBinSortAtomonlyGhost, 
           NP_FULL | NP_BIN | NP_NEWTON | NP_SORT |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOMONLY | NP_GHOST)

#else

#ifndef CAC_NPAIR_FULL_BIN_SORT_ATOMONLY_GHOST_H
#define CAC_NPAIR_FULL_BIN_SORT_ATOMONLY_GHOST_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinSortAtomonlyGhost : public NPair {
 public:
  NPairFullBinSortAtomonlyGhost(class CAC *);
  ~NPairFullBinSortAtomonlyGhost(); 
  void build(class NeighList *);

};

}

#endif
#endif
