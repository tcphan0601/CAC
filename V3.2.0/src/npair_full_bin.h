#ifdef NPAIR_CLASS

NPairStyle(full/bin, 
           NPairFullBin, 
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_GAUSS)

#else

#ifndef CAC_NPAIR_FULL_BIN_H
#define CAC_NPAIR_FULL_BIN_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBin : public NPair {
 public:
  NPairFullBin(class CAC *);
  ~NPairFullBin() {}
  void build(class NeighList *);
};

}

#endif
#endif
