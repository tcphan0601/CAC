#ifdef NPAIR_CLASS

NPairStyle(half/bin/newtoff, 
           NPairHalfBinNewtoff, 
           NP_HALF | NP_BIN | NP_NEWTOFF |
           NP_ORTHO | NP_TRI | NP_ATOM | NP_GAUSS)

#else

#ifndef CAC_NPAIR_HALF_BIN_NEWTOFF_H
#define CAC_NPAIR_HALF_BIN_NEWTOFF_H

#include "npair.h"

namespace CAC_NS {

class NPairHalfBinNewtoff : public NPair {
 public:
  NPairHalfBinNewtoff(class CAC *);
  ~NPairHalfBinNewtoff() {}
  void build(class NeighList *);
};

}

#endif
#endif
