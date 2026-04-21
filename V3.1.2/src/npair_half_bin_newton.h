#ifdef NPAIR_CLASS

NPairStyle(half/bin/newton, 
           NPairHalfBinNewton, 
           NP_HALF | NP_BIN | NP_NEWTON | 
           NP_ORTHO | NP_ATOM | NP_GAUSS)
 
#else

#ifndef CAC_NPAIR_HALF_BIN_NEWTON_H
#define CAC_NPAIR_HALF_BIN_NEWTON_H

#include "npair.h"

namespace CAC_NS {

class NPairHalfBinNewton : public NPair {
 public:
  NPairHalfBinNewton(class CAC *);
  ~NPairHalfBinNewton() {}
  void build(class NeighList *);
};

}

#endif
#endif
