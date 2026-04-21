#ifdef NPAIR_CLASS

NPairStyle(halffull/newton, 
           NPairHalffullNewton, 
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_ATOM | NP_GAUSS)

#else

#ifndef CAC_NPAIR_HALFFULL_NEWTON_H
#define CAC_NPAIR_HALFFULL_NEWTON_H

#include "npair.h"

namespace CAC_NS {

class NPairHalffullNewton : public NPair {
 public:
  NPairHalffullNewton(class CAC *);
  ~NPairHalffullNewton() {}
  void build(class NeighList *);
};

}

#endif
#endif

