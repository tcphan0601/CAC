#ifdef NPAIR_CLASS

NPairStyle(halffull/newton/intg,
           NPairHalffullNewtonIntg,
           NP_HALFFULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTG)

#else

#ifndef CAC_NPAIR_HALFFULL_NEWTON_INTG_H
#define CAC_NPAIR_HALFFULL_NEWTON_INTG_H

#include "npair.h"

namespace CAC_NS {

class NPairHalffullNewtonIntg : public NPair {
 public:
  NPairHalffullNewtonIntg(class CAC *);
  ~NPairHalffullNewtonIntg() {}
  void build(class NeighList *);
};

}

#endif
#endif

