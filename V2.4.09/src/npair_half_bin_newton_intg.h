#ifdef NPAIR_CLASS

NPairStyle(half/bin/newton/intg,
           NPairHalfBinNewtonIntg,
           NP_HALF | NP_BIN | NP_NEWTON | NP_ORTHO |
           NP_TRI | NP_ATOM | NP_ELEM | NP_INTG)

#else

#ifndef CAC_NPAIR_HALF_BIN_NEWTON_INTG_H
#define CAC_NPAIR_HALF_BIN_NEWTON_INTG_H

#include "npair.h"

namespace CAC_NS {

class NPairHalfBinNewtonIntg : public NPair {
 public:
  NPairHalfBinNewtonIntg(class CAC *);
  ~NPairHalfBinNewtonIntg() {}
  void build(class NeighList *);
};

}

#endif
#endif
