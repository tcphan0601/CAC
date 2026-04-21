#ifdef NPAIR_CLASS

NPairStyle(half/bin/newton/tri, 
           NPairHalfBinNewtonTri, 
           NP_HALF | NP_BIN | NP_NEWTON | 
           NP_TRI | NP_ATOM | NP_GAUSS)

#else

#ifndef CAC_NPAIR_HALF_BIN_NEWTON_TRI_H
#define CAC_NPAIR_HALF_BIN_NEWTON_TRI_H

#include "npair.h"

namespace CAC_NS {

class NPairHalfBinNewtonTri : public NPair {
 public:
  NPairHalfBinNewtonTri(class CAC *);
  ~NPairHalfBinNewtonTri() {}
  void build(class NeighList *);
};

}

#endif
#endif
