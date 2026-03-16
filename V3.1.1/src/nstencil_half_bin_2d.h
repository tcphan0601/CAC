#ifdef NSTENCIL_CLASS

NStencilStyle(half/bin/2d, 
              NStencilHalfBin2d, 
              NS_HALF | NS_BIN | NS_2D | NS_ORTHO);

#else

#ifndef CAC_NSTENCIL_HALF_BIN_2D_H
#define CAC_NSTENCIL_HALF_BIN_2D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilHalfBin2d : public NStencil {
 public:
  NStencilHalfBin2d(class CAC *);
  ~NStencilHalfBin2d() {}
  void create();
};

}

#endif
#endif
