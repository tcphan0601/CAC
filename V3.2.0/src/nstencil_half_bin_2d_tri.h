#ifdef NSTENCIL_CLASS

NStencilStyle(half/bin/2d/tri, 
              NStencilHalfBin2dTri, 
              NS_HALF | NS_BIN | NS_2D | NS_TRI);

#else

#ifndef CAC_NSTENCIL_HALF_BIN_2D_TRI_H
#define CAC_NSTENCIL_HALF_BIN_2D_TRI_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilHalfBin2dTri : public NStencil {
 public:
  NStencilHalfBin2dTri(class CAC *);
  ~NStencilHalfBin2dTri() {}
  void create();
};

}

#endif
#endif
