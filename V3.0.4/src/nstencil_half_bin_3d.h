#ifdef NSTENCIL_CLASS

NStencilStyle(half/bin/3d, 
              NStencilHalfBin3d, 
              NS_HALF | NS_BIN | NS_3D | NS_ORTHO);

#else

#ifndef CAC_NSTENCIL_HALF_BIN_3D_H
#define CAC_NSTENCIL_HALF_BIN_3D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilHalfBin3d : public NStencil {
 public:
  NStencilHalfBin3d(class CAC *);
  ~NStencilHalfBin3d() {}
  void create();
};

}

#endif
#endif
