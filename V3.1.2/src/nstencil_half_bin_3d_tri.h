#ifdef NSTENCIL_CLASS

NStencilStyle(half/bin/3d/tri, 
              NStencilHalfBin3dTri, 
              NS_HALF | NS_BIN | NS_3D | NS_TRI);

#else

#ifndef CAC_NSTENCIL_HALF_BIN_3D_TRI_H
#define CAC_NSTENCIL_HALF_BIN_3D_TRI_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilHalfBin3dTri : public NStencil {
 public:
  NStencilHalfBin3dTri(class CAC *);
  ~NStencilHalfBin3dTri() {}
  void create();
};

}

#endif
#endif
