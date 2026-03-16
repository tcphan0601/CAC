#ifdef NSTENCIL_CLASS

NStencilStyle(full/bin/3d,
              NStencilFullBin3d,
              NS_FULL | NS_BIN | NS_3D)

#else

#ifndef CAC_NSTENCIL_FULL_BIN_3D_H
#define CAC_NSTENCIL_FULL_BIN_3D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilFullBin3d : public NStencil {
 public:
  NStencilFullBin3d(class CAC *);
  ~NStencilFullBin3d() {}
  void create();
};

}

#endif
#endif
