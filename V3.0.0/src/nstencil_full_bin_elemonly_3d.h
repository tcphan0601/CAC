#ifdef NSTENCIL_CLASS

NStencilStyle(full/bin/elemonly/3d, 
              NStencilFullBinElemonly3d, 
              NS_FULL | NS_BIN | NS_3D |
              NS_ORTHO | NS_TRI | NS_ELEMONLY);

#else

#ifndef CAC_NSTENCIL_FULL_BIN_ELEMONLY_3D_H
#define CAC_NSTENCIL_FULL_BIN_ELEMONLY_3D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilFullBinElemonly3d : public NStencil {
 public:
  NStencilFullBinElemonly3d(class CAC *);
  ~NStencilFullBinElemonly3d() {}
  void create();
  void create_setup();
};

}

#endif
#endif
