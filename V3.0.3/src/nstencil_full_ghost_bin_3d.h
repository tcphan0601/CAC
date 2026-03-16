#ifdef NSTENCIL_CLASS

NStencilStyle(full/ghost/bin/3d, 
              NStencilFullGhostBin3d, 
              NS_FULL | NS_GHOST | NS_BIN | NS_3D |
              NS_ORTHO | NS_TRI);

#else

#ifndef CAC_NSTENCIL_FULL_GHOST_BIN_3D_H
#define CAC_NSTENCIL_FULL_GHOST_BIN_3D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilFullGhostBin3d : public NStencil {
 public:
  NStencilFullGhostBin3d(class CAC *);
  ~NStencilFullGhostBin3d() {}
  void create();
};

}

#endif
#endif
