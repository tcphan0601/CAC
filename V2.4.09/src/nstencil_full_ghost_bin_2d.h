#ifdef NSTENCIL_CLASS

NStencilStyle(full/ghost/bin/2d,
              NStencilFullGhostBin2d,
              NS_FULL | NS_GHOST | NS_BIN | NS_2D |
              NS_NEWTON | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef CAC_NSTENCIL_FULL_GHOST_BIN_2D_H
#define CAC_NSTENCIL_FULL_GHOST_BIN_2D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilFullGhostBin2d : public NStencil {
 public:
  NStencilFullGhostBin2d(class CAC *);
  ~NStencilFullGhostBin2d() {}
  void create();
};

}

#endif
#endif

