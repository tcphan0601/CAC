#ifdef NSTENCIL_CLASS

NStencilStyle(full/bin/2d, 
              NStencilFullBin2d, 
              NS_FULL | NS_BIN | NS_2D | NS_ORTHO | NS_TRI);

#else

#ifndef CAC_NSTENCIL_FULL_BIN_2D_H
#define CAC_NSTENCIL_FULL_BIN_2D_H

#include "nstencil.h"

namespace CAC_NS {

class NStencilFullBin2d : public NStencil {
 public:
  NStencilFullBin2d(class CAC *);
  ~NStencilFullBin2d() {}
  void create();
};

}

#endif
#endif
