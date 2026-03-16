#ifdef REGION_CLASS

RegionStyle(block,RegBlock)

#else

#ifndef CAC_REGION_BLOCK_H
#define CAC_REGION_BLOCK_H

#include "region.h"

namespace CAC_NS {

class RegBlock : public Region {

  public: 
   RegBlock(class CAC *, int, char **);
   ~RegBlock();
   int inside(double *);

  private:
   double xlo,xhi,ylo,yhi,zlo,zhi;
};

}

#endif
#endif
