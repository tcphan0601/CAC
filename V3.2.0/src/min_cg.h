#ifdef MINIMIZE_CLASS

MinimizeStyle(cg, MinCG)

#else

#ifndef CAC_MIN_CG_H
#define CAC_MIN_CG_H

#include "min_linesearch.h"

namespace CAC_NS {

class MinCG : public MinLineSearch {
 public:
  MinCG(class CAC *);
  int iterate(int);
};

}

#endif
#endif
