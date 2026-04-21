#ifdef MINIMIZE_CLASS

MinimizeStyle(sd, MinSD)

#else

#ifndef CAC_MIN_SD_H
#define CAC_MIN_SD_H

#include "min_linesearch.h"

namespace CAC_NS {

class MinSD : public MinLineSearch {
 public:
  MinSD(class CAC *);
  int iterate(int);
};

}

#endif
#endif
