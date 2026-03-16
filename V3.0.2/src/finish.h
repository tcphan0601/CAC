#ifndef CAC_FINISH_H
#define CAC_FINISH_H

#include "pointers.h"

namespace CAC_NS {

class Finish : protected Pointers {
 public:
  Finish(class CAC *);
  void end(int);

 private:
  void stats(int, double *, double *, double *, double *, int, int *);
};

}

#endif
/*  ERROR/WARNING messages:

 */
