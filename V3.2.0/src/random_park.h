#ifndef CAC_RANPARK_H
#define CAC_RANPARK_H

#include "pointers.h"

namespace CAC_NS {

class RanPark : protected Pointers {
 public:
  RanPark(class CAC *, int);
  double uniform();
  double gaussian();
  void reset(int);
  void reset(int, double *);
  int state();

 private:
  int seed, save;
  double second;
};

}

#endif

/*  ERROR/WARNING messages:

E: Invalid seed for Park random # generator

The initial seed for this random number generator must be a positive
integer.

 */
