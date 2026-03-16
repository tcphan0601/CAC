#ifndef CAC_RANMARS_H
#define CAC_RANMARS_H

#include "pointers.h"

namespace CAC_NS {

class RanMars : protected Pointers {
 public:
  RanMars(class CAC *, int);
  ~RanMars();
  double uniform();
  double gaussian();

 private:
  int save;
  double second;
  double *u;
  int i97,j97;
  double c,cd,cm;
};

}

#endif

