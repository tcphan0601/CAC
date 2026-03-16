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
  double gaussian(double mu, double sigma);
  double rayleigh(double sigma);
  double besselexp(double theta, double alpha, double cp);
  void select_subset(bigint, int, int *, int *);


 private:
  int save;
  double second;
  double *u;
  int i97, j97;
  double c, cd, cm;
};

}

#endif

