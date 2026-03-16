#ifdef REGION_CLASS

RegionStyle(sphere,RegSphere)

#else

#ifndef CAC_REGION_SPHERE_H
#define CAC_REGION_SPHERE_H

#include "region.h"

namespace CAC_NS {

class RegSphere : public Region {
 public:
  RegSphere(class CAC *, int, char **);
  ~RegSphere();
  void init();
  int inside(double *);

 private:
  double xc,yc,zc;
  double radius;
  int rstyle,rvar;
  char *rstr;

  void variable_check();
};

}

#endif
#endif
