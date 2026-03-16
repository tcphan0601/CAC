#ifdef REGION_CLASS

RegionStyle(cylinder,RegCylinder)

#else

#ifndef CAC_REGION_CYLINDER_H
#define CAC_REGION_CYLINDER_H

#include "region.h"

namespace CAC_NS {

class RegCylinder : public Region {

 public:
  RegCylinder(class CAC *, int, char **);
  ~RegCylinder();
  void init();
  int inside(double *);

 private:
  char axis;
  double c1,c2;
  double radius,rsq;
  double lo,hi;

};

}

#endif
#endif
