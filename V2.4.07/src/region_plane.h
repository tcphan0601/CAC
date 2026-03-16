#ifdef REGION_CLASS

RegionStyle(plane,RegPlane)

#else

#ifndef CAC_REGION_PLANE_H
#define CAC_REGION_PLANE_H

#include "region.h"

namespace CAC_NS {

class RegPlane : public Region {
 public:
  RegPlane(class CAC *, int, char **);
  ~RegPlane();
  int inside(double, double, double);
  //int surface_interior(double *, double);
  //int surface_exterior(double *, double);

 private:
  double xp,yp,zp;
  double normal[3];
};

}

#endif
#endif
