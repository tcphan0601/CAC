#ifdef REGION_CLASS

RegionStyle(intersect,RegIntersect)

#else

#ifndef CAC_REGION_INTERSECT_H
#define CAC_REGION_INTERSECT_H

#include "region.h"

namespace CAC_NS {

class RegIntersect : public Region {
 public:
  RegIntersect(class CAC *, int, char **);
  ~RegIntersect();
  void init();
  int inside(double, double, double);
  //int surface_interior(double *, double);
  //int surface_exterior(double *, double);
  //void shape_update();
  void pretransform();
  //void set_velocity();
  //void length_restart_string(int&);
  //void write_restart(FILE *);
  //int restart(char *, int&);
  //void reset_vel();

 private:
  char **idsub;
};

}

#endif
#endif

