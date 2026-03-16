#ifdef REGION_CLASS

RegionStyle(union, RegUnion)

#else

#ifndef CAC_REGION_UNION_H
#define CAC_REGION_UNION_H

#include "region.h"

namespace CAC_NS {

class RegUnion : public Region {
 public:
  RegUnion(class CAC *, int, char **);
  ~RegUnion();
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

