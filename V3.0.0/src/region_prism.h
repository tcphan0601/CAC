#ifdef REGION_CLASS

RegionStyle(prism, RegPrism)

#else

#ifndef CAC_REGION_PRISM_H
#define CAC_REGION_PRISM_H

#include "region.h"

namespace CAC_NS {

class RegPrism : public Region {
  friend class CreateBox;

 public:
  RegPrism(class CAC *, int, char **);
  ~RegPrism();
  int inside(double, double, double);
  //int surface_interior(double *, double);
  //int surface_exterior(double *, double);

 private:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;
  double h[3][3], hinv[3][3];
  int dimension;
  double a[3], b[3], c[3];       // edge vectors of region
  double clo[3], chi[3];        // opposite corners of prism
  double face[6][3];           // unit normals of 6 prism faces
  double corners[8][3];        // 8 corner pts of prism
  int tri[12][3];              // 3 corner pts of 12 triangles (2 per face)

  //void find_nearest(double *, double &, double &, double &);
  //int inside_tri(double *, double *, double *, double *, double *);
  //double closest(double *, double *, double *, double);
};

}

#endif
#endif

