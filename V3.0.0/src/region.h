#ifndef CAC_REGION_H
#define CAC_REGION_H

#include "pointers.h"

namespace CAC_NS {

class Region : protected Pointers {
 public:
  char *id, *style;
  int interior;                               // 1 for interior, 0 for exterior
  int scaleflag;                              // 1 for lattice, 0 for box
  double xscale, yscale, zscale;              // scale factors for box/lattice units
  double extent_xlo, extent_xhi;              // bounding box on region
  double extent_ylo, extent_yhi;
  double extent_zlo, extent_zhi;  
  int bboxflag;                               // 1 if bounding box is computable
  int varshape;                               // 1 if region shape changes over time
  int dynamic;                                // 1 if position/orient changes over time
  int moveflag, rotateflag;                    // 1 if position/orientation changes
  int openflag;			                          // 1 if any face is open
  int open_faces[6];		                      // flags for which faces are open


  // contact = particle near region surface

  struct Contact {
    double r;                                 // distance between particle & surf, r > 0.0
    double delx, dely, delz;                  // vector from surface pt to particle
    double radius;                            // curvature of region at contact point
    int iwall;                                // unique id of wall for storing shear history
    int varflag;                              // 1 if wall can be variable-controlled

  };
  Contact *contact;                           // list of contacts
  int cmax;                                   // max # of contacts possible with region
  int tmax;                                   // max # of touching contacts possible

  double dx, dy, dz, theta;      // current displacement and orientation
  Region(class CAC *, int, char **);
  virtual ~Region();
  virtual void init();
  int dynamic_check();

  // called by other classes to check point versus region
  
  void prematch();
  int match(double, double, double);
  int match(double *);

  // implemented by each region, not called by other classes
  
  virtual int inside(double, double, double) = 0;
  //virtual int surface_interior(double *, double) = 0;
  //virtual int surface_exterior(double *, double) = 0;
  virtual void pretransform();

 protected:
  void options(int, char **);
  void forward_transform(double &, double &, double &);
  double point[3], runit[3];

 private:
  char *xstr, *ystr, *zstr, *tstr;
  int xvar, yvar, zvar, tvar;
  double axis[3];

  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, double);
};
}
#endif
