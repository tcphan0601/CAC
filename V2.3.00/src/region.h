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

  int openflag;			    // 1 if any face is open
  int open_faces[6];		    // flags for which faces are open


  // contact = particle near region surface

  struct Contact {
    double r;                                  // distance between particle & surf, r > 0.0
    double delx, dely, delz;                   // vector from surface pt to particle
  };
  Contact *contact;                           // list of contacts
  int cmax;                                   // max # of contacts possible with region

  Region(class CAC *, int, char **);
  virtual ~Region();
  virtual void init();
  void prematch();
  int match(double *);

  virtual int inside(double *) = 0;

 protected:
  void options(int, char **);

 private:
  int moveflag,rotateflag;                    // 1 if position/orientation changes

  double point[3],axis[3],runit[3];
  char *xstr, *ystr, *zstr, *tstr;
  int xvar,yvar,zvar,tvar;
  double dx,dy,dz,theta;
};
}
#endif
