#ifdef FIX_CLASS

FixStyle(aveforce, FixAveForce)

#else

#ifndef CAC_FIX_AVEFORCE_H
#define CAC_FIX_AVEFORCE_H

#include "fix.h"

namespace CAC_NS {

class FixAveForce : public Fix {
 public:
  FixAveForce(class CAC *, int, char **);
  ~FixAveForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  //void post_force_respa(int, int, int);
  //void min_post_force(int);
  double compute_vector(int);

 private:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  char *idregion;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  int iregion;
  double foriginal_all[4];
  //int nlevels_respa, ilevel_respa;
};

}

#endif
#endif
