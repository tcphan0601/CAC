#ifdef FIX_CLASS

FixStyle(addforce, FixAddForce)

#else

#ifndef CAC_FIX_ADDFORCE_H
#define CAC_FIX_ADDFORCE_H

#include "fix.h"

namespace CAC_NS {

class FixAddForce : public Fix {
 public:
  FixAddForce(class CAC *, int, char **);
  ~FixAddForce();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double xvalue, yvalue, zvalue;
  int varflag, iregion;
  char *xstr, *ystr, *zstr;
  char *idregion;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal[4], foriginal_all[4];
  int force_flag;

  int maxatom;
  double **sforce;
};

}

#endif
#endif

