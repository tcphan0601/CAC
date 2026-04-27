#ifdef FIX_CLASS

FixStyle(nve, FixNVE)

#else

#ifndef CAC_FIX_NVE_H
#define CAC_FIX_NVE_H

#include "fix.h"

namespace CAC_NS {

class FixNVE : public Fix {

 public:
  FixNVE(class CAC *, int, char **);
  virtual ~FixNVE() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void reset_dt();

 protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;


};

}

#endif

#endif
