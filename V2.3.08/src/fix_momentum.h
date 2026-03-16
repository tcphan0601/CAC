#ifdef FIX_CLASS

FixStyle(momentum,FixMomentum)

#else

#ifndef CAC_FIX_MOMENTUM_H
#define CAC_FIX_MOMENTUM_H

#include "fix.h"

namespace CAC_NS {

class FixMomentum : public Fix {
 public:
  FixMomentum(class CAC *, int, char **);
  int setmask();
  void init();
  void end_of_step();
  double compute_vector(int);

 protected:
  int linear,angular,rescale;
  int xflag,yflag,zflag;
  double masstotal, lcm[6];
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Fix momentum group has no atoms

Self-explanatory.

*/
