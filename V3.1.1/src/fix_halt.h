#ifdef FIX_CLASS

FixStyle(halt, FixHalt)

#else

#ifndef CAC_FIX_HALT_H
#define CAC_FIX_HALT_H

#include "fix.h"

namespace CAC_NS {

class FixHalt : public Fix {
 public:
  FixHalt(class CAC *, int, char **);
  ~FixHalt();
  int setmask();
  void init();
  void end_of_step();
  void min_post_force(int);
  void post_run();

 private:
  int attribute, operation, eflag, msgflag, ivar;
  bigint nextstep, thisstep;
  double value, tratio;
  char *idvar;

  //double bondmax();
  double tlimit();
};

}

#endif
#endif

