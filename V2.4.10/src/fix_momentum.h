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

 protected:
  int linear,angular,rescale;
  int xflag,yflag,zflag;
  double masstotal;
};

}

#endif
#endif
