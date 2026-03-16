#ifdef FIX_CLASS

FixStyle(enforce2d, FixEnforce2D)

#else

#ifndef CAC_FIX_ENFORCE2D_H
#define CAC_FIX_ENFORCE2D_H

#include "fix.h"

namespace CAC_NS {

class FixEnforce2D : public Fix {
 public:
  FixEnforce2D(class CAC *, int, char **);
  ~FixEnforce2D();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);

 private:
  int nfixlist;
  class Fix **flist;
};

}

#endif
#endif

