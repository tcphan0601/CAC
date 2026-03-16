#ifdef FIX_CLASS

FixStyle(print, FixPrint)

#else

#ifndef CAC_FIX_PRINT_H
#define CAC_FIX_PRINT_H

#include <stdio.h>
#include "fix.h"

namespace CAC_NS {

class FixPrint : public Fix {
 public:
  FixPrint(class CAC *, int, char **);
  ~FixPrint();
  int setmask();
  void setup(int);
  void end_of_step();
  //void setup(int);

 private:
  int me, screenflag;
  FILE *fp;
  char *string, *copy, *work;
  int maxcopy, maxwork;
};

}

#endif
#endif
