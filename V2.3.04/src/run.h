#ifdef COMMAND_CLASS

CommandStyle(run,Run)

#else

#ifndef CAC_RUN_H
#define CAC_RUN_H

#include "pointers.h"

namespace CAC_NS {

class Run : protected Pointers {
 public:
  Run(class CAC *);
  void command(int, char **);
};

}

#endif
#endif
