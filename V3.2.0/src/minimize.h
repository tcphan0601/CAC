#ifdef COMMAND_CLASS

CommandStyle(minimize,Minimize)

#else

#ifndef CAC_MINIMIZE_H
#define CAC_MINIMIZE_H

#include "pointers.h"

namespace CAC_NS {

class Minimize : protected Pointers {
 public:
  Minimize(class CAC *);
  void command(int, char **);
};

}

#endif
#endif
