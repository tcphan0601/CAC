#ifdef COMMAND_CLASS

CommandStyle(unwrap, Unwrap)

#else

#ifndef CAC_UNWRAP_H
#define CAC_UNWRAP_H

#include "pointers.h"

namespace CAC_NS {

class Unwrap : protected Pointers {
 public:
  Unwrap(class CAC *);
  ~Unwrap();
  void command(int, char **);

};

}

#endif
#endif

