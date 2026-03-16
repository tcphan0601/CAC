#ifdef COMMAND_CLASS

CommandStyle(scale,Scale)

#else

#ifndef CAC_SCALE_H
#define CAC_SCALE_H

#include "pointers.h"

namespace CAC_NS {

class Scale : protected Pointers {
 public:
  Scale(class CAC *);
  ~Scale();
  void command(int, char **);

};

}

#endif
#endif

