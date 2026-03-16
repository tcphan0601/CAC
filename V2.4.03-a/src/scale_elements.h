#ifdef COMMAND_CLASS

CommandStyle(scale_elements,ScaleElements)

#else

#ifndef CAC_SCALEELEMENTS_H
#define CAC_SCALEELEMENTS_H

#include "pointers.h"

namespace CAC_NS {

class ScaleElements : protected Pointers {
 public:
  ScaleElements(class CAC *);
  ~ScaleElements();
  void command(int, char **);

};

}

#endif
#endif

