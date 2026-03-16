#ifdef COMMAND_CLASS

CommandStyle(create_box,CreateBox)

#else

#ifndef CAC_CREATE_BOX_H
#define CAC_CREATE_BOX_H

#include "pointers.h"

namespace CAC_NS {

class CreateBox : protected Pointers {
 public:
  CreateBox(class CAC *);
  void command(int, char **);
};

}

#endif
#endif

