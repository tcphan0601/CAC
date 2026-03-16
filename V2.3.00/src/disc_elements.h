#ifdef COMMAND_CLASS

CommandStyle(disc_elements,DiscElements)

#else

#ifndef CAC_DISC_ELEMENTS_H
#define CAC_DISC_ELEMENTS_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class DiscElements : protected Pointers {
 public:
  DiscElements(class CAC *);
  void command(int, char **);

 private:
  int *disc_flag_list,compress_flag;
  int me;

  void list_disc_region(int, char **);
  void list_disc_group(int, char **);
  void options(int, char **);
};

}

#endif
#endif

