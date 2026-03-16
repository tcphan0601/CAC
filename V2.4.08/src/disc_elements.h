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
  int disc_style;           // 0 discretize to atoms
                            // 1 discretize to smaller size elements
  int origin_shape_id;
  int split_cell,split_dim;
  int split_style;
  int split_dir,split_axis,larger_side;
  int *disc_flag_list;
  int migrate_flag,compress_flag;
  int me;

  void list_disc_region(int, char **);
  void list_disc_group(int, char **);
  void options(int, char **);
};

}

#endif
#endif

