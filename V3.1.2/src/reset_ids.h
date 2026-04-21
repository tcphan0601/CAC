#ifdef COMMAND_CLASS

CommandStyle(reset_ids, ResetIDs)

#else

#ifndef CAC_RESET_IDS_H
#define CAC_RESET_IDS_H

#include "pointers.h"

namespace CAC_NS {

class ResetIDs : protected Pointers {
 public:
  int order_flag;
  ResetIDs(class CAC *);
  void command(int, char **);

 private:
};

}

#endif
#endif
