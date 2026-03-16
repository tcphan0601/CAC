#ifdef COMMAND_CLASS

CommandStyle(write_dump_subelem_debug,WriteDumpSubelemDebug)

#else

#ifndef CAC_WRITE_DUMP_SUBELEM_DEBUG_H
#define CAC_WRITE_DUMP_SUBELEM_DEBUG_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDumpSubelemDebug : protected Pointers {
 public:
  WriteDumpSubelemDebug(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int id,ielem;
  FILE *fp;

  void header();
};

}

#endif
#endif


