#ifdef COMMAND_CLASS

CommandStyle(write_dump_ghost_debug,WriteDumpGhostDebug)

#else

#ifndef CAC_WRITE_DUMP_GHOST_DEBUG_H
#define CAC_WRITE_DUMP_GHOST_DEBUG_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDumpGhostDebug : protected Pointers {
 public:
  WriteDumpGhostDebug(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;

  void header();
};

}

#endif
#endif


