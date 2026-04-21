#ifdef COMMAND_CLASS

CommandStyle(write_data_ghost_debug,WriteDataGhostDebug)

#else

#ifndef CAC_WRITE_DATA_GHOST_DEBUG_H
#define CAC_WRITE_DATA_GHOST_DEBUG_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDataGhostDebug : protected Pointers {
 public:
  WriteDataGhostDebug(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int max_tag_atom;
  int max_tag_elem;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void header();
  void atoms();
  void elements();
  void nodes();
};

}

#endif
#endif


