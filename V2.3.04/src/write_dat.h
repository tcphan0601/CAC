#ifdef COMMAND_CLASS

CommandStyle(write_dat,WriteDat)

#else

#ifndef CAC_WRITE_DAT_H
#define CAC_WRITE_DAT_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDat : protected Pointers {
 public:
  WriteDat(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void atoms();
  void node_connect();
  void nodes();
};

}

#endif
#endif


