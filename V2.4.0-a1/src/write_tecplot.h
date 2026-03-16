#ifdef COMMAND_CLASS

CommandStyle(write_tecplot,WriteTecplot)

#else

#ifndef CAC_WRITE_TECPLOT_H
#define CAC_WRITE_TECPLOT_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteTecplot : protected Pointers {
 public:
  WriteTecplot(class CAC *);
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


