#ifdef COMMAND_CLASS

CommandStyle(write_tecplot,WriteTecplot)

#else

#ifndef CAC_WRITE_TECPLOT_H
#define CAC_WRITE_TECPLOT_H

#include "TECIO.h"
#include <stdio.h>
#include <stdlib.h>
#include "pointers.h"

namespace CAC_NS {

class WriteTecplot : protected Pointers {
 public:
  WriteTecplot(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int datapackingflag;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  // parameters for TECIO

  int debug,zonetype;
  double soltime;
  int fileformat; // 0 == PLT, 1 == SZPLT, 2 == ASCII
  int dummy,dummy1;

  // functions to output each section

  void atoms();
  void node_connect();
  void nodes();
};

}

#endif
#endif


