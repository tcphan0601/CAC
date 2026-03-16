#ifdef COMMAND_CLASS

CommandStyle(write_data,WriteData)

#else

#ifndef CAC_WRITE_DATA_H
#define CAC_WRITE_DATA_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteData : protected Pointers {
 public:
  WriteData(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void header();
  void type_arrays();
  void atoms();
  void elements();
  void nodes();
  void atom_velocities();
  void node_velocities();
};

}

#endif
#endif


