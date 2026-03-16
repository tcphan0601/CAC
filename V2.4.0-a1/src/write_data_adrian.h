#ifdef COMMAND_CLASS

CommandStyle(write_data_adrian,WriteDataAdrian)

#else

#ifndef CAC_WRITE_DATA_ADRIAN_H
#define CAC_WRITE_DATA_ADRIAN_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDataAdrian : protected Pointers {
 public:
  WriteDataAdrian(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void header();
  void type_arrays();
  //void atoms();
  void elements();
  void nodes();
  //void atom_velocities();
  //void node_velocities();
};

}

#endif
#endif


