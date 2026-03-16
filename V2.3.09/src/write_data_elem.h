#ifdef COMMAND_CLASS

CommandStyle(write_data_elem,WriteDataElem)

#else

#ifndef CAC_WRITE_DATA_ELEM_H
#define CAC_WRITE_DATA_ELEM_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDataElem : protected Pointers {
 public:
  WriteDataElem(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int nodeflag;
//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void header();
  void type_arrays();
  void atoms();
  void nodes(tagint);
  void elements(tagint);
  void atom_velocities();
  void elem_velocities(tagint);
  void node_velocities(tagint);
};

}

#endif
#endif


