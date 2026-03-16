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
  int nnodes_local;
  int velocityflag;
  int elementflag;
  int atomflag;

//  int pairflag;
//  int coeffflag;
  FILE *fp;

  void header();
  void type_arrays();
  void atoms();
  void elements();
  void nodes(int);
  void atom_velocities();
  void node_velocities(int);
  void element_clusters();
};

}

#endif
#endif


