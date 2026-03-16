#ifdef COMMAND_CLASS

CommandStyle(write_data_adrian, WriteDataAdrian)

#else

#ifndef CAC_WRITE_DATA_ADRIAN_H
#define CAC_WRITE_DATA_ADRIAN_H

#include "pointers.h"
#include <stdio.h>

namespace CAC_NS {

class WriteDataAdrian : protected Pointers {
public:
  WriteDataAdrian(class CAC *);
  void command(int, char **);
  void write(char *);

private:
  int me, nprocs;
  //  int pairflag;
  //  int coeffflag;
  FILE *fp;

  void header();
  void atoms();
  void elements();
  // void atom_velocities();
  // void node_velocities();
};

} // namespace CAC_NS

#endif
#endif
