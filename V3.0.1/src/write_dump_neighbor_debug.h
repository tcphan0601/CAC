#ifdef COMMAND_CLASS

CommandStyle(write_dump_neighbor_debug,WriteDumpNeighborDebug)

#else

#ifndef CAC_WRITE_DUMP_NEIGHBOR_DEBUG_H
#define CAC_WRITE_DUMP_NEIGHBOR_DEBUG_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class WriteDumpNeighborDebug : protected Pointers {
 public:
  WriteDumpNeighborDebug(class CAC *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int id,ielem,igcell,ibasis,allflag;
  int ii;
  int iatom;
  int singleflag;
  FILE *fp;

  class NeighList *list;        // copy pointer from pair neighbor list
  class NeighList *listhalf;        // copy pointer from pair neighbor list
  class NeighList *listfull;        // copy pointer from pair neighbor list

  void header();
};

}

#endif
#endif


