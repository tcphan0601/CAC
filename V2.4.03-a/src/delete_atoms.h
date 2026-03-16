#ifdef COMMAND_CLASS

CommandStyle(delete_atoms,DeleteAtoms)

#else

#ifndef CAC_DELETE_ATOMS_H
#define CAC_DELETE_ATOMS_H

#include "pointers.h"
#include <map>

namespace CAC_NS {

class DeleteAtoms : protected Pointers {
 public:
  DeleteAtoms(class CAC *);
  void command(int, char **);

 private:
  int *del_flag_list;
  int allflag,compress_flag;
//  std::map<tagint,int> *hash;

  void list_delete_group(int, char **);
  void list_delete_region(int, char **);
//  void delete_overlap(int, char **);
//  void delete_porosity(int, char **);
  void options(int, char **);

//  inline int sbmask(int j) const {
//    return j >> SBBITS & 3;
//  }
};

}

#endif
#endif

