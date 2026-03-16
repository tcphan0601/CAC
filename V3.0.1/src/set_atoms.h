#ifdef COMMAND_CLASS

CommandStyle(set_atoms, SetAtoms)

#else

#ifndef CAC_SET_ATOMS_H
#define CAC_SET_ATOMS_H

#include "pointers.h"

namespace CAC_NS {

class SetAtoms : protected Pointers {
 public:
  SetAtoms(class CAC *cac) : Pointers(cac) {};
  void command(int, char **);

 private:
  char *id;
  int *select;
  int style, ivalue, newtype, count;
  int ximage, yimage, zimage, ximageflag, yimageflag, zimageflag;
  bigint nsubset;
  double dvalue, xvalue, yvalue, zvalue, fraction;

  int varflag, varflag1, varflag2, varflag3;
  int ivar1, ivar2, ivar3;
  double *vec1, *vec2, *vec3;


  void selection(int);
  void set(int);
  void setrandom(int);
  void varparse(char *, int);
};

}

#endif
#endif

