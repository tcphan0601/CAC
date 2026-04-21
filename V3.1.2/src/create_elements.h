#ifdef COMMAND_CLASS

CommandStyle(create_elements, CreateElements)

#else

#ifndef CAC_CREATE_ELEMENTS_H
#define CAC_CREATE_ELEMENTS_H

#include "pointers.h"

namespace CAC_NS {

class CreateElements : protected Pointers {
 public:
  CreateElements(class CAC *);
  void command(int, char **);
  int **ctypelist;
  double **xtemplist;

 private:
  int netype, nctype, style, iregion, nbasis, nrandom, seed;
  int *basisctype;
  double xone[3], quatone[4];
  int remapflag;
  int offsetplane, offset1, offset2, offsetdir1, offsetdir2;

  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  char *xstr_copy, *ystr_copy, *zstr_copy;

  int triclinic;
  double sublo[3], subhi[3];   // epsilon-extended proc sub-box for adding atoms

  void add_single();
  //void add_random();
  void add_lattice();
  int vartest(double *);        // evaluate a variable with new atom position
};

}

#endif
#endif

