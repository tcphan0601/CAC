#ifdef COMMAND_CLASS

CommandStyle(create_atoms,CreateAtoms)

#else

#ifndef CAC_CREATE_ATOMS_H
#define CAC_CREATE_ATOMS_H

#include "pointers.h"

namespace CAC_NS {

class CreateAtoms : protected Pointers {
 public:
  CreateAtoms(class CAC *);
  void command(int, char **);

 private:
  int ntype,style,nregion,nbasis,nrandom,seed;
  int *basistype;
  double xone[3],quatone[4];
  int remapflag;

  int varflag,vvar,xvar,yvar,zvar;
  char *vstr,*xstr,*ystr,*zstr;
  char *xstr_copy,*ystr_copy,*zstr_copy;

  int triclinic;
  double sublo[3],subhi[3];   // epsilon-extended proc sub-box for adding atoms

  void add_single();
  void add_random();
  void add_lattice();
  int vartest(double *);        // evaluate a variable with new atom position
};

}

#endif
#endif

