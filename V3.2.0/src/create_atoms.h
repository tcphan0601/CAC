#ifdef COMMAND_CLASS

CommandStyle(create_atoms, CreateAtoms)

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
  int me, nprocs;
  int ntype, style, iregion, nbasis, nrandom, seed;
  int *basistype;
  double xone[3], quatone[4];
  int remapflag;
  int wholecellflag;
  int subsetflag;
  int gtag;
  bigint nsubset;
  double subsetfrac;

  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  char *xstr_copy, *ystr_copy, *zstr_copy;

  int ilo, ihi, jlo, jhi, klo, khi;

  int nlatt;                  // number of owned lattice sites
  int nlatt_overflow;         // 1 if local nlatt exceeds a 32-bit int

  int *flag;                  // flag subset of particles to insert on lattice
  int *next;

  class RanMars *ranlatt;

  int triclinic;
  double sublo[3], subhi[3];   // epsilon-extended proc sub-box for adding atoms

  void add_single();
  void add_random();
  void add_lattice();
  void loop_lattice(int);
  int vartest(double *);        // evaluate a variable with new atom position
};

}

#endif
#endif

