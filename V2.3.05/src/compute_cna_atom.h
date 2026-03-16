#ifdef COMPUTE_CLASS

ComputeStyle(cna/atom,ComputeCNAAtom)

#else

#ifndef CAC_COMPUTE_CNA_ATOM_H
#define CAC_COMPUTE_CNA_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputeCNAAtom : public Compute {
 public:
  ComputeCNAAtom(class CAC *, int, char **);
  ~ComputeCNAAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int namax; // max # of atoms locally
  int nemax; // max # of elements locally
  double cutsq;
  class NeighList *list;
  int **anearest;
  int **anearesttype;     // ATOM = atom neighbor
                          // INTPL = interpolated atom neighbor
  int ***ianearest;
  int ***ianearesttype;   // ATOM = atom neighbor
                          // INTPL = interpolated atom neighbor
  bigint *annearest;
  int **iannearest;
  double *apattern;
  double **iapattern;

  void grow_atom(int);
  void grow_intpl(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Compute cna/atom requires a pair style be defined

Self-explanatory.

E: Compute cna/atom cutoff is longer than pairwise cutoff

Self-explanatory.

W: Compute cna/atom cutoff may be too large to find ghost atom neighbors

The neighbor cutoff used may not encompass enough ghost atoms
to perform this operation correctly.

W: More than one compute cna/atom defined

It is not efficient to use compute cna/atom  more than once.

W: Too many neighbors in CNA for %d atoms

More than the maximum # of neighbors was found multiple times.  This
was unexpected.

W: Too many common neighbors in CNA %d times

More than the maximum # of neighbors was found multiple times.  This
was unexpected.

*/
