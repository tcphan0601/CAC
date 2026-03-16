#ifdef COMPUTE_CLASS

ComputeStyle(pe,ComputePE)

#else

#ifndef CAC_COMPUTE_PE_H
#define CAC_COMPUTE_PE_H

#include "compute.h"

namespace CAC_NS {

class ComputePE : public Compute {
 public:
  ComputePE(class CAC *, int, char **);
  ~ComputePE() {}
  void init() {}
  double compute_scalar();

 private:
  int pairflag,fixflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Compute pe must use group all

Energies computed by potentials (pair, bond, etc) are computed on all
atoms.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
