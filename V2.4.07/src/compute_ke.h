#ifdef COMPUTE_CLASS

ComputeStyle(ke,ComputeKE)

#else

#ifndef CAC_COMPUTE_KE_H
#define CAC_COMPUTE_KE_H

#include "compute.h"

namespace CAC_NS {

class ComputeKE : public Compute {
 public:
  ComputeKE(class CAC *, int, char **);
  void init();
  double compute_scalar();

 private:
  double pfactor;
  int pairflag,fixflag;
  int accurate_flag;
  int mass_style;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pe must use group all

*/
