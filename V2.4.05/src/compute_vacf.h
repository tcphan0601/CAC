#ifdef COMPUTE_CLASS

ComputeStyle(vacf,ComputeVACF)

#else

#ifndef CAC_COMPUTE_VACF_H
#define CAC_COMPUTE_VACF_H

#include "compute.h"

namespace CAC_NS {

class ComputeVACF : public Compute {
 public:
  ComputeVACF(class CAC *, int, char **);
  ~ComputeVACF();
  void init();
  virtual void compute_vector();
  void set_atom_arrays(int);
  void set_elem_arrays(int);

 protected:
  bigint nvacf;
  int interpolate_flag;
  char *id_fix;
  class FixStore *fix;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute vacf fix ID

Self-explanatory.

*/
