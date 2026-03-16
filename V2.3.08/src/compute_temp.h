#ifdef COMPUTE_CLASS

ComputeStyle(temp,ComputeTemp)

#else

#ifndef CAC_COMPUTE_TEMP_H
#define CAC_COMPUTE_TEMP_H

#include "compute.h"

namespace CAC_NS {

class ComputeTemp : public Compute {
 public:
  ComputeTemp(class CAC *, int, char **);
  virtual ~ComputeTemp();
  void init();
  virtual double compute_scalar();
  virtual void compute_vector();

 protected:
  double tfactor;
  int    sumflag;

  virtual void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
