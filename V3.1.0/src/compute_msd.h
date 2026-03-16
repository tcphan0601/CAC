#ifdef COMPUTE_CLASS

ComputeStyle(msd, ComputeMSD)

#else

#ifndef CAC_COMPUTE_MSD_H
#define CAC_COMPUTE_MSD_H

#include "compute.h"

namespace CAC_NS {

class ComputeMSD : public Compute {
 public:
  ComputeMSD(class CAC *, int, char **);
  virtual ~ComputeMSD();
  void init();
  virtual void compute_vector();
  void set_atom_arrays(int);
  void set_elem_arrays(int);

 protected:
  int comflag;   // comflag = 1 if reference moves with center of mass
  int avflag;    // avflag = 1 if using average position as reference
  int naverage;  // number of samples for average position
  bigint nmsd;
  double masstotal;
  char *id_fix;
  class FixStore *fix;
};

}

#endif
#endif

/*  ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Could not find compute msd fix ID

Self-explanatory.

 */
