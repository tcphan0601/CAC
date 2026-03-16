#ifdef COMPUTE_CLASS

ComputeStyle(temp, ComputeTemp)

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
  int accurate_flag;
  int mass_style;

  virtual void dof_compute();
};

}

#endif
#endif

