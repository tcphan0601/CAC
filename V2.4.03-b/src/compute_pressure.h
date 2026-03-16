#ifdef COMPUTE_CLASS

ComputeStyle(pressure,ComputePressure)

#else

#ifndef CAC_COMPUTE_PRESSURE_H
#define CAC_COMPUTE_PRESSURE_H

#include "compute.h"

namespace CAC_NS {

class ComputePressure : public Compute {
 public:
  ComputePressure(class CAC *, int, char **);
  virtual ~ComputePressure();
  virtual void init();
  virtual double compute_scalar();
  virtual void compute_vector();
  void reset_extra_compute_fix(const char *);

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double **vptr;
  //Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag,pairflag;
  int fixflag;

  void virial_compute(int, int);
};

}

#endif
#endif
