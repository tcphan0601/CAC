#ifdef COMPUTE_CLASS

ComputeStyle(stress/mech, ComputeStressMech)

#else

#ifndef CAC_COMPUTE_STRESS_MECH_H
#define CAC_COMPUTE_STRESS_MECH_H

#include "compute.h"

namespace CAC_NS {

class ComputeStressMech : public Compute {
 public:
  ComputeStressMech(class CAC *, int, char **);
  ~ComputeStressMech();
  void init();
  void compute_peratom();
  int pack_atom_reverse_comm(int, int, double *);
  void unpack_atom_reverse_comm(int, int *, double *);
  int pack_elem_reverse_comm(int, int, double *);
  void unpack_elem_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  double a;     // volume stress element edge 
  int style;
  int atom_nmax, elem_nmax;
  int maxnpe;
  double **atom_stress, ****node_stress;
};

}

#endif
#endif
