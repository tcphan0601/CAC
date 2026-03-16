#ifdef COMPUTE_CLASS

ComputeStyle(stress/atom,ComputeStressAtom)

#else

#ifndef CAC_COMPUTE_STRESS_ATOM_H
#define CAC_COMPUTE_STRESS_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputeStressAtom : public Compute {
 public:
  ComputeStressAtom(class CAC *, int, char **);
  ~ComputeStressAtom();
  void init();
  void compute_peratom();
  int pack_atom_reverse_comm(int, int, double *);
  void unpack_atom_reverse_comm(int, int *, double *);
  int pack_elem_reverse_comm(int, int, double *);
  void unpack_elem_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int keflag,pairflag;
  int fixflag,biasflag;
  char *id_temp;

  int atom_nmax,elem_nmax;
  int max_npe;
  double **atom_stress,***node_stress;
};

}

#endif
#endif
