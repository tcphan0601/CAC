#ifdef COMPUTE_CLASS

ComputeStyle(centroid/stress/atom, ComputeCentroidStressAtom)

#else

#ifndef CAC_COMPUTE_CENTROID_STRESS_ATOM_H
#define CAC_COMPUTE_CENTROID_STRESS_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputeCentroidStressAtom : public Compute {
 public:
  ComputeCentroidStressAtom(class CAC *, int, char **);
  ~ComputeCentroidStressAtom();
  void init();
  void compute_peratom();
  int pack_atom_reverse_comm(int, int, double *);
  void unpack_atom_reverse_comm(int, int *, double *);
  int pack_elem_reverse_comm(int, int, double *);
  void unpack_elem_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int keflag, pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int kspaceflag, fixflag, biasflag;
  //Compute *temperature;
  char *id_temp;

  int atom_nmax, elem_nmax;
  int maxnpe, maxapc;
  double **atom_stress, ****node_stress;
};

}

#endif
#endif

