#ifdef COMPUTE_CLASS

ComputeStyle(pe/atom,ComputePEAtom)

#else

#ifndef CAC_COMPUTE_PE_ATOM_H
#define CAC_COMPUTE_PE_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputePEAtom : public Compute {
 public:
  ComputePEAtom(class CAC *, int, char **);
  ~ComputePEAtom();
  void init() {}
  void compute_peratom();
  int pack_atom_reverse_comm(int, int, double *);
  void unpack_atom_reverse_comm(int, int *, double *);
  int pack_elem_reverse_comm(int, int, double *);
  void unpack_elem_reverse_comm(int, int *, double *);

  double memory_usage();

 private:
  int pairflag,fixflag;//,bondflag,angleflag,dihedralflag,improperflag;
  //int kspaceflag,fixflag;
  int atom_nmax;
  int elem_nmax;
  int max_npe;
  double *atom_energy;
  double **node_energy;
};

}

#endif
#endif

