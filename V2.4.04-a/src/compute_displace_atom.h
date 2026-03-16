#ifdef COMPUTE_CLASS

ComputeStyle(displace/atom,ComputeDisplaceAtom)

#else

#ifndef CAC_COMPUTE_DISPLACE_ATOM_H
#define CAC_COMPUTE_DISPLACE_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputeDisplaceAtom : public Compute {
 public:
  ComputeDisplaceAtom(class CAC *, int, char **);
  ~ComputeDisplaceAtom();
  void init();
  void compute_peratom();
  void set_atom_arrays(int);
  void set_elem_arrays(int);
  void pass_on_atom_arrays(int, int, int);
  void pass_on_elem_arrays(int, int, int*);

  //void refresh();
  double memory_usage();

 private:
  int namax,nemax;
  double **atomdisplace;
  double ***nodedisplace;
  char *id_fix;
  class FixStore *fix;

  //int refreshflag,ivar,nvmax;    // refresh option is enabled
  //char *rvar;                    // for incremental dumps
  //double *varatom;
};

}

#endif
#endif
