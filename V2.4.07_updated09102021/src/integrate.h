#ifndef CAC_INTEGRATE_H
#define CAC_INTEGRATE_H

#include "pointers.h"

namespace CAC_NS {

class Integrate : protected Pointers {
 public:
  Integrate(class CAC *, int, char **);
  virtual ~Integrate();
  virtual void init();
  virtual void setup(int flag=1) = 0;
//  virtual void setup_minimal(int) = 0;
  virtual void run(int) = 0;
  virtual void cleanup() {}
  virtual void reset_dt() {}
  virtual bigint memory_usage() {return 0;}

 protected:
  int eflag,vflag;                  // flags for energy/virial computation
  int virial_style;                 // compute virial explicitly or implicitly
  int external_force_clear;         // clear forces locally or externally

  int nelist_global,nelist_atom;    // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_global;     // lists of PE,virial Computes
  class Compute **elist_atom;
  class Compute **vlist_global;
  class Compute **vlist_atom;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

  void ev_setup();
  void ev_set(bigint);
};

}

#endif

/* ERROR/WARNING messages:

*/
