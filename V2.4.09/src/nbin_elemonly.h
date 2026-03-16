#ifdef NBIN_CLASS

NBinStyle(elemonly,
          NBinElemonly,
          NB_STANDARD | NB_ELEMONLY)

#else

#ifndef CAC_NBIN_ELEMONLY_H
#define CAC_NBIN_ELEMONLY_H

#include "nbin.h"

namespace CAC_NS {

class NBinElemonly : public NBin {
 public:
  NBinElemonly(class CAC *);
  ~NBinElemonly() {}
  void setup_bins(int);
  void bin_all();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
