#ifdef NBIN_CLASS

NBinStyle(standard,
          NBinStandard,
          0)

#else

#ifndef CAC_NBIN_STANDARD_H
#define CAC_NBIN_STANDARD_H

#include "nbin.h"

namespace CAC_NS {

class NBinStandard : public NBin {
 public:
  NBinStandard(class CAC *);
  ~NBinStandard() {}
  void setup_bins(int);
  void bin_all();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
