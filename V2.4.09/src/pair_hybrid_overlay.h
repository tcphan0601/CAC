#ifdef PAIR_CLASS

PairStyle(hybrid/overlay,PairHybridOverlay)

#else

#ifndef CAC_PAIR_HYBRID_OVERLAY_H
#define CAC_PAIR_HYBRID_OVERLAY_H

#include "pair_hybrid.h"

namespace CAC_NS {

class PairHybridOverlay : public PairHybrid {
 public:
  PairHybridOverlay(class CAC *);
  virtual ~PairHybridOverlay() {}
  void coeff(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

*/
