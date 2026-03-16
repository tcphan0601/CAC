#ifdef NPAIR_CLASS

NPairStyle(copy,
           NPairCopy,
           NP_COPY)

#else

#ifndef CAC_NPAIR_COPY_H
#define CAC_NPAIR_COPY_H

#include "npair.h"

namespace CAC_NS {

class NPairCopy : public NPair {
 public:
  NPairCopy(class CAC *);
  ~NPairCopy() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
