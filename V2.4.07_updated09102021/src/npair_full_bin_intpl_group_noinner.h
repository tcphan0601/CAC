#ifdef NPAIR_CLASS

NPairStyle(full/bin/intpl/group/noinner,
           NPairFullBinIntplGroupNoInner,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTPL |
           NP_GROUP | NP_NOINNER)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTPL_GROUP_NOINNER_H
#define CAC_NPAIR_FULL_BIN_INTPL_GROUP_NOINNER_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntplGroupNoInner : public NPair {
 public:
  NPairFullBinIntplGroupNoInner(class CAC *);
  ~NPairFullBinIntplGroupNoInner() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
