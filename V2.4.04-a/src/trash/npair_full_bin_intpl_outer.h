#ifdef NPAIR_CLASS

NPairStyle(full/bin/intpl/outer,
           NPairFullBinIntplOuter,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTPL_OUTER)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTPL_OUTER_H
#define CAC_NPAIR_FULL_BIN_INTPL_OUTER_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntplOuter : public NPair {
 public:
  NPairFullBinIntplOuter(class CAC *);
  ~NPairFullBinIntplOuter() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
