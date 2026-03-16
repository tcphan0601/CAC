#ifdef NPAIR_CLASS

NPairStyle(full/bin/intpl,
           NPairFullBinIntpl,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ATOM | NP_ELEM | 
           NP_INTPL)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTPL_H
#define CAC_NPAIR_FULL_BIN_INTPL_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntpl : public NPair {
 public:
  NPairFullBinIntpl(class CAC *);
  ~NPairFullBinIntpl() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
