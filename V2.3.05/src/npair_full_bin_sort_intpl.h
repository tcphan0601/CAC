#ifdef NPAIR_CLASS

NPairStyle(full/bin/sort/intpl,
           NPairFullBinSortIntpl,
           NP_FULL | NP_BIN | NP_NEWTON | NP_SORT |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTPL)

#else

#ifndef CAC_NPAIR_FULL_BIN_SORT_INTPL_H
#define CAC_NPAIR_FULL_BIN_SORT_INTPL_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinSortIntpl : public NPair {
 public:
  NPairFullBinSortIntpl(class CAC *);
  ~NPairFullBinSortIntpl() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
