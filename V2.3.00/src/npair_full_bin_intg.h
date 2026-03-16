#ifdef NPAIR_CLASS

NPairStyle(full/bin/intg,
           NPairFullBinIntg,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ATOM | NP_ELEM | 
           NP_INTG)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTG_H
#define CAC_NPAIR_FULL_BIN_INTG_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntg : public NPair {
 public:
  NPairFullBinIntg(class CAC *);
  ~NPairFullBinIntg() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
