#ifdef NPAIR_CLASS

NPairStyle(full/bin/intg/dbl,
           NPairFullBinIntgDouble,
           NP_FULL | NP_BIN | NP_NEWTON | NP_DOUBLE | NP_THREEBODY |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTG)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTG_DOUBLE_H
#define CAC_NPAIR_FULL_BIN_INTG_DOUBLE_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntgDouble : public NPair {
 public:
  NPairFullBinIntgDouble(class CAC *);
  ~NPairFullBinIntgDouble();
  void build(class NeighList *);
};

}

#endif
#endif
