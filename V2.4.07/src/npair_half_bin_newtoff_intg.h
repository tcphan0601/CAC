#ifdef NPAIR_CLASS

NPairStyle(half/bin/newtoff/intg,
           NPairHalfBinNewtoffIntg,
           NP_HALF | NP_BIN | NP_NEWTOFF | NP_ORTHO |
           NP_TRI | NP_ATOM | NP_ELEM | NP_INTG)

#else

#ifndef CAC_NPAIR_HALF_BIN_NEWTOFF_INTG_H
#define CAC_NPAIR_HALF_BIN_NEWTOFF_INTG_H

#include "npair.h"

namespace CAC_NS {

class NPairHalfBinNewtoffIntg : public NPair {
 public:
  NPairHalfBinNewtoffIntg(class CAC *);
  ~NPairHalfBinNewtoffIntg() {}
  void build(class NeighList *);
};

}

#endif
#endif
