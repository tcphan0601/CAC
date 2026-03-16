#ifdef NPAIR_CLASS

NPairStyle(full/bin/intg/threebody,
           NPairFullBinIntgThreebody,
           NP_FULL | NP_BIN | NP_NEWTON | NP_THREEBODY |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_INTG)

#else

#ifndef CAC_NPAIR_FULL_BIN_INTG_THREEBODY_H
#define CAC_NPAIR_FULL_BIN_INTG_THREEBODY_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinIntgThreebody : public NPair {
 public:
  NPairFullBinIntgThreebody(class CAC *);
  ~NPairFullBinIntgThreebody() {}
  void build(class NeighList *);
};

}

#endif
#endif
