#ifdef NPAIR_CLASS

NPairStyle(full/bin/elem,
           NPairFullBinElem,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ELEM | NP_GROUP)

#else

#ifndef CAC_NPAIR_FULL_BIN_ELEM_H
#define CAC_NPAIR_FULL_BIN_ELEM_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinElem : public NPair {
 public:
  NPairFullBinElem(class CAC *);
  ~NPairFullBinElem() {}
  void build(class NeighList *);
};

}

#endif
#endif
