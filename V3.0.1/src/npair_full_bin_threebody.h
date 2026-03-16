#ifdef NPAIR_CLASS

NPairStyle(full/bin/threebody, 
           NPairFullBinThreebody, 
           NP_FULL | NP_BIN | NP_NEWTON | NP_THREEBODY |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_GAUSS)

#else

#ifndef CAC_NPAIR_FULL_BIN_THREEBODY_H
#define CAC_NPAIR_FULL_BIN_THREEBODY_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinThreebody : public NPair {
 public:
  NPairFullBinThreebody(class CAC *);
  ~NPairFullBinThreebody() {}
  void build(class NeighList *);
};

}

#endif
#endif
