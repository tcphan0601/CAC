#ifdef NPAIR_CLASS

NPairStyle(full/bin/vatom, 
           NPairFullBinVatom, 
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_VATOM | NP_SORT |
           NP_GROUP | NP_NOINNER)

#else

#ifndef CAC_NPAIR_FULL_BIN_VATOM_H
#define CAC_NPAIR_FULL_BIN_VATOM_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinVatom : public NPair {
 public:
  NPairFullBinVatom(class CAC *);
  ~NPairFullBinVatom();
  void build(class NeighList *);
};

}

#endif
#endif

/*  ERROR/WARNING messages:

 */
