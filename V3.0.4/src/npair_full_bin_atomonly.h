#ifdef NPAIR_CLASS

NPairStyle(full/bin/atomonly, 
           NPairFullBinAtomonly, 
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOMONLY)

#else

#ifndef CAC_NPAIR_FULL_BIN_ATOMONLY_H
#define CAC_NPAIR_FULL_BIN_ATOMONLY_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinAtomonly : public NPair {
 public:
  NPairFullBinAtomonly(class CAC *);
  ~NPairFullBinAtomonly() {}
  void build(class NeighList *);
};

}

#endif
#endif
