#ifdef NPAIR_CLASS

NPairStyle(full/bin/atom, 
           NPairFullBinAtom, 
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM)

#else

#ifndef CAC_NPAIR_FULL_BIN_ATOM_H
#define CAC_NPAIR_FULL_BIN_ATOM_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinAtom : public NPair {
 public:
  NPairFullBinAtom(class CAC *);
  ~NPairFullBinAtom() {}
  void build(class NeighList *);
};

}

#endif
#endif
