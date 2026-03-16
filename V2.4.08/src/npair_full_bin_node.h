#ifdef NPAIR_CLASS

NPairStyle(full/bin/node,
           NPairFullBinNode,
           NP_FULL | NP_BIN | NP_NEWTON | 
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOM | NP_ELEM | NP_NODE)

#else

#ifndef CAC_NPAIR_FULL_BIN_NODE_H
#define CAC_NPAIR_FULL_BIN_NODE_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinNode : public NPair {
 public:
  NPairFullBinNode(class CAC *);
  ~NPairFullBinNode() {}
  void build(class NeighList *);
};

}

#endif
#endif
