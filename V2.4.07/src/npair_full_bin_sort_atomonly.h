#ifdef NPAIR_CLASS

NPairStyle(full/bin/sort/atomonly,
           NPairFullBinSortAtomonly,
           NP_FULL | NP_BIN | NP_NEWTON | NP_SORT |
           NP_NEWTOFF | NP_ORTHO | NP_TRI |
           NP_ATOMONLY)

#else

#ifndef CAC_NPAIR_FULL_BIN_SORT_ATOMONLY_H
#define CAC_NPAIR_FULL_BIN_SORT_ATOMONLY_H

#include "npair.h"

namespace CAC_NS {

class NPairFullBinSortAtomonly : public NPair {
 public:
  NPairFullBinSortAtomonly(class CAC *);
  ~NPairFullBinSortAtomonly(); 
  void build(class NeighList *);
 private:
  int *ibucket;
  double *rbucket;
  int maxbucket;

};

}

#endif
#endif
