#ifdef PAIR_CLASS

PairStyle(eam/fs, PairEAMFS)

#else

#ifndef CAC_PAIR_EAM_FS_H
#define CAC_PAIR_EAM_FS_H

#include "pair_eam.h"

namespace CAC_NS {

// need virtual public b/c of how eam/fs/opt inherits from it

class PairEAMFS : virtual public PairEAM {
 public:
  PairEAMFS(class CAC *);
  virtual ~PairEAMFS() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif

