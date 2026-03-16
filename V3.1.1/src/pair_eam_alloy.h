#ifdef PAIR_CLASS

PairStyle(eam/alloy, PairEAMAlloy)

#else

#ifndef CAC_PAIR_EAM_ALLOY_H
#define CAC_PAIR_EAM_ALLOY_H

#include "pair_eam.h"

namespace CAC_NS {

// need virtual public b/c of how eam/alloy/opt inherits from it

class PairEAMAlloy : virtual public PairEAM {
 public:
  PairEAMAlloy(class CAC *);
  virtual ~PairEAMAlloy() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif

/*  ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: No matching element in EAM potential file

The EAM potential file does not contain elements that match the
requested elements.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect element names in EAM potential file

The element names in the EAM file do not match those requested.

 */
