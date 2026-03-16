#ifdef FIX_CLASS

FixStyle(thermal/conductivity,FixThermalConductivity)

#else

#ifndef CAC_FIX_THERMAL_CONDUCTIVITY_H
#define CAC_FIX_THERMAL_CONDUCTIVITY_H

#include "fix.h"

namespace CAC_NS {

class FixThermalConductivity : public Fix {
 public:
  FixThermalConductivity(class CAC *, int, char **);
  ~FixThermalConductivity();
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();

 private:
  int me;
  int edim,nbin,periodicity;
  int nswap;
  double prd,boxlo,boxhi;
  double slablo_lo,slablo_hi,slabhi_lo,slabhi_hi;
  double e_exchange;

  int nlo,nhi;
  int (*index_lo)[2],(*index_hi)[2];
  double *ke_lo,*ke_hi;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Fix thermal/conductivity swap value must be positive

Self-explanatory.

W: Fix thermal/conductivity comes before fix ave/spatial

The order of these 2 fixes in your input script is such that fix
thermal/conductivity comes first.  If you are using fix ave/spatial to
measure the temperature profile induced by fix viscosity, then this
may cause a glitch in the profile since you are averaging immediately
after swaps have occurred.  Flipping the order of the 2 fixes
typically helps.

*/
