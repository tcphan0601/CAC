#ifdef FIX_CLASS

FixStyle(temp/rescale, FixTempRescale)

#else

#ifndef CAC_FIX_TEMP_RESCALE_H
#define CAC_FIX_TEMP_RESCALE_H

#include "fix.h"

namespace CAC_NS {

class FixTempRescale : public Fix {
 public:
  FixTempRescale(class CAC *, int, char **);
  virtual ~FixTempRescale();
  int setmask();
  void init();
  virtual void end_of_step();
  int modify_param(int, char **);
  void reset_target(double);
  double compute_scalar();
  virtual void *extract(const char *, int &);

 protected:
  int which;
  double t_start, t_stop, t_window, t_target;
  double fraction, energy, efactor;
  int tstyle, tvar;
  char *tstr;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}

#endif
#endif

/*  ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running CAC to see the offending line.

E: Variable name for fix temp/rescale does not exist

Self-explanatory.

E: Variable for fix temp/rescale is invalid style

Only equal-style variables can be used.

E: Temperature ID for fix temp/rescale does not exist

Self-explanatory.

E: Computed temperature for fix temp/rescale cannot be 0.0

Cannot rescale the temperature to a new value if the current
temperature is 0.0.

E: Fix temp/rescale variable returned negative temperature

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Group for fix_modify temp != fix group

The fix_modify command is specifying a temperature computation that
computes a temperature on a different group of atoms than the fix
itself operates on.  This is probably not what you want to do.

 */
