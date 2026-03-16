#ifdef FIX_CLASS

FixStyle(temp/press, FixTempPress)

#else

#ifndef CAC_FIX_TEMP_PRESS_H
#define CAC_FIX_TEMP_PRESS_H

#include "fix.h"

namespace CAC_NS {

class FixTempPress : public Fix {
 public:
  FixTempPress(class CAC *, int, char **);
  ~FixTempPress();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double press, boltz, mean_press, press_fluct_amplitude;
  double t_start, t_stop, t_period, t_target;
  double ucell_volume_user;
  double ucell_volume_style;
  class RanMars *random;
  int seed;
  int zero_mean_flag;
  int style;
  int direction_flag[3];
  char *tstr;
  char *idregion;
  int tvar, tstyle;
  double pnorm[3];
  double net_force[7], net_force_all[7];
  int force_flag;
  
  int nfreq;
  double *ffreq, *vfreq;

  //int maxatom;
  //double **sforce;
  
  virtual void read_file(char*);
  void compute_target();
};

}

#endif
#endif

