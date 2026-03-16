#ifdef FIX_CLASS

FixStyle(temp/press,FixTempPress)

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
  double tvalue,pvalue,boltz,arandom,pvalue_c,adynamic;
  int varflag,iregion,seed,rand_flag,const_flag,dynamic_flag;
  int shear_direc,srand_flag;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  double pnorm[3];
  double net_force[7],net_force_all[7];
  int force_flag;
  
  int nfreq;
  double *ffreq,*vfreq;

  //int maxatom;
  //double **sforce;
  
  virtual void read_file(char*);
};

}

#endif
#endif

