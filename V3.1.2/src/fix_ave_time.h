#ifdef FIX_CLASS

FixStyle(ave/time,FixAveTime)

#else

#ifndef CAC_FIX_AVE_TIME_H
#define CAC_FIX_AVE_TIME_H

#include "fix.h"

namespace CAC_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class CAC *, int, char **);
  ~FixAveTime();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid,nvalid_last;
  int *which,*argindex,*value2index,*offcol;
  int *varlen;               // 1 if value is from variable-length compute
  char **ids;
  FILE *fp;
  int nrows;
  int any_variable_length;
  int all_variable_length;

  int ave,nwindow,startstep,mode;
  int noff,overwrite;
  int *offlist;
  char *format,*format_user;
  char *title1,*title2,*title3;
  long filepos;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;
  double *column;
  double **array;
  double **array_total;
  double ***array_list;

  int column_length(int);
  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, int, char **);
  void allocate_arrays();
  bigint nextvalid();
};

}

#endif
#endif
