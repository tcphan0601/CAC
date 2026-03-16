#ifdef COMMAND_CLASS

CommandStyle(balance,Balance)

#else

#ifndef CAC_BALANCE_H
#define CAC_BALANCE_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class Balance : protected Pointers {
 public:
  class RCB *rcb;
  int outflag;                    // 1 for output of balance results to file

  Balance(class CAC *);
  ~Balance();
  void command(int, char **);
  int options(int, int, char **);
  double imbalance_factor(double &, double &);
  void shift_setup(char *, int, double);
  int shift();
  int *bisection(int sortflag = 0);
  void dumpout(bigint);

 private:
  int me,nprocs;

  double thresh;                                    // threshold to perform LB
  int style;                                        // style of LB
  int xflag,yflag,zflag;                            // xyz LB flags
  double *user_xsplit,*user_ysplit,*user_zsplit;    // params for xyz LB

  int nitermax;              // params for shift LB
  double stopthresh;
  char bstr[4];

  int shift_allocate;        // 1 if SHIFT vectors have been allocated
  int ndim;                  // length of balance string bstr
  int *bdim;                 // XYZ for each character in bstr
  double *onecost;           // work vector of counts in one dim
  double *allcost;           // counts for slices in one dim
  double *sum;               // cumulative count for slices in one dim
  double *target;            // target sum for slices in one dim
  double *lo,*hi;            // lo/hi split coords that bound each target
  double *losum,*hisum;      // cumulative counts at lo/hi coords
  int rho;                   // 0 for geometric recursion
                             // 1 for density weighted recursion

  double *proccost;          // particle cost per processor
  double *allproccost;       // proccost summed across procs

 
  double user_element_weight;

  FILE *fp;                  // balance output file
  int firststep;

  double imbalance_splits();
  void shift_setup_static(char *);
  void tally(int, int, double *);
  int adjust(int, double *);
  int binary_search(double, int, double *);
#ifdef BALANCE_DEBUG
  void debug_shift_output(int, int, int, double *);
#endif
};

}

#endif
#endif

