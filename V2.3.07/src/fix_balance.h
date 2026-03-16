#ifdef FIX_CLASS

FixStyle(balance,FixBalance)

#else

#ifndef CAC_FIX_BALANCE_H
#define CAC_FIX_BALANCE_H

#include <stdio.h>
#include "fix.h"

namespace CAC_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class CAC *, int, char **);
  ~FixBalance();
  int setmask();
  void post_constructor();
  void init();
  void setup(int);
  void setup_pre_exchange();
  void pre_exchange();
  void pre_neighbor();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int nevery,lbstyle,nitermax;
  int screenflag;
  double thresh,stopthresh;
  char bstr[4];

  double imbnow;                // current imbalance factor
  double imbprev;               // imbalance factor before last rebalancing
  double imbfinal;              // imbalance factor after last rebalancing
  double maxloadperproc;        // max load on any processor
  double minloadperproc;        // min load on any processor
  int itercount;                // iteration count of last call to Balance
  int pending;
  bigint lastbalance;           // last timestep balancing was attempted

  class Balance *balance;
  class IrregularComm *irregular_comm;

  void rebalance();
};

}

#endif
#endif
