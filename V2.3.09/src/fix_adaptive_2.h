#ifdef FIX_CLASS

FixStyle(adaptive2,FixAdaptive2)

#else

#ifndef CAC_FIX_ADAPTIVE_2_H
#define CAC_FIX_ADAPTIVE_2_H

#include "fix.h"

namespace CAC_NS {

class FixAdaptive2 : public Fix {
 public:
  FixAdaptive2(class CAC *, int, char **);
  ~FixAdaptive2();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_exchange();
  void setup_pre_exchange();
  void post_integrate();
  //double compute_vector(int);

  double memory_usage();

 private:

  void grow_splitlist();
  int **splitlist;
  int maxsplit;

  class NeighList *list;

  int pending;
  int perfect_crystal;                        // compute value for perfect crystal

  char *id_compute;
  int icompute;
  class Compute *compute;



};

}

#endif
#endif
