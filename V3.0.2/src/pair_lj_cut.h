#ifdef PAIR_CLASS

PairStyle(lj/cut, PairLJCut)

#else

#ifndef CAC_PAIR_LJ_CUT_H
#define CAC_PAIR_LJ_CUT_H

#include "pair.h"

namespace CAC_NS {

class PairLJCut : public Pair {

 public:

  PairLJCut(class CAC * );
  virtual ~PairLJCut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_data_all(FILE *);
  double init_one(int, int);
  void init_style();

 protected:

  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;

  virtual void allocate();

};

}

#endif

#endif
