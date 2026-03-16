#ifdef PAIR_CLASS

PairStyle(gauss/cut, PairGaussCut)

#else

#ifndef CAC_PAIR_GAUSS_CUT_H
#define CAC_PAIR_GAUSS_CUT_H

#include "pair.h"

namespace CAC_NS {

class PairGaussCut : public Pair {

 public:
  PairGaussCut(class CAC *);
  ~PairGaussCut();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void write_data(FILE *);
  virtual void write_data_all(FILE *);
  virtual double init_one(int, int);

  virtual double memory_usage();
 protected:
  double cut_global;
  double **cut;
  double **hgauss, **sigmah, **rmh;
  double **pgauss, **offset;

  void allocate();
};

} 

#endif

#endif
