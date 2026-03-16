#ifdef PAIR_CLASS

PairStyle(lj/gauss/cut, PairLJGaussCut)

#else

#ifndef CAC_PAIR_LJ_GAUSS_CUT_H
#define CAC_PAIR_LJ_GAUSS_CUT_H

#include "pair.h"

namespace CAC_NS {

class PairLJGaussCut : public Pair {

 public:
  PairLJGaussCut(class CAC *);
  virtual ~PairLJGaussCut();
  virtual void compute(int, int);
  virtual void pair_a2a(int, int, double *);
  virtual void pair_a2ia(int, int, int, double *, double *);
  virtual void pair_ia2ia(int, int, double *, int, int, double *, double *);
  void settings(int, char **);
  void coeff(int, char **);
  void write_data_all(FILE *);
  double init_one(int, int);
  void init_style();
  double memory_usage();

 protected:
  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;

  double h,sigmah,rmh;


  double pgauss;
  double nemax;
  virtual void allocate();
};

} // namespace CAC_NS

#endif

#endif
