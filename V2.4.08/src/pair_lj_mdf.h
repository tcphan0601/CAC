#ifdef PAIR_CLASS

PairStyle(lj/mdf, PairLJMDF)

#else

#ifndef CAC_PAIR_LJ_MDF_H
#define CAC_PAIR_LJ_MDF_H

#include "pair.h"

namespace CAC_NS {

class PairLJMDF : public Pair {

public:
  PairLJMDF(class CAC *);
  virtual ~PairLJMDF();
  virtual void compute(int, int);
  virtual void pair_a2a(int, int, double *);
  virtual void pair_a2ia(int, int, int, double *, double *);
  virtual void pair_ia2ia(int, int, double *, int, int, double *, double *);
  void settings(int, char **);
  void coeff(int, char **);
  void write_data_all(FILE *);
  double init_one(int, int);
  // void init_style();

protected:
  double cut_global, cut_inner_global;
  double **cut, **cut_inner, **cut_inner_sq;
  double **aparm, **bparm;
  double **lj1, **lj2, **lj3, **lj4;

  double nemax;
  virtual void allocate();
};

} // namespace CAC_NS

#endif

#endif
