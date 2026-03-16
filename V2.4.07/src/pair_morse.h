#ifdef PAIR_CLASS

PairStyle(morse, PairMorse)

#else

#ifndef CAC_PAIR_MORSE_H
#define CAC_PAIR_MORSE_H

#include "pair.h"

namespace CAC_NS {

class PairMorse : public Pair {

public:
  PairMorse(class CAC *);
  virtual ~PairMorse();
  virtual void compute(int, int);
  virtual void pair_a2a(int, int, double *);
  virtual void pair_a2ia(int, int, int, double *, double *);
  virtual void pair_ia2ia(int, int, double *, int, int, double *, double *);
  void settings(int, char **);
  void coeff(int, char **);
  void write_data_all(FILE *);
  double init_one(int, int);
  void init_style();

protected:
  // double alpha;
  double cut_global;
  double **cut;
  double **d0, **alpha, **r0;
  double **morse1;
  double **offset;

  virtual void allocate();
};
} // namespace CAC_NS

#endif

#endif
