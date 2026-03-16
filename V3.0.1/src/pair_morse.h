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

  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_data(FILE *);
  void write_data_all(FILE *);

 protected:
  double cut_global;
  double **cut;
  double **d0, **alpha, **r0;
  double **morse1;
  double **offset;

  virtual void allocate();
};

} 

#endif

#endif
