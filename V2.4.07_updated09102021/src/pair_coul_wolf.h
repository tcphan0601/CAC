#ifdef PAIR_CLASS

PairStyle(coul/wolf, PairCoulWolf)

#else

#ifndef CAC_PAIR_COUL_WOLF_H
#define CAC_PAIR_COUL_WOLF_H

#include "pair.h"

namespace CAC_NS {

class PairCoulWolf : public Pair {

public:
  PairCoulWolf(class CAC *);
  virtual ~PairCoulWolf();
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
  double alpha;
  double cut_global;
  double **cut;
  double **q_one_matrix, **q_two_matrix;
  double *cut_respa;

  double nemax;
  virtual void allocate();
};
} // namespace CAC_NS

#endif

#endif
