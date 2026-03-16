#ifdef PAIR_CLASS

PairStyle(buck, PairBuck)

#else

#ifndef CAC_PAIR_BUCK_H
#define CAC_PAIR_BUCK_H

#include "pair.h"

namespace CAC_NS {

class PairBuck : public Pair {
 public:
  PairBuck(class CAC *);
  virtual ~PairBuck();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  //void write_restart(FILE *);
  //void read_restart(FILE *);
  //void write_restart_settings(FILE *);
  //void read_restart_settings(FILE *);
  //void write_data(FILE *);
  //void write_data_all(FILE *);
  //double single(int, int, int, int, double, double, double, double &);
  //void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut;
  double **a, **rho, **c;
  double **rhoinv, **buck1, **buck2, **offset;

  virtual void allocate();
};

}

#endif
#endif
