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
  double cut_global, cut_inner_global;
  double **cut, **cut_inner, **cut_inner_sq;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4;

  void allocate();
};

}

#endif
#endif

