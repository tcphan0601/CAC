#ifdef PAIR_CLASS

PairStyle(coul/wolf,PairCoulWolf)

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
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  //void write_restart(FILE *);
  //void read_restart(FILE *);
  //void write_restart_settings(FILE *);
  //void read_restart_settings(FILE *);
  //double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_coul, cut_coulsq, alf;

  void allocate();
};

}

#endif
#endif
