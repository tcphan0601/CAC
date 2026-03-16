#ifdef PAIR_CLASS

PairStyle(sw,PairSW)

#else

#ifndef CAC_PAIR_SW_H
#define CAC_PAIR_SW_H

#include "pair.h"

namespace CAC_NS {

class PairSW : public Pair {
 public:
  PairSW(class CAC *);
  virtual ~PairSW();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();

  struct Param {
    double epsilon,sigma;
    double littlea,lambda,gamma,costheta;
    double biga,bigb;
    double powerp,powerq;
    double tol;
    double cut,cutsq;
    double sigma_gamma,lambda_epsilon,lambda_epsilon2;
    double c1,c2,c3,c4,c5,c6;
    int ielement,jelement,kelement;
  };

  double memory_usage();

 protected:
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction
  int maxshort_a;               // size of short atom neighbor list array
  int maxshort_ia;              // size of short interpolated atom neighbor list array
  int *neighshort_a;            // short atom neighbor list array
  int *neighshort_ia;           // short interpolated atom neighbor list array - parent element index
  int *neighshort_ia_index;     // short interpolated atom neighbor list array - index within elements

  // short nia neighbor list
  // list of nia index that need to looped through in the current timestep
  
  int nniamax;
  int *shortnia_list;          // list of shortnia: 0 = nia index
  int *shortnia_flag;

  int nemax;                    // allocated size of per-element arrays
  double ***intgf;              // accumulated force at integration points

  virtual void allocate();
  void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double , double , 
      double, double, double, double, double *, double *, int, double &);
};

}

#endif
#endif

