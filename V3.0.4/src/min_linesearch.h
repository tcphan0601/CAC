#ifndef CAC_MIN_LSRCH_H
#define CAC_MIN_LSRCH_H

#include "min.h"

namespace CAC_NS {

class MinLineSearch : public Min {
 public:
  MinLineSearch(class CAC *);
  ~MinLineSearch();
  void init();
  void setup_style();
  void reset_vectors();

 protected:
  // vectors needed by linesearch minimizers
  // combined locally to one array

  double *x0;                 // coords at start of linesearch
  double *g;                  // old gradient vector
  double *h;                  // search direction vector

  // arrays for atoms/nodes allocated and stored by fix_minimize
  // pointer copy from fix_minimize
  double *atomx0, *nodex0;
  double *atomg, *nodeg;
  double *atomh, *nodeh;

  double *gextra;             // g,h for extra global dof, x0 is stored by fix
  double *hextra;

  double **x0extra_atom;      // x0,g,h for extra per-atom dof
  double **gextra_atom;
  double **hextra_atom;

  typedef int (MinLineSearch::*FnPtr)(double, double &);
  FnPtr linemin;
  int linemin_backtrack(double, double &);
  int linemin_quadratic(double, double &);
  int linemin_forcezero(double, double &);

  double alpha_step(double, int);
  double compute_dir_deriv(double &);

  virtual void copy_vectors();
  virtual void copy_force();
};

}

#endif
