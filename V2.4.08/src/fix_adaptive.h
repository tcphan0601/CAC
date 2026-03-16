#ifdef FIX_CLASS

FixStyle(adaptive,FixAdaptive)

#else

#ifndef CAC_FIX_ADAPTIVE_H
#define CAC_FIX_ADAPTIVE_H

#include "fix.h"

namespace CAC_NS {

class FixAdaptive : public Fix {
 public:
  FixAdaptive(class CAC *, int, char **);
  ~FixAdaptive();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_exchange();
  void setup_pre_exchange();
  void post_integrate();
  //double compute_vector(int);

  double memory_usage();

 private:

  void grow_splitlist();
  int **splitlist;
  int maxsplit;
  int nsplit;
  int nsplits_total;

  class NeighList *list;

  int pending;
  int perfect_crystal;                        // compute value for perfect crystal
  double tol,cut_atom_neigh;
  double eig_frac_thres;
  double split_width;
  double thres;                               // threshold for checking dislocation from atomisitc region

  char *id_structure;
  int istructure;
  class Compute *structure;

  //char *id_stress;
  //int istress;
  //class Compute *stress;

  // debug
  int debug;
  void write_neighbor(int);
  void write_surface_intpl(int);
  void write_split_element(int, int);
  void write_split_atom(int, int, int *); 

};

}

#endif
#endif
