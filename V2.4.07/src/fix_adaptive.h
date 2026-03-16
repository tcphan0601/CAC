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
  void setup(int);
  void pre_exchange();
  void setup_pre_exchange();
  void post_integrate();
  //double compute_vector(int);
  //
  
  void grow_atom_arrays(int);
  void copy_atom_arrays(int, int, int);
  void set_atom_arrays(int);
  void pass_on_atom_arrays(int, int, int);
  int pack_atom_exchange(int, double *);
  int unpack_atom_exchange(int, double *);

  void grow_elem_arrays(int);
  void copy_elem_arrays(int, int, int);
  void set_elem_arrays(int);
  void pass_on_elem_arrays(int, int, int *);
  int pack_elem_exchange(int, double *);
  int unpack_elem_exchange(int, double *);

  double memory_usage();

 private:
  //int nsplitgroups,ntestgroups;
  //int **test_ids,**split_ids;
  int split_dim;                
  double blo[3],bhi[3];
  //double *lo,*hi;
  int nplanes;
  int *splitflag;                             // 0 if has not been splitted
  int *splitlist;
                                              // 1 if has been splitted
  int *disorder_flag;                         // disorder flag for each plane in local proc;
  int *disorder_flagall;                      // disorder flag for each plane in system
 
  //char *id_store;
  //class FixStore *fix;

  double *atom_plane_id,*elem_plane_id;       // storing plane ids for atoms and elements
                                              // must be double since elem_plane_id might be *.5
                                              // not checking if < 0
  int perfect_crystal;                        // compute value for perfect crystal

  char *id_compute;
  int icompute;
  class Compute *compute;

  int pending;
  double l;

  int find_plane(double *);

};

}

#endif
#endif
