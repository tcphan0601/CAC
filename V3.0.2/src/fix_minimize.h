#ifdef FIX_CLASS

FixStyle(MINIMIZE,FixMinimize)

#else

#ifndef CAC_FIX_MINIMIZE_H
#define CAC_FIX_MINIMIZE_H

#include "fix.h"

namespace CAC_NS {

class FixMinimize : public Fix {
  friend class MinLineSearch;

 public:
  FixMinimize(class CAC *, int, char **);
  virtual ~FixMinimize();
  int setmask();
  virtual void init() {}

  double memory_usage();
  virtual void grow_atom_arrays(int);
  virtual void copy_atom_arrays(int, int, int);
  virtual int pack_atom_exchange(int, double *);
  virtual int unpack_atom_exchange(int, double *);
  virtual void grow_elem_arrays(int);
  virtual void copy_elem_arrays(int, int, int);
  virtual int pack_elem_exchange(int, double *);
  virtual int unpack_elem_exchange(int, double *);

  virtual void add_vector(int);
  double *request_vector(int, int);
  void store_box();
  void reset_coords();

 protected:
  int nvector;
  int *peratom;
  double **atomvectors;
  double **nodevectors;
  double boxlo[3], boxhi[3];

  void box_swap();
};

}

#endif
#endif
