#ifdef FIX_CLASS

FixStyle(STORE,FixStore)

#else

#ifndef CAC_FIX_STORE_H
#define CAC_FIX_STORE_H

#include <cstdio>
#include "fix.h"

namespace CAC_NS {

class FixStore : public Fix {
 public:
  int nrow,ncol;         // size of global data array
  int nvalues;           // number of per-atom values
  double *vstore;        // vector storage for GLOBAL or PERATOM
  double **astore;       // array storage for GLOBAL or PERATOM
  double ***a3store;     // 3D array storage for PERATOM (element)
  int disable;        // 1 if operations (except grow) are currently disabled

  FixStore(class CAC *, int, char **);
  ~FixStore();
  int setmask();
  void reset_global(int, int);

  //void write_restart(FILE *);
  //void restart(char *);

  void grow_atom_arrays(int);
  void copy_atom_arrays(int, int, int);
  int pack_atom_exchange(int, double *);
  int unpack_atom_exchange(int, double *);
  void grow_elem_arrays(int);
  void copy_elem_arrays(int, int, int);
  int pack_elem_exchange(int, double *);
  int unpack_elem_exchange(int, double *);

  //int pack_restart(int, double *);
  //void unpack_restart(int, int);
  //int size_restart(int);
  //int maxsize_restart();

  double memory_usage();

 private:
  int flavor;                   // GLOBAL or PERATOM
  int vecflag;                  // 1 if ncol=1 or nvalues=1

  double *rbuf;                 // restart buffer for GLOBAL vec/array
};

}

#endif
#endif

