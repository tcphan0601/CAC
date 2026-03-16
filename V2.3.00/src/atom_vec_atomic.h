#ifdef ATOM_CLASS

AtomStyle(atomic,AtomVecAtomic)

#else

#ifndef CAC_ATOM_VEC_ATOMIC_H
#define CAC_ATOM_VEC_ATOMIC_H

#include "atom_vec.h"

namespace CAC_NS {

class AtomVecAtomic : public AtomVec {

 public:
  AtomVecAtomic(class CAC *);
  virtual ~AtomVecAtomic() {}
  void grow(int);
  void copy(int, int, int);
  void data_atom(double *, imageint, char **);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);
  virtual void pack_data(double **);
  virtual void write_data(FILE *, int, double **);
  virtual void pack_dat(double **);
  virtual void write_dat(FILE *, int, double **);
  virtual void create_atom(double *,int,int);

  bigint memory_usage();
 
 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;

};

}

#endif
#endif
