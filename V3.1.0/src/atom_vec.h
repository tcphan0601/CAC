#ifndef CAC_ATOM_VEC_H
#define CAC_ATOM_VEC_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class AtomVec : protected Pointers {
 public:
  int mass_type;                       // 1 if per-type masses
  int forceclearflag;                  // 1 if has forceclear() method
  int comm_x_only;                     // 1 if only exchange x in forward comm
  int comm_f_only;                     // 1 if only exchange f in reverse comm
  int size_forward;                    // # of values per atom in comm
  int size_reverse;                    // # in reverse comm
  int size_border;                     // # in border comm
  int size_velocity;                   // # of velocity based quantities
  int size_data_atom;                  // number of values in Atom line in data file
  int size_data_atom_V2;               // number of values in Atom line in data file
  int size_tecplot_atom;               // number of values in Discrete Atom line in .dat file for tecplot
  int size_data_vel;                   // number of values in Velocity line
  int size_data_vel_V2;                // number of values in Velocity line
  int size_data_bonus;                 // number of values in Bonus line
  int xcol_data;                       // column (1-N) where x is in Atom line
  int xcol_data_V2;                       // column (1-N) where x is in Atom line
  int nargcopy;          // copy of command-line args for atom_style command
  char **argcopy;        // used when AtomVec is realloced (restart, replicate)

  AtomVec(class CAC *);
  virtual ~AtomVec();
  void store_args(int, char **);
  virtual void process_args(int, char **);
  virtual void init();

  virtual void grow(int) = 0;
  bigint roundup(bigint);
  virtual void copy(int, int, int) = 0;
  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;
  virtual int pack_border_vel(int, int *, double *, int, int *) = 0;
  virtual void unpack_border(int, int, double *) = 0;
  virtual void unpack_border_vel(int, int, double *) = 0;
  virtual int pack_border(int, int *, double *, int, int *) = 0;
  virtual int pack_comm(int, int *, double *, int, int *) = 0;
  virtual int pack_comm_vel(int, int *, double *, int, int *) = 0;
  virtual void unpack_comm(int, int, double *) = 0;
  virtual void unpack_comm_vel(int, int, double *) = 0;
  virtual int pack_reverse(int, int, double *) = 0;
  virtual void unpack_reverse(int, int *, double *) = 0;
  virtual void data_atom(double *, int, imageint, char **) = 0;
  virtual void data_atom_V2(double *, int, imageint, char **) = 0;
  virtual void data_atom_strain(double *, double *, imageint, char **) = 0;
  virtual void data_vel(int, char **);
  virtual void pack_data(double **) = 0;
  virtual void pack_data_elem(double **, int) = 0;
  virtual void pack_tecplot(double **) = 0;
  virtual void pack_tecplot_binary(double *, int) = 0;
  virtual void pack_vel(double **);
  virtual void write_data(FILE *, int, double **) = 0;
  virtual void write_tecplot(FILE *, int, double **) = 0;
  virtual void write_vel(FILE *, int, double **);
  virtual void create_atom(double *, int, tagint, int gtag = -1) = 0;
  virtual bigint memory_usage() = 0;

 protected:
  int nmax;                             // local copy of atom->nmax
  int deform_vremap;                    // local copy of domain properties
  int deform_groupbit;
  double *h_rate;

  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // this avoids aliasing issues by having 2 pointers (double, int)
  //   to same buf memory
  // constructor for 32-bit int prevents compiler
  //   from possibly calling the double constructor when passed an int
  // copy to a double *buf:
  //   buf[m++] = ubuf(foo).d, where foo is a 32-bit or 64-bit int
  // copy from a double *buf:
  //   foo = (int) ubuf(buf[m++]).i;, where (int) or (tagint) match foo
  //   the cast prevents compiler warnings about possible truncation
  //
  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };
  
  void grow_nmax();
};

}

#endif
