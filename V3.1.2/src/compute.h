#ifndef CAC_COMPUTE_H
#define CAC_COMPUTE_H

#include "pointers.h"

namespace CAC_NS {

class Compute : protected Pointers {

 public:
  enum {
    INVOKED_NONE    = 0,
    INVOKED_SCALAR  = 1<<0,
    INVOKED_VECTOR  = 1<<1,
    INVOKED_ARRAY   = 1<<2,
    INVOKED_PERATOM = 1<<3,
    INVOKED_LOCAL   = 1<<4,
    INVOKED_PERGRID = 1<<5,
  };

  static int instance_total;     // # of Compute classes ever instantiated

  char *id, *style;
  int igroup, groupbit;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double **array;           // computed global array
  double *vector_atom;      // computed per-atom vector
  double ***vector_node;     // computed per-node vector
  double ***vector_vatom;      // computed per-ucell vector
  double **array_atom;      // computed per-atom array
  double ****array_node;     // computed per-node array
  double ****array_vatom;      // computed per-ucell-atom array
  double scalar_local;      // computed local scapar
  double *vector_local;     // computed local vector
  double **array_local;     // computed local array

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int group_flag;
  int array_flag;           // 0/1 if compute_array() function exists
  int size_vector;          // length of global vector
  int size_array_rows;      // rows in global array
  int size_array_cols;      // columns in global array
  int size_vector_variable;      // 1 if vec length is unknown in advance
  int size_array_rows_variable;  // 1 if array rows is unknown in advance

  int peratom_flag;         // 0/1 if compute_peratom() function exists
  int size_peratom_cols;    // 0 = vector, N = columns in peratom array

  int local_flag;           // 0 = no compute local
                            // 1 = compute scalar_local() exists
                            // 2 = compute vector_local() exists
                            // 3 = compute array_local() exists

  int size_vector_local;          // rows in local vector
  int size_array_local_rows;      // rows in local vector
  int size_array_local_cols;      // columns in local array
  int size_vector_local_variable;      // 1 if local vec length is unknown in advance
  int size_array_local_rows_variable;  // 1 if local array rows is unknown in advance



  int extscalar;            // 0/1 if global scalar is intensive/extensive
  int extvector;            // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;             // list of 0/1 int/ext for each vec component
  int extarray;             // 0/1 if global array is all intensive/extensive

  int preneighflag;        // to be set by caller if invoke build_one() before
                           // neigh list is re-build on the same step
  int tempflag;       // 1 if Compute can be used as temperature
                      // must have both compute_scalar, compute_vector
  int pressflag;      // 1 if Compute can be used as pressure (uses virial)
                      // must have both compute_scalar, compute_vector
  int pressatomflag;  // 1 if Compute calculates per-atom virial
                      // 2 if Compute calculates per-atom centroid virial
  int mechatomflag;   // 1 if Compute calculates per-atom mechanical stress
  int ghostskinflag;  // 1 if Compute requires additional ghost cutoff
  int peflag;         // 1 if Compute calculates PE (uses Force energies)
  int peatomflag;     // 1 if Compute calculates per-atom PE
  int create_attribute;    // 1 if compute stores attributes that need
                           // setting when a new atom is created

  int tempbias;       // 0/1 if Compute temp includes self/extra bias

  int timeflag;       // 1 if Compute stores list of timesteps it's called on
  int ntime;          // # of entries in time list
  int maxtime;        // max # of entries time list can hold
  bigint *tlist;      // list of timesteps the Compute is called on

  int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
  bigint invoked_scalar;  // last timestep on which compute_scalar() was invoked
  bigint invoked_vector;  // ditto for compute_vector()
  bigint invoked_array;   // ditto for compute_array()
  bigint invoked_peratom;   // ditto for compute_peratom()
  bigint invoked_scalar_local;  // ditto for compute_array_local()
  bigint invoked_vector_local;  // ditto for compute_vector_local()
  bigint invoked_array_local;   // ditto for compute_array_local()

  double ghostskin;       // additional ghost cutoff
  double dof;         // degrees-of-freedom for temperature

  int comm_atom_forward;         // size of forward communication (0 if none)
  int comm_atom_reverse;         // size of reverse communication (0 if none)
  int comm_elem_forward;         // size of forward communication (0 if none)
  int comm_elem_reverse;         // size of reverse communication (0 if none)
  int dynamic_group_allow;  // 1 if can be used with dynamic group, else 0

  unsigned int datamask;
  unsigned int datamask_ext;

  int cudable;              // 1 if compute is CUDA-enabled

  Compute(class CAC *, int, char **);
  virtual ~Compute();
  void reset_extra_dof();
  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_peratom() {} // may also compute pernode quantities
  virtual double compute_scalar_local() {return 0.0;}
  virtual void compute_vector_local() {}
  virtual void compute_array_local() {}
  virtual void set_atom_arrays(int) {}
  virtual void set_elem_arrays(int) {}
  virtual void pass_on_atom_arrays(int, int, int) {}
  virtual void pass_on_elem_arrays(int, int, int *) {}

  virtual int pack_atom_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_atom_forward_comm(int, int, double *) {}
  virtual int pack_atom_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_atom_reverse_comm(int, int *, double *) {}
  virtual int pack_elem_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_elem_forward_comm(int, int, double *) {}
  virtual int pack_elem_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_elem_reverse_comm(int, int *, double *) {}

  virtual void reset_extra_compute_fix(const char *);
  virtual void refresh() {}

  void addstep(bigint);
  int matchstep(bigint);
  void clearstep();
  virtual double memory_usage() {return 0.0;}

 protected:
  int instance_me;             // which Compute class instantiation I am

  double natoms_temp;          // # of atoms used for temperature calculation
  int extra_dof;               // extra DOF for temperature computes
  int fix_dof;                 // DOF due to fixes
  int dynamic;                 // recount atoms for temperature computes
  int dynamic_user;            // user request for temp compute to be dynamic
  int thermoflag;              // 1 if include fix PE for PE computes

  double vbias[3];             // stored velocity bias for one atom
  double **vbiasall;           // stored velocity bias for all atoms
  int maxbias;                 // size of vbiasall array

  

};

}

#endif
