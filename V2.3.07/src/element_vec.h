#ifndef CAC_ELEMENT_VEC_H
#define CAC_ELEMENT_VEC_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class ElementVec : protected Pointers {
 public:
  int size_forward;                    // # of values per atom in comm
  int size_reverse;                    // # in reverse comm
  int size_border;                     // # in border comm
  int size_velocity;                   // # of velocity based quantities
  int size_data_element;               // number of values in Element line
  int size_data_node;                  // number of values in Node line in data file
  int size_dat_node;                   // number of values in Coarse Grain line in dat file
  int size_dat_node_connect;           // number of values in Node Connectivity line in dat file
  int size_data_vel;                   // number of values in Node Velocity line
  int size_data_bonus;                 // number of values in Bonus line
  int xcol_data;                       // column (1-N) where x is in Element line
  int nargcopy;          // copy of command-line args for element_style command
  char **argcopy;        // used when ElementVec is realloced (restart,replicate)

  int *integration_setflag;
  int *interpolate_setflag;

  ElementVec(class CAC *);
  virtual ~ElementVec();
  void store_args(int, char **);
  void update_apc(int);
  virtual void process_args(int, char **);
  virtual void init();
  virtual void add_atoms(int, char **);


  virtual void update_npe() = 0;
  virtual int element2atom(int) = 0;
  virtual int element2element(int, int, tagint) = 0;
  virtual void check_type() = 0;
  virtual void check_element_size() = 0;
  virtual void set_interpolate(const char *, int) = 0;
  virtual void set_integration(const char *, int) = 0;
  virtual void setup_integration_point(int) = 0;
  virtual void set_intg_user(char *, int) = 0;
  virtual void set_nintg_user(char *, int) = 0;
  virtual void update_center_coord() = 0;
  virtual void update_node_coord(int, int) = 0;
  virtual void update_node_coord() = 0;
  virtual void setup_sub_element() = 0;
  virtual void setup_sub_element(int) = 0;
  virtual void grow_etype_arrays(int) = 0;
  virtual void grow(int) = 0;
  virtual void copy(int, int, int) = 0;
  virtual void create_element(double *, double **, int, int, int) = 0;
  virtual void create_element(double *, int, int, int) = 0;
  virtual bigint create_element_template(int, int, int, 
      int, int, int, int, int *, double **, int *) = 0; 
  virtual void update_slip_plane() = 0;
  virtual void interpolate(double *, double ***, int, int, int) = 0;

  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;
  virtual void unpack_border(int, int, double *) = 0;
  virtual int pack_border(int, int *, double *, int, int *) = 0;
  virtual int pack_comm(int, int *, double *, int, int *) = 0;
  virtual void unpack_comm(int, int, double *) = 0;
  virtual int pack_reverse(int, int, double *) = 0;
  virtual void unpack_reverse(int, int *, double *) = 0;
  virtual void data_element(double *, imageint, char **) = 0;
  virtual void data_node(int, char **, int, double *) = 0;
  virtual void data_vel(int, char **) = 0;
  virtual void pack_element_data(double **) = 0;
  virtual void write_element_data(FILE *, int, double **) = 0;
  virtual void pack_node_data(double **) = 0;
  virtual void pack_node_dat(double **) = 0;
  virtual void pack_node_connect_dat(double **) = 0;
  virtual void write_node_data(FILE *, int, double **) = 0;
  virtual void write_node_dat(FILE *, int, double **) = 0;
  virtual void write_node_connect_dat(FILE *, int, double **) = 0;
  virtual void pack_vel(double **) = 0;
  virtual void write_vel(FILE *, int, double **) = 0;
  virtual void write_interpolate(FILE *) = 0;
  virtual void write_integration(FILE *) = 0;
  virtual bigint memory_usage() = 0;

 protected:
  
  int nmax;                             // local copy of element->nmax
  int nmaxintg;                        // local copy of element->nmax_intg
  bigint nmaxintpl;                    // local copy of element->nmax_intpl
  int npe;                              // local copy of element->npe
  int apc;                              // local copy of element->apn
  int nsubelem;                         // local copy of element->nsubelem
  int subelemflag;

   // union data struct for packing 32-bit and 64-bit ints into double bufs
  // this avoids aliasing issues by having 2 pointers (double,int)
  //   to same buf memory
  // constructor for 32-bit int prevents compiler
  //   from possibly calling the double constructor when passed an int
  // copy to a double *buf:
  //   buf[m++] = ubuf(foo).d, where foo is a 32-bit or 64-bit int
  // copy from a double *buf:
  //   foo = (int) ubuf(buf[m++]).i;, where (int) or (tagint) match foo
  //   the cast prevents compiler warnings about possible truncation
  
  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };
  
  void grow_nmax();
  void grow_nmaxintg();
  void grow_nmaxintpl();

};

}

#endif
