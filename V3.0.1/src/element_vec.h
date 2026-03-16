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

  // for write_data_elem command
  
  int size_data_elem;
  int size_data_elem_vel;                   // number of values in Velocity line

  // for read_data/write_data command
  
  int size_data_element;               // number of values in Elements line
  int size_data_node;                  // number of values in Nodes line
  int size_data_vel;                   // number of values in Node Velocities line
  int size_data_element_V2;               // number of values in Elements line
  int size_data_node_V2;                  // number of values in Nodes line
  int size_data_vel_V2;                   // number of values in Node Velocities line


  // for write_tecplot command
  
  int size_tecplot_node;                   // number of values in Coarse Grain line in write_tecplot file
  int size_tecplot_node_connect;           // number of values in Node Connectivity line in write_tecplot file
  

  int size_data_bonus;                 // number of values in Bonus line

  int xcol_data;                       // column (1-N) where x is in Element line
  int xcol_data_V2;                       // column (1-N) where x is in Element line
  int nargcopy;          // copy of command-line args for element_style command
  char **argcopy;        // used when ElementVec is realloced (restart, replicate)

  int *element_type_setflag;

  // data to specify virtual atoms and gaussian cells

  int **ncells;     // number of cells along each direction
  int **ngcells;    // number of gaussian cells along each direction


  ElementVec(class CAC *);
  virtual ~ElementVec();
  void store_args(int, char **);
  void update_apc(int);
  virtual void process_args(int, char **);
  virtual void init();
  virtual void add_requested_etype(int);
  virtual void request_new_etype(int *);
  virtual int find_etype(int *);

  virtual void set_data_size() = 0;
  virtual void add_atoms(int, char **) = 0;
  virtual int element2atom(int) = 0;
  virtual int element2element(int, int, tagint) = 0;
  virtual void hex2hex(int, int, int, int) = 0;
  virtual void hex2multihex(int, int, int, int, int) = 0;
  virtual void hex2hexquad(int, int, int, int) = 0;
  virtual void hex2wedge(int, int, int, int, int) = 0;
  virtual void hex2pyrtet(int, int, int, int, int) = 0;
  virtual void wedge2pyrtet(int, int, int, int) = 0;
  virtual void wedge2tet(int, int, int, int) = 0;
  virtual void pyr2tet(int, int, int, int) = 0;
  virtual void set_element_types(const char *, int) = 0;
  virtual void setup_integration_point(int) = 0;
  virtual void update_center_coord() = 0;
  virtual void update_node_coord(int, int) = 0;
  virtual void update_node_coord() = 0;
  virtual void setup_sub_element(int) = 0;
  virtual void grow_etype_arrays(int) = 0;
  virtual void grow(int) = 0;
  virtual void copy(int, int, int) = 0;
  virtual void create_element(double ***, int, int *, int) = 0;
  virtual void create_element(double *, int, int *, int) = 0;
  virtual void create_pass_on_element(int, int *, int, int *, tagint) = 0;
  virtual void create_pass_on_atom(int, int, int, int, tagint) = 0;
  virtual void compute_element_geometry(int) = 0;
  virtual void interpolate(double *, double ****, int, int, int, int) = 0;
  virtual void interpolate(double *, int ****, int, int, int, int) = 0;
  virtual double interpolate(double ****, int, int, int, int) = 0;
  virtual double interpolate(double ***, int, int, int) = 0;
  virtual double interpolate(int ****, int, int, int, int) = 0;
  virtual double interpolate(int ***, int, int, int) = 0;

  // functions for comm
  
  virtual int pack_exchange(int, double *) = 0;
  virtual int unpack_exchange(double *) = 0;
  virtual void unpack_border(int, int, double *) = 0;
  virtual int pack_border(int, int *, double *, int, int *) = 0;
  virtual int pack_comm(int, int *, double *, int, int *) = 0;
  virtual void unpack_comm(int, int, double *) = 0;
  virtual int pack_reverse(int, int, double *) = 0;
  virtual void unpack_reverse(int, int *, double *) = 0;

  // functions for atomic_strain commands

  virtual void data_node_strain(int, char **) = 0;

  // functions for read_data commands
 
  virtual void data_element(double *, imageint, char **) = 0;
  virtual void data_node(int, char **, int, int, double *) = 0;
  virtual void data_node_reference(int, char **, int, int, double *, int) = 0;
  virtual void data_vel(int, char **) = 0;
  virtual void data_element_V2(double *, imageint, char **) = 0;
  virtual void data_node_V2(int, char **, int, double *) = 0;
  virtual void data_node_reference_V2(int, char **, int, double *) = 0;
  virtual void data_vel_V2(int, char **) = 0;


  // functions for write_data command
  
  virtual void pack_element_data(double **) = 0;
  virtual void write_element_data(FILE *, int, double **) = 0;
  virtual void pack_node_data(double **) = 0;
  virtual void write_node_data(FILE *, int, double **) = 0;
  virtual void pack_vel_data(double **) = 0;
  virtual void write_vel_data(FILE *, int, double **) = 0;
  virtual void write_element_types(FILE *) = 0;

  // functions for write_data_elem command
  
  virtual void pack_element_data_elem(double **, tagint, int) = 0;
  virtual void pack_element_vel_data_elem(double **, tagint) = 0;
  virtual void pack_node_data_elem(double **, tagint, int) = 0;
  virtual void pack_node_vel_data_elem(double **, tagint) = 0;
  virtual void write_data_elem(FILE *, int, double **) = 0;
  virtual void write_vel_data_elem(FILE *, int, double **) = 0;

  // functions for write_tecplot command
  
  virtual void pack_node_tecplot(double **, int) = 0;
  virtual void pack_node_tecplot_binary(double *, int, int) = 0;
  virtual void write_node_tecplot(FILE *, int, double **) = 0;
  virtual void pack_node_connect_tecplot(int *, int **) = 0;
  virtual void pack_node_connect_tecplot(int *, int ***) = 0;
  virtual void write_node_connect_tecplot(FILE *, int, int *) = 0;
  

  virtual bigint memory_usage() = 0;

 protected:
  
  int nmax;                            // local copy of element->nmax
  int nmaxnode;                        // local copy of element->nmaxnode
  int nmaxgcell;                        // local copy of element->nmaxgcell
  bigint nmaxucell;                    // local copy of element->nmaxucell
  int *npe;                            // local copy of element->npe
  int *apc;                            // local copy of element->apc
  double *nodal_weight;                // local copy of element->nodal_weight
  int *element_shape_ids;              // local copy of element->element_shape_ids
  char **element_shape_names;            // local copy of element->element_shape_names
  int *requested_etype;
  int nrequested_etype;
  int maxrequested_etype;

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
