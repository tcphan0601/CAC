#ifdef ELEMENT_CLASS

ElementStyle(rhb,ElementVecRHB)

#else

#ifndef CAC_ELEMENT_VEC_RHB_H
#define CAC_ELEMENT_VEC_RHB_H

#include "element_vec.h"

namespace CAC_NS {

class ElementVecRHB : public ElementVec {

 public:
  ElementVecRHB(class CAC *);
  ~ElementVecRHB();
  void grow(int);
  void copy(int, int, int);
  void set_interpolate(const char *, int);
  void set_integration(const char *, int);
  //int get_size(int, int);
  void setup_integration_point(int);
  void grow_etype_arrays(int);
  void set_intg_user(char *, int);
  void set_nintg_user(char *, int);
  void check_element_size();
  void update_center_coord();
  void update_node_coord(int, int);
  void update_node_coord();
  void setup_sub_element();  
  void setup_sub_element(int);
  int element2atom(int);
  int element2element(int, int, tagint);
  void check_type();
  int split_element(int, int, double, int, int, double *);
  int split_element(int **, int);
  int check_split_elem(int, int, int, int &, int &);
  void add_requested_etype();
  void request_new_etype(int *);
  int find_etype(int *);
  virtual void create_element(double *, double **, int, int, int);
  virtual void create_element(double *, int, int, int);
  virtual void create_pass_on_element(int, int *, int, int, tagint);
  virtual void create_pass_on_atom(int, int, int, tagint);
  virtual bigint create_element_template(int, int, int, 
      int, int, int, int, int *, double **, int *); 
  virtual void add_atoms(int, char **);
  virtual void update_slip_plane();
  virtual void interpolate(double *, double ***, int, int, int);
  bigint memory_usage();
  virtual void update_npe();

  // functions for comm
 
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);

  // functions for read_data commands

  void data_element(double *, imageint, char **);
  void data_node(int, char **, int, double *);
  void data_vel(int, char **);

  // functions for write_data command

  virtual void pack_element_data(double **);
  virtual void write_element_data(FILE *, int, double **);
  virtual void pack_node_data(double **);
  virtual void write_node_data(FILE *, int, double **);
  virtual void pack_vel_data(double **);
  virtual void write_vel_data(FILE *, int, double **);
  virtual void write_interpolate(FILE *);
  virtual void write_integration(FILE *);

  // functions for write_data_elem command
  
  virtual void pack_element_data_elem(double **, tagint);
  virtual void pack_element_vel_data_elem(double **, tagint);
  virtual void pack_node_data_elem(double **, tagint);
  virtual void pack_node_vel_data_elem(double **, tagint);
  virtual void write_data_elem(FILE *, int, double **);
  virtual void write_vel_data_elem(FILE *, int, double **);

  // functions for write_dat command
 
  virtual void pack_node_dat(double **);
  virtual void pack_node_connect_dat(double **);
  virtual void write_node_dat(FILE *, int, double **);
  virtual void write_node_connect_dat(FILE *, int, double **);


 
 protected:
  
  // copy pointers from element class
  
  tagint *tag,**nodetag;
  int *etype,*ctype,*mask,**nodemask;
  imageint *image;
  double **x,***nodex,***nodev,***nodef,***slip_plane,**cell_size;
  double ***shape_array,**shape_array_center_subelem;
  double ***weighted_shape_array;
  double **weight,*elem_size,*initial_size;
  int **i2ia,**ia2i,**i2n,**n2i,**n2ia,***ias2ia,*nintpl,*nintg;
  int ***surface_intpl,**nsurface_intpl,***edge_intpl,**nedge_intpl;
  int **natom_subelem;

  // data to specify interpolated atoms and integration points

  int **naae;
  int **niae;
  int ***intg;

  int node_set_3D[3][2][4];       // set of nodes for each face of element (+-x, +-y, +-z)   
  int node_set_2D[2][2][2];       // ditto for 2D

  void setup_interpolate_element(int);


};

}

#endif
#endif
