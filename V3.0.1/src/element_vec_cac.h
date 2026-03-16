#ifdef ELEMENT_CLASS

ElementStyle(cac, ElementVecCAC)

#else

#ifndef CAC_ELEMENT_VEC_CAC_H
#define CAC_ELEMENT_VEC_CAC_H

#include "element_vec.h"

namespace CAC_NS {

class ElementVecCAC : public ElementVec {

 public:
  ElementVecCAC(class CAC *);
  ~ElementVecCAC();
  void grow(int);
  void copy(int, int, int);
  void set_element_types(const char *, int);
  void set_data_size();
  void setup_integration_point(int);
  void grow_etype_arrays(int);
  void update_center_coord();
  void update_node_coord(int, int);
  void update_node_coord();
  void setup_sub_element(int);
  int element2atom(int);
  int element2element(int, int, tagint);
  void hex2hex(int, int, int, int);
  void hex2multihex(int, int, int, int, int);
  void hex2hexquad(int, int, int, int);
  void hex2wedge(int, int, int, int, int);
  void wedge2pyrtet(int, int, int, int);
  void wedge2tet(int, int, int, int);
  void pyr2tet(int, int, int, int);
  void hex2pyrtet(int, int, int, int, int);
  virtual void create_element(double ***, int, int *, int);
  virtual void create_element(double *, int, int *, int);
  virtual void create_pass_on_element(int, int *, int, int *, tagint);
  virtual void create_pass_on_atom(int, int, int, int, tagint);
  virtual void add_atoms(int, char **);
  virtual void compute_element_geometry(int);
  virtual void interpolate(double *, double ****, int, int, int, int);
  virtual void interpolate(double *, int ****, int, int, int, int);
  virtual double interpolate(double ****, int, int, int, int);
  virtual double interpolate(double ***, int, int, int);
  virtual double interpolate(int ****, int, int, int, int);
  virtual double interpolate(int ***, int, int, int);

  bigint memory_usage();

  // functions for comm
 
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);

  // functions for atomic_strain commands

  void data_node_strain(int, char **);

  // functions for read_data commands

  void data_element(double *, imageint, char **);
  void data_node(int, char **, int, int, double *);
  void data_node_reference(int, char **, int, int, double *, int);
  void data_vel(int, char **);
  void data_element_V2(double *, imageint, char **);
  void data_node_V2(int, char **, int, double *);
  void data_node_reference_V2(int, char **, int, double *);
  void data_vel_V2(int, char **);

  // functions for write_data command

  virtual void pack_element_data(double **);
  virtual void write_element_data(FILE *, int, double **);
  virtual void pack_node_data(double **);
  virtual void write_node_data(FILE *, int, double **);
  virtual void pack_vel_data(double **);
  virtual void write_vel_data(FILE *, int, double **);
  virtual void write_element_types(FILE *);

  // functions for write_data_elem command
  
  virtual void pack_element_data_elem(double **, tagint, int);
  virtual void pack_element_vel_data_elem(double **, tagint);
  virtual void pack_node_data_elem(double **, tagint, int);
  virtual void pack_node_vel_data_elem(double **, tagint);
  virtual void write_data_elem(FILE *, int, double **);
  virtual void write_vel_data_elem(FILE *, int, double **);

  // functions for write_tecplot command
 
  virtual void pack_node_tecplot(double **, int);
  virtual void pack_node_tecplot_binary(double *, int, int);
  virtual void pack_node_connect_tecplot(int *, int **);
  virtual void pack_node_connect_tecplot(int *, int ***);
  virtual void write_node_tecplot(FILE *, int, double **);
  virtual void write_node_connect_tecplot(FILE *, int, int *);

 
 protected:
  
  // copy pointers from element class
  
  tagint *tag;
  int *etype, **ctype, *mask, ***nodemask;
  imageint *image;
  double **x, ****nodex, ****nodex_current, ****nodev, ****nodef, ****gaussf;
  double ***surface_plane, **element_box2xi, **element_box2xi_scale, **cell_size;
  double ***shape_array, ***shape_array_center_subelem, ****shape_array_corner_subelem;
  double ***weighted_shape_array;
  double **weight, *subelem_size, **initial_box_size, **element_bound_box;
  int **g2u, **u2g, **g2n, **is_outer, **n2g, **n2u, ***us2u, *nucell, *ngcell;
  int ***surface_ucell, **nsurface_ucell, ***edge_ucell, **nedge_ucell;
  int *nsubelem;                        // local copy of element->nsubelem
  int **nucell_subelem, *subsplit;


  void setup_interpolate_element(int);

  // debug pointers
  
  int **debug_gcell_type, **debug_gcell_inactive;

};

}

#endif
#endif
