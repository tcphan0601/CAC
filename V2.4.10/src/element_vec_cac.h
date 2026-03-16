#ifdef ELEMENT_CLASS

ElementStyle(cac,ElementVecCAC)

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
  void hex2hexquad(int, int, int, int);
  void hex2wedge(int, int, int, int, int);
  void wedge2pyrtet(int, int, int, int);
  void wedge2tet(int, int, int, int);
  void pyr2tet(int, int, int, int);
  void hex2pyrtet(int, int, int, int, int);
  //int split_element(int, int, double, int, int, double *);
  int split_element(int **, int);
  int check_split_element(int, int, int, double, int *, int &, int);
  int check_split_element(int, double [3], double [3], double, double, int *, int &);
  virtual void create_element(double **, int, int, int);
  virtual void create_element(double *, int, int, int);
  virtual void create_pass_on_element(int, int *, int, int, tagint);
  virtual void create_pass_on_atom(int, int, int, tagint);
  virtual void add_atoms(int, char **);
  virtual void update_surface_plane();
  virtual void interpolate(double *, double ***, int, int, int);
  virtual double interpolate(double ***, int, int, int);
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
  void data_node(int, char **, int, double *);
  void data_node_reference(int, char **, int, double *);
  void data_vel(int, char **);

  // functions for write_data command

  virtual void pack_element_data(double **);
  virtual void write_element_data(FILE *, int, double **);
  virtual void pack_element_cluster_data(double **);
  virtual void write_element_cluster_data(FILE *, int, double **);
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
 
  virtual void pack_node_tecplot(double **);
  virtual void pack_node_tecplot_binary(double *, int);
  virtual void pack_node_connect_tecplot(int *);
  virtual void write_node_tecplot(FILE *, int, double **);
  virtual void write_node_connect_tecplot(FILE *, int, int *);

  // functions for write_data_adrian command
  
  virtual void pack_element_data_adrian(double **);
  virtual void write_element_data_adrian(FILE *, int, double **);


 
 protected:
  
  // copy pointers from element class
  
  tagint *tag,**nodetag;
  int *etype,*ctype,*mask,**nodemask;
  imageint *image;
  double **x,***nodex,***nodex_current,***nodev,***nodef,***surface_plane,**cell_size;
  double ***shape_array,***shape_array_center_subelem,****shape_array_corner_subelem;
  double ***weighted_shape_array;
  double **weight,*subelem_size,**initial_box_size,**element_bound_box;
  int **i2ia,**ia2i,**i2n,**is_outer,**n2i,**n2ia,***ias2ia,*nintpl,*nintg;
  int ***surface_intpl,**nsurface_intpl,***edge_intpl,**nedge_intpl;
  int *nsubelem;                        // local copy of element->nsubelem
  int **natom_subelem,*subsplit;
  int **element_clusters;
  int element_cluster_flag;

  int node_set_3D[3][2][4];       // set of nodes for each face of element (+-x, +-y, +-z)   
  int node_set_2D[2][2][2];       // ditto for 2D

  void setup_interpolate_element(int);

  // debug pointers
  
  int **debug_intg_type,**debug_intg_inactive;

};

}

#endif
#endif
