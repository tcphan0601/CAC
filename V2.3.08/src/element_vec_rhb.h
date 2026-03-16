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
  void data_element(double *, imageint, char **);
  void data_node(int, char **);
  void data_vel(int, char **);

  void set_interpolate(const char *, int);
  void set_integration(const char *, int);
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

  virtual void update_npe();
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);
  virtual void pack_element_data(double **);
  virtual void write_element_data(FILE *, int, double **);
  virtual void pack_node_data(double **);
  virtual void pack_node_dat(double **);
  virtual void pack_node_connect_dat(double **);
  virtual void write_node_data(FILE *, int, double **);
  virtual void write_node_dat(FILE *, int, double **);
  virtual void write_node_connect_dat(FILE *, int, double **);
  virtual void pack_vel(double **);
  virtual void write_vel(FILE *, int, double **);
  virtual void write_interpolate(FILE *);
  virtual void write_integration(FILE *);
  virtual void create_element(double *, double **, int, int, int);
  virtual void create_element(double *, int, int, int);
  virtual bigint create_element_template(int, int, int, 
      int, int, int, int, int *, double **, int *); 
  virtual void add_atoms(int, char **);
  virtual void update_slip_plane();
  virtual void interpolate(double *, double ***, int, int, int);
  virtual void update_ele_surf_area();
  virtual void update_ele_rand_amp(double, int);
  virtual void update_surface_normal_vector();
  bigint memory_usage();
 
 protected:
  
  // copy pointers from element class
  
  tagint *tag,**nodetag;
  int *etype,*ctype,*mask,**nodemask;
  imageint *image;
  double **x,***nodex,***nodev,***nodef,***slip_plane,**ele_surf_area,***surface_normvec;
  double ***shape_array,**shape_array_center_subelem;
  double ***weighted_shape_array;
  double **weight,*elem_size,*initial_size,*ele_rand_amp;
  int **i2ia,**ia2i,**i2n,**n2i,***ias2ia,*nintpl,*nintg;
  int **natom_subelem;
  int **elemlink;
  double linkcut;

  // data to specify interpolated atoms and integration points

  int **naae;
  int **niae;
  int ***intg;

  void setup_interpolate_element(int);
};

}

#endif
#endif
