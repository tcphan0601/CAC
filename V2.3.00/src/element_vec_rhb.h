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
  void check_type();

  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
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
  virtual void update_slip_plane();
  virtual void interpolate(double *, double ***, int, int, int);
  bigint memory_usage();
 
 protected:
  
  // copy pointers from element class
  
  tagint *tag,**nodetag;
  int *etype,*ctype,*mask,**nodemask;
  imageint *image;
  double **x,***nodex,***nodev,***nodef,***slip_plane;
  double ***shape_array,**shape_array_center_subelem;
  double **weight,*elem_size,*initial_size;
  int **i2ia,**i2n,**n2i,***ias2ia,*nintpl,*nintg;
  int **natom_subelem;

  // data to specify interpolated atoms and integration points

  int **naae;
  int **niae;
  int ***intg;
  int *integration_setflag;
  int *interpolate_setflag;

  void setup_integration_point(int);
  void setup_interpolate_element(int);
};

}

#endif
#endif
