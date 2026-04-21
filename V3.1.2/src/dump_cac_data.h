#ifdef DUMP_CLASS

DumpStyle(cac/data, DumpCACData)

#else

#ifndef CAC_DUMP_CAC_DATA_H
#define CAC_DUMP_CAC_DATA_H

#include "dump.h"

namespace CAC_NS {

class DumpCACData : public Dump {
 public:
  DumpCACData(CAC *, int, char **);
  ~DumpCACData();

 protected:
  int nevery;
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int reference_flag;
  int precision_flag;
  int velocity_flag;
  char *id_fix;
  class FixStore *fix;


  void set_atom_arrays(int);
  void set_elem_arrays(int);
  void pass_on_atom_arrays(int, int, int);
  void pass_on_elem_arrays(int, int, int * );

  void init_style();
  void write_header(bigint);
  void write_atom_header(bigint);
  void write_elem_header(bigint);
  void write_node_header(bigint);
  void pack_atom(tagint *);
  void pack_elem(tagint *);
  void pack_node(tagint *);
  int convert_atom_string(int, int, double *);
  int convert_node_string(int, int, double *);
  int convert_elem_string(int, int, double *);
  void write_atom_data(int, double *);
  void write_node_data(int, double *);
  void write_elem_data(int, double *);
  int count_atoms();
  int count_elements();
  int count_nodes();

  // functions & parameters for reference configuration
  
  // domain values at reference configuration
  double boxxlo_ref, boxxhi_ref;
  double boxylo_ref, boxyhi_ref;
  double boxzlo_ref, boxzhi_ref;
  double boxxy_ref, boxxz_ref, boxyz_ref;
  double reference_triclinic;
  double external_file_reference;
  int *atom_read_flag;
  int ***node_read_flag;
  FILE *ref_fp;

  void read_file();
  void header();
  void atoms();
  void nodes();
  void parse_keyword(int);
  char *line, *keyword, *buffer;

  typedef void (DumpCACData:: * FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  void header_item(bigint);
  void header_item_triclinic(bigint);
  void header_item_reference(bigint);
  void header_item_reference_triclinic(bigint);

  typedef void (DumpCACData:: * FnPtrWrite)(int, double *);
  FnPtrWrite write_atom_choice;
  FnPtrWrite write_elem_choice;
  FnPtrWrite write_node_choice;
  void write_atom_lines_noimage(int, double *);
  void write_elem_lines_noimage(int, double *);
  void write_node_lines_noimage(int, double *);
  void write_string(int, double *);
  void write_atom_lines_reference(int, double *);
  void write_elem_lines_reference(int, double *);
  void write_node_lines_reference(int, double *);

  typedef void (DumpCACData:: * FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  FnPtrPack pack_elem_choice;               // ptr to pack functions
  FnPtrPack pack_node_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);
  void pack_elem_noscale_noimage(tagint *);
  void pack_node_noscale_noimage(tagint *);
  void pack_atom_noscale_reference(tagint *);
  void pack_elem_noscale_reference(tagint *);
  void pack_node_noscale_reference(tagint *);
 

  typedef int (DumpCACData:: * FnPtrConvert)(int, int, double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  FnPtrConvert convert_elem_choice;          // ptr to convert data functions
  FnPtrConvert convert_node_choice;          // ptr to convert data functions
  int convert_atom_noimage(int, int, double *);
  int convert_elem_noimage(int, int, double *);
  int convert_node_noimage(int, int, double *);
  int convert_atom_reference(int, int, double *);
  int convert_node_reference(int, int, double *);
  int convert_elem_reference(int, int, double *);

};
}

#endif
#endif
