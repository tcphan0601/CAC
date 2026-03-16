#ifdef DUMP_CLASS

DumpStyle(tecplot, DumpTecplot)

#else

#ifndef CAC_DUMP_TECPLOT_H
#define CAC_DUMP_TECPLOT_H

#include "dump.h"

namespace CAC_NS {

class DumpTecplot : public Dump {
 public:
  DumpTecplot(CAC *, int, char **);
  ~DumpTecplot();
  void write();
  void openfile();
  
 protected:
  int nevery;
  class Compute *stress;
  class Compute *displace;
  class Compute *pe;
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int npe_connect;               // number of nodes in the connectivity section (8 for 3D and 4 for 2D)
  int stress_flag;
  int velocity_flag;
  int force_flag;
  int displace_flag;
  int pe_flag;
  int owner_proc_flag;
  int average_flag;              // 1 if average nodal values of elements in the same unit cell
  int **nodecell_ids;
  int ***node_ids;
  int maxelem;

  // parameters for TECIO
  
  char *title, *values;
  int *valuelocation;
  int debug, zonetype;
  double soltime;
  int fileformat; // 0 == PLT, 1 == SZPLT, 2 == ASCII
  int dummy, dummy1;
  int ival;            // value counter
  int valueflag;

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
  int count_elements();
  int count_nodes(); 
  void comm_buf_tecio(int);
  void grow_buf(int);

  typedef void (DumpTecplot:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  FnPtrHeader atom_header_choice;                       // ptr to write header functions
  FnPtrHeader elem_header_choice;                       // ptr to write header functions
  FnPtrHeader node_header_choice;                       // ptr to write header functions
  void header_item(bigint);
  void atom_header_item(bigint);
  void node_header_item(bigint);
  void elem_header_item(bigint);

  typedef void (DumpTecplot:: * FnPtrWrite)(int, double *);
  FnPtrWrite write_atom_choice;
  FnPtrWrite write_elem_choice;
  FnPtrWrite write_node_choice;
  void write_string(int, double *);
  void write_lines_noimage(int, double *);
  void write_elem_lines_noimage(int, double *);


  typedef void (DumpTecplot:: * FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  FnPtrPack pack_node_choice;               // ptr to pack functions
  FnPtrPack pack_elem_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);
  void pack_atom_noscale_noimage_binary(tagint *);
  void pack_elem_noscale_noimage(tagint *);
  void pack_elem_noscale_noimage_binary(tagint *);
  void pack_node_noscale_noimage(tagint *);
  void pack_node_noscale_noimage_binary(tagint *);

  typedef int (DumpTecplot:: * FnPtrConvert)(int, int, double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  FnPtrConvert convert_elem_choice;          // ptr to convert data functions
  FnPtrConvert convert_node_choice;          // ptr to convert data functions
  int convert_atom_noimage(int, int, double *);
  int convert_node_noimage(int, int, double *);
  int convert_elem_noimage(int, int, double *);

};
}

#endif
#endif
