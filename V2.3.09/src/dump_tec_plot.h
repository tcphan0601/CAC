#ifdef DUMP_CLASS

DumpStyle(tecplot, DumpTecPlot)

#else

#ifndef CAC_DUMP_TEC_PLOT_H
#define CAC_DUMP_TEC_PLOT_H

#include "dump.h"

namespace CAC_NS {

class DumpTecPlot : public Dump {
 public:
  DumpTecPlot(CAC *, int,char**);
  ~DumpTecPlot();
  int pack_elem_forward_comm(int, int *, double *, int, int *);
  void unpack_elem_forward_comm(int, int, double *);

 protected:
  int nevery;
  class Compute *stress;
  class Compute *displace;
  int nclusters;
  int maxclusters;
  int **clusters;
  int apc;
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int stress_flag;
  int force_flag;
  int displace_flag;
  int link_flag;
  double cutlink;
  tagint next_tag;

  void init_style();
  void setup_pre_write();
  void write_header(bigint);
  void write_atom_header(bigint);
  void write_elem_header(bigint);
  void write_node_header(bigint);
  void pack_atom(tagint *);
  void pack_elem(tagint *);
  void pack_node(tagint *);
  int convert_atom_string(int,int,double *);
  int convert_node_string(int,int,double *);
  int convert_elem_string(int,int,double *);
  void write_atom_data(int, double *);
  void write_node_data(int, double *);
  void write_elem_data(int, double *);
  int count_elements();
  int count_nodes(); 

  typedef void (DumpTecPlot:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  FnPtrHeader atom_header_choice;                       // ptr to write header functions
  FnPtrHeader elem_header_choice;                       // ptr to write header functions
  FnPtrHeader node_header_choice;                       // ptr to write header functions
  void header_item(bigint);
  void atom_header_item(bigint);
  void node_header_item(bigint);
  void elem_header_item(bigint);

  typedef void (DumpTecPlot::*FnPtrWrite)(int,double *);
  FnPtrWrite write_atom_choice;
  FnPtrWrite write_elem_choice;
  FnPtrWrite write_node_choice;
  void write_string(int, double *);
  void write_lines_noimage(int, double *);
  void write_elem_lines_noimage(int, double *);


  typedef void (DumpTecPlot::*FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  FnPtrPack pack_node_choice;               // ptr to pack functions
  FnPtrPack pack_elem_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);
  void pack_elem_noscale_noimage(tagint *);
  void pack_node_noscale_noimage(tagint *);

  typedef int (DumpTecPlot::*FnPtrConvert)(int,int,double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  FnPtrConvert convert_elem_choice;          // ptr to convert data functions
  FnPtrConvert convert_node_choice;          // ptr to convert data functions
  int convert_atom_noimage(int,int,double *);
  int convert_node_noimage(int,int,double *);
  int convert_elem_noimage(int,int,double *);

};
}

#endif
#endif
