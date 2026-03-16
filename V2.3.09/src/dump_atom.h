#ifdef DUMP_CLASS

DumpStyle(atom, DumpAtom)

#else

#ifndef CAC_DUMP_ATOM_H
#define CAC_DUMP_ATOM_H

#include "dump.h"

namespace CAC_NS {

class DumpAtom : public Dump {
 public:
  DumpAtom(CAC *, int,char**);
  ~DumpAtom();

 protected:
  int nevery;
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int stress_flag;
  int force_flag;
  int displace_flag;
  class Compute *stress;
  class Compute *displace;
 
  void init_style();
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

  typedef void (DumpAtom::*FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  void header_item(bigint);

  typedef void (DumpAtom::*FnPtrWrite)(int,double *);
  FnPtrWrite write_atom_choice;
  FnPtrWrite write_elem_choice;
  void write_string(int, double *);
  void write_atom_lines_noimage(int, double *);
  void write_elem_lines_noimage(int, double *);


  typedef void (DumpAtom::*FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  FnPtrPack pack_elem_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);
  void pack_elem_noscale_noimage(tagint *);
 

  typedef int (DumpAtom::*FnPtrConvert)(int,int,double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  FnPtrConvert convert_elem_choice;          // ptr to convert data functions
  int convert_atom_noimage(int,int,double *);
  int convert_elem_noimage(int,int,double *);

};
}

#endif
#endif
