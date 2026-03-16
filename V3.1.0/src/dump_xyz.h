#ifdef DUMP_CLASS

DumpStyle(xyz, DumpXYZ)

#else

#ifndef CAC_DUMP_XYZ_H
#define CAC_DUMP_XYZ_H

#include "dump.h"

namespace CAC_NS {

class DumpXYZ : public Dump {
 public:
  DumpXYZ(CAC *, int, char **);
  ~DumpXYZ();

 protected:
  int ntypes;
  char **typenames;
  int dump_element_style;
  int wrap_flag;

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
  int modify_param(int, char **);


  typedef void (DumpXYZ:: * FnPtrWrite)(int, double *);
  FnPtrWrite write_atom_choice;
  void write_string(int, double *);
  void write_atom_lines(int, double *);

};
}

#endif
#endif
