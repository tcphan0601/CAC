#ifdef DUMP_CLASS

DumpStyle(stress/mech, DumpStressMech)

#else

#ifndef CAC_DUMP_STRESS_MECH_H
#define CAC_DUMP_STRESS_MECH_H

#include "dump.h"

namespace CAC_NS {

class DumpStressMech : public Dump {
 public:
  DumpStressMech(CAC *, int, char **);

 protected:
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no

  class Compute *compute;        // ptr to compute stress/mech class

  void init_style();
  void write_header(bigint);
  void write_atom_header(bigint);
  void write_elem_header(bigint) {}
  void write_node_header(bigint) {}
  void pack_atom(tagint *);
  void pack_elem(tagint *) {}
  void pack_node(tagint *) {}

  int count_atoms();
  int count_elements() { return 0; }
  int count_nodes() { return 0; }
  int convert_atom_string(int, int, double *);
  int convert_node_string(int, int, double *) { return 0; }
  int convert_elem_string(int, int, double *) { return 0; }
  void write_atom_data(int, int, double *);
  void write_node_data(int, int, double *) {}
  void write_elem_data(int, int, double *) {}

  typedef void (DumpStressMech:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                // ptr to write header functions
  void header_item(bigint);

  typedef void (DumpStressMech:: * FnPtrWrite)(int, int, double *);
  FnPtrWrite write_atom_choice;
  void write_string(int, int, double *);
  void write_surfaces(int, int, double *);

  typedef void (DumpStressMech:: * FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  void pack_surfaces(tagint *);

  typedef int (DumpStressMech:: * FnPtrConvert)(int, int, double *);
  FnPtrConvert convert_atom_choice;         // ptr to convert data functions
  int convert_surfaces(int, int, double *);

};
}

#endif
#endif
