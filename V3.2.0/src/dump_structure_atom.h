#ifdef DUMP_CLASS

DumpStyle(structure/atom, DumpStructureAtom)

#else

#ifndef CAC_DUMP_STRUCTURE_ATOM_H
#define CAC_DUMP_STRUCTURE_ATOM_H

#include "dump.h"

namespace CAC_NS {

class DumpStructureAtom : public Dump {
 public:
  DumpStructureAtom(CAC *, int, char **);

 protected:
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int allflag;

  int structure_style;
  char *id_compute;
  class Compute *struct_compute; // ptr to compute structure/atom class to perform structure calculation
  int nkeep;                     // number of structure to keep for output
  double keep[6];                // arrays to store the structure types to keep for output

  void init_style();
  void write_header(bigint);
  void write_atom_header(bigint);
  void write_elem_header(bigint);
  void write_node_header(bigint);
  void pack_atom(tagint *);
  void pack_elem(tagint *);
  void pack_node(tagint *);

  int count_atoms();
  int count_elements() { return 0; }
  int count_nodes() { return 0; }
  int convert_atom_string(int, int, double *);
  int convert_node_string(int, int, double *) { return 0; }
  int convert_elem_string(int, int, double *) { return 0; }
  void write_atom_data(int, int, double *);
  void write_node_data(int, int, double *) {}
  void write_elem_data(int, int, double *) {}

  typedef void (DumpStructureAtom:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  void header_item(bigint);
  void header_item_triclinic(bigint);

  typedef void (DumpStructureAtom:: * FnPtrWrite)(int, int, double *);
  FnPtrWrite write_atom_choice;
  void write_string(int, int, double *);
  void write_atom_lines_noimage(int, int, double *);

  typedef void (DumpStructureAtom:: * FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);

  typedef int (DumpStructureAtom:: * FnPtrConvert)(int, int, double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  int convert_atom_noimage(int, int, double *);

};
}

#endif
#endif
