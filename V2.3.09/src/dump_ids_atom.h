#ifdef DUMP_CLASS

DumpStyle(ids/atom, DumpIDSAtom)

#else

#ifndef CAC_DUMP_IDS_ATOM_H
#define CAC_DUMP_IDS_ATOM_H

#include "dump.h"

namespace CAC_NS {

class DumpIDSAtom : public Dump {
 public:
  DumpIDSAtom(CAC *, int,char**);
  ~DumpIDSAtom();

 protected:
  int scale_flag;                // 1 if element coords are scaled, 0 if no
  int image_flag;                // 1 if append box count to element coords, 0 if no
  int allflag;                   // 1 if output all atoms
  char *id_ids;

  class Compute *compute_ids;    // ptr to compute ids/atom class to perform Identify IDS Structure calculation
  int nkeep;                     // number of structure to keep for output
  double keep[8];                   // arrays to store the structure types to keep for output

  void init_style();
  void write_header(bigint);
  void write_atom_header(bigint);
  void write_elem_header(bigint);
  void write_node_header(bigint);
  void pack_atom(tagint *);
  void pack_elem(tagint *);
  void pack_node(tagint *);

  int count_atoms();
  int count_intpl_atoms();
  int count_elements();
  int count_nodes();
  int convert_atom_string(int, int, double *);
  int convert_node_string(int, int, double *);
  int convert_elem_string(int, int, double *);
  void write_atom_data(int, double *);
  void write_node_data(int, double *);
  void write_elem_data(int, double *);

  typedef void (DumpIDSAtom:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  void header_item(bigint);

  typedef void (DumpIDSAtom::*FnPtrWrite)(int,double *);
  FnPtrWrite write_atom_choice;
  void write_string(int, double *);
  void write_atom_lines_noimage(int, double *);

  typedef void (DumpIDSAtom::*FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  void pack_atom_noscale_noimage(tagint *);

  typedef int (DumpIDSAtom::*FnPtrConvert)(int,int,double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  int convert_atom_noimage(int,int,double *);


};
}

#endif
#endif
