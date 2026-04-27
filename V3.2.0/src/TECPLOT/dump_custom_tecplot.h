#ifdef DUMP_CLASS

DumpStyle(custom/tecplot, DumpCustomTecplot)

#else

#ifndef CAC_DUMP_CUSTOM_TECPLOT_H
#define CAC_DUMP_CUSTOM_TECPLOT_H

#include "dump_custom.h"

namespace CAC_NS {

// Replaces dump tecplot. Outputs atoms + CG element FE zones using the binary
// TECIO API (or ASCII .dat), with user-specified columns instead of fixed flags.
//
// Syntax:
//   dump ID group custom/tecplot N file [format plt|szplt|ascii] [average yes|no] col1 col2 ...
//
// Defaults: format detected from filename extension (.plt/.szplt/.dat); average yes.
// Atom zone  : discrete real atoms, columns applied directly.
// Element zone: CG nodes with connectivity; compute columns read from vatom storage.

class DumpCustomTecplot : public DumpCustom {
 public:
  DumpCustomTecplot(class CAC *, int, char **);
  ~DumpCustomTecplot();

  void write() override;
  void openfile() override;

 protected:
  int average_flag;   // 1 = average nodal values over bases; 0 = per-basis nodes
  int fileformat;     // 0=PLT, 1=SZPLT, 2=ASCII
  int npe_connect;    // 8 (3D) or 4 (2D) — connectivity entries per element

  // TECIO binary state
  char *tec_title;
  char *tec_values;   // space-separated column name string
  int  *valuelocation; // 1=nodal for all vars (size ncols)
  int   dummy, dummy1;
  double soltime;
  int   debug_tec;
  int   zonetype;

  // element connectivity buffers (reused across writes)
  int maxelem_buf;
  int **nodecell_ids;   // [maxelem][npe_connect] — averaged connectivity
  int ***node_ids;      // [maxelem][maxapc][npe_connect] — per-basis connectivity
  // ibuf and maxibuf are inherited from Dump base class

  void init_style() override;

  // atom zone (real atoms only — no vatoms)
  void pack_atom(tagint *) override;
  int  count_atoms() override;

  // node zone (FE nodes of CG elements)
  void pack_node(tagint *) override;
  int  count_nodes() override;

  // element connectivity zone
  void pack_elem(tagint *) override;
  int  count_elements() override;

  // header / write stubs (only used in ASCII path via base Dump::write())
  void write_header(bigint) override;
  void write_atom_header(bigint) override;
  void write_node_header(bigint) override;
  void write_elem_header(bigint) override {}
  int  convert_atom_string(int, int, double *) override;
  int  convert_node_string(int, int, double *) override;
  int  convert_elem_string(int, int, double *) override;
  void write_atom_data(int, double *) override;
  void write_node_data(int, double *) override;
  void write_elem_data(int, double *) override;

  // TECIO helpers
  void comm_buf_tecio(int nme);
  void grow_buf(int nme);
  void write_box_zone();

  // per-column binary pack helpers
  void pack_atom_col(int k, tagint *ids);   // column k, real atoms
  void pack_node_col(int k, tagint *ids);   // column k, FE nodes

  // value accessor for FE nodes
  double get_node_val(int i, int j, int inode, const ColDef &cd);
};

} // namespace CAC_NS

#endif
#endif
