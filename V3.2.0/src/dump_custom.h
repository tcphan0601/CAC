#ifdef DUMP_CLASS

DumpStyle(custom, DumpCustom)

#else


#ifndef CAC_DUMP_CUSTOM_H
#define CAC_DUMP_CUSTOM_H

#include "dump.h"

namespace CAC_NS {

class DumpCustom : public Dump {
 public:
  DumpCustom(class CAC *, int, char **);
  ~DumpCustom();
  void write();
  void openfile();

 private:
  enum DumpStyle  { 
    PLT = 0, 
    SZPLT = 1, 
    ASCII = 2, 
    FULLMAP = 3, 
    SURFACE = 4, 
    GAUSSPOINT = 5, 
    NODE = 6, 
    CENTER = 7 
    };
  enum DumpFormat { LAMMPS, TECPLOT };

  int dump_format;  // DumpFormat enum: default LAMMPS
  int dump_style;   // DumpStyle enum; default FULLMAP
  int average_flag;
  int wrap_flag;
  int precision_flag;
  int current_col;
  enum ColType {
    COL_ID, 
    COL_TYPE,
    COL_X,  COL_Y,  COL_Z,
    COL_VX, COL_VY, COL_VZ,
    COL_FX, COL_FY, COL_FZ,
    COL_IX, COL_IY, COL_IZ,
    COL_XU, COL_YU, COL_ZU,
    COL_XU_TRI, COL_YU_TRI, COL_ZU_TRI,
    COL_PROC, // MPI rank that owns this atom/vatom
    COL_PROCP1, // MPI rank + 1 that owns this atom/vatom
    COL_MASS,
    COL_COMPUTE_SCALAR,   // c_ID   (scalar per-atom compute)
    COL_COMPUTE_ARRAY,    // c_ID[N] (array per-atom compute, 0-based col)
    COL_FIX_SCALAR,       // f_ID   (scalar per-atom fix)
    COL_FIX_ARRAY,        // f_ID[N] (array per-atom fix, 0-based col)
    // debug columns — virtual atoms only; real atoms always output -1
    COL_EID,        // element global tag
    COL_ETYPE,      // element type index
    COL_IBASIS,     // basis index within element
    COL_IGAUSS,     // gauss-point index (GAUSSPOINT style only; else -1)
    COL_INODE,      // node index (NODE style only; else -1)
    COL_IUCELL,     // ucell index (FULLMAP/OUTER styles only; else -1)

  };

  struct ColDef {
    ColType type;
    class Compute *compute;  // non-null for COL_COMPUTE_*
    class Fix     *fix;      // non-null for COL_FIX_*
    int col;                 // 0-based column index for COL_COMPUTE_ARRAY / COL_FIX_ARRAY
    int vtype;               //
    char *col_label;
  };

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

  int nevery;
  int ncols;
  ColDef  *coldefs;

  // referenced computes (c_ID / c_ID[N] columns)
  int ncompute;
  char           **id_compute;
  class Compute  **computes;

  // referenced fixes (f_ID / f_ID[N] columns)
  int nfix;
  char      **id_fix;
  class Fix **fixes;


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
  void write_atom_data(int, int, double *);
  void write_node_data(int, int, double *);
  void write_elem_data(int, int, double *);
  int count_elements();
  int count_nodes(); 
  int count_atoms();

  // ---- value accessors ----
  // Value for real atom i, column cd.
  double get_atom_val(int i, const ColDef &cd);
  // Value for virtual atom (element i, basis ibasis, ucell iucell).
  //   igcell   = gauss-cell index for this point (-1 if not a gauss point).
  //   inode    = node index for this point (-1 if not a node).
  //   vatom_id = virtual atom ID (for COL_ID).
  //   coord    = interpolated position [3] (for COL_X/Y/Z).
  double get_vatom_val(int i, int ibasis, int iucell, int igcell, int inode,
                       tagint vatom_id, const double *coord, const ColDef &cd);
  // value accessor for FE nodes
  double get_node_val(int i, int j, int inode, const ColDef &cd);


  typedef void (DumpCustom:: *FnPtrHeader) (bigint);
  FnPtrHeader header_choice;                       // ptr to write header functions
  void header_item_tecplot(bigint); 
  void header_item_lammps(bigint);
  void header_item_lammps_triclinic(bigint);

  FnPtrHeader atom_header_choice;                       // ptr to write header functions
  void atom_header_item_tecplot(bigint);

  FnPtrHeader elem_header_choice;                       // ptr to write header functions
  void node_header_item_tecplot(bigint);

  FnPtrHeader node_header_choice;                       // ptr to write header functions
  void empty_header(bigint) {}

  typedef void (DumpCustom:: * FnPtrWrite)(int, int,  double *);
  FnPtrWrite write_atom_choice;
  FnPtrWrite write_elem_choice;
  FnPtrWrite write_node_choice;
  void write_string(int, int, double *);
  void write_lines(int, int, double *);
  void write_lines_node_connect(int, int, double *);
  void write_empty(int, int, double *) {}


  typedef void (DumpCustom:: * FnPtrPack)(tagint *);
  FnPtrPack pack_atom_choice;               // ptr to pack functions
  FnPtrPack pack_node_choice;               // ptr to pack functions
  FnPtrPack pack_elem_choice;               // ptr to pack functions
  void pack_atom_lammps(tagint *);
  void pack_atom_tecplot(tagint *);
  void pack_elem_tecplot(tagint *);
  void pack_node_tecplot(tagint *);
  void pack_empty(tagint *) {}

  typedef int (DumpCustom:: * FnPtrConvert)(int, int, double *);
  FnPtrConvert convert_atom_choice;          // ptr to convert data functions
  FnPtrConvert convert_elem_choice;          // ptr to convert data functions
  FnPtrConvert convert_node_choice;          // ptr to convert data functions
  int convert_string(int, int, double *);
  int convert_string_node_connect(int, int, double *);
  int convert_empty(int, int, double *) { return 0; }

#ifdef TECIO_ENABLED        
  void pack_atom_tecplot_binary(tagint *);
  void pack_elem_tecplot_binary(tagint *);
  void pack_node_tecplot_binary(tagint *);
  void comm_buf_tecio(int);
#endif

  void grow_buf(int);
  // Parse column keywords from arg[istart..narg-1]; skip anything already consumed.
  void parse_columns(int istart, int narg, char **arg);
  void write_box_zone();
};

} // namespace CAC_NS

#endif
#endif
