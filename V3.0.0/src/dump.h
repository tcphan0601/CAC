#ifndef CAC_DUMP_H
#define CAC_DUMP_H

#include <mpi.h>
#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  char *filename;            // user-specified file
  int igroup, groupbit;       // group that Dump is performed on

  int first_flag;            // 0 if no initial dump, 1 if yes initial dump
  int clearstep;             // 1 if dump invokes computes, 0 if not
  int create_attribute;          // 1 if dump stores attributes that need

  int comm_atom_forward;         // size of forward communication (0 if none)
  int comm_atom_reverse;         // size of reverse communication (0 if none)
  int comm_elem_forward;         // size of forward communication (0 if none)
  int comm_elem_reverse;         // size of reverse communication (0 if none)

  static Dump *dumpptr;         // holds a ptr to Dump currently being used

  Dump(class CAC *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write();

  virtual int pack_atom_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_atom_forward_comm(int, int, double *) {}
  virtual int pack_atom_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_atom_reverse_comm(int, int *, double *) {}
  virtual int pack_elem_forward_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_elem_forward_comm(int, int, double *) {}
  virtual int pack_elem_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_elem_reverse_comm(int, int *, double *) {}


  void modify_params(int, char **);
  virtual bigint memory_usage();

  virtual void set_atom_arrays(int) {}
  virtual void set_elem_arrays(int) {}
  virtual void pass_on_atom_arrays(int, int, int) {}
  virtual void pass_on_elem_arrays(int, int, int *) {}


 protected:
  int me, nprocs;             // proc info

  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 

  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  char *multiname;           // filename with % converted to cluster ID
  MPI_Comm clustercomm;      // MPI communicator within my cluster of procs

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int sort_flag;             // 1 if sorted output
  int append_flag;           // 1 if open file in append mode, 0 if not
  int buffer_allow;          // 1 if style allows for buffer_flag, 0 if not
  int buffer_flag;           // 1 if buffer output as one big string, 0 if not
  int padflag;               // timestep padding in filename
  int pbcflag;               // 1 if remap dumped atoms via PBC, 0 if not
  int singlefile_opened;     // 1 = one big file, already opened, else 0
  int sortcol;               // 0 to sort on ID, 1-N on columns
  int sortcolm1;             // sortcol - 1
  int sortorder;             // ASCEND or DESCEND

  char boundstr[9];          // encoding of boundary flags
  char *columns;
  char *id_stress;
  char *id_displace;
  char *id_pe;
  char *format;
  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format_elem;              // format string for the file write
  char *format_atom;              // format string for the file write
  char *format_node;              // format string for the file write
  char *format_stress;
  char *format_velocity;
  char *format_force;
  char *format_displace;
  char *format_pe;
  char *format_owner_proc;
  char *format_element;      // format string for element info
  char *format_line_user;    // user-specified format strings
  char *format_float_user;
  char *format_int_user;
  char *format_bigint_user;
  char **format_column_user;


  FILE *fp;                  // file to write dump to
  int atom_size_one;              // # of quantities for one atom
  int node_size_one;              // # of quantities for one node 
  int elem_size_one;              // # of quantities for one element in node connectivity section
  int max_size_one;     
  int nucell_me;
  int natom_me;                   // # of atoms in this dump from me
  int nelem_me;                  // # of elements in this dump fror me
  int nnode_me;                  // # of nodes in this dump fror me
  int nselem_me;                 // # of chars in string output in node connectivity section from me
  int nsnode_me;                 // # of chars in string output in node section from me
  int nsatom_me;                  // # of chars in string output in atom section from me

  double boxxlo, boxxhi;      // local copies of domain values
  double boxylo, boxyhi;      // lo/hi are bounding box for triclinic
  double boxzlo, boxzhi;
  double boxxy, boxxz, boxyz;

  int maxfiles;              // max number of files created, -1 == infinite
  int numfiles;              // number of files in names list
  int fileidx;               // index of file in names list
  char **nameslist;          // list of history file names

  bigint natom_total;             // total # of per-atom lines in snapshot
  bigint nelem_total;            // total # of element lines in snapshot
  bigint nnode_total;            // total # of node lines in snapshot
  int reorderflag;           // 1 if OK to reorder instead of sort
  int ntotal_reorder;        // # of atoms that must be in snapshot
  int nme_reorder;           // # of atoms I must own in snapshot
  tagint idlo;               // lowest ID I own when reordering
  tagint max_atom_tag;

  int maxbuf;                // size of buf
  double *buf;               // memory for atom/element/node quantities
  int maxsbuf;               // size of sbuf
  char *sbuf;                // memory for atom/element/node quantities in string format
  int maxibuf;
  int *ibuf;                 // memory for atom/element/node quantities  


  int maxids;                // size of ids
  int maxsort;               // size of bufsort, idsort, index
  int maxproc;               // size of proclist
  tagint *ids;               // list of atom IDs, if sorting on IDs
  double *bufsort;
  tagint *idsort;
  int *index, *proclist;
  
  virtual void init_style() = 0;
  virtual void openfile();
  virtual int modify_param(int, char **) {return 0;}
  virtual void write_header(bigint) = 0;
  virtual void write_atom_header(bigint) = 0;
  virtual void write_elem_header(bigint) = 0;
  virtual void write_node_header(bigint) = 0;

  // count number of local lines for each section
  
  virtual int count_atoms();
  virtual int count_elements();
  virtual int count_nodes();

  virtual void pack_atom(tagint *) = 0;
  virtual void pack_elem(tagint *) = 0;
  virtual void pack_node(tagint *) = 0;

  virtual int convert_atom_string(int, int, double *) {return 0;}
  virtual int convert_elem_string(int, int, double *) {return 0;}
  virtual int convert_node_string(int, int, double *) {return 0;}

  virtual void write_atom_data(int, double *) = 0;
  virtual void write_node_data(int, double *) = 0;
  virtual void write_elem_data(int, double *) = 0;


};
}

#endif
