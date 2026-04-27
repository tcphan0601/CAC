#ifdef COMMAND_CLASS

CommandStyle(atomic_strain, AtomicStrain)

#else

#ifndef CAC_ATOMIC_STRAIN_H
#define CAC_ATOMIC_STRAIN_H

#include "stdio.h"
#include "pointers.h"

namespace CAC_NS {

class AtomicStrain : protected Pointers {
 public:
  AtomicStrain(class CAC *);
  ~AtomicStrain();
  void command(int, char **);

 private:
  int me, nprocs;
  char *line, *keyword, *buffer, *style;
  FILE *fp;
  char **arg;
  int narg, maxarg;

  int *atomValidity, ***nodeValidity;
  double **atomD, *atomVonMisesStrain, *atomVolStrain, **atomF, **atomQ, **atomW, **atomStretchTensor, **atomStrainTensor, *atomNSD;
  double ****nodeD, ***nodeVonMisesStrain, ***nodeVolStrain, ****nodeF, ****nodeQ, ****nodeW, ****nodeStretchTensor, ****nodeStrainTensor, ***nodeNSD;

  class NeighList *list;

  bigint natoms;
  bigint nelements;
  bigint nnodes;
  int natypes;
  int netypes;

  MPI_Status status;
  MPI_Request request;

  class Balance *balance;

  double boxlo[3], boxhi[3];
  double xy, xz, yz;
  double boxlo_current[3], boxhi_current[3];
  double xy_current, xz_current, yz_current;

  int triclinic;

  int noutputs;
  int *output_format_list;
  int *output_configuration_list;
  char **output_filename_list;
  int element_info_flag;
  int deformation_gradient_flag;
  int displacement_flag;
  int rotation_quat_flag;
  int rotation_vect_flag;
  int nsd_flag;
  int strain_tensor_flag;
  int stretch_tensor_flag;
  int wrap_flag;
  int average_flag;
  int igroup, groupbit;
  char *region_id;
  int ntimestep;
  int no_elem_flag, no_atom_flag;
  double cut;

  // read input functions
  void read_input(char *);
  void read_ref_file(char *);
  void open_input(char *);
  void header();
  void header_ref();
  void parse_keyword(int);
  void skip_lines(bigint);
  int style_match(const char *, const char *);
  void atoms();
  void atoms_ref();
  void elements();
  void nodes();
  void nodes_ref();
  void element_types();

  // write output functions/parameters

  int size_one;
  double *buf;
  int *ibuf;
  int maxbuf, maxibuf;
 
  void write_output(char *, int, int);
  void write_header_tecplot(char *, int);
  void write_header_atom(char *, int, bigint);
  void write_tecplot_atoms(int, int);
  void write_tecplot_element_connectivity(int);
  void write_tecplot_nodes(int, int);
  void open_output(char *);

  void pack_atom_all(int, int);
  void pack_atom_tecplot(int, int);
  void pack_node_tecplot(int, int);
  void pack_element_connectivity();
  void write_lines_tecplot(int, double *);
  void write_lines_atom(int, double *);
#ifdef TECPLOT_ENABLE
  void comm_buf_tecio(int);
#endif

  void grow_buf(int);
  void grow_ibuf(int);

  // parameters for TECIO
  
  int *valuelocation;
  int debug, zonetype;
  double soltime;
  int fileformat; // 0 == PLT, 1 == SZPLT, 2 == ASCII, 3 == ATOM
  int dummy, dummy1;
  int ival;            // value counter
  int valueflag;

  // strain calculation functions
  
  void compute();
  int compute_strain(int, int *, int *, double *, double *, 
      double *, double *, double *, double *, double *, 
      double *, double *, double *);

};
}

#endif
#endif

