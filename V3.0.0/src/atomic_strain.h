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
  double *atomVonMisesStrain, *atomVolStrain, **atomF, **atomQ, **atomStretchTensor, **atomStrainTensor, *atomNSD;
  double ***nodeVonMisesStrain, ***nodeVolStrain, ****nodeF, ****nodeQ, ****nodeStretchTensor, ****nodeStrainTensor, ***nodeNSD;

  double **atom_displacement;
  double ****node_displacement;

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

  int output_format;
  int output_configuration;
  int element_info_flag;
  int deformation_gradient_flag;
  int rotation_flag;
  int nsd_flag;
  int strain_tensor_flag;
  int stretch_tensor_flag;
  int wrap_flag;
  int average_flag;
  double cut;

  // read input functions
  void read_input(char *);
  void open_input(char *);
  void header();
  void parse_keyword(int);
  void skip_lines(bigint);
  int style_match(const char *, const char *);
  void atoms();
  void elements();
  void nodes();
  void element_types();

  // write output functions/parameters

  int size_one;
  double *buf;
  int *ibuf;
  int maxbuf, maxibuf;
 
  void write_output(char *);
  void write_header_tecplot(char *);
  void write_header_atom(char *, bigint);
  void write_tecplot_atoms();
  void write_tecplot_element_connectivity();
  void write_tecplot_nodes();
  void open_output(char *);

  void pack_atom_all();
  void pack_atom_tecplot();
  void pack_node_tecplot();
  void pack_element_connectivity();
  void write_lines_tecplot(int, double *);
  void write_lines_atom(int, double *);
  void comm_buf_tecio(int);

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
      double *, double *);

};
}

#endif
#endif

