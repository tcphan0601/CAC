#ifdef PAIR_CLASS

PairStyle(meam, PairMEAM)

#else

#ifndef CAC_PAIR_MEAM_H
#define CAC_PAIR_MEAM_H

#include "pair.h"

namespace CAC_NS {

class PairMEAM : public Pair {
 public:
  PairMEAM(class CAC *);
  ~PairMEAM();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  double init_one_inner(int, int);
  virtual void *extract(const char *, int &);

  int pack_atom_forward_comm(int, int *, double *, int, int *);
  void unpack_atom_forward_comm(int, int, double *);
  int pack_atom_reverse_comm(int, int, double *);
  void unpack_atom_reverse_comm(int, int *, double *);
  int pack_elem_forward_comm(int, int *, double *, int, int *);
  void unpack_elem_forward_comm(int, int, double *);
  int pack_elem_reverse_comm(int, int, double *);
  void unpack_elem_reverse_comm(int, int *, double *);

  double memory_usage();

 private:
  class MEAM *meam_inst;
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique chemical elements (not to be confused with FEs)
  char **elements;              // names of unique chemical elements (not to be confused with FEs)
  double *mass;                 // mass of each element

  int *map;                     // mapping from atom types (1-indexed) to elements (1-indexed)
  double **scale;               // scaling factor for adapt



  void allocate();
  void read_files(char *, char *);
  void neigh_strip(int, int *, int *, int **);
};

}

#endif
#endif


