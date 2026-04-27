#ifdef COMPUTE_CLASS

ComputeStyle(ptm/atom, ComputePTMAtom)

#else

#ifndef CAC_COMPUTE_PTM_ATOM_H
#define CAC_COMPUTE_PTM_ATOM_H

#include "compute.h"
#include <cstdint>

// per-particle output column indices
#define PTM_COL_TYPE   0
#define PTM_COL_RMSD   1
#define PTM_COL_IDIST  2  // interatomic distance
#define PTM_COL_QUATW  3
#define PTM_COL_QUATX  4
#define PTM_COL_QUATY  5
#define PTM_COL_QUATZ  6
#define PTM_NCOLS      7

namespace CAC_NS {

class ComputePTMAtom : public Compute {
 public:
  ComputePTMAtom(class CAC *, int, char **);
  ~ComputePTMAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

  virtual int pack_atom_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_atom_forward_comm(int, int, double *);
  virtual int pack_elem_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_elem_forward_comm(int, int, double *);

 private:
  int namax, nemax;
  int32_t input_flags;   // bitmask: which structure types to check
  double rmsd_threshold; // RMSD cutoff (INFINITY = no cutoff)
  int sortnum;           // neighbor list sort depth (set from input_flags)
  int inner_flag;        // 0 = outer ucells only (default), 1 = all ucells

  class NeighList *list;

  // per-particle output: [natom][PTM_NCOLS] and [nelem][apc][ucell][PTM_NCOLS]
  double **atompattern;
  double ****vatompattern;

  void grow_atom(int);
  void grow_vatom(int, int, int);
};

}

#endif
#endif

/*  ERROR/WARNING messages:

E: Illegal compute ptm/atom command

Syntax: compute ID group ptm/atom <struct-string> <rmsd_threshold> [inner yes/no]
  struct-string: hyphen-separated list of fcc/hcp/bcc/ico/sc/dcub/dhex/graphene,
                 or preset keyword "default" (fcc-hcp-bcc-ico) or "all".
  rmsd_threshold: RMSD cutoff (0 = no threshold = accept any match).
  inner: yes = compute PTM for all virtual atoms including interior ucells;
         no (default) = outer (surface) ucells only.

E: Compute ptm/atom requires a pair style be defined

Self-explanatory.

E: PTM global initialization failed

PTM library failed to initialize its graph data structures.

W: More than one compute ptm/atom defined

It is not efficient to use compute ptm/atom more than once.

*/
