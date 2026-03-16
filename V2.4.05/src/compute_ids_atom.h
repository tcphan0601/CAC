#ifdef COMPUTE_CLASS

ComputeStyle(ids/atom,ComputeIDSAtom)

#else

#ifndef CAC_COMPUTE_IDS_ATOM_H
#define CAC_COMPUTE_IDS_ATOM_H

#include "compute.h"

namespace CAC_NS {

class ComputeIDSAtom : public Compute {
 public:
  ComputeIDSAtom(class CAC *, int, char **);
  ~ComputeIDSAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int namax; // max # of atoms locally
  int nemax; // max # of elements locally

  class NeighList *list;

  double *apattern;
  double **iapattern;

  double rscale;
  double rcutsq;          // cutoff for CNA structure calculation
  double cutsq;           // cutoff to find second nearest neighbor of ghost

  int *ibucket;
  int *indexbucket;
  double *rbucket;
  int maxbucket;

  void grow_atom(int);
  void grow_intpl(int, int);
  void grow_bucket();
};

}

#endif
#endif

