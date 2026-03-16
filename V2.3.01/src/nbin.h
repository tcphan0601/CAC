#ifndef CAC_NBIN_H
#define CAC_NBIN_H

#include "pointers.h"

namespace CAC_NS {

class NBin : protected Pointers {
 public:
  int istyle;                      // 1-N index into binnames
  bigint last_bin;                 // last timestep objs were binned

  int nbinx,nbiny,nbinz;           // # of global bins
  int mbins;                       // # of local bins and offset on this proc
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int *atombinhead;                // index of first atom in each bin
  int *elembinhead;                // index of first element in each bin
  int *atombins;                   // index of next atom in same bin
  int *elembins;                   // index of next element in same bin
  int *atom2bin;               // bin assignment for each atom (local+ghost)
  int *elem2bin;               // bin assignment for each element (local+ghost)

  double cutoff_custom;        // cutoff set by requestor

  NBin(class CAC *);
  ~NBin();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  virtual void bin_setup(int,int);
  bigint memory_usage();

  virtual void setup_bins(int) = 0;
  virtual void bin_all() = 0;
  int coord2bin(double *);

 protected:

  // data from Neighbor class

  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int binsizeflag;
  double binsize_user;
  double *bboxlo,*bboxhi;

  // data common to all NBin variants

  int dimension;
  int triclinic;

  int maxbin;                       // size of binhead array
  int maxatom;                      // size of per-atom bins array
  int maxelem;                      // size of per-elem bins array

};

}

#endif

/* ERROR/WARNING messages:

*/
