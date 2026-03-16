#ifndef CAC_NBIN_H
#define CAC_NBIN_H

#include "pointers.h"

namespace CAC_NS {

class NBin : protected Pointers {
 public:
  int istyle;                      // 1-N index into binnames
  bigint last_bin;                 // last timestep objs were binned

  // atom bin info
 
  int natombinx,natombiny,natombinz;           // # of global atom bins
  int matombins;                               // # of local bins and offset on this proc
  int matombinx,matombiny,matombinz;
  int matombinxlo,matombinylo,matombinzlo;
  double atombinsizex,atombinsizey,atombinsizez;  // atom bin sizes and inverse sizes
  double atombininvx,atombininvy,atombininvz;

  int *atombinhead;                // index of first atom in each bin
  int *atombins;                   // index of next atom in same bin
  int *atom2bin;                   // bin assignment for each atom (local+ghost)

  // element bin info
 
  int nelembinx,nelembiny,nelembinz;           // # of global element bins
  int melembins;                       // # of local bins and offset on this proc
  int melembinx,melembiny,melembinz;
  int melembinxlo,melembinylo,melembinzlo;
  double elembinsizex,elembinsizey,elembinsizez;  // element bin sizes and inverse sizes
  double elembininvx,elembininvy,elembininvz;

  int *elembinhead;            // index of first element in each bin
  int *elembins;               // index of next element in same bin
  int *elem2bin;               // bin assignment for each element (local+ghost)

  double cutoff_custom;        // cutoff set by requestor
  int atomonly;

  NBin(class CAC *);
  ~NBin();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  virtual void bin_setup(int, int);
  bigint memory_usage();

  virtual void setup_bins(int) = 0;
  virtual void bin_all() = 0;
  virtual int coord2atombin(double, double, double);
  virtual int coord2elembin(double, double, double);

 protected:

  // data from Neighbor class

  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int atombinsizeflag;
  int elembinsizeflag;
  double atombinsize_user;
  double elembinsize_user;
  double *bboxlo,*bboxhi;

  // data common to all NBin variants

  int dimension;
  int triclinic;

  int maxatombin;                   // size of atombinhead array
  int maxelembin;                   // size of elembinhead array
  int maxatom;                      // size of per-atom bins array
  int maxelem;                      // size of per-elem bins array

};

}

#endif

/* ERROR/WARNING messages:

*/
