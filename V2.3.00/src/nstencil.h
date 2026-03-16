#ifndef CAC_NSTENCIL_H
#define CAC_NSTENCIL_H

#include "pointers.h"

namespace CAC_NS {

class NStencil : protected Pointers {
 public:
  int istyle;                      // 1-N index into binnames
  class NBin *nb;                  // ptr to NBin instance I depend on
  bigint last_stencil;             // last timestep stencil was created

  int object_type;             // ATOM or ELEMENT
  int nstencil;                    // # of bins in stencil
  int *stencil;                    // list of bin offsets
  int **stencilxyz;                // bin offsets in xyz dims
  int *nstencil_multi;             // # bins in each type-based multi stencil
  int **stencil_multi;             // list of bin offsets in each stencil
  double **distsq_multi;           // sq distances to bins in each stencil

  double cutoff_custom;            // cutoff set by requestor

  NStencil(class CAC *);
  virtual ~NStencil();
  void post_constructor(class NeighRequest *);
  void copy_neighbor_info();
  virtual void create_setup();
  bigint memory_usage();

  virtual void create() = 0;

  inline int get_maxstencil() {return maxstencil;}

 protected:

  // data from Neighbor class

  int neighstyle;
  double cutneighmax;
  double cutneighmaxsq;
  double *cuttypesq;

  // data from NBin class

  int mbinx,mbiny,mbinz;
  double binsizex,binsizey,binsizez;
  double bininvx,bininvy,bininvz;

  // data common to all NStencil variants

  int xyzflag;                     // 1 if stencilxyz is allocated
  int maxstencil;                  // max size of stencil
  int maxstencil_multi;            // max sizes of stencils
  int sx,sy,sz;                    // extent of stencil in each dim

  int dimension;

  // methods for all NStencil variants

  void copy_bin_info();                     // copy info from NBin class
  double bin_distance(int, int, int);       // distance between bin corners
};

}

#endif

/* ERROR/WARNING messages:

*/
