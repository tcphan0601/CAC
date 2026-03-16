#ifndef CAC_NSTENCIL_H
#define CAC_NSTENCIL_H

#include "pointers.h"

namespace CAC_NS {

class NStencil : protected Pointers {
 public:
  int istyle;                      // 1-N index into stencilnames
  class NBin *nb;                  // ptr to NBin instance I depend on
  bigint last_stencil;             // last timestep stencil was created

  // atom stencil
  
  int natomstencil;                // # of bins in atomstencil
  int *atomstencil;                // list of bin offsets in atomstencil
  int **atomstencilxyz;                // bin offsets in xyz dims
  int *natomstencil_multi;             // # bins in each type-based multi stencil
  int **atomstencil_multi;             // list of bin offsets in each stencil
  double **atomdistsq_multi;           // sq distances to bins in each stencil

  // element stencil
  
  int nelemstencil;                // # of bins in elemstencil
  int *elemstencil;                // list of bin offsets in elemstencil
  int **elemstencilxyz;                // bin offsets in xyz dims
  int *nelemstencil_multi;             // # bins in each type-based multi stencil
  int **elemstencil_multi;             // list of bin offsets in each stencil
  double **elemdistsq_multi;           // sq distances to bins in each stencil

  double cutoff_custom;            // cutoff set by requestor
  int atomonly;

  NStencil(class CAC *);
  virtual ~NStencil();
  void post_constructor(class NeighRequest *);
  void copy_neighbor_info();
  virtual void create_setup();
  bigint memory_usage();

  virtual void create() = 0;

  inline int get_maxatomstencil() {return maxatomstencil;}
  inline int get_maxelemstencil() {return maxelemstencil;}

 protected:

  // data from Neighbor class

  int neighstyle;
  double cutneighmaxatom;
  double cutneighmaxatomsq;
  double cutelemx, cutelemy, cutelemz;

  double *cuttypesq;

  // data from NBin class

  int matombinx, matombiny, matombinz;
  double atombinsizex, atombinsizey, atombinsizez;
  double atombininvx, atombininvy, atombininvz;

  int melembinx, melembiny, melembinz;
  double elembinsizex, elembinsizey, elembinsizez;
  double elembininvx, elembininvy, elembininvz;

  // data common to all NStencil variants

  int xyzflag;                     // 1 if stencilxyz is allocated
  int maxatomstencil;              // max size of atom stencil
  int maxatomstencil_multi;        // max sizes of atom stencils
  int atomsx, atomsy, atomsz;        // extent of atom stencil in each dim

  int maxelemstencil;              // max size of element stencil
  int maxelemstencil_multi;        // max sizes of element stencils
  int elemsx, elemsy, elemsz;        // extent of element stencil in each dim

  int dimension;

  // methods for all NStencil variants

  void copy_bin_info();                     // copy info from NBin class
  double atombin_distance(int, int, int);   // distance between atom bin corners
  double elembin_distance(int, int);        // distance between element bin in a direction
};

}

#endif
