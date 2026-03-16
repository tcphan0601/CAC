#ifndef CAC_NPAIR_H
#define CAC_NPAIR_H

#include "pointers.h"
#include "my_page.h"

namespace CAC_NS {

class NPair : protected Pointers {
 public:
  int istyle;                   // 1-N index into pairnames
  class NBin *nb;               // ptr to NBin instance I depend on

  class NStencil *ns;           // ptr to NStencil instance I depend on

  bigint last_build;            // last timestep build performed

  double cutoff_custom;         // cutoff set by requestor
  int group;                    // if requested group
  int groupbit;                 // group to find neighbor set by requestor
  int force_rebuild;            // 1 if need to build regardless of last re-neighboring
  int custom_domain;            // 1 if not using domain size from domain class

  NPair(class CAC *);
  virtual ~NPair();
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  void build_setup();
  virtual void build(class NeighList *) = 0;

 protected:
  double **mycutneigh;           // per-type cutoffs when user specified
  double **mycutneighsq;         // per-type cutoffs when user specified

  // data from Neighbor class

  int exclude;
  double skin;
  double **cutneighsq;
  double **cutneighghostsq;
  double **cutneigh;

  double *bboxlo,*bboxhi;
  //double **cutneighghostsq;

  // filter 
  
  int atom_filter;
  int atom_filter_value;
  int atom_filter_preneighflag;
  class Compute *atom_filter_compute;

  // exclusion data from Neighbor class
  
  int nex_type;                    // # of entries in type exclusion list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  //int *special_flag;

  // data from NBin class
  
  int natombinx,natombiny,natombinz;
  int matombins;
  int matombinx,matombiny,matombinz;
  int matombinxlo,matombinylo,matombinzlo;
  double atombininvx,atombininvy,atombininvz;
  int *atom2bin,*atombins;
  int *atombinhead;
  
  int nelembinx,nelembiny,nelembinz;
  int melembins;
  int melembinx,melembiny,melembinz;
  int melembinxlo,melembinylo,melembinzlo;
  double elembininvx,elembininvy,elembininvz;
  int *elem2bin,*elembins;
  int *elembinhead;

  // data from NStencil class

  int natomstencil;
  int *atomstencil;
  int **atomstencilxyz;
  int *natomstencil_multi;
  int **atomstencil_multi;
  double **atomdistsq_multi;

  int nelemstencil;
  int *elemstencil;
  int **elemstencilxyz;
  int *nelemstencil_multi;
  int **elemstencil_multi;
  double **elemdistsq_multi;

  // methods for all NPair variants

  virtual void copy_bin_info();
  virtual void copy_stencil_info();

  int exclusion(int, int, int, int) const;

  int coord2atombin(double *, int &, int &, int &);
  int coord2elembin(double *, int &, int &, int &);

  // find_special: determine if atom j is in special list of atom i
  // if it is not, return 0
  // if it is and special flag is 0 (both coeffs are 0.0), return -1
  // if it is and special flag is 1 (both coeffs are 1.0), return 0
  // if it is and special flag is 2 (otherwise), return 1,2,3
  //   for which level of neighbor it is (and which coeff it maps to)
/*
  inline int find_special(const tagint *list, const int *nspecial,
                          const tagint tag) const {
    const int n1 = nspecial[0];
    const int n2 = nspecial[1];
    const int n3 = nspecial[2];

    for (int i = 0; i < n3; i++) {
      if (list[i] == tag) {
        if (i < n1) {
          if (special_flag[1] == 0) return -1;
          else if (special_flag[1] == 1) return 0;
          else return 1;
        } else if (i < n2) {
          if (special_flag[2] == 0) return -1;
          else if (special_flag[2] == 1) return 0;
          else return 2;
        } else {
          if (special_flag[3] == 0) return -1;
          else if (special_flag[3] == 1) return 0;
          else return 3;
        }
      }
    }
    return 0;
  };*/
};

}

#endif

