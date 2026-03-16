#ifdef COMMAND_CLASS

CommandStyle(coarse_graining,CoarseGraining)

#else

#ifndef CAC_COARSE_GRAINING_H
#define CAC_COARSE_GRAINING_H

#include <stdio.h>
#include "pointers.h"

namespace CAC_NS {

class CoarseGraining : protected Pointers {
 public:
  CoarseGraining(class CAC *);
  void command(int, char **);

 private:
  int compute_style;
  int structure_type;
  class Compute *compute;
  char *id_compute;
  int user_compute;
  double cut,cutsq;

  class AtomVec *avec;
  class ElementVec *evec;

  int ngrains;
  int grain_tag_enable;
  int *grain_tag;
  int *seed_ids;
  int **seed_neigh_ids;
  int grain_etype;               // etype for each grain
  int grain_ctype;               // ctype for each grain
  double a1[3],a2[3],a3[3];
  double tol,mag_tol;;           // tol = angle tolerant
                                 // mag_tol = magnitude tolerant
  double *pattern;
  
  int nbasis;
  double **basis;
  int inum;
  int *ilist,*numneigh,**firstneigh;

  int ncellx,ncelly,ncellz;
  int maxneigh_user;
  int maxneigh;
  
  int ****atom_list;
  double ***nodecoord;

  int *select_flag;
  int centroid_count;
  double **centroid;
  int style;
  
  int migrate_flag,compress_flag,wrapx,wrapy,wrapz;
  int me;

  
  void options(int, char **);
  void coarse_grain(int, int, double *, double *, 
      double *, int, int, int, int first_flag=0);
  int collect_atoms(int, double *, double *, double *, 
      int, int, int);
  int align(int, int, double *, double);
  void add_basis(double, double, double);
  int check_neighbors(int, int);
  int is_neighbor(int, int);
  void read_centroid(char *);


  // debug function and variables
/*  
  FILE *fp;
  
  int index;
  int grain_test;
  void header();
*/
};

}

#endif
#endif

