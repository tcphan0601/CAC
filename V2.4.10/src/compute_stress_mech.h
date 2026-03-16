#ifdef COMPUTE_CLASS

ComputeStyle(stress/mech,ComputeStressMech)

#else

#ifndef CAC_COMPUTE_STRESS_MECH_H
#define CAC_COMPUTE_STRESS_MECH_H

#include "compute.h"
#include "my_page.h"
#include "nbin.h"
#include "nstencil.h"

namespace CAC_NS {

class ComputeStressMech : public Compute {
 public:
  ComputeStressMech(class CAC *, int, char **);
  ~ComputeStressMech();
  void init();
  void compute_array();
  double compute_scalar();
  double memory_usage();

 private:
  double lo[3],hi[3];               // bound of box to calculate mechanical stress
  double delta[3];                  // size of surface and distance between surfaces
  int nsurface[3];                        // number of surface in each direction (number of cells)
  int surface_arrange_style;        // 1 = traditional (default)
                                    // add more to include more styles
  int triclinic;

  double cutneighsq[3];
  
  int surface_direction;            // 0, 1, or 2 = x, y, or z 
  double **surface;                 // information about surfaces 0 1 2 -> x y z position of surface center
                                    // 3 4 5 -> stress in x y z direction
  double surface_area;              // area of a surface
  double xsurface[3];               // position of surface within cell, must be between 0 and 1 

  // local surface neighbor list 
  // surface above neighbors pages (ab)

  MyPage<int> *a_ab_page;           // pages of surface to atom neighbor indices
  MyPage<int> *ia_ab_page;          // pages of surface to element neighbor indices
  MyPage<int> *ia_index_ab_page;    // pages of surface to interpolated atom

  int *numneigh_a_ab;               // # of J atom neighbors to surface I
  int *numneigh_ia_ab;              // # of J interpolated atom neighbors to surface I
  int **firstneigh_a_ab;            // ptr to array (size numneigh_a) of J atom tag of each I surface
  int **firstneigh_ia_ab;           // ptr to array (size numneigh_ia) of J interpolated atom's parent element tag of each I surface
  int **firstneigh_ia_index_ab;     // ptr to array (size numneigh_ia) of J interpolated atom index within its parent element of each I surface

                                    // neighbor's parent element indices
  
  // surface below neighbors pages (bl)
  
  MyPage<int> *a_bl_page;           // pages of surface to atom neighbor indices
  MyPage<int> *ia_bl_page;          // pages of surface to element neighbor indices
  MyPage<int> *ia_index_bl_page;    // pages of surface to interpolated atom
                                    // neighbor's parent element indices

  int *numneigh_a_bl;               // # of J atom neighbors to surface I
  int *numneigh_ia_bl;              // # of J interpolated atom neighbors to surface I
  int **firstneigh_a_bl;            // ptr to array (size numneigh_a) of J atom tag of each I surface
  int **firstneigh_ia_bl;           // ptr to array (size numneigh_ia) of J interpolated atom's parent element tag of each I surface
  int **firstneigh_ia_index_bl;     // ptr to array (size numneigh_ia) of J interpolated atom index within its parent element of each I surface

  int pgsize;                     // size of each page
  int oneatom;                    // max size for one atom
  int oneatom_user;
  int snum;                       // number of surfaces

  int is_crossing(double *, double *, double *);
  void setup_surfaces();
  void allocate();
  void surface_neighbor();
  
};

}

#endif
#endif
