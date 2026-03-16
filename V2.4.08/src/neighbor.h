#ifndef CAC_NEIGHBOR_H
#define CAC_NEIGHBOR_H

#include "pointers.h"
#include <map>

namespace CAC_NS {

class Neighbor : protected Pointers {
  friend class ComputeStressMech;
 public:
  enum{NSQ,BIN,MULTI};            
  int style;                       // 0,1,2 = nsq, bin, multi
  int every;                       // build every this many steps
  int delay;                       // delay build for this many steps
  int dist_check;                  // 0 = always build, 1 = only if 1/2 dist
  int ago;                         // how many steps ago neighboring occurred
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  int build_once;                  // 1 if only build lists once per run

  double skin;                     // skin distance
  double cutneighmin;              // min neighbor cutoff for all type pairs
  double cutneighmax;              // max neighbor cutoff for all type pairs
  double cutneighmaxsq;            // cutneighmax squared
  double **cutneighsq;             // neighbor cutneigh sq for each type pair
  double **cutneigh;               // neighbor cutneigh for each type pair
  double **cutneighghostsq;        // cutneigh sq for each ghost type pair
  double *cuttype;                 // for each type, max neigh cut w/ others
  double *cuttypesq;               // cuttype squared

  // cutoffs for element neighboring
  double cut_all;
  double cut_elem;
  double cut_elemsq;
  double cut_subelem;
  double cut_subelemsq;

  int atombinsizeflag;             // user-chosen bin size for atoms
  int elembinsizeflag;             // user-chosen bin size for elements
  double atombinsize_user;         // set externally by some accelerator pkgs
  double elembinsize_user;         // set externally by some accelerator pkgs

  bigint ncalls;                   // # of times build has been called
  bigint ndanger;                  // # of dangerous builds
  bigint lastcall;                 // timestep of last neighbor::build() call

  // geometry and static info, used by other Neigh classes

  double *bboxlo,*bboxhi;          // ptrs to full domain bounding box
                                   // different for orthog vs triclinic
  double *zeroes;                  // vector of zeroes for shear history init

  // exclusion info, used by NeighPair

  int exclude;                     // 0 if no type/group exclusions, 1 if yes

  int nex_type;                    // # of entries in type exclusion list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  // special info, used by NeighPair

  int special_flag[4];             // flags for 1-2, 1-3, 1-4 neighbors

  // cluster setting, used by NeighTopo

  int cluster_check;               // 1 if check bond/angle/etc satisfies minimg

  // pairwise neighbor lists and corresponding requests

  int nlist;                           // # of pairwise neighbor lists
  int nrequest;                        // # of requests, same as nlist
  int old_nrequest;                    // # of requests for previous run

  class NeighList **lists;
  class NeighRequest **requests;       // from Pair,Fix,Compute,Command classes
  class NeighRequest **old_requests;   // copy of requests to compare to


  // public methods

  Neighbor(class CAC *);
  virtual ~Neighbor();
  virtual void init();
  int request(void *, int instance=0);
  int decide();                     // decide whether to build or not
  virtual int check_distance();     // check max distance moved since last build
  void setup_bins();                // setup bins based on box and cutoff
  void build_bins();                 
  virtual void build(int topoflag=1);  // build all perpetual neighbor lists
  void build_one(class NeighList *list, int preflag=0);
                                    // create a one-time pairwise neigh list
  void set(int, char **);           // set neighbor style and skin distance
  void reset_timestep(bigint);      // reset of timestep counter
  void modify_params(int, char**);  // modify params that control builds
  void exclusion_group_group_delete(int, int);  // rm a group-group exclusion
  int exclude_setting();            // return exclude value to accelerator pkg

  class NeighRequest *find_request(void *);  // find a neighbor request

  bigint memory_usage();

 protected:
  int me,nprocs;
  int firsttime;                   // flag for calling init_styles() only once

  int dimension;                   // 2/3 for 2d/3d
  int triclinic;                   // 0 if domain is orthog, 1 if triclinic
  int newton_pair;                 // 0 if newton off for pairwise, 1 if on

  int must_check;                  // 1 if must check other classes to reneigh
//  int restart_check;               // 1 if restart enabled, 0 if no
  int fix_check;                   // # of fixes that induce reneigh
  int *fixchecklist;               // which fixes to check

  bigint last_setup_bins;          // step of last neighbor::setup_bins() call

  double triggersq;                // trigger = build when atom moves this dist

  double **atomxhold;                  // atom coords at last neighbor build
  double ***nodexhold;                 // node coords at last neighbor build
  int atommaxhold;                     // size of atomxhold array
  int nodemaxhold;                     // size of nodexhold array

  int boxcheck;                        // 1 if need to store box size
  double boxlo_hold[3],boxhi_hold[3];  // box size at last neighbor build
  double corners_hold[8][3];           // box corners at last neighbor build
  double (*corners)[3];                // ptr to 8 corners of triclinic box


  int old_style,old_triclinic;     // previous run info
  int old_pgsize,old_oneatom;      // used to avoid re-creating neigh lists

  int nstencil_perpetual;         // # of perpetual NeighStencil classes
  int npair_perpetual;            // #x of perpetual NeighPair classes
  int *slist;                     // indices of them in neigh_stencil
  int *plist;                     // indices of them in neigh_pair
  int maxex_type;                  // max # in exclusion type list
  int maxex_group;                 // max # in exclusion group list
  int maxrequest;                  // max size of NeighRequest list
  int maxwt;                       // max weighting factor applied + 1

  // info for other Neigh classes: NBin,NStencil,NPair,NTopo

  int nbin,nstencil;
  int nbclass,nsclass,npclass;

  typedef class NBin *(*BinCreator)(class CAC *);
  BinCreator *binclass;
  char **binnames;
  int *binmasks;
  class NBin **neigh_bin;

  typedef class NStencil *(*StencilCreator)(class CAC *);
  StencilCreator *stencilclass;
  char **stencilnames;
  int *stencilmasks;
  class NStencil **neigh_stencil;

  typedef class NPair *(*PairCreator)(class CAC *);
  PairCreator *pairclass;
  char **pairnames;
  int *pairmasks;
  class NPair **neigh_pair;

  // internal methods
  // including creator methods for Nbin,Nstencil,Npair instances

  void init_styles();
  int init_pair();

  void morph_other();
//  void morph_skip();
//  void morph_granular();
//  void morph_halffull();
  void morph_copy();

  void print_pairwise_info();
  void requests_new2old();

  int choose_bin(class NeighRequest *);
  int choose_stencil(class NeighRequest *);
  int choose_pair(class NeighRequest *);

  template <typename T> static NBin *bin_creator(class CAC *);
  template <typename T> static NStencil *stencil_creator(class CAC *);
  template <typename T> static NPair *pair_creator(class CAC *);

  // dummy functions provided by NeighborKokkos, called in init()
  // otherwise NeighborKokkos would have to overwrite init()

  int copymode;
};

namespace NeighConst {

  // masks for bin methods

  static const int NB_STANDARD = 1<<0;
  static const int NB_ELEMONLY = 1<<1;
  
  // masks for stencil methods
  
  static const int NS_BIN     = 1<<0;
  static const int NS_MULTI   = 1<<1;
  static const int NS_HALF    = 1<<2;
  static const int NS_FULL    = 1<<3;
  static const int NS_2D      = 1<<4;
  static const int NS_3D      = 1<<5;
  static const int NS_NEWTON  = 1<<6;
  static const int NS_NEWTOFF = 1<<7;
  static const int NS_ORTHO   = 1<<8;
  static const int NS_TRI     = 1<<9;
  static const int NS_GHOST   = 1<<10;
  static const int NS_ELEMONLY   = 1<<11;

  // masks for neigh pair methods
  
  static const int NP_NSQ           = 1<<0;
  static const int NP_BIN           = 1<<1;
  static const int NP_MULTI         = 1<<2;

  // half or full list
  
  static const int NP_HALF          = 1<<3;
  static const int NP_FULL          = 1<<4;

  // orthogonal or triclinic box

  static const int NP_ORTHO         = 1<<5;
  static const int NP_TRI           = 1<<6;

  // contain neighbor list for atoms, elements, integration points, and/or interpolated atoms (or outers)
  
  static const int NP_ATOM          = 1<<7;
  static const int NP_ATOMONLY      = 1<<8;
  static const int NP_ELEM          = 1<<9;
  static const int NP_NODE          = 1<<10;
  static const int NP_INTG          = 1<<11;
  static const int NP_INTPL         = 1<<12;
  static const int NP_INTPL_OUTER   = 1<<13;

  // newton setting
  
  static const int NP_NEWTON        = 1<<14;
  static const int NP_NEWTOFF       = 1<<15;

  // neighborlist for ghost
  
  static const int NP_GHOST         = 1<<16;

  // threebody flag: neighborlist for interpolated atoms neighbor

  static const int NP_THREEBODY     = 1<<17;

  // a copy from another list
  
  static const int NP_COPY          = 1<<18;

  // perform neighbor sorting

  static const int NP_SORT          = 1<<19;

  // neighbor list only for specific group
  
  static const int NP_GROUP         = 1<<20;

  // neighbor list exclude inner interpolated atoms
  
  static const int NP_NOINNER       = 1<<21;
}

}

#endif

