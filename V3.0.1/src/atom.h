#ifndef CAC_ATOM_H
#define CAC_ATOM_H

#include "pointers.h"
#include <map>
#include <string>

namespace CAC_NS {

class Atom : protected Pointers {
 public:
  char *atom_style;
  class AtomVec *avec;

  //atom counts

  bigint natoms;                     //total # of atoms in systems, could be 0
  //natoms may not be current if atoms lost

  int nlocal, nghost;                 // # of owned and ghost atoms on this proc
  int nmax;                          // max # of owned+ghost in arrays on this proc
  int tag_enable;                    // 0/1 if atom ID tags are defined
  int ngrains;                       // # of crystal grains

  int ntypes;

  int firstgroup;               // store atoms in this group first, -1 if unset
  int nfirst;                   // # of atoms in first group on this proc
  char *firstgroupname;         // group-ID to store first, NULL if unset

  // per-atom arrays
  // customize by adding new array

  tagint *tag;
  int *type, *mask;
  imageint *image;
  double **x, **v, **f;
  int *grain_tag;               // grain id starts from 1, 2, ...
                                // grain id = 0 means atoms don't belong to any grains

  // custom arrays used by fix property/atom

  int **ivector;
  double **dvector;
  char **iname, **dname;
  int nivector, ndvector;

  // Peridynamics scale factor, used by dump cfg

  double pdscale;

  // extra peratom info in restart file destined for fix & diag

  double **extra;

  // per-type arrays

  double *mass;
  int *mass_setflag;

  // callback ptrs for atom arrays managed by fix classes
  // list of fixes to call during exchange, border, or restart

  int nextra_grow, nextra_restart, nextra_border;  // # of callbacks of each type
  int *extra_grow, *extra_restart, *extra_border;  // index of fix to callback to
  int nextra_grow_max, nextra_restart_max;        // size of callback lists
  int nextra_border_max;
  int nextra_store;

  int map_style;                  // style of atom map: 0=none, 1=array, 2=hash
  int map_user;                   // user selected style = same 0, 1, 2
  tagint map_tag_max;             // max atom ID that map() is setup for

  // spatial sorting of atoms

  int sortfreq;             // sort atoms every this many steps, 0 = off
  bigint nextsort;          // next timestep to sort on
  double userbinsize;       // requested sort bin size

  // indices of atoms with same ID

  int *sametag;              // sametag[I] = next atom with same ID, -1 if no more

  // AtomVec factory types and map

  typedef AtomVec *( * AtomVecCreator)(CAC *);
  typedef std::map<std::string, AtomVecCreator> AtomVecCreatorMap;
  AtomVecCreatorMap *avec_map;
 
  // functions

  Atom(class CAC *);
  ~Atom();

  void add_callback(int);
  void delete_callback(const char *, int);
  void update_callback(int);

  void create_avec(const char *, int, char **);
  virtual class AtomVec *new_avec(const char *);
  void init();
  void setup();

  void modify_params(int, char **);
  void tag_check();
  void tag_extend();
  tagint maxtag();
  
  virtual void allocate_type_arrays();

  void data_atoms(int, char *, tagint, int, int, double *);
  void data_atoms_reference(int, char *, tagint, int, int, double *, int);
  void data_vels(int, char *, tagint);

  void data_atoms_V2(int, char *, tagint, int, int, double *);
  void data_atoms_reference_V2(int, char *, tagint, int, int, double *);
  void data_vels_V2(int, char *, tagint);

  void data_fix_compute_dump_variable(int, int);
  void data_pass_on_fix_compute_dump_variable(int, int, int);
  void set_mass(const char *, int, const char *, int);
  void set_mass(const char *, int, int, double);
  void set_mass(const char *, int, int, char **);
  void set_mass(double *);
  void check_mass(const char *, int);
  virtual void sort();

  inline int * get_map_array() {return map_array;};
  inline int get_map_size() {return map_tag_max+1;};

  bigint memory_usage();
  int memcheck(const char *);

  // functions for global to local ID mapping
  // map lookup function inlined for efficiency
  // return -1 if no map defined

  inline int map(tagint global) {
    if (map_style == 1) return map_array[global];
    else if (map_style == 2) return map_find_hash(global);
    else return -1;
  };

  void map_init(int check = 1);
  void map_clear();
  void map_set();
  void map_one(tagint, int);
  int map_style_set();
  void map_delete();
  int map_find_hash(tagint);

  // for atomic_strain command

  int atom_strain_flag;
  double **x_current;
  void data_atoms_strain(int, char *);

 protected:

  // global to local ID mapping

  int *map_array;       // direct map via array that holds map_tag_max
  int map_maxarray;     // allocated size of map_array (1 larger than this)

  struct HashElem {     // hashed map
    tagint global;      // key to search on = global ID
    int local;          // value associated with key = local index
    int next;           // next entry in this bucket, -1 if last
  };
  int map_nhash;        // # of entries hash table can hold
  int map_nused;        // # of actual entries in hash table
  int map_free;         // ptr to 1st unused entry in hash table
  int map_nbucket;      // # of hash buckets
  int *map_bucket;      // ptr to 1st entry in each bucket
  HashElem *map_hash;   // hash table

  int max_same;         // allocated size of sametag

  // spatial sorting of atoms

  int nbins;                      // # of sorting bins
  int nbinx, nbiny, nbinz;          // bins in each dimension
  int maxbin;                     // max # of bins
  int maxnext;                    // max size of next, permute
  int *binhead;                   // 1st atom in each bin
  int *next;                      // next atom in bin
  int *permute;                   // permutation vector
  double bininvx, bininvy, bininvz; // inverse actual bin sizes
  double bboxlo[3], bboxhi[3];     // bounding box of my sub-domain

  int memlength;                  // allocated size of memstr
  char *memstr;                   // string of array names already counted

  void setup_sort_bins();

 private:
  template <typename T> static AtomVec *avec_creator(CAC *);
};

}

#endif
