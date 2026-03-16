#ifndef CAC_ELEMENT_H
#define CAC_ELEMENT_H

#include "pointers.h"
#include <map>
#include <string>


namespace CAC_NS {

class Element : protected Pointers {
 public: 
  char *element_style;
  class ElementVec *evec;

  // element counts
  
  bigint nelements;
  int nlocal,nghost;     // number of local elements and ghost elements
  int netypes;           // number of element types

  // element data & per-element arrays
  
  double **x;
  tagint *tag;
  int    *ctype;
  int    *etype;
  int    *mask;
  int     nsurface;
  int     nedge;
  //int    **elemlink;
  //int    linkcut;
  //int    linkstyle;     // 1 = link by neighbor
                        // 2 = link by tag id
  imageint *image;
  double *elem_size;    // size of each element
  double *initial_size; // initial size of each of element

  double **cell_size;   // size of elements along each local axes
  double ***slip_plane; // unit vectors of local axes in global coordinates

  double max_size_array[3],max_size;       // maximum element cubical size (x, y, or z) (for ghost elements)
  double max_diag_size;  // maximum element diagonal size (for neighboring)
  int nmax;              // maximum # of owned+ghost elements that the allocated array can store
  int tag_enable;        // 0/1 if element ID tags are defined
  int map_style;        // style of element map: 0 = none, 1 = array, 2 = hash
  int map_user;         // user selected style = same 0,1,2
  tagint map_tag_max;
  int *sametag;

  // node data & per-node arrays

  bigint nnodes; 
  double ***nodex, ***nodev, ***nodef;
  tagint **nodetag;       // use only for node connectivity in .dat file for tecplot visualization
  int **nodemask;
  int **n2ia;             // mapping from node index to interpolation point index in element
  int **n2i;              // mapping from node index to integration point index in element
  int npe;                // number of nodes per element dependent on the dimension of the simulations
  int apc;                // number of element per clusters;

  // interpolated atoms data & per-intpl arrays

  bigint nintpls;
  int max_nintpl;                 // store the max # of interpolation points in all types of element
  double ***shape_array;          // shape function inside each type of elements for all interpolated atoms
  int *nintpl;                    // number of interpolated atoms in an element for each type of element
  int ***surface_intpl;           // indices of interpolated atoms on the surfaces for each type of element
  int **nsurface_intpl;           // number of interpolated atoms on the surfaces for each type of element
  int max_surface_intpl;          // max number of interpolated atoms on a surface
  int ***edge_intpl;              // indices of interpolated atoms on the edges for each type of element
  int **nedge_intpl;              // number of interpolated atoms on the edges for each type of element
  int max_edge_intpl;             // max number of interpolated atoms on an edge

  bigint nmaxintpl;               // maximum # of owned+ghost interpolated atoms

  // integration points data & per-intg arrays

  double **weight;           // weight of each integration point
  double ***weighted_shape_array; // weighted shape array for force distribution from integration points to nodes.
  int **i2ia;                // mapping from integration point index to interpolated atom index in element
  int **ia2i;                // mapping from interpolated atom index to intergration point index in element 
                             // -1 if not an integration point
  int **i2n;                 // mapping from integration point index to node index in element
                             // -1 if not a node
  int max_nintg;             // store the max # of integration points in all types of element
  int *nintg;                // number of integration points in an element for each type of element
  int nmaxintg;              // maximum # of owned+ghost integration points
  int intg_input_style;      // 0/1/2 for not set, calculate, user-specify respectively
  double **wgauss,**xgauss;
  double node_weight;                // node weight

  // sub-element data & per-sub-element arrays

  int **natom_subelem;                 // # of interpolated atoms inside the sub element
  int esplit;                          // number of sub elements in each direction
  int nsubelem;                        // number of sub-elements in each element
  double **shape_array_center_subelem; // shape function of the center of the sub elements
  int ***ias2ia;                       // mapping from interpolated atom index in sub element to index in element

  // maximum data for each element

  int maxintpl;
  int maxintg;
  int maxsubelem;
  int maxsubintpl;
  double maxelemchange;
  
  // extra peratom info in restart file destined for fix & diag

  double **extra;


  // callback ptrs for atom arrays managed by fix classes
  // list of fixes to call during exchange, border, or restart

  int nextra_grow,nextra_restart,nextra_border;  // # of callbacks of each type
  int *extra_grow,*extra_restart,*extra_border;  // index of fix to callback to
  int nextra_grow_max,nextra_restart_max;        // size of callback lists
  int nextra_border_max;
  int nextra_store;

  // ElementVec factory types and map

  typedef ElementVec *(*ElementVecCreator)(CAC *);
  typedef std::map<std::string,ElementVecCreator> ElementVecCreatorMap;
  ElementVecCreatorMap *evec_map;
 
  // functions

  Element(class CAC *);
  ~Element();

  void add_callback(int);
  void delete_callback(const char *, int);
  void update_callback(int);

  void create_evec(const char *, int, char **);
  class ElementVec *new_evec(const char *);

  void init();
  void setup() {}
  void grow(int);
  void tag_check();
  void tag_extend();

  void data_elements(int, char *, tagint, int, int, int, double *);
  void data_nodes(int, char *, tagint, int, double *);
  void data_vels(int, char *, tagint);
  void data_fix_compute_variable(int, int);
  void data_pass_on_fix_compute_variable(int, int, int *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_border(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void modify_params(int, char **);
  void modify_elements(int, char **);
  void add_etype(int, char **);
  int element2atom(int); 
  void reset_node_tags(int);
  void check_node_coords(int);
  void count_intpl();

  bigint memory_usage();
  int memcheck(const char *);

  // map functions

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

 protected:

  // per-element arrays

  int *map_array;
  int map_maxarray;

  struct HashElem {     // hashed map
    tagint global;      // key to search on = global ID
    int local;          // value associated with key = local index
    int next;           // next entry in this bucket, -1 if last
  };    
  int map_nhash;
  int map_nused;
  int map_free;
  int map_nbucket;
  int *map_bucket;
  HashElem *map_hash;

  int max_same;
 
  int memlength;                  // allocated size of memstr
  char *memstr;                   // string of array names already counted

 private:
  template <typename T> static ElementVec *evec_creator(CAC *);

};
};

#endif
