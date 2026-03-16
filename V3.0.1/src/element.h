#ifndef CAC_ELEMENT_H
#define CAC_ELEMENT_H

#include "pointers.h"
#include <map>
#include <string>


namespace CAC_NS {

class Element : protected Pointers {
 public: 
  char *element_style;
  char **element_shape_list;  // list of available element shapes

  double ***mass_matrices;             // list of consistent mass matrices for each element shape
  double ***inversed_mass_matrices;    // list of inversed consistent mass matrices for each element shape
  int ***surface_node_list;   // list of nodes for each surface of element
  int **surface_num_node;     // number of node on each element surface
  int mass_style;             // 0 = lumped mass matrix
                              // 1 = consistent mass matrix
  int num_element_shapes;     // number of available element shapes
  int inner_rigid_flag;       // whether inner gaussian cells feel force from atoms
  class ElementVec *evec;

  // element shapes
  
  enum{QUADRILATERAL, TRIANGLE, HEXAHEDRON, TETRAHEDRON, OCTAHEDRON, PYRAMID, WEDGE};

  // mass style
  
  enum{LUMPED, CONSISTENT};

  // weight scale style

  enum{NOSCALE, UNIFORM, NONUNIFORM};

  // element counts
  
  bigint nelements;
  int nlocal, nghost;     // number of local elements and ghost elements
  int netypes;           // number of element types

  // element data & per-element arrays
  
  double **x;           // element center position
  tagint *tag;
  int *grain_tag;
  int **ctype;        // chemical type of elements (for pair force calculation)
  int *etype;        // type of elements (for shape and interpolation)
  int *mask;
  int *nsurface;     // number of surfaces for each element type
  int *nedge;        // number of edges for each element type
  int *element_shape_ids;     // ids of element shape for each element type
  char **element_shape_names;  // names of element shape for each element type

  int nodex_lamda_flag;
  int x_lamda_flag;

  imageint *image;
  double *subelem_size;    // size of each element diagonally
  double **initial_box_size; // initial size of each of element
                             // for triclinic box, box size is measured in lamda coords

  double **cell_size;      // size of elements along each local axes
  double ***surface_plane; // unit vectors of local axes in global coordinates
  double **element_box2xi;  // box to natural coordinate transformation matrix for each element
  double **element_box2xi_scale;

  // element bound box
  // orthogonal: in box coords
  // triclinic: in lamda coords
  
  double max_element_bound_box_size[3];  // maximum element bound box size locally (for ghost cutoff and element bins)
  double local_element_bound_box[6];     // local bound box for all my local elements (for ghost cutoff)
  double **element_bound_box;    // bounding box for each element
                                 // 0, 1, 2 = x y z lower bound
                                 // 3, 4, 5 = x y z upper bound
  double max_diag_size;  // maximum element diagonal size (for neighboring)
  int nmax;              // maximum # of owned+ghost elements that the allocated array can store
  int tag_enable;        // 0/1 if element ID tags are defined
  int map_style;         // style of element map: 0 = none, 1 = array, 2 = hash
  int map_user;          // user selected style = same 0, 1, 2
  tagint map_tag_max;
  int *sametag;

  // node data & per-node arrays

  bigint nnodes;    // # of nodes 
  bigint nncells;   // # of nodecells 

  double ****nodex, ****nodev, ****nodef, ****gaussf;    // [ielem][ibasis][inode][dim]
  //tagint **nodetag;       // use only for node connectivity in .dat file for tecplot visualization
                          // also use to store old atom IDs when using coarse_graining command
  int ***nodemask;
  int **n2u;             // mapping from node index to interpolation point index in element
  int **n2g;              // mapping from node index to gaussian cell index in element
  int maxnpe;            // maximum number of nodes for each element;
  int *npe;               // number of nodes per element for each element type
  int maxapc;            // maximum number of atoms per unit cell
  int *apc;               // number of atoms per unit cells for each element type
  double *nodal_weight;   // nodal weight for each element type

  // unit cells data & per-ucell arrays

  bigint nucells;   // # of unit cells
  bigint nvatoms;   // # of virtual atoms (including nodes)
  int max_nucell;                 // store the max # of interpolation points in all types of element
  double ***shape_array;          // shape function inside each type of elements for all virtual atoms
  int *nucell;                    // number of virtual atoms in an element for each type of element
  int ***surface_ucell;           // indices of virtual atoms on the surfaces for each type of element
  int **nsurface_ucell;           // number of virtual atoms on the surfaces for each type of element
  int max_outer_ucell;            // max number of virtual atoms on the element boundary
  int max_surface_ucell;          // max number of virtual atoms on a surface
  int max_surface;                // max number of surfaces on an element
  int ***edge_ucell;              // indices of virtual atoms on the edges for each type of element
  int **nedge_ucell;              // number of virtual atoms on the edges for each type of element
  int **u2g;                // mapping from virtual atom index to intergration point index in element 
  int max_edge_ucell;             // max number of virtual atoms on an edge
  int max_edge;                   // max number of edges on an element
  int **is_outer;                 // = 0 if interpolated point is inner
                                  // = index of outer otherwise


  // gaussian cells data & per-gcell arrays

  double **weight;           // weight of each gaussian cell
  double weight_scale[4];     // weight scale of each gaussian cell type (node, edge, surface, inner)
  int weight_scale_flag;    // if weight scale of each gaussian cell type (node, edge, surface, inner) is set
  double ***weighted_shape_array; // weighted shape array for force distribution from gaussian cells to nodes.
  int **g2u;                // mapping from gaussian cell index to virtual atom index in element
                             // -1 if not an gaussian cell
  int **g2n;                 // mapping from gaussian cell index to node index in element
                             // -1 if not a node
  int max_ngcell;             // store the max # of gaussian cells in all types of element
  int *ngcell;                // number of gaussian cells in an element for each type of element
  double **wgauss, **xgauss;

  // sub-element data & per-sub-element arrays

  int **nucell_subelem;                    // # of virtual atoms inside the sub element
  int *subsplit;                          // sub element split division for each element type
  double subelem_size_factor;             // scale up factor for subelem size to account for element large deformation error
  int *nsubelem;                          // number of sub-elements in each element type
  double ***shape_array_center_subelem;   // shape function of the center of the sub elements for each type
  int ***us2u;                          // mapping from virtual atom index in sub element to index in element

  // maximum data for each element

  int maxucell;
  int maxgcell;
  int maxsubelem;
  int maxsubucell;
  double maxelemchange;
  
  // extra peratom info in restart file destined for fix & diag

  double **extra;


  // callback ptrs for atom arrays managed by fix classes
  // list of fixes to call during exchange, border, or restart

  int nextra_grow, nextra_restart, nextra_border;  // # of callbacks of each type
  int *extra_grow, *extra_restart, *extra_border;  // index of fix to callback to
  int nextra_grow_max, nextra_restart_max;        // size of callback lists
  int nextra_border_max;
  int nextra_store;

  // ElementVec factory types and map

  typedef ElementVec *( * ElementVecCreator)(CAC *);
  typedef std::map<std::string, ElementVecCreator> ElementVecCreatorMap;
  ElementVecCreatorMap *evec_map;

  int node_set_3D[3][2][4];       // set of nodes for each face of element (+-x, +-y, +-z)   
  int node_set_2D[2][2][2];       // ditto for 2D

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
  void compress_etype();
  void update_subelem_size();
  void update_element_bound_box(int check = 0);
  void tag_extend();
  int find_element_shape_id(char *);
  void destroy_etype_arrays();

  void data_elements(int, char *, tagint, int, int, double *);
  void data_nodes(int, char *, tagint, int, int, double *);
  void data_nodes_reference(int, char *, tagint, int, int, double *, int);
  void data_vels(int, char *, tagint);

  void data_elements_V2(int, char *, tagint, int, int, int, double *);
  void data_nodes_V2(int, char *, tagint, int, double *);
  void data_nodes_reference_V2(int, char *, tagint, int, double *);
  void data_vels_V2(int, char *, tagint);

  void data_fix_compute_dump_variable(int, int);
  void data_pass_on_fix_compute_dump_variable(int, int, int *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_border(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void modify_params(int, char **);
  void modify_elements(int, char **);
  void add_etype(int, char **);
  int element2atom(int); 
  void set_node_connectivities(int, int **);
  void set_node_connectivities(int, int ***);
  void check_node_coords(int);
  bigint count_vatom();
  bigint count_ucell();
  int count_nodes(int);
  int count_polygons(int);
  int count_node_cells(int);
  void compute_nodef();
  int node_connectivity(int *, int, int *);
  void box2natural(double *, double *, int);


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

  // for atom_strain command

  double ****nodex_current;
  void data_nodes_strain(int, char *);

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
  void define_mass_matrices();
  void define_inversed_mass_matrices();
  void define_surface_node_list();

  // debug 
  
 public:
  int **debug_gcell_type; // gaussian cell type 
                         // 0 = inner
                         // 1 = surface
                         // 2 = edge
                         // 3 = node
                         
  int debug_mode;        // 1 = on, 0 = off
  int **debug_gcell_inactive;

  void debug_element_gcell_off(int, char **);
};
};

#endif
