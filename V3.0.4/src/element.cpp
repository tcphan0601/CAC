#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "element.h"
#include "style_element.h"
#include "element_vec.h"
#include "region.h"
#include "group.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "input.h"
#include "output.h"
#include "variable.h"
#include "compute.h"
#include "dump.h"
#include "memory.h"
#include "modify.h"
#include "fix.h"
#include "error.h"
#include "universe.h"
#include "math_const.h"

#include "math_extra.h"

using namespace CAC_NS;
using namespace MathConst;
using namespace MathExtra;

#define NSHAPES       7     // increase NSHAPES after adding new element shapes
#define DELTA         1
#define DELTA_MEMSTR  1024
#define EPSILON       1.0e-6
#define EXTRA         1000
#define BIG           1e30

Element::Element(CAC *cac) : Pointers(cac)
{
  inner_rigid_flag = 0;
  num_element_shapes = NSHAPES;
  memory->create(element_shape_list, num_element_shapes, 256, "element:element_shape_list");
  strcpy(element_shape_list[QUADRILATERAL], "Quad");
  strcpy(element_shape_list[TRIANGLE], "Tri");
  strcpy(element_shape_list[HEXAHEDRON], "Hex");
  strcpy(element_shape_list[TETRAHEDRON], "Tet");
  strcpy(element_shape_list[OCTAHEDRON], "Oct");
  strcpy(element_shape_list[PYRAMID], "Pyr");
  strcpy(element_shape_list[WEDGE], "Wedge");

  memory->create(mass_matrices, num_element_shapes, 8, 8, "element:mass_matrices");
  memory->create(inversed_mass_matrices, num_element_shapes, 8, 8, "element:inversed_mass_matrices");
  memory->create(surface_node_list, num_element_shapes, 6, 4, "element:surface_node_list");
  memory->create(surface_num_node, num_element_shapes, 6, "element:surface_num_node");
  memory->create(nsurface, num_element_shapes, "element:nsurface");
  memory->create(nedge, num_element_shapes, "element:nedge");

  nodal_force_style = CONSISTENT;
  define_mass_matrices();
  define_inversed_mass_matrices();
  define_surface_node_list();

  nsurface[QUADRILATERAL] = 1;
  nsurface[TRIANGLE] = 1;
  nsurface[HEXAHEDRON] = 6;
  nsurface[TETRAHEDRON] = 4;
  //nsurface[OCTAHEDRON] = 8;
  nsurface[PYRAMID] = 5;
  nsurface[WEDGE] = 5;
  nedge[QUADRILATERAL] = 4;
  nedge[TRIANGLE] = 3;
  nedge[HEXAHEDRON] = 12;
  nedge[TETRAHEDRON] = 6;
  //nedge[OCTAHEDRON] = 12;
  nedge[PYRAMID] = 8;
  nedge[WEDGE] = 9;

  // default parameter for each element
  // can be modified through element_modify command

  max_outer_ucell = 2400;
  max_surface_ucell = 400;
  max_surface = 8;
  max_edge_ucell = 20;
  max_edge = 12;
  maxucell = 8000;  
  maxgcell = 64;
  maxsubucell = 500;
  maxsubelem = 125;
  maxelemchange = 2.0;
  subelem_size_factor = 1.1;

  nodex_lamda_flag = x_lamda_flag = 0;

  nmax = 0;
  max_ngcell = max_nucell = 0;
  nelements = nnodes = nucells = nvatoms = 0;
  netypes = 0; 
  nlocal = 0;
  nghost = 0;

  // per-element arrays

  nodex = NULL;
  nodev = NULL;
  nodef = NULL;
  gaussf = NULL;
  x = NULL;
  image = NULL;
  surface_plane = NULL;
  element_box2xi = NULL;
  element_box2xi_scale = NULL;
  cell_size = NULL;
  tag = NULL;
  grain_tag = NULL;
  ctype = NULL;
  etype = NULL;
  mask = NULL;
  nodemask = NULL;
  subelem_size = NULL;
  initial_box_size = NULL;
  element_bound_box = NULL;

  // netype-length arrays

  apc = NULL;
  npe = NULL;
  element_shape_ids = NULL;
  element_shape_names = NULL;
  nucell = NULL;
  surface_ucell = NULL;
  nsurface_ucell = NULL;
  edge_ucell = NULL;
  nedge_ucell = NULL;
  ngcell = NULL;
  us2u = NULL;
  g2u = NULL; 
  u2g = NULL;
  g2n = NULL;
  is_outer = NULL;
  n2g = NULL;
  n2u = NULL;
  nucell_subelem = NULL;
  weight = NULL;
  weight_scale_flag = NOSCALE;
  weight_scale[0] = weight_scale[1] = -1;
  weight_scale[2] = weight_scale[3] = -1;
  shape_array = NULL; 
  weighted_shape_array = NULL;
  nodal_weight = NULL;

  nsubelem = NULL;
  shape_array_center_subelem = NULL;
  subsplit = NULL;

  max_element_bound_box_size[0] = 0.0;
  max_element_bound_box_size[1] = 0.0;
  max_element_bound_box_size[2] = 0.0;
  local_element_bound_box[0] = BIG;
  local_element_bound_box[1] = BIG;
  local_element_bound_box[2] = BIG;
  local_element_bound_box[3] = -BIG;
  local_element_bound_box[4] = -BIG;
  local_element_bound_box[5] = -BIG;

  max_same = 0;

  // element map

  map_array = NULL;
  map_bucket = NULL;
  map_hash = NULL;
  sametag = NULL;
  map_style = map_user = 0;
  map_tag_max = -1;
  map_maxarray = map_nhash = -1; 

  tag_enable = 1;

  // gauss gaussian cells and weights for each element shape
  // modify for more gaussian cells 

  memory->create(xgauss, 6, 6, "element:xgauss");
  memory->create(wgauss, 6, 6, "element:wgauss");

  xgauss[1][1] = 0.0;
  wgauss[1][1] = 2.0;

  xgauss[2][1] = -sqrt(1.0/3.0);
  xgauss[2][2] =  sqrt(1.0/3.0);
  wgauss[2][1] = 1.0;
  wgauss[2][2] = 1.0;

  xgauss[3][1] = -sqrt(3.0/5.0);
  xgauss[3][2] = 0.0;
  xgauss[3][3] =  sqrt(3.0/5.0);
  wgauss[3][1] = 5.0/9.0;
  wgauss[3][2] = 8.0/9.0;
  wgauss[3][3] = 5.0/9.0;

  xgauss[4][1] = -sqrt(3.0/7.0+2.0/7.0 * sqrt(1.2));
  xgauss[4][1] = -sqrt(3.0/7.0-2.0/7.0 * sqrt(1.2));
  xgauss[4][1] =  sqrt(3.0/7.0-2.0/7.0 * sqrt(1.2));
  xgauss[4][1] =  sqrt(3.0/7.0+2.0/7.0 * sqrt(1.2));
  wgauss[4][1] = (18.0-sqrt(30.0))/36.0;
  wgauss[4][2] = (18.0+sqrt(30.0))/36.0;
  wgauss[4][3] = (18.0+sqrt(30.0))/36.0;
  wgauss[4][4] = (18.0-sqrt(30.0))/36.0;

  xgauss[5][1] = -sqrt(5.0+2.0 * sqrt(10.0/7.0))/3.0;
  xgauss[5][2] = -sqrt(5.0-2.0 * sqrt(10.0/7.0))/3.0;
  xgauss[5][3] = 0;
  xgauss[5][4] =  sqrt(5.0-2.0 * sqrt(10.0/7.0))/3.0;
  xgauss[5][5] =  sqrt(5.0+2.0 * sqrt(10.0/7.0))/3.0;
  wgauss[5][1] = (322.0-13.0 * sqrt(70.0))/900.0;
  wgauss[5][2] = (322.0+13.0 * sqrt(70.0))/900.0;
  wgauss[5][3] = 128.0/225.0;
  wgauss[5][4] = (322.0+13.0 * sqrt(70.0))/900.0;
  wgauss[5][5] = (322.0-13.0 * sqrt(70.0))/900.0;

  // callback lists & extra restart info

  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = NULL;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = NULL;

  element_style = NULL;
  evec = NULL;

  evec_map = new ElementVecCreatorMap();

  nodex_current = NULL;

  node_set_3D[0][0][0] = 0;
  node_set_3D[0][0][1] = 3;
  node_set_3D[0][0][2] = 7;
  node_set_3D[0][0][3] = 4;
  node_set_3D[0][1][0] = 1;
  node_set_3D[0][1][1] = 2;
  node_set_3D[0][1][2] = 6;
  node_set_3D[0][1][3] = 5;

  node_set_3D[1][0][0] = 0;
  node_set_3D[1][0][1] = 1;
  node_set_3D[1][0][2] = 5;
  node_set_3D[1][0][3] = 4;
  node_set_3D[1][1][0] = 3;
  node_set_3D[1][1][1] = 2;
  node_set_3D[1][1][2] = 6;
  node_set_3D[1][1][3] = 7;

  node_set_3D[2][0][0] = 0;
  node_set_3D[2][0][1] = 1;
  node_set_3D[2][0][2] = 2;
  node_set_3D[2][0][3] = 3;
  node_set_3D[2][1][0] = 4;
  node_set_3D[2][1][1] = 5;
  node_set_3D[2][1][2] = 6;
  node_set_3D[2][1][3] = 7;

  node_set_2D[0][0][0] = 0;
  node_set_2D[0][0][1] = 3;
  node_set_2D[0][1][0] = 1;
  node_set_2D[0][1][1] = 2;

  node_set_2D[1][0][0] = 0;
  node_set_2D[1][0][1] = 1;
  node_set_2D[1][1][0] = 3;
  node_set_2D[1][1][1] = 2;


#define ELEMENT_CLASS
#define ElementStyle(key, Class) \
  (*evec_map)[#key] = &evec_creator<Class>;
#include "style_element.h"
#undef ElementStyle
#undef ELEMENT_CLASS

  // debug
  
  debug_mode = 0;
  debug_gcell_type = NULL;
  debug_gcell_inactive = NULL;
}

/*  ----------------------------------------------------------------------  */

Element::~Element()
{
  memory->destroy(xgauss);
  memory->destroy(wgauss);
  memory->destroy(element_shape_list);
  memory->destroy(mass_matrices);
  memory->destroy(inversed_mass_matrices);
  memory->destroy(surface_node_list);
  memory->destroy(surface_num_node);
  memory->destroy(nsurface);
  memory->destroy(nedge);

  delete [] element_style;
  delete evec;

  // delete per-element arrays

  memory->destroy(nodex);
  memory->destroy(nodex_current);
  memory->destroy(nodev);
  memory->destroy(nodef);
  memory->destroy(gaussf);
  memory->destroy(x);
  memory->destroy(image);
  memory->destroy(surface_plane);
  memory->destroy(element_box2xi);
  memory->destroy(element_box2xi_scale);
  memory->destroy(cell_size);
  memory->destroy(tag);
  memory->destroy(mask);
  memory->destroy(nodemask);
  memory->destroy(ctype);
  memory->destroy(etype);
  memory->destroy(subelem_size);
  memory->destroy(initial_box_size);
  memory->destroy(element_bound_box);

  // delete per-etype arrays

  destroy_etype_arrays();
  
  // delete mapping data structures

  map_delete();

}

/*  ----------------------------------------------------------------------
   create an ElementVec style
   called from cac.cpp, input script
   -------------------------------------------------------------------------  */

void Element::create_evec(const char *style, int narg, char **arg)
{

  delete [] element_style;
  if (evec) delete evec;
  element_style = NULL;
  evec = NULL;

  // create instance of ElementVec
  // use grow() to initialize element-based arrays to length 1
  // so that x[0][0] and nodex[0][0][0] can always be referenced even if proc has no elements

  evec = new_evec(style);
  evec->store_args(narg, arg);

  evec->process_args(narg, arg); 
  evec->grow(1); 

  int n = strlen(style) + 1;
  element_style = new char[n];
  strcpy(element_style, style);

}

/*  ----------------------------------------------------------------------
   generate an ElementVec class
   -------------------------------------------------------------------------  */

ElementVec *Element::new_evec(const char *style)
{
  if (evec_map->find(style) != evec_map->end()) {
    ElementVecCreator evec_creator = (*evec_map)[style];
    return evec_creator(cac);
  } else error->all(FLERR, "Unknown element style");
  return NULL;
}

/*  ----------------------------------------------------------------------
   one instance per ElementVec style in style_element.h
   -------------------------------------------------------------------------  */

  template <typename T>
ElementVec *Element::evec_creator(CAC *cac)
{
  return new T(cac);
}

/*  ----------------------------------------------------------------------
   register a callback to a fix so it can manage elem-based arrays
   happens when fix is created
   flag = 0 for grow, 1 for restart, 2 for border comm
   -------------------------------------------------------------------------  */

void Element::add_callback(int flag)
{
  int ifix;

  // find the fix
  // if find NULL ptr:
  //   it's this one, since it is being replaced and has just been deleted
  //   at this point in re-creation
  // if don't find NULL ptr:
  //   i is set to nfix = new one currently being added at end of list

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix] == NULL) break;

  // add callback to lists, reallocating if necessary

  if (flag == 0) {
    if (nextra_grow == nextra_grow_max) {
      nextra_grow_max += DELTA;
      memory->grow(extra_grow, nextra_grow_max, "element:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
  } else if (flag == 1) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
      memory->grow(extra_restart, nextra_restart_max, "element:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
  } else if (flag == 2) {
    if (nextra_border == nextra_border_max) {
      nextra_border_max += DELTA;
      memory->grow(extra_border, nextra_border_max, "element:extra_border");
    }
    extra_border[nextra_border] = ifix;
    nextra_border++;
  }
}

/*  ----------------------------------------------------------------------
   unregister a callback to a fix
   happens when fix is deleted, called by its destructor
   flag = 0 for grow, 1 for restart, 2 for border
   -------------------------------------------------------------------------  */

void Element::delete_callback(const char *id, int flag)
{
  if (id == NULL) return;

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(id, modify->fix[ifix]->id) == 0) break;

  // compact the list of callbacks

  if (flag == 0) {
    int match;
    for (match = 0; match < nextra_grow; match++)
      if (extra_grow[match] == ifix) break;
    for (int i = match; i < nextra_grow-1; i++)
      extra_grow[i] = extra_grow[i+1];
    nextra_grow--;

  } else if (flag == 1) {
    int match;
    for (match = 0; match < nextra_restart; match++)
      if (extra_restart[match] == ifix) break;
    for (int i = match; i < nextra_restart-1; i++)
      extra_restart[i] = extra_restart[i+1];
    nextra_restart--;

  } else if (flag == 2) {
    int match;
    for (match = 0; match < nextra_border; match++)
      if (extra_border[match] == ifix) break;
    for (int i = match; i < nextra_border-1; i++)
      extra_border[i] = extra_border[i+1];
    nextra_border--;
  }
}

/*  ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
   -------------------------------------------------------------------------  */

void Element::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix) extra_border[i]--;
}

/*  ----------------------------------------------------------------------  */

void Element::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy(extra);
    extra = NULL;
    nextra_store = 0;
  }

  // init ElemVec

  evec->init();

  // call update_element_bound_box to initialize initial_size
  // if triclinic convert node coord to lamda cord 
  // (no need to convert center x coord to lamda)

  if (domain->triclinic)
    domain->nodex2lamda(element->nlocal, element->nodex);

  update_element_bound_box();

  if (domain->triclinic) 
    domain->lamda2nodex(element->nlocal, element->nodex);


}

/*  ----------------------------------------------------------------------
   unpack n lines from Element section of data file
   call style-specific routine to parse line
   -------------------------------------------------------------------------  */

void Element::data_elements(int n, char *buf, tagint id_offset, 
    int etype_offset, int shiftflag, double *shift)
{
  int m, xptr, iptr;
  imageint imagedata;
  double xdata[3], lamda[3];
  double *coord;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_element && nwords != evec->size_data_element + 3)
    error->all(FLERR, "Incorrect element format in data file");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  char **values = new char*[nwords];

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3], subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != Comm::LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }
  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // xptr = which word in line starts xyz coords of element
  // iptr = which word in line starts ix, iy, iz image flags

  xptr = evec->xcol_data - 1;
  int imageflag = 0;
  if (nwords > evec->size_data_element) imageflag = 1;
  if (imageflag) iptr = nwords - 3; 

  // loop over lines of element data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap element into simulation box

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect element format in data file");
    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL, " \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR, "Incorrect element format in data file");
    }

    if (imageflag)
      imagedata = ((imageint) (atoi(values[iptr]) + IMGMAX) & IMGMASK) |
        (((imageint) (atoi(values[iptr+1]) + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (atoi(values[iptr+2]) + IMGMAX) & IMGMASK) << IMG2BITS);
    else imagedata = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    if (shiftflag) {
      xdata[0] += shift[0];
      xdata[1] += shift[1];
      xdata[2] += shift[2];
    }

    domain->remap(xdata, imagedata);
    if (triclinic) {
      domain->x2lamda(xdata, lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      evec->data_element(xdata, imagedata, values);

      if (id_offset) tag[nlocal-1] += id_offset;
      if (etype_offset) {
        etype[nlocal-1] += etype_offset;
        if (etype[nlocal-1] > netypes)
          error->one(FLERR, "Invalid element type in Elements section of data file");
      }
    }
    buf = next+1; 
  }
  delete [] values;
}

/*  ----------------------------------------------------------------------
   unpack n lines from Node Velocity section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   -------------------------------------------------------------------------  */

void Element::data_vels(int n, char *buf, tagint id_offset)
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_vel)
    error->all(FLERR, "Incorrect node velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of node velocities
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL, " \t\n\r\f");

    tagdata = ATOTAGINT(values[0]) + id_offset;
    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Node Velocities section of data file");
    if ((m = map(tagdata)) >= 0) evec->data_vel(m, &values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/*  ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   -------------------------------------------------------------------------  */

void Element::data_nodes(int n, char *buf, tagint id_offset, int ctype_offset, int shiftflag, double *shift) 
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node)
    error->all(FLERR, "Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect node format in data file");
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL, " \t\n\r\f");
      if (values[j] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
    }

    tagdata = ATOTAGINT(values[0]) + id_offset;

    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Nodes section of data file");
    if ((m = map(tagdata)) >= 0)
      evec->data_node(m, &values[1], ctype_offset, shiftflag, shift);

    buf = next + 1; 
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_nodes_reference(int n, char *buf, tagint id_offset, 
    int ctype_offset, int shiftflag, double *shift, int flag) 
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node + 3)
    error->all(FLERR, "Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect node format in data file");
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL, " \t\n\r\f");
      if (values[j] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
    }

    tagdata = ATOTAGINT(values[0]) + id_offset;

    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Nodes section of data file");
    if ((m = map(tagdata)) >= 0)
      evec->data_node_reference(m, &values[1], ctype_offset, shiftflag, shift, flag);
    buf = next + 1; 
  }

  delete [] values;
}


/*  ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   -------------------------------------------------------------------------  */

void Element::data_nodes_strain(int n, char *buf) 
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node + 3)
    error->all(FLERR, "Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect node format in data file");
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL, " \t\n\r\f");
      if (values[j] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
    }

    tagdata = ATOTAGINT(values[0]);

    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Nodes section of data file");
    if ((m = map(tagdata)) >= 0) 
      evec->data_node_strain(m, &values[1]);

    buf = next + 1; 
  }

  delete [] values;
}


/* --------------------------------------------------------------------------------------------------------------
  return # of bytes of allocated memory
  ---------------------------------------------------------------------------------------------------------------- */

bigint Element::memory_usage()
{
  // per-element arrays

  memlength = DELTA_MEMSTR;
  memory->create(memstr, memlength, "element:memstr");
  memstr[0] = '\0';
  bigint bytes = evec->memory_usage();

  memory->destroy(memstr);

  // element map

  bytes += max_same * sizeof(int);
  if (map_style == 1)
    bytes += memory->usage(map_array, map_maxarray);
  else if (map_style == 2) {
    bytes += map_nbucket * sizeof(int);
    bytes += map_nhash * sizeof(HashElem);
  }

  // per-type arrays

  int n = netypes + 1;
  bytes += memory->usage(npe, n);
  bytes += memory->usage(apc, n);
  bytes += memory->usage(element_shape_ids, n);
  bytes += memory->usage(element_shape_names, n, 256);
  bytes += memory->usage(weight, n, maxgcell);
  bytes += memory->usage(g2u, n, maxgcell);
  bytes += memory->usage(u2g, n, maxucell);
  bytes += memory->usage(g2n, n, maxgcell);
  bytes += memory->usage(is_outer, n, maxucell);
  bytes += memory->usage(n2g, n, maxnpe);
  bytes += memory->usage(n2u, n, maxnpe);
  bytes += memory->usage(nucell, n);
  bytes += memory->usage(surface_ucell, n, max_surface, max_surface_ucell);
  bytes += memory->usage(nsurface_ucell, n, max_surface);
  bytes += memory->usage(edge_ucell, n, max_edge, max_edge_ucell);
  bytes += memory->usage(nedge_ucell, n, max_edge);
  bytes += memory->usage(shape_array, n, maxucell, maxnpe);
  bytes += memory->usage(weighted_shape_array, n, maxgcell, maxnpe);
  bytes += memory->usage(shape_array_center_subelem, n, maxsubelem, maxnpe);
  bytes += memory->usage(nsubelem, n);
  bytes += memory->usage(subsplit, n);
  bytes += memory->usage(nucell_subelem, n, maxsubelem);
  bytes += memory->usage(us2u, n, maxsubelem, maxsubucell);

  return bytes;
}

/*  ----------------------------------------------------------------------
   accumulate per-element vec names in memstr, padded by spaces
   return 1 if padded str is not already in memlist, else 0
   -------------------------------------------------------------------------  */

int Element::memcheck(const char *str)
{
  int n = strlen(str) + 3;
  char *padded = new char[n];
  strcpy(padded, " ");
  strcat(padded, str);
  strcat(padded, " ");

  if (strstr(memstr, padded)) {
    delete [] padded;
    return 0;
  }

  if (strlen(memstr) + n >= memlength) {
    memlength += DELTA_MEMSTR;
    memory->grow(memstr, memlength, "element:memstr");
  }

  strcat(memstr, padded);
  delete [] padded;
  return 1;
}


/*  -----------------------------------
   Modify element parameters
   Called from element_modify command
   ------------------------------------ */

void Element::modify_params(int narg, char **arg)
{
  if (domain->box_exist) error->all(FLERR, "Element_modify command after simulation box is defined");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "id") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      if (domain->box_exist)
        error->all(FLERR, 
            "Element_modify id command after simulation box is defined");
      if (strcmp(arg[iarg+1], "yes") == 0) tag_enable = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) tag_enable = 0;
      else error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "map") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      if (domain->box_exist)
        error->all(FLERR, 
            "Element_modify map command after simulation box is defined");
      if (strcmp(arg[iarg+1], "array") == 0) map_user = 1;
      else if (strcmp(arg[iarg+1], "hash") == 0) map_user = 2;
      else if (strcmp(arg[iarg+1], "yes") == 0) map_user = 3;
      else error->all(FLERR, "Illegal element_modify command");
      map_style = map_user;
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/ucell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      maxucell = universe->inumeric(FLERR, arg[iarg+1]);
      if (maxucell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/gcell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      maxgcell = universe->inumeric(FLERR, arg[iarg+1]);
      if (maxgcell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/sub/ucell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      maxsubucell = universe->inumeric(FLERR, arg[iarg+1]);
      if (maxsubucell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/sub/elem") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      maxsubelem = universe->inumeric(FLERR, arg[iarg+1]);
      if (maxsubelem <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/elem/change") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      maxelemchange = universe->inumeric(FLERR, arg[iarg+1]);
      if (maxelemchange <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/surface/ucell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      max_surface_ucell = universe->inumeric(FLERR, arg[iarg+1]);
      if (max_surface_ucell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/outer/ucell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      max_outer_ucell = universe->inumeric(FLERR, arg[iarg+1]);
      if (max_outer_ucell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/edge/ucell") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      max_edge_ucell = universe->inumeric(FLERR, arg[iarg+1]);
      if (max_edge_ucell <= 0) error->all(FLERR, "Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "nodal/force/style") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      if (strcmp(arg[iarg+1], "lumped") == 0) nodal_force_style = LUMPED;
      else if (strcmp(arg[iarg+1], "consistent") == 0) nodal_force_style = CONSISTENT;
      else if (strcmp(arg[iarg+1], "node/only") == 0) nodal_force_style = NODEONLY;
      else error->all(FLERR, "Illegal element_modify command, invalid mass matrix type");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rigid") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal element_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) inner_rigid_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) inner_rigid_flag = 0;
      else error->all(FLERR, "Illegal element_modify command, invalid rigid option");
      iarg += 2;
    } else error->all(FLERR, "Illegal element_modify command");
  }
}

/*  ----------------------------------------------------------------------
   add element type, call from input script as add_etype command 
   -------------------------------------------------------------------------  */

void Element::add_etype(int narg, char **arg)
{
  if (narg != 9 && narg != 8) error->all(FLERR, "Illegal add_etype command");
  int itype;

  char *new_type_string;
  new_type_string = new char[256];
  int offset = 0;
  if (narg == 8) {
    itype = netypes+1;
    offset += sprintf(&new_type_string[offset], "%d ", itype);
  } else {
    itype = universe->inumeric(FLERR, arg[0]);
    if (itype <= netypes) 
      error->all(FLERR, "Illegal element type in add_etype command: etype already been set");
  }

  for (int i = 0; i < narg; i++)
    offset += sprintf(&new_type_string[offset], "%s ", arg[i]);

  // init to make sure subelements are setup

  init();  
  evec->grow_etype_arrays(itype);
  evec->set_element_types(new_type_string, 0);
  delete [] new_type_string;
}

/*  ----------------------------------------------------------------------
   assign node ids for node connectivity in dump tecplot command
   assign 1 id for each node cell
   called during dump tecplot, write_tecplot
   -------------------------------------------------------------------------  */

void Element::set_node_connectivities(int igroup, int **nodeid)
{

  // maxtag = # of nodes I own in group
  // maxtag_sum = # of total nodes on procs <= me in group

  int groupbit = group->bitmask[igroup];
  int maxtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) maxtag += npe[etype[i]]; 

  int maxtag_sum;

  MPI_Scan(&maxtag, &maxtag_sum, 1, MPI_CAC_TAGINT, MPI_SUM, world);

  // itag = 1st tag that my nodes in group should use

  tagint itag = maxtag_sum - maxtag + 1;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int j = 0; j < npe[etype[i]]; j++)
        nodeid[i][j] = itag++;
}

/*  ----------------------------------------------------------------------
   assign node ids for node connectivity in dump tecplot command
   assign unique id for each atoms in the node cell
   called during dump tecplot, write_tecplot
   -------------------------------------------------------------------------  */

void Element::set_node_connectivities(int igroup, int ***nodeid)
{

  // maxtag = # of nodes I own in group
  // maxtag_sum = # of total nodes on procs <= me in group

  int groupbit = group->bitmask[igroup];
  int maxtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) maxtag += npe[etype[i]] * apc[etype[i]]; 

  tagint maxtag_sum;

  MPI_Scan(&maxtag, &maxtag_sum, 1, MPI_CAC_TAGINT, MPI_SUM, world);

  // itag = 1st tag that my nodes in group should use

  int itag = maxtag_sum - maxtag + 1;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int j = 0; j < apc[etype[i]]; j++)
        for (int k = 0; k < npe[etype[i]]; k++)
          nodeid[i][j][k] = itag++;
}


/*  ----------------------------------------------------------------------
   add unique tags to any elements with tag = 0
   new tags are grouped by proc and start after max current tag
   called after creating new elements
   error if new tags will exceed MAXTAGINT
   -------------------------------------------------------------------------  */

void Element::tag_extend()
{
  // maxtag_all = max tag for all elements

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag, tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_CAC_TAGINT, MPI_MAX, world);

  // notag = # of elements I own with no tag (tag = 0)
  // notag_sum = # of total elements on procs <= me with no tag

  bigint notag = 0;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) notag++;

  bigint notag_total;
  MPI_Allreduce(&notag, &notag_total, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (notag_total >= MAXTAGINT)
    error->all(FLERR, "New element IDs exceed maximum allowed ID");

  bigint notag_sum;
  MPI_Scan(&notag, &notag_sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);

  // itag = 1st new tag that my untagged elements should use

  tagint itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) tag[i] = itag++;

}

/*  ----------------------------------------------------------------------
   check if node coords match up with center coords
   update_flag = 1: update node coords if not match 
   = 0: throw error if not match
   -------------------------------------------------------------------------  */

void Element::check_node_coords(int update_flag) 
{
  if (update_flag)
    evec->update_node_coord();
  else {
    double coord[3];
    int inpe;
    for (int i = 0; i < nlocal; i++) {
      inpe = npe[etype[i]];
      coord[0] = coord[1] = coord[2] = 0.0;
      for (int j = 0; j < apc[etype[i]]; j++) 
        for (int k = 0; k < npe[etype[i]]; k++) {
          coord[0] += nodex[i][j][k][0];
          coord[1] += nodex[i][j][k][1];
          coord[2] += nodex[i][j][k][2];
        }
      coord[0] /= npe[etype[i]] * apc[etype[i]];
      coord[1] /= npe[etype[i]] * apc[etype[i]];
      coord[2] /= npe[etype[i]] * apc[etype[i]];
      domain->remap(coord);
      if (fabs(coord[0]-x[i][0]) > EPSILON ||
          fabs(coord[1]-x[i][1]) > EPSILON ||
          fabs(coord[2]-x[i][2]) > EPSILON) {
        error->one(FLERR, "Node coords does not match with element center coords");
      }
    }
  }
}

/*  ----------------------------------------------------------------------
    check that element IDs are valid
    error if any element ID < 0 or element ID = MAXTAGINT
    if any element ID > 0, error if any element ID == 0
    if any element ID > 0, error if tag_enable = 0
    if all element IDs = 0, tag_enable must be 0
    if max element IDs < nelements, must be duplicates
    OK if max element IDs > nelements
NOTE: not fully checking that element IDs are unique
-------------------------------------------------------------------------  */

void Element::tag_check()
{
  tagint min = MAXTAGINT;
  tagint max = 0;

  for (int i = 0; i < nlocal; i++) {
    min = MIN(min, tag[i]);
    max = MAX(max, tag[i]);
  }

  tagint minall, maxall;
  MPI_Allreduce(&min, &minall, 1, MPI_CAC_TAGINT, MPI_MIN, world);
  MPI_Allreduce(&max, &maxall, 1, MPI_CAC_TAGINT, MPI_MAX, world);

  if (minall < 0) error->all(FLERR, "One or more element IDs is negative");
  if (maxall >= MAXTAGINT) error->all(FLERR, "One or more element IDs is too big");
  if (maxall > 0 && minall == 0)
    error->all(FLERR, "One or more element IDs is zero");
  if (maxall > 0 && tag_enable == 0)
    error->all(FLERR, "Non-zero element IDs with element_modify id = no");
  if (maxall == 0 && nelements && tag_enable)
    error->all(FLERR, "All element IDs = 0 but element_modify id = yes");
  if (tag_enable && maxall < nelements)
    error->all(FLERR, "Duplicate element IDs exist");
}

/*  ----------------------------------------------------------------------
    count number of virtual atoms
    -------------------------------------------------------------------------  */

bigint Element::count_vatom() 
{
  bigint n = 0;
  for (int i = 0; i < nlocal; i++)
    n += nucell[etype[i]] * apc[etype[i]];
  MPI_Allreduce(&n, &nvatoms, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  return n;
}


/*  ----------------------------------------------------------------------
    count number of unit cells
    -------------------------------------------------------------------------  */

bigint Element::count_ucell() 
{
  bigint n = 0;
  for (int i = 0; i < nlocal; i++)
    n += nucell[etype[i]];
  MPI_Allreduce(&n, &nucells, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  return n;
}

/*  ----------------------------------------------------------------------
    count number of nodes
    if flag = 1, sum up total number of nodes across procs
    -------------------------------------------------------------------------  */

int Element::count_node_cells(int flag) 
{
  int n = 0;
  int ntotal;
  for (int i = 0; i < nlocal; i++)
    n += npe[etype[i]];
  if (flag) {
    bigint bn = n;
    MPI_Allreduce(&bn, &nncells, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  }
  return n;
}

/*  ----------------------------------------------------------------------
    count number of FE polygons
    if flag = 1, sum up total number of polygons across procs
    -------------------------------------------------------------------------  */

int Element::count_polygons(int flag) 
{
  int n = 0;
  int ntotal;
  for (int i = 0; i < nlocal; i++)
    n += apc[etype[i]];
  if (flag) 
    MPI_Allreduce(&n, &ntotal, 1, MPI_INT, MPI_SUM, world);
  else ntotal = n;
  
  return ntotal;
}

/*  ----------------------------------------------------------------------
    count number of nodes
    if flag = 1, sum up total number of nodes across procs
    -------------------------------------------------------------------------  */

int Element::count_nodes(int flag) 
{
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    n += npe[etype[i]] * apc[etype[i]];
  if (flag) {
    bigint bn = n;
    MPI_Allreduce(&bn, &nnodes, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  }
  return n;
}


/*  ----------------------------------------------------------------------
    init per-element fix/compute/dump/variable values for newly created elements
    called from create_elements, read_data, read_dump, 
    fixes, computes, dumps, variables may or may not exist when called
    -------------------------------------------------------------------------  */

void Element::data_fix_compute_dump_variable(int nprev, int nnew)
{
  for (int m = 0; m < modify->nfix; m++) {
    Fix *fix = modify->fix[m];
    if (fix->create_attribute)
      for (int i = nprev; i < nnew; i++)
        fix->set_elem_arrays(i);
  }

  for (int m = 0; m < modify->ncompute; m++) {
    Compute *compute = modify->compute[m];
    if (compute->create_attribute)
      for (int i = nprev; i < nnew; i++)
        compute->set_elem_arrays(i);
  }

  for (int m = 0; m < output->ndump; m++) {
    Dump *dump = output->dump[m];
    if (dump->create_attribute)
      for (int i = nprev; i < nnew; i++)
        dump->set_elem_arrays(i);
  }


  for (int i = nprev; i < nnew; i++)
    input->variable->set_elem_arrays(i);
}

/*  ----------------------------------------------------------------------
    init per-element fix/compute/dump/variable values for a newly created element that is passed on from an old element
    called from fix_adaptive
    fixes, computes, dumps, variables may or may not exist when called
    -------------------------------------------------------------------------  */

void Element::data_pass_on_fix_compute_dump_variable(int i, int i_old, int *indexlist)
{
  for (int m = 0; m < modify->nfix; m++) {
    Fix *fix = modify->fix[m];
    if (fix->create_attribute)
      fix->pass_on_elem_arrays(i, i_old, indexlist);
  }

  for (int m = 0; m < modify->ncompute; m++) {
    Compute *compute = modify->compute[m];
    if (compute->create_attribute)
      compute->pass_on_elem_arrays(i, i_old, indexlist);
  }

  for (int m = 0; m < output->ndump; m++) {
    Dump *dump = output->dump[m];
    if (dump->create_attribute)
      dump->pass_on_elem_arrays(i, i_old, indexlist);
  }

  input->variable->pass_on_elem_arrays(i, i_old, indexlist);
}

/* ----------------------------------------
   update element bound box for acquiring ghost atoms/elements
   must be called in comm->setup()
   if check = 1, check for element irregular deformation
    caller must be careful in triclinic case where coords can be in box or lamda coord

------------------------------------------ */

void Element::update_element_bound_box(int check)
{
  if (nelements == 0) return;
  int *periodicity = domain->periodicity;

  double element_box_size;
  int i, j, k, dim;
  double maxcoord, mincoord;
  double ***inodex;
  int dimension = domain->dimension;
  int itype, inpe, iapc;

  max_element_bound_box_size[0] = 0.0;
  max_element_bound_box_size[1] = 0.0;
  max_element_bound_box_size[2] = 0.0;
  local_element_bound_box[0] = BIG;
  local_element_bound_box[1] = BIG;
  local_element_bound_box[2] = BIG;
  local_element_bound_box[3] = -BIG;
  local_element_bound_box[4] = -BIG;
  local_element_bound_box[5] = -BIG;
  for (i = 0; i < nlocal + nghost; i++) {
    itype = etype[i];
    inpe = npe[itype];
    iapc = apc[itype];
    inodex = nodex[i];

    // determine element bound box and max bound box sizes

    for (dim = 0; dim < dimension; dim++) {
      maxcoord = -BIG;
      mincoord = BIG;
      for (j = 0; j < iapc; j++) 
        for (k = 0; k < inpe; k++) {
          maxcoord = MAX(maxcoord, inodex[j][k][dim]);
          mincoord = MIN(mincoord, inodex[j][k][dim]);
        }
      element_bound_box[i][dim] = mincoord;
      element_bound_box[i][dim+3] = maxcoord;

      // skip the rest for ghost elements

      if (i >= nlocal) continue;

      element_box_size = maxcoord - mincoord;
      max_element_bound_box_size[dim] = MAX(element_box_size, max_element_bound_box_size[dim]);
      local_element_bound_box[dim] = MIN(mincoord, local_element_bound_box[dim]);
      local_element_bound_box[dim+3] = MAX(maxcoord, local_element_bound_box[dim+3]);

      // check element deformation

      if (check) {
        if (initial_box_size[i][dim] < 0.0) initial_box_size[i][dim] = element_box_size;
        else if (element_box_size/initial_box_size[i][dim] > maxelemchange || 
            initial_box_size[i][dim]/element_box_size > maxelemchange) {
          char *errstr = new char[500];
          sprintf(errstr, "Simulation unstable! Element ID " TAGINT_FORMAT " deformed by more than %3.1f times initial_size = %g, current size = %g", 
              tag[i], maxelemchange, initial_box_size[i][dim], element_box_size);
          error->one(FLERR, errstr);
          delete [] errstr;
        }
      }
    }
  }
}

/* ----------------------------------------
   update element subelem_size for building neighbor list
   must be called before building any neighbor lists
   ------------------------------------------ */

void Element::update_subelem_size()
{
  double delx, dely, delz, rsq;
  double ***inodex;
  double ucellbound[6], maxsize_ucell;
  int itype, iapc, inpe;

  // calculate subelem size
  // add new element shapes here

  for (int i = 0; i < nlocal + nghost; i++) {
    subelem_size[i] = 0.0;
    inodex = nodex[i];
    itype = etype[i];
    iapc = apc[itype];
    inpe = npe[itype];
    if (element_shape_ids[itype] == Element::QUADRILATERAL) {
      for (int j = 0; j < 2; j++) {
        delx = dely = 0.0;
        for (int k = 0; k < iapc; k++) {
          delx += inodex[k][j][0] - inodex[k][j+2][0];
          dely += inodex[k][j][1] - inodex[k][j+2][1];
        }
        delx /= iapc;
        dely /= iapc;
        rsq = delx * delx + dely * dely;
        subelem_size[i] = MAX(subelem_size[i], rsq); 
      }
    } else if (element_shape_ids[itype] == Element::TRIANGLE) {

      // use 2 vectors 0->1, 0->2 to form a parallelogram and find max diag

      delx = dely = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] + inodex[k][2][0] - 3.0 * inodex[k][0][0];
        dely += inodex[k][1][1] + inodex[k][2][1] - 3.0 * inodex[k][0][1];
      }
      delx /= iapc;
      dely /= iapc;
      rsq = delx * delx + dely * dely;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] - inodex[k][2][0];
        dely += inodex[k][1][1] - inodex[k][2][1];
      }
      delx /= iapc;
      dely /= iapc;
      rsq = delx * delx + dely * dely;
      subelem_size[i] = MAX(subelem_size[i], rsq); 
      subelem_size[i] *= 1.21;

    } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {
      for (int j = 0; j < 4; j++) {
        delx = dely = delz = 0.0;
        for (int k = 0; k < iapc; k++) 
          if (j < 2) {
            delx += inodex[k][j][0] - inodex[k][j+6][0];
            dely += inodex[k][j][1] - inodex[k][j+6][1];
            delz += inodex[k][j][2] - inodex[k][j+6][2];
          } else {
            delx += inodex[k][j][0] - inodex[k][j+2][0];
            dely += inodex[k][j][1] - inodex[k][j+2][1];
            delz += inodex[k][j][2] - inodex[k][j+2][2];
          }
        delx /= iapc;
        dely /= iapc;
        delz /= iapc;

        rsq = delx * delx + dely * dely + delz * delz;
        subelem_size[i] = MAX(subelem_size[i], rsq); 
      }

    } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

      // use 3 vectors 0->1, 0->2, 0->3 to form a rhombohedron and find max diag

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] + inodex[k][2][0] + inodex[k][3][0] - 4.0 * inodex[k][0][0];
        dely += inodex[k][1][1] + inodex[k][2][1] + inodex[k][3][1] - 4.0 * inodex[k][0][1];
        delz += inodex[k][1][2] + inodex[k][2][2] + inodex[k][3][2] - 4.0 * inodex[k][0][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] + inodex[k][2][0] - inodex[k][3][0] - inodex[k][0][0];
        dely += inodex[k][1][1] + inodex[k][2][1] - inodex[k][3][1] - inodex[k][0][1];
        delz += inodex[k][1][2] + inodex[k][2][2] - inodex[k][3][2] - inodex[k][0][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] - inodex[k][2][0] + inodex[k][3][0] - inodex[k][0][0];
        dely += inodex[k][1][1] - inodex[k][2][1] + inodex[k][3][1] - inodex[k][0][1];
        delz += inodex[k][1][2] - inodex[k][2][2] + inodex[k][3][2] - inodex[k][0][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += - inodex[k][1][0] + inodex[k][2][0] + inodex[k][3][0] - inodex[k][0][0];
        dely += - inodex[k][1][1] + inodex[k][2][1] + inodex[k][3][1] - inodex[k][0][1];
        delz += - inodex[k][1][2] + inodex[k][2][2] + inodex[k][3][2] - inodex[k][0][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 
      subelem_size[i] *= 1.21;

    } else if (element_shape_ids[itype] == Element::PYRAMID) {

      // offset the 4 corners of the base by the difference between the center of the base and 
      // the peak to form a hexahedron and then calculate the maximum diagonal distance of the hexahedron

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][4][0] + 0.75 * inodex[k][0][0] 
          - 0.25 * inodex[k][1][0] - 1.25 * inodex[k][2][0] - 0.25 * inodex[k][3][0];
        dely += inodex[k][4][1] + 0.75 * inodex[k][0][1] 
          - 0.25 * inodex[k][1][1] - 1.25 * inodex[k][2][1] - 0.25 * inodex[k][3][1];
        delz += inodex[k][4][2] + 0.75 * inodex[k][0][2] 
          - 0.25 * inodex[k][1][2] - 1.25 * inodex[k][2][2] - 0.25 * inodex[k][3][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][4][0] - 0.25 * inodex[k][0][0] 
          + 0.75 * inodex[k][1][0] - 1.25 * inodex[k][2][0] - 1.25 * inodex[k][3][0];
        dely += inodex[k][4][1] - 0.25 * inodex[k][0][1] 
          + 0.75 * inodex[k][1][1] - 1.25 * inodex[k][2][1] - 1.25 * inodex[k][3][1];
        delz += inodex[k][4][2] - 0.25 * inodex[k][0][2] 
          + 0.75 * inodex[k][1][2] - 1.25 * inodex[k][2][2] - 1.25 * inodex[k][3][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][4][0] - 1.25 * inodex[k][0][0] 
          - 0.25 * inodex[k][1][0] + 0.75 * inodex[k][2][0] - 0.25 * inodex[k][3][0];
        dely += inodex[k][4][1] - 1.25 * inodex[k][0][1] 
          - 0.25 * inodex[k][1][1] + 0.75 * inodex[k][2][1] - 0.25 * inodex[k][3][1];
        delz += inodex[k][4][2] - 1.25 * inodex[k][0][2] 
          - 0.25 * inodex[k][1][2] + 0.75 * inodex[k][2][2] - 0.25 * inodex[k][3][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][4][0] - 0.25 * inodex[k][0][0] 
          - 1.25 * inodex[k][1][0] - 0.25 * inodex[k][2][0] + 0.75 * inodex[k][3][0];
        dely += inodex[k][4][1] - 0.25 * inodex[k][0][1] 
          - 1.25 * inodex[k][1][1] - 0.25 * inodex[k][2][1] + 0.75 * inodex[k][3][1];
        delz += inodex[k][4][2] - 0.25 * inodex[k][0][2] 
          - 1.25 * inodex[k][1][2] - 0.25 * inodex[k][2][2] + 0.75 * inodex[k][3][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      subelem_size[i] *= 1.21;
    } else if (element_shape_ids[itype] == Element::WEDGE) {

      // sum 2 vectors 0->1, 0->2 to get 7th corner
      // sum 2 vectors 3->4, 3->5 to get 8th corner
      // check the max diag of the formed hexahedron

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] - inodex[k][5][0];
        dely += inodex[k][1][1] - inodex[k][5][1];
        delz += inodex[k][1][2] - inodex[k][5][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][2][0] - inodex[k][4][0];
        dely += inodex[k][2][1] - inodex[k][4][1];
        delz += inodex[k][2][2] - inodex[k][4][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][1][0] + inodex[k][2][0] - 2.0 * inodex[k][0][0] - inodex[k][3][0];
        delx += inodex[k][1][1] + inodex[k][2][1] - 2.0 * inodex[k][0][1] - inodex[k][3][1];
        delx += inodex[k][1][2] + inodex[k][2][2] - 2.0 * inodex[k][0][2] - inodex[k][3][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 

      delx = dely = delz = 0.0;
      for (int k = 0; k < iapc; k++) {
        delx += inodex[k][4][0] + inodex[k][5][0] - 2.0 * inodex[k][3][0] - inodex[k][0][0];
        delx += inodex[k][4][1] + inodex[k][5][1] - 2.0 * inodex[k][3][1] - inodex[k][0][1];
        delx += inodex[k][4][2] + inodex[k][5][2] - 2.0 * inodex[k][3][2] - inodex[k][0][2];
      }
      delx /= iapc;
      dely /= iapc;
      delz /= iapc;
      rsq = delx * delx + dely * dely + delz * delz;
      subelem_size[i] = MAX(subelem_size[i], rsq); 
      subelem_size[i] *= 1.21;
    }

    maxsize_ucell = -BIG;
    for (int j = 0; j < inpe; j++) {
      ucellbound[0] = ucellbound[1] = ucellbound[2] = BIG;
      ucellbound[3] = ucellbound[4] = ucellbound[5] = -BIG;
      for (int k = 0; k < iapc; k++) 
        for (int dim = 0; dim < 3; dim++) {
          ucellbound[dim] = MIN(ucellbound[dim], inodex[k][j][dim]);
          ucellbound[dim+3] = MAX(ucellbound[dim+3], inodex[k][j][dim]);
        }
      delx = ucellbound[3] - ucellbound[0];
      dely = ucellbound[4] - ucellbound[1];
      delz = ucellbound[5] - ucellbound[2];
      rsq = delx * delx + dely * dely + delz * delz;
      maxsize_ucell = MAX(sqrt(rsq), maxsize_ucell);
    }
    subelem_size[i] = (sqrt(subelem_size[i]) / subsplit[itype] / 2.0 + maxsize_ucell) * subelem_size_factor;
  }
}

/* ----------------------------------------
   define inversed consistent mass matrix for each element shape 
   ------------------------------------------ */

void Element::define_inversed_mass_matrices()
{
  // QUADRILATERAL element

  inversed_mass_matrices[QUADRILATERAL][0][0] = 16.0;
  inversed_mass_matrices[QUADRILATERAL][0][1] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][0][2] = 4.0;
  inversed_mass_matrices[QUADRILATERAL][0][3] = -8.0;

  inversed_mass_matrices[QUADRILATERAL][1][0] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][1][1] = 16.0;
  inversed_mass_matrices[QUADRILATERAL][1][2] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][1][3] = 4.0;

  inversed_mass_matrices[QUADRILATERAL][2][0] = 4.0;
  inversed_mass_matrices[QUADRILATERAL][2][1] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][2][2] = 16.0;
  inversed_mass_matrices[QUADRILATERAL][2][3] = -8.0;

  inversed_mass_matrices[QUADRILATERAL][3][0] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][3][1] = 4.0;
  inversed_mass_matrices[QUADRILATERAL][3][2] = -8.0;
  inversed_mass_matrices[QUADRILATERAL][3][3] = 16.0;

  // TRIANGLE element

  inversed_mass_matrices[TRIANGLE][0][0] = 9.0;
  inversed_mass_matrices[TRIANGLE][0][1] = -3.0;
  inversed_mass_matrices[TRIANGLE][0][2] = -3.0;

  inversed_mass_matrices[TRIANGLE][1][0] = -3.0;
  inversed_mass_matrices[TRIANGLE][1][1] = 9.0;
  inversed_mass_matrices[TRIANGLE][1][2] = -3.0;

  inversed_mass_matrices[TRIANGLE][2][0] = -3.0;
  inversed_mass_matrices[TRIANGLE][2][1] = -3.0;
  inversed_mass_matrices[TRIANGLE][2][2] = 9.0;

  // HEXAHEDRON element

  inversed_mass_matrices[HEXAHEDRON][0][0] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][0][1] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][0][2] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][0][3] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][0][4] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][0][5] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][0][6] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][0][7] = 16.0;

  inversed_mass_matrices[HEXAHEDRON][1][0] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][1][1] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][1][2] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][1][3] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][1][4] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][1][5] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][1][6] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][1][7] = -8.0;

  inversed_mass_matrices[HEXAHEDRON][2][0] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][2][1] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][2][2] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][2][3] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][2][4] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][2][5] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][2][6] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][2][7] = 16.0;

  inversed_mass_matrices[HEXAHEDRON][3][0] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][3][1] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][3][2] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][3][3] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][3][4] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][3][5] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][3][6] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][3][7] = -32.0;

  inversed_mass_matrices[HEXAHEDRON][4][0] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][4][1] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][4][2] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][4][3] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][4][4] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][4][5] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][4][6] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][4][7] = -32.0;

  inversed_mass_matrices[HEXAHEDRON][5][0] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][5][1] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][5][2] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][5][3] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][5][4] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][5][5] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][5][6] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][5][7] = 16.0;

  inversed_mass_matrices[HEXAHEDRON][6][0] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][6][1] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][6][2] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][6][3] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][6][4] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][6][5] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][6][6] = 64.0;
  inversed_mass_matrices[HEXAHEDRON][6][7] = -32.0;

  inversed_mass_matrices[HEXAHEDRON][7][0] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][7][1] = -8.0;
  inversed_mass_matrices[HEXAHEDRON][7][2] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][7][3] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][7][4] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][7][5] = 16.0;
  inversed_mass_matrices[HEXAHEDRON][7][6] = -32.0;
  inversed_mass_matrices[HEXAHEDRON][7][7] = 64.0;

  // TETRAHEDRON element

  inversed_mass_matrices[TETRAHEDRON][0][0] = 16.0;
  inversed_mass_matrices[TETRAHEDRON][0][1] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][0][2] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][0][3] = -4.0;

  inversed_mass_matrices[TETRAHEDRON][1][0] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][1][1] = 16.0;
  inversed_mass_matrices[TETRAHEDRON][1][2] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][1][3] = -4.0;

  inversed_mass_matrices[TETRAHEDRON][2][0] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][2][1] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][2][2] = 16.0;
  inversed_mass_matrices[TETRAHEDRON][2][3] = -4.0;

  inversed_mass_matrices[TETRAHEDRON][3][0] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][3][1] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][3][2] = -4.0;
  inversed_mass_matrices[TETRAHEDRON][3][3] = 16.0;

  // PYRAMID element

  inversed_mass_matrices[PYRAMID][0][0] = 83.0/3.0;
  inversed_mass_matrices[PYRAMID][0][1] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][0][2] = 23.0/3.0;
  inversed_mass_matrices[PYRAMID][0][3] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][0][4] = -4.0;

  inversed_mass_matrices[PYRAMID][1][0] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][1][1] = 83.0/3.0;
  inversed_mass_matrices[PYRAMID][1][2] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][1][3] = 23.0/3.0;
  inversed_mass_matrices[PYRAMID][1][4] = -4.0;

  inversed_mass_matrices[PYRAMID][2][0] = 23.0/3.0;
  inversed_mass_matrices[PYRAMID][2][1] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][2][2] = 83.0/3.0;
  inversed_mass_matrices[PYRAMID][2][3] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][2][4] = -4.0;

  inversed_mass_matrices[PYRAMID][3][0] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][3][1] = 23.0/3.0;
  inversed_mass_matrices[PYRAMID][3][2] = -37.0/3.0;
  inversed_mass_matrices[PYRAMID][3][3] = 83.0/3.0;
  inversed_mass_matrices[PYRAMID][3][4] = -4.0;

  inversed_mass_matrices[PYRAMID][4][0] = -4.0;
  inversed_mass_matrices[PYRAMID][4][1] = -4.0;
  inversed_mass_matrices[PYRAMID][4][2] = -4.0;
  inversed_mass_matrices[PYRAMID][4][3] = -4.0;
  inversed_mass_matrices[PYRAMID][4][4] = 12.0;

  // WEDGE element

  inversed_mass_matrices[WEDGE][0][0] = 36.0;
  inversed_mass_matrices[WEDGE][0][1] = -12.0;
  inversed_mass_matrices[WEDGE][0][2] = -12.0;
  inversed_mass_matrices[WEDGE][0][3] = -18.0;
  inversed_mass_matrices[WEDGE][0][4] = 6.0;
  inversed_mass_matrices[WEDGE][0][5] = 6.0;

  inversed_mass_matrices[WEDGE][1][0] = -12.0;
  inversed_mass_matrices[WEDGE][1][1] = 36.0; 
  inversed_mass_matrices[WEDGE][1][2] = -12.0;
  inversed_mass_matrices[WEDGE][1][3] = 6.0;
  inversed_mass_matrices[WEDGE][1][4] = -18.0;
  inversed_mass_matrices[WEDGE][1][5] = 6.0;

  inversed_mass_matrices[WEDGE][2][0] = -12.0;
  inversed_mass_matrices[WEDGE][2][1] = -12.0;
  inversed_mass_matrices[WEDGE][2][2] = 36.0;
  inversed_mass_matrices[WEDGE][2][3] = 6.0;
  inversed_mass_matrices[WEDGE][2][4] = 6.0;
  inversed_mass_matrices[WEDGE][2][5] = -18.0;

  inversed_mass_matrices[WEDGE][3][0] = -18.0;
  inversed_mass_matrices[WEDGE][3][1] = 6.0;
  inversed_mass_matrices[WEDGE][3][2] = 6.0;
  inversed_mass_matrices[WEDGE][3][3] = 36.0;
  inversed_mass_matrices[WEDGE][3][4] = -12.0;
  inversed_mass_matrices[WEDGE][3][5] = -12.0;

  inversed_mass_matrices[WEDGE][4][0] = 6.0;
  inversed_mass_matrices[WEDGE][4][1] = -18.0;
  inversed_mass_matrices[WEDGE][4][2] = 6.0;
  inversed_mass_matrices[WEDGE][4][3] = -12.0;
  inversed_mass_matrices[WEDGE][4][4] = 36.0;
  inversed_mass_matrices[WEDGE][4][5] = -12.0;

  inversed_mass_matrices[WEDGE][5][0] = 6.0;
  inversed_mass_matrices[WEDGE][5][1] = 6.0;
  inversed_mass_matrices[WEDGE][5][2] = -18.0;
  inversed_mass_matrices[WEDGE][5][3] = -12.0;
  inversed_mass_matrices[WEDGE][5][4] = -12.0;
  inversed_mass_matrices[WEDGE][5][5] = 36.0;
}

/* ----------------------------------------
   define consistent mass matrix for each element shape 
   ------------------------------------------ */

void Element::define_mass_matrices()
{
  // QUADRILATERAL element

  mass_matrices[QUADRILATERAL][0][0] = 1.0 / 9.0;
  mass_matrices[QUADRILATERAL][0][1] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][0][2] = 1.0 / 36.0;
  mass_matrices[QUADRILATERAL][0][3] = 1.0 / 18.0;

  mass_matrices[QUADRILATERAL][1][0] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][1][1] = 1.0 / 9.0;
  mass_matrices[QUADRILATERAL][1][2] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][1][3] = 1.0 / 36.0;

  mass_matrices[QUADRILATERAL][2][0] = 1.0 / 36.0;
  mass_matrices[QUADRILATERAL][2][1] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][2][2] = 1.0 / 9.0;
  mass_matrices[QUADRILATERAL][2][3] = 1.0 / 18.0;

  mass_matrices[QUADRILATERAL][3][0] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][3][1] = 1.0 / 36.0;
  mass_matrices[QUADRILATERAL][3][2] = 1.0 / 18.0;
  mass_matrices[QUADRILATERAL][3][3] = 1.0 / 9.0;

  // TRIANGLE element

  mass_matrices[TRIANGLE][0][0] = 1.0 / 6.0;
  mass_matrices[TRIANGLE][0][1] = 1.0 / 12.0;
  mass_matrices[TRIANGLE][0][2] = 1.0 / 12.0;

  mass_matrices[TRIANGLE][1][0] = 1.0 / 12.0;
  mass_matrices[TRIANGLE][1][1] = 1.0 / 6.0;
  mass_matrices[TRIANGLE][1][2] = 1.0 / 12.0;

  mass_matrices[TRIANGLE][2][0] = 1.0 / 12.0;
  mass_matrices[TRIANGLE][2][1] = 1.0 / 12.0;
  mass_matrices[TRIANGLE][2][2] = 1.0 / 6.0;

  // HEXAHEDRON element

  mass_matrices[HEXAHEDRON][0][0] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][0][1] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][0][2] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][0][3] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][0][4] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][0][5] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][0][6] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][0][7] = 1.0 / 108.0;

  mass_matrices[HEXAHEDRON][1][0] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][1][1] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][1][2] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][1][3] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][1][4] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][1][5] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][1][6] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][1][7] = 1.0 / 216.0;

  mass_matrices[HEXAHEDRON][2][0] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][2][1] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][2][2] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][2][3] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][2][4] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][2][5] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][2][6] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][2][7] = 1.0 / 108.0;

  mass_matrices[HEXAHEDRON][3][0] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][3][1] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][3][2] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][3][3] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][3][4] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][3][5] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][3][6] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][3][7] = 1.0 / 54.0;

  mass_matrices[HEXAHEDRON][4][0] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][4][1] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][4][2] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][4][3] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][4][4] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][4][5] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][4][6] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][4][7] = 1.0 / 54.0;

  mass_matrices[HEXAHEDRON][5][0] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][5][1] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][5][2] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][5][3] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][5][4] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][5][5] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][5][6] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][5][7] = 1.0 / 108.0;

  mass_matrices[HEXAHEDRON][6][0] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][6][1] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][6][2] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][6][3] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][6][4] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][6][5] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][6][6] = 1.0 / 27.0;
  mass_matrices[HEXAHEDRON][6][7] = 1.0 / 54.0;

  mass_matrices[HEXAHEDRON][7][0] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][7][1] = 1.0 / 216.0;
  mass_matrices[HEXAHEDRON][7][2] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][7][3] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][7][4] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][7][5] = 1.0 / 108.0;
  mass_matrices[HEXAHEDRON][7][6] = 1.0 / 54.0;
  mass_matrices[HEXAHEDRON][7][7] = 1.0 / 27.0;

  // TETRAHEDRON element

  mass_matrices[TETRAHEDRON][0][0] = 1.0 / 10.0;
  mass_matrices[TETRAHEDRON][0][1] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][0][2] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][0][3] = 1.0 / 20.0;

  mass_matrices[TETRAHEDRON][1][0] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][1][1] = 1.0 / 10.0;
  mass_matrices[TETRAHEDRON][1][2] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][1][3] = 1.0 / 20.0;

  mass_matrices[TETRAHEDRON][2][0] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][2][1] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][2][2] = 1.0 / 10.0;
  mass_matrices[TETRAHEDRON][2][3] = 1.0 / 20.0;

  mass_matrices[TETRAHEDRON][3][0] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][3][1] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][3][2] = 1.0 / 20.0;
  mass_matrices[TETRAHEDRON][3][3] = 1.0 / 10.0;

  // PYRAMID element

  mass_matrices[PYRAMID][0][0] = 1.0 / 15.0;
  mass_matrices[PYRAMID][0][1] = 1.0 / 30.0;
  mass_matrices[PYRAMID][0][2] = 1.0 / 60.0;
  mass_matrices[PYRAMID][0][3] = 1.0 / 30.0;
  mass_matrices[PYRAMID][0][4] = 3.0 / 80.0;

  mass_matrices[PYRAMID][1][0] = 1.0 / 30.0;
  mass_matrices[PYRAMID][1][1] = 1.0 / 15.0;
  mass_matrices[PYRAMID][1][2] = 1.0 / 30.0;
  mass_matrices[PYRAMID][1][3] = 1.0 / 60.0;
  mass_matrices[PYRAMID][1][4] = 3.0 / 80.0;

  mass_matrices[PYRAMID][2][0] = 1.0 / 60.0;
  mass_matrices[PYRAMID][2][1] = 1.0 / 30.0;
  mass_matrices[PYRAMID][2][2] = 1.0 / 15.0;
  mass_matrices[PYRAMID][2][3] = 1.0 / 30.0;
  mass_matrices[PYRAMID][2][4] = 3.0 / 80.0;

  mass_matrices[PYRAMID][3][0] = 1.0 / 30.0;
  mass_matrices[PYRAMID][3][1] = 1.0 / 60.0;
  mass_matrices[PYRAMID][3][2] = 1.0 / 30.0;
  mass_matrices[PYRAMID][3][3] = 1.0 / 15.0;
  mass_matrices[PYRAMID][3][4] = 3.0 / 80.0;

  mass_matrices[PYRAMID][4][0] = 3.0 / 80.0;
  mass_matrices[PYRAMID][4][1] = 3.0 / 80.0;
  mass_matrices[PYRAMID][4][2] = 3.0 / 80.0;
  mass_matrices[PYRAMID][4][3] = 3.0 / 80.0;
  mass_matrices[PYRAMID][4][4] = 1.0 / 10.0;

  // WEDGE element

  mass_matrices[WEDGE][0][0] = 1.0 / 18.0;
  mass_matrices[WEDGE][0][1] = 1.0 / 36.0;
  mass_matrices[WEDGE][0][2] = 1.0 / 36.0;
  mass_matrices[WEDGE][0][3] = 1.0 / 36.0;
  mass_matrices[WEDGE][0][4] = 1.0 / 72.0;
  mass_matrices[WEDGE][0][5] = 1.0 / 72.0;

  mass_matrices[WEDGE][1][0] = 1.0 / 36.0;
  mass_matrices[WEDGE][1][1] = 1.0 / 18.0;
  mass_matrices[WEDGE][1][2] = 1.0 / 36.0;
  mass_matrices[WEDGE][1][3] = 1.0 / 72.0;
  mass_matrices[WEDGE][1][4] = 1.0 / 36.0;
  mass_matrices[WEDGE][1][5] = 1.0 / 72.0;

  mass_matrices[WEDGE][2][0] = 1.0 / 36.0;
  mass_matrices[WEDGE][2][1] = 1.0 / 36.0;
  mass_matrices[WEDGE][2][2] = 1.0 / 18.0;
  mass_matrices[WEDGE][2][3] = 1.0 / 72.0;
  mass_matrices[WEDGE][2][4] = 1.0 / 72.0;
  mass_matrices[WEDGE][2][5] = 1.0 / 36.0;

  mass_matrices[WEDGE][3][0] = 1.0 / 36.0;
  mass_matrices[WEDGE][3][1] = 1.0 / 72.0;
  mass_matrices[WEDGE][3][2] = 1.0 / 72.0;
  mass_matrices[WEDGE][3][3] = 1.0 / 18.0;
  mass_matrices[WEDGE][3][4] = 1.0 / 36.0;
  mass_matrices[WEDGE][3][5] = 1.0 / 36.0;

  mass_matrices[WEDGE][4][0] = 1.0 / 72.0;
  mass_matrices[WEDGE][4][1] = 1.0 / 36.0;
  mass_matrices[WEDGE][4][2] = 1.0 / 72.0;
  mass_matrices[WEDGE][4][3] = 1.0 / 36.0;
  mass_matrices[WEDGE][4][4] = 1.0 / 18.0;
  mass_matrices[WEDGE][4][5] = 1.0 / 36.0;

  mass_matrices[WEDGE][5][0] = 1.0 / 72.0;
  mass_matrices[WEDGE][5][1] = 1.0 / 72.0;
  mass_matrices[WEDGE][5][2] = 1.0 / 36.0;
  mass_matrices[WEDGE][5][3] = 1.0 / 36.0;
  mass_matrices[WEDGE][5][4] = 1.0 / 36.0;
  mass_matrices[WEDGE][5][5] = 1.0 / 18.0;

}

/*  ----------------------------------------------------------------------
    determine nodef for all local element nodes
    should be called after reverse_comm 
    Consistent mass style: compute nodef from force_columns and inversed mass matrix
    Lumped mass style: tally force columns to nodef directly
    -------------------------------------------------------------------------  */

void Element::compute_nodef()
{
  int ietype, inpe, iapc, ingcell;

  if (nodal_force_style == CONSISTENT) {
    double **force_columns;
    memory->create(force_columns, maxnpe, 3, "element:force_columns");

    for (int i = 0; i < nlocal; i++) {
      ietype = etype[i];
      inpe = npe[ietype];
      iapc = apc[ietype];
      ingcell = ngcell[ietype];
      double **inversed_mass_matrix = inversed_mass_matrices[element_shape_ids[ietype]];
      for (int ibasis = 0; ibasis < iapc; ibasis++) {
        memset(&force_columns[0][0], 0, sizeof(double) * 3 * maxnpe);
        for (int igcell = 0; igcell < ingcell; igcell++) 
          for (int node = 0; node < inpe; node++) {
            force_columns[node][0] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][0];
            force_columns[node][1] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][1];
            force_columns[node][2] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][2];
          }
        for (int inode = 0; inode < inpe; inode++)
          for (int j = 0; j < inpe; j++) {
            nodef[i][ibasis][inode][0] += inversed_mass_matrix[inode][j] * force_columns[j][0];
            nodef[i][ibasis][inode][1] += inversed_mass_matrix[inode][j] * force_columns[j][1];
            nodef[i][ibasis][inode][2] += inversed_mass_matrix[inode][j] * force_columns[j][2];
          }
      } 
    } 
    memory->destroy(force_columns);
  } else if (nodal_force_style == LUMPED) {
    for (int i = 0; i < nlocal; i++) {
      ietype = etype[i];
      inpe = npe[ietype];
      iapc = apc[ietype];
      ingcell = ngcell[ietype];
      for (int ibasis = 0; ibasis < iapc; ibasis++) 
        for (int igcell = 0; igcell < ingcell; igcell++) 
          for (int node = 0; node < inpe; node++) {
            nodef[i][ibasis][node][0] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][0];
            nodef[i][ibasis][node][1] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][1];
            nodef[i][ibasis][node][2] += 
              weighted_shape_array[ietype][igcell][node] * gaussf[i][ibasis][igcell][2];
          }
    }
  } else if (nodal_force_style == NODEONLY) {
    for (int i = 0; i < nlocal; i++) {
      ietype = etype[i];
      inpe = npe[ietype];
      iapc = apc[ietype];
      int *in2g = n2g[ietype];
      for (int ibasis = 0; ibasis < iapc; ibasis++) 
        for (int inode = 0; inode < inpe; inode++) {
          int igcell = in2g[inode];
          nodef[i][ibasis][inode][0] = gaussf[i][ibasis][igcell][0];
          nodef[i][ibasis][inode][1] = gaussf[i][ibasis][igcell][1];
          nodef[i][ibasis][inode][2] = gaussf[i][ibasis][igcell][2];
        }
    }
  } else {

  }
}

/*  ---------------------------------------------
    Find element shape id from element shape name
    --------------------------------------------- */

int Element::find_element_shape_id(char *element_shape_name) 
{
  for (int i = 0; i < num_element_shapes; i++) {
    if (strcmp(element_shape_name, element_shape_list[i]) == 0)
      return i;
  }
  return -1;
}

/*  ---------------------------------------------
    Define surface node list
    In consistent with the index of surface in surface_ucell in evec->set_interpolate_element() function
    Node indices are placed in the order following right-hand rule pointing outward
    --------------------------------------------- */

void Element::define_surface_node_list()
{
  // HEXAHEDRON element

  // -x 

  surface_node_list[HEXAHEDRON][0][0] = 0;
  surface_node_list[HEXAHEDRON][0][1] = 4;
  surface_node_list[HEXAHEDRON][0][2] = 7;
  surface_node_list[HEXAHEDRON][0][3] = 3;
  surface_num_node[HEXAHEDRON][0] = 4;

  // +x 

  surface_node_list[HEXAHEDRON][1][0] = 1;
  surface_node_list[HEXAHEDRON][1][1] = 2;
  surface_node_list[HEXAHEDRON][1][2] = 6;
  surface_node_list[HEXAHEDRON][1][3] = 5;
  surface_num_node[HEXAHEDRON][1] = 4;

  // -y 

  surface_node_list[HEXAHEDRON][2][0] = 0;
  surface_node_list[HEXAHEDRON][2][1] = 1;
  surface_node_list[HEXAHEDRON][2][2] = 5;
  surface_node_list[HEXAHEDRON][2][3] = 4;
  surface_num_node[HEXAHEDRON][2] = 4;

  // +y 

  surface_node_list[HEXAHEDRON][3][0] = 3;
  surface_node_list[HEXAHEDRON][3][1] = 7;
  surface_node_list[HEXAHEDRON][3][2] = 6;
  surface_node_list[HEXAHEDRON][3][3] = 2;
  surface_num_node[HEXAHEDRON][3] = 4;

  // -z 

  surface_node_list[HEXAHEDRON][4][0] = 0;
  surface_node_list[HEXAHEDRON][4][1] = 3;
  surface_node_list[HEXAHEDRON][4][2] = 2;
  surface_node_list[HEXAHEDRON][4][3] = 1;
  surface_num_node[HEXAHEDRON][4] = 4;

  // +z 

  surface_node_list[HEXAHEDRON][5][0] = 4;
  surface_node_list[HEXAHEDRON][5][1] = 5;
  surface_node_list[HEXAHEDRON][5][2] = 6;
  surface_node_list[HEXAHEDRON][5][3] = 7;
  surface_num_node[HEXAHEDRON][5] = 4;

  // TETRAHEDRON element

  // xy plane 

  surface_node_list[TETRAHEDRON][0][0] = 0;
  surface_node_list[TETRAHEDRON][0][1] = 2;
  surface_node_list[TETRAHEDRON][0][2] = 1;
  surface_num_node[TETRAHEDRON][0] = 3;

  // yz plane 

  surface_node_list[TETRAHEDRON][1][0] = 0;
  surface_node_list[TETRAHEDRON][1][1] = 3;
  surface_node_list[TETRAHEDRON][1][2] = 2;
  surface_num_node[TETRAHEDRON][1] = 3;

  // xz plane 

  surface_node_list[TETRAHEDRON][2][0] = 0;
  surface_node_list[TETRAHEDRON][2][1] = 1;
  surface_node_list[TETRAHEDRON][2][2] = 3;
  surface_num_node[TETRAHEDRON][2] = 3;

  // x + y + z = 1 plane 

  surface_node_list[TETRAHEDRON][3][0] = 1;
  surface_node_list[TETRAHEDRON][3][1] = 2;
  surface_node_list[TETRAHEDRON][3][2] = 3;
  surface_num_node[TETRAHEDRON][3] = 3;

  // PYRAMID element

  // xy plane 

  surface_node_list[PYRAMID][0][0] = 3;
  surface_node_list[PYRAMID][0][1] = 2;
  surface_node_list[PYRAMID][0][2] = 1;
  surface_node_list[PYRAMID][0][3] = 0;
  surface_num_node[PYRAMID][0] = 4;

  // -x plane 

  surface_node_list[PYRAMID][1][0] = 0;
  surface_node_list[PYRAMID][1][1] = 4;
  surface_node_list[PYRAMID][1][2] = 3;
  surface_num_node[PYRAMID][1] = 3;

  // -y plane 

  surface_node_list[PYRAMID][2][0] = 0;
  surface_node_list[PYRAMID][2][1] = 1;
  surface_node_list[PYRAMID][2][2] = 4;
  surface_num_node[PYRAMID][2] = 3;

  // +x plane 

  surface_node_list[PYRAMID][3][0] = 1;
  surface_node_list[PYRAMID][3][1] = 2;
  surface_node_list[PYRAMID][3][2] = 4;
  surface_num_node[PYRAMID][3] = 3;

  // +y plane 

  surface_node_list[PYRAMID][4][0] = 3;
  surface_node_list[PYRAMID][4][1] = 4;
  surface_node_list[PYRAMID][4][2] = 2;
  surface_num_node[PYRAMID][4] = 3;

  // WEDGE element

  // +yz plane 

  surface_node_list[WEDGE][0][0] = 3;
  surface_node_list[WEDGE][0][1] = 4;
  surface_node_list[WEDGE][0][2] = 5;
  surface_num_node[WEDGE][0] = 3;

  // -yz plane 

  surface_node_list[WEDGE][1][0] = 0;
  surface_node_list[WEDGE][1][1] = 2;
  surface_node_list[WEDGE][1][2] = 1;
  surface_num_node[WEDGE][1] = 3;

  // xy plane 

  surface_node_list[WEDGE][2][0] = 0;
  surface_node_list[WEDGE][2][1] = 1;
  surface_node_list[WEDGE][2][2] = 4;
  surface_node_list[WEDGE][2][3] = 3;
  surface_num_node[WEDGE][2] = 4;

  // xz plane 

  surface_node_list[WEDGE][3][0] = 0;
  surface_node_list[WEDGE][3][1] = 3;
  surface_node_list[WEDGE][3][2] = 5;
  surface_node_list[WEDGE][3][3] = 2;
  surface_num_node[WEDGE][3] = 4;

  // y + z = 1 plane 

  surface_node_list[WEDGE][4][0] = 1;
  surface_node_list[WEDGE][4][1] = 2;
  surface_node_list[WEDGE][4][2] = 5;
  surface_node_list[WEDGE][4][3] = 4;
  surface_num_node[WEDGE][4] = 4;

}

/*  ---------------------------------------------
    Specify node connectivity for each element shape
    --------------------------------------------- */

int Element::node_connectivity(int *node_connect, int shape_id, int *node_ids)
{
  int m = 0;
  if (shape_id == Element::QUADRILATERAL) {
    for (int j = 0; j < 4; j++) 
      node_connect[m++] = node_ids[j];
    if (domain->dimension == 3) {
      for (int j = 0; j < 4; j++) 
        node_connect[m++] = node_ids[j];
    }

  } else if (shape_id == Element::TRIANGLE) {

    // node connectivity scheme: 1 2 3 3

    for (int j = 0; j < 3; j++) 
      node_connect[m++] = node_ids[j];
    node_connect[m++] = node_ids[2];

  } else if (shape_id == Element::HEXAHEDRON) {

    for (int j = 0; j < 8; j++) 
      node_connect[m++] = node_ids[j];

  } else if (shape_id == Element::PYRAMID) {

    // node connectivity scheme: 1 2 3 4 5 5 5 5

    for (int j = 0; j < 4; j++) 
      node_connect[m++] = node_ids[j];
    for (int j = 4; j < 8; j++)
      node_connect[m++] = node_ids[4];

  } else if (shape_id == Element::TETRAHEDRON) {

    // node connectivity scheme: 1 2 3 3 4 4 4 4

    for (int j = 0; j < 3; j++) 
      node_connect[m++] = node_ids[j];
    node_connect[m++] = node_ids[2];
    for (int j = 4; j < 8; j++) 
      node_connect[m++] = node_ids[3];
  } else if (shape_id == Element::WEDGE) {

    // node connectivity scheme: 1 2 3 3 4 5 6 6

    for (int j = 0; j < 3; j++) 
      node_connect[m++] = node_ids[j];
    for (int j = 3; j < 7; j++) 
      node_connect[m++] = node_ids[j-1];
    node_connect[m++] = node_ids[5];
  }
  return m;
}
/*  ---------------------------------------------
    concert box coords to element I natural coords (xi)
    xi = H * (coord - x_center)
    --------------------------------------------- */

void Element::box2natural(double *coord, double *xi, int i)
{
  double delta[3];
  delta[0] = coord[0] - x[i][0];
  delta[1] = coord[1] - x[i][1];
  delta[2] = coord[2] - x[i][2];
  double *h = element_box2xi[i]; 
  xi[0] = h[0] * delta[0] + h[1] * delta[1] + h[2] * delta[2];
  xi[1] = h[3] * delta[0] + h[4] * delta[1] + h[5] * delta[2];
  xi[2] = h[6] * delta[0] + h[7] * delta[1] + h[8] * delta[2];
}

/*  ------------------------------------------------------------------------ */

void Element::compress_etype()
{
  int *etype_flag = new int[netypes+1];
  int *etype_flag_all = new int[netypes+1];
  int netypes_new;
  int netypes_old = netypes;
  int **old_type_info;
  memory->create(old_type_info, netypes + 1, 8, "element:old_type_info");

  for (int i = 0; i <= netypes; i++)
    etype_flag[i] = 0;

  for (int i = 0; i < nlocal; i++) 
    etype_flag[etype[i]] = 1;
  MPI_Allreduce(etype_flag, etype_flag_all, netypes + 1, MPI_INT, MPI_MAX, world);

  for (int i = 1; i <= netypes; i++)
    netypes_new += etype_flag_all[i];

  // store old etype info 

  for (int i = 1; i <= netypes; i++) {
    old_type_info[i][0] = evec->ncells[i][0];
    old_type_info[i][1] = evec->ncells[i][1];
    old_type_info[i][2] = evec->ncells[i][2];
    old_type_info[i][3] = evec->ngcells[i][0];
    old_type_info[i][4] = evec->ngcells[i][1];
    old_type_info[i][5] = evec->ngcells[i][2];
    old_type_info[i][6] = element_shape_ids[i];
    old_type_info[i][7] = apc[i];
  }

  // destroy old etype data structure

  evec->grow_etype_arrays(0);

  // add new etypes
  // request new etypes

  for (int i = 1; i <= netypes_old; i++) {
    if (etype_flag[i]) {
      evec->request_new_etype(old_type_info[i]);
    }
  }

  evec->add_requested_etype(1);

  // reassign etype to elements (including ghost elements)

  for (int i = 0; i < nlocal+nghost; i++) {
    int ietype = evec->find_etype(old_type_info[etype[i]]);
    if (ietype)
      etype[i] = ietype;
    else {
      error->one(FLERR, "Cannot find matching etype");

    }
  }

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen, "Deleted %d etypes:\n", netypes_old-netypes_new);
      for (int i = 1; i <= netypes_old; i++)
        if (!etype_flag_all[i])
          fprintf(screen, "  %d %s %d %d %d %d %d %d %d\n", i
              , element_shape_list[old_type_info[i][6]], old_type_info[i][7]
              , old_type_info[i][0], old_type_info[i][1], old_type_info[i][2]
              , old_type_info[i][3], old_type_info[i][4], old_type_info[i][5]);
      fprintf(screen, "New etypes:\n");
      for (int i = 1; i <= netypes; i++)
        fprintf(screen, "  %d %s %d %d %d %d %d %d %d\n", i
            , element_shape_list[element_shape_ids[i]], apc[i]
            , evec->ncells[i][0], evec->ncells[i][1], evec->ncells[i][2]
            , evec->ngcells[i][0], evec->ngcells[i][1], evec->ngcells[i][2]);
    }
    if (logfile) {
      fprintf(logfile, "Deleted %d etypes:\n", netypes_old-netypes_new);
      for (int i = 1; i <= netypes_old; i++)
        if (!etype_flag_all[i])
          fprintf(logfile, "  %d %s %d %d %d %d %d %d %d\n", i
              , element_shape_list[old_type_info[i][6]], old_type_info[i][7]
              , old_type_info[i][0], old_type_info[i][1], old_type_info[i][2]
              , old_type_info[i][3], old_type_info[i][4], old_type_info[i][5]);
      fprintf(logfile, "New etypes:\n");
      for (int i = 1; i <= netypes; i++)
        fprintf(logfile, "  %d %s %d %d %d %d %d %d %d\n", i
            , element_shape_list[element_shape_ids[i]], apc[i]
            , evec->ncells[i][0], evec->ncells[i][1], evec->ncells[i][2]
            , evec->ngcells[i][0], evec->ngcells[i][1], evec->ngcells[i][2]);

    }
  }
  memory->destroy(old_type_info);
  delete [] etype_flag;
  delete [] etype_flag_all;
}

/*  ------------------------------------------------------------------------- */

void Element::destroy_etype_arrays()
{
  memory->destroy(npe);
  memory->destroy(apc);
  memory->destroy(nodal_weight);
  memory->destroy(element_shape_ids);
  memory->destroy(element_shape_names);
  memory->destroy(shape_array);
  memory->destroy(weighted_shape_array);
  memory->destroy(nucell);
  memory->destroy(surface_ucell);
  memory->destroy(nsurface_ucell);
  memory->destroy(edge_ucell);
  memory->destroy(nedge_ucell);
  memory->destroy(ngcell);
  memory->destroy(weight);
  memory->destroy(g2u);
  memory->destroy(u2g);
  memory->destroy(g2n);
  memory->destroy(is_outer);
  memory->destroy(n2g);
  memory->destroy(n2u);
  memory->destroy(us2u);
  memory->destroy(nucell_subelem);
  memory->destroy(shape_array_center_subelem);
  memory->destroy(nsubelem);
  memory->destroy(subsplit);

  memory->destroy(debug_gcell_type);
  memory->destroy(debug_gcell_inactive);

}

/*  ------------------------------------------------------------------------- */
// debug command
void Element::debug_element_gcell_off(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "debug_element_gcell command before simulation box is defined");
  if (narg < 2) error->all(FLERR, "Invalid debug_element_gcell command");
  debug_mode = 1;
  int ietype = universe->inumeric(FLERR, arg[0]); // ietype = 0 means applying to all element types
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "inner") == 0) {
      if (ietype)
        debug_gcell_inactive[ietype][0] = 1;
      else 
        for (int itype = 1; itype <= netypes; itype++)
          debug_gcell_inactive[itype][0] = 1;
    } else if (strcmp(arg[iarg], "surface") == 0) {
      if (ietype)
        debug_gcell_inactive[ietype][1] = 1;
      else 
        for (int itype = 1; itype <= netypes; itype++)
          debug_gcell_inactive[itype][1] = 1;
    } else if (strcmp(arg[iarg], "edge") == 0) {
      if (ietype)
        debug_gcell_inactive[ietype][2] = 1;
      else 
        for (int itype = 1; itype <= netypes; itype++)
          debug_gcell_inactive[itype][2] = 1;
    } else if (strcmp(arg[iarg], "node") == 0) {
      if (ietype)
        debug_gcell_inactive[ietype][3] = 1;
      else 
        for (int itype = 1; itype <= netypes; itype++)
          debug_gcell_inactive[itype][3] = 1;
    } else error->all(FLERR, "Invalid debug_element_gcell command");
    iarg++;
  }
}

/* ----------------------------------------------------------------------
   unpack n lines from Element section of data file
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_elements_V2(int n, char *buf, tagint id_offset, int ctype_offset, 
    int etype_offset, int shiftflag, double *shift)

{
  int m, xptr, iptr;
  imageint imagedata;
  double xdata[3], lamda[3];
  double *coord;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_element_V2 && nwords != evec->size_data_element_V2 + 3)
    error->all(FLERR, "Incorrect element format in data file");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  char **values = new char*[nwords];

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3], subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != Comm::LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }
  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // xptr = which word in line starts xyz coords of element
  // iptr = which word in line starts ix, iy, iz image flags

  xptr = evec->xcol_data_V2 - 1;
  int imageflag = 0;
  if (nwords > evec->size_data_element_V2) imageflag = 1;
  if (imageflag) iptr = nwords - 3; 

  // loop over lines of element data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap element into simulation box

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect element format in data file");
    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL, " \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR, "Incorrect element format in data file");
    }

    if (imageflag)
      imagedata = ((imageint) (atoi(values[iptr]) + IMGMAX) & IMGMASK) |
        (((imageint) (atoi(values[iptr+1]) + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (atoi(values[iptr+2]) + IMGMAX) & IMGMASK) << IMG2BITS);
    else imagedata = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    if (shiftflag) {
      xdata[0] += shift[0];
      xdata[1] += shift[1];
      xdata[2] += shift[2];
    }

    domain->remap(xdata, imagedata);
    if (triclinic) {
      domain->x2lamda(xdata, lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      evec->data_element_V2(xdata, imagedata, values);

      if (id_offset) tag[nlocal-1] += id_offset;
      if (ctype_offset) {
        ctype[nlocal-1][0] += ctype_offset;
        if (ctype[nlocal-1][0] > atom->ntypes)
          error->one(FLERR, "Invalid atom type in Elements section of data file");
      }
      if (etype_offset) {
        etype[nlocal-1] += etype_offset;
        if (etype[nlocal-1] > netypes)
          error->one(FLERR, "Invalid element type in Elements section of data file");
      }
    }
    buf = next+1; 
  }
  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack n lines from Node Velocity section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_vels_V2(int n, char *buf, tagint id_offset)
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_vel_V2)
    error->all(FLERR, "Incorrect node velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of node velocities
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL, " \t\n\r\f");

    tagdata = ATOTAGINT(values[0]) + id_offset;
    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Node Velocities section of data file");
    if ((m = map(tagdata)) >= 0) evec->data_vel_V2(m, &values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_nodes_V2(int n, char *buf, tagint id_offset, int shiftflag, double *shift) 
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node_V2)
    error->all(FLERR, "Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect node format in data file");
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL, " \t\n\r\f");
      if (values[j] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
    }

    tagdata = ATOTAGINT(values[0]) + id_offset;

    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Nodes section of data file");
    if ((m = map(tagdata)) >= 0)
      evec->data_node_V2(m, &values[1], shiftflag, shift);

    buf = next + 1; 
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_nodes_reference_V2(int n, char *buf, tagint id_offset, int shiftflag, double *shift) 
{
  int j, m;
  tagint tagdata;
  char *next;

  next = strchr(buf, '\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node_V2 + 3)
    error->all(FLERR, "Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');

    values[0] = strtok(buf, " \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR, "Incorrect node format in data file");
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL, " \t\n\r\f");
      if (values[j] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
    }

    tagdata = ATOTAGINT(values[0]) + id_offset;

    if (tagdata <= id_offset || tagdata > map_tag_max)
      error->one(FLERR, "Invalid element ID in Nodes section of data file");
    if ((m = map(tagdata)) >= 0)
      evec->data_node_reference_V2(m, &values[1], shiftflag, shift);

    buf = next + 1; 
  }

  delete [] values;
}

