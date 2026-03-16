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
#include "variable.h"
#include "compute.h"
#include "memory.h"
#include "modify.h"
#include "fix.h"
#include "error.h"
#include "universe.h"
#include "math_const.h"

using namespace CAC_NS;
using namespace MathConst;

#define DELTA 1
#define DELTA_MEMSTR 1024
#define EPSILON 1.0e-6
#define EXTRA 1000


enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};

Element::Element(CAC *cac) : Pointers(cac)
{
  apc = 1;
  npe = 8;

  // default parameter for each element

  maxintpl = 4000;  
  maxintg = 300;
  maxsubintpl = 500;
  maxsubelem = 125;
  maxelemchange = 2.0;

  esplit = 2; 

  nmax = 0;
  nmaxintg = 0;
  nmaxintpl = 0;
  max_nintg = max_nintpl = 0;
  nelements = nnodes = nintpls = 0;
  nintpls = 0;
  nctypes = 0;
  netypes = 0; 
  nlocal = 0;
  nghost = 0;

  // per-element arrays

  nodex = NULL;
  nodev = NULL;
  nodef = NULL;
  nodetag = NULL;
  x = NULL;
  image = NULL;
  slip_plane = NULL;
  tag = NULL;
  ctype = NULL;
  etype = NULL;
  mask = NULL;
  nodemask = NULL;
  elem_size = NULL;
  initial_size = NULL;

  //netype-length arrays

  nintpl = NULL;
  nintg = NULL;
  ias2ia = NULL;
  i2ia = NULL; 
  ia2i = NULL;
  i2n = NULL;
  n2i = NULL;
  natom_subelem = NULL;
  weight = NULL;
  shape_array = NULL; 
  weighted_shape_array = NULL;
  shape_array_center_subelem = NULL;

  max_size = 0.0;
  max_diag_size = 0.0;
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
  intg_input_style = 0;

  // gauss integration points and weights
  // modify for more integration points

  memory->create(xgauss,6,6,"element:xgauss");
  memory->create(wgauss,6,6,"element:wgauss");

  xgauss[1][1] = 0.0;

  xgauss[2][1] = -sqrt(1.0/3.0);
  xgauss[2][2] =  sqrt(1.0/3.0);

  xgauss[3][1] = -sqrt(3.0/5.0);
  xgauss[3][2] = 0.0;
  xgauss[3][3] =  sqrt(3.0/5.0);

  xgauss[4][1] = -sqrt(3.0/7.0+2.0/7.0*sqrt(1.2));
  xgauss[4][1] = -sqrt(3.0/7.0-2.0/7.0*sqrt(1.2));
  xgauss[4][1] =  sqrt(3.0/7.0-2.0/7.0*sqrt(1.2));
  xgauss[4][1] =  sqrt(3.0/7.0+2.0/7.0*sqrt(1.2));

  xgauss[5][1] = -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;
  xgauss[5][2] = -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
  xgauss[5][3] = 0;
  xgauss[5][4] =  sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0;
  xgauss[5][5] =  sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0;

  wgauss[1][1] = 2.0;

  wgauss[2][1] = 1.0;
  wgauss[2][2] = 1.0;

  wgauss[3][1] = 5.0/9.0;
  wgauss[3][2] = 8.0/9.0;
  wgauss[3][3] = 5.0/9.0;

  wgauss[4][1] = (18.0-sqrt(30.0))/36.0;
  wgauss[4][2] = (18.0+sqrt(30.0))/36.0;
  wgauss[4][3] = (18.0+sqrt(30.0))/36.0;
  wgauss[4][4] = (18.0-sqrt(30.0))/36.0;

  wgauss[5][1] = (322.0-13.0*sqrt(70.0))/900.0;
  wgauss[5][2] = (322.0+13.0*sqrt(70.0))/900.0;
  wgauss[5][3] = 128.0/225.0;
  wgauss[5][4] = (322.0+13.0*sqrt(70.0))/900.0;
  wgauss[5][5] = (322.0-13.0*sqrt(70.0))/900.0;

  // callback lists & extra restart info
 
  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = NULL;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = NULL;

  element_style = NULL;
  evec = NULL;

  evec_map = new ElementVecCreatorMap();

#define ELEMENT_CLASS
#define ElementStyle(key,Class) \
  (*evec_map)[#key] = &evec_creator<Class>;
#include "style_element.h"
#undef ElementStyle
#undef ELEMENT_CLASS

}

/* ---------------------------------------------------------------------- */

Element::~Element()
{
  memory->destroy(xgauss);
  memory->destroy(wgauss);

  delete [] element_style;
  delete evec;

  // delete per-element arrays

  memory->destroy(nodex);
  memory->destroy(nodev);
  memory->destroy(nodef);
  memory->destroy(nodetag);
  memory->destroy(x);
  memory->destroy(image);
  memory->destroy(slip_plane);
  memory->destroy(tag);
  memory->destroy(mask);
  memory->destroy(nodemask);
  memory->destroy(ctype);
  memory->destroy(etype);
  memory->destroy(elem_size);
  memory->destroy(initial_size);

  // delete per-etype arrays

  memory->destroy(shape_array);
  memory->destroy(weighted_shape_array);
  memory->destroy(nintpl);
  memory->destroy(nintg);
  memory->destroy(weight);
  memory->destroy(i2ia);
  memory->destroy(ia2i);
  memory->destroy(i2n);
  memory->destroy(n2i);
  memory->destroy(ias2ia);
  memory->destroy(natom_subelem);
  memory->destroy(shape_array_center_subelem);

  // delete mapping data structures

  map_delete();
}

/* ----------------------------------------------------------------------
   create an ElementVec style
   called from cac.cpp, input script
   ------------------------------------------------------------------------- */

void Element::create_evec(const char *style, int narg, char **arg)
{

  delete [] element_style;
  if (evec) delete evec;
  element_style = NULL;
  evec = NULL;

  // create instance of ElementVec
  // use grow() to initialize element-based arrays to length 1
  // so that x[0][0] and nodex[0][0][0] can always be referenced even if proc has no atoms

  evec = new_evec(style);
  evec->store_args(narg,arg);
  evec->process_args(narg,arg); 
  evec->grow(1); 
}

/* ----------------------------------------------------------------------
   generate an ElementVec class
   ------------------------------------------------------------------------- */

ElementVec *Element::new_evec(const char *style)
{
  if (evec_map->find(style) != evec_map->end()) {
    ElementVecCreator evec_creator = (*evec_map)[style];
    return evec_creator(cac);
  } else error->all(FLERR,"Unknown element style");
  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per ElementVec style in style_element.h
   ------------------------------------------------------------------------- */

  template <typename T>
ElementVec *Element::evec_creator(CAC *cac)
{
  return new T(cac);
}

/* ----------------------------------------------------------------------
   register a callback to a fix so it can manage elem-based arrays
   happens when fix is created
   flag = 0 for grow, 1 for restart, 2 for border comm
------------------------------------------------------------------------- */

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
      memory->grow(extra_grow,nextra_grow_max,"atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
  } else if (flag == 1) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
      memory->grow(extra_restart,nextra_restart_max,"atom:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
  } else if (flag == 2) {
    if (nextra_border == nextra_border_max) {
      nextra_border_max += DELTA;
      memory->grow(extra_border,nextra_border_max,"atom:extra_border");
    }
    extra_border[nextra_border] = ifix;
    nextra_border++;
  }
}

/* ----------------------------------------------------------------------
   unregister a callback to a fix
   happens when fix is deleted, called by its destructor
   flag = 0 for grow, 1 for restart, 2 for border
------------------------------------------------------------------------- */

void Element::delete_callback(const char *id, int flag)
{
  if (id == NULL) return;

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(id,modify->fix[ifix]->id) == 0) break;

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

/* ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
------------------------------------------------------------------------- */

void Element::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix) extra_border[i]--;
}

/* ---------------------------------------------------------------------- */

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
}

/* ----------------------------------------------------------------------
   unpack n lines from Element section of data file
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_elements(int n, char *buf)
{
  int m,xptr,iptr;
  imageint imagedata;
  double xdata[3],lamda[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
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

  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
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


  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    values[0] = strtok(buf," \t\n\r\f");

    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL," \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR,"Incorrect element format in data file");
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

    domain->remap(xdata,imagedata);
    if (triclinic) {
      domain->x2lamda(xdata,lamda);
      coord = lamda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) 
      evec->data_element(xdata,imagedata,values);

    buf = next+1; 
  }
  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack n lines from Node Velocity section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_vels(int n, char *buf)
{
  int j,m;
  tagint tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_vel)
    error->all(FLERR,"Incorrect node velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of node velocities
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]);
    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid element ID in Node Velocities section of data file");
    if ((m = map(tagdata)) >= 0) evec->data_vel(m,&values[1]);

    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Node section of data file
   check that element IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
   ------------------------------------------------------------------------- */

void Element::data_nodes(int n, char *buf) {
  int j,m;
  tagint tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != evec->size_data_node)
    error->all(FLERR,"Incorrect node format in data file");

  char **values = new char*[nwords];

  // loop over lines of node 
  // tokenize the line into values
  // if I own element tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]);

    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid element ID in Node Velocities section of data file");
    if ((m = map(tagdata)) >= 0) evec->data_node(m,&values[1]);
    buf = next + 1; 
  }

  delete [] values;
}

/*--------------------------------------------------------------------------------------------------------------
  return # of bytes of allocated memory
  ----------------------------------------------------------------------------------------------------------------*/

bigint Element::memory_usage()
{
  // per-element arrays

  memlength = DELTA_MEMSTR;
  memory->create(memstr,memlength,"element:memstr");
  memstr[0] = '\0';
  bigint bytes = evec->memory_usage();

  memory->destroy(memstr);

  // element map

  bytes += max_same*sizeof(int);
  if (map_style == 1)
    bytes += memory->usage(map_array,map_maxarray);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }

  // per-type arrays

  int n = netypes + 1;
  bytes += memory->usage(weight,n,maxintg);
  bytes += memory->usage(i2ia,n,maxintg);
  bytes += memory->usage(ia2i,n,maxintpl);
  bytes += memory->usage(i2n,n,maxintg);
  bytes += memory->usage(n2i,n,npe);
  bytes += memory->usage(shape_array,n,maxintpl,npe);
  bytes += memory->usage(weighted_shape_array,n,maxintg,npe);
  bytes += memory->usage(ias2ia,n,nsubelem,maxsubintpl);

  return bytes;
}

/* ----------------------------------------------------------------------
   accumulate per-element vec names in memstr, padded by spaces
   return 1 if padded str is not already in memlist, else 0
   ------------------------------------------------------------------------- */

int Element::memcheck(const char *str)
{
  int n = strlen(str) + 3;
  char *padded = new char[n];
  strcpy(padded," ");
  strcat(padded,str);
  strcat(padded," ");

  if (strstr(memstr,padded)) {
    delete [] padded;
    return 0;
  }

  if (strlen(memstr) + n >= memlength) {
    memlength += DELTA_MEMSTR;
    memory->grow(memstr,memlength,"element:memstr");
  }

  strcat(memstr,padded);
  delete [] padded;
  return 1;
}


/* -----------------------------------
   Modify parameters
   Called from element_modify command
   ------------------------------------*/

void Element::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"esplit") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      esplit = universe->inumeric(FLERR,arg[iarg+1]);
      evec->setup_sub_element();
      if (esplit <= 0) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxintpl") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      maxintpl = universe->inumeric(FLERR,arg[iarg+1]);
      if (maxintpl <= 0) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxintg") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      maxintg = universe->inumeric(FLERR,arg[iarg+1]);
      if (maxintg <= 0) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxsubintpl") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      maxsubintpl = universe->inumeric(FLERR,arg[iarg+1]);
      if (maxsubintpl <= 0) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxelemchange") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      maxelemchange = universe->inumeric(FLERR,arg[iarg+1]);
      if (maxelemchange <= 0) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"apc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal element_modify command");
      if (domain->box_exist) error->all(FLERR,"Cannot change apc after box exist");
      apc = universe->inumeric(FLERR,arg[iarg+1]);
      if (apc < 1) error->all(FLERR,"Illegal element_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal element_modify command");
  } 
}

/* -----------------------------------
   Modify element types
   Called from modify_elements command
   ------------------------------------*/

void Element::modify_elements(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal modify_elements command");

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find modify_elements group ID");
  int groupbit = group->bitmask[igroup];

  // check which type to modify

  int *attribute;
  int value = universe->inumeric(FLERR,arg[2]);
  if (strcmp(arg[1],"etype") == 0) {
    attribute = etype;
    if (value < 0 || value > netypes) 
      error->all(FLERR,"Invalid element type in modify_elements command");
  } else if (strcmp(arg[1],"ctype") == 0) { 
    attribute = ctype;
    if (value < 0 || value > nctypes) 
      error->all(FLERR,"Invalid element atomic type in modify_elements command");
  } else error->all(FLERR,"Illegal modify_elements command");

  // modify elements

  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) 
      attribute[i] = value;

  count_intpl();

}

/* ----------------------------------------------------------------------
   add element type, call from input script as add_etype command 
   ------------------------------------------------------------------------- */

void Element::add_etype(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Illegal add_etype command");
  int itype = universe->inumeric(FLERR,arg[0]);
  if (itype <= netypes) error->all(FLERR,"Illegal element type in add_etype command");
  char *sintg,*sintpl;
  sintg = new char[100];
  sintpl = new char[100];
  sprintf(sintg,"%d ",itype);
  sprintf(sintpl,"%d ",itype);

  int iarg = 1;

  // set flag
  // -1 not set yet
  //  0 being set
  //  1 already set

  int intpl_flag = -1;
  int intg_flag = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"interpolate") == 0) {
      if (intpl_flag != -1) error->all(FLERR,"Illegal add_etype command");
      if (intg_flag == 0) intg_flag = 1;
      intpl_flag = 0;
      iarg++;
    } else if(strcmp(arg[iarg],"integration") == 0) {
      if (intg_flag != -1) error->all(FLERR,"Illegal add_etype command");
      if (intpl_flag == 0) intpl_flag = 1;
      intg_flag = 0;
      iarg++;
    } else {
      if (intg_flag == 0) sprintf(sintg+strlen(sintg),"%s ",arg[iarg]);
      else if (intpl_flag == 0) sprintf(sintpl+strlen(sintpl),"%s ",arg[iarg]);
      else error->all(FLERR,"Illegal add_etype command");
      iarg++;
    }
  }

  if (intg_flag == -1 || intpl_flag == -1) error->all(FLERR,"Illegal add_etype command");

  // init to make sure subelements are setup
  init();  
  evec->grow_etype_arrays(itype);
  evec->set_interpolate(sintpl,0);
  evec->set_integration(sintg,0);
  delete [] sintg;
  delete [] sintpl;

}

/* ----------------------------------------------------------------------
   reset node tags for node connectivity in dump cac command
   called during dump cac
   ------------------------------------------------------------------------- */

void Element::reset_node_tags(int groupbit)
{

  // maxtag = # of nodes I own in group
  // maxtag_sum = # of total nodes on procs <= me in group

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) maxtag += npe; 
  tagint maxtag_sum;

  MPI_Scan(&maxtag,&maxtag_sum,1,MPI_CAC_TAGINT,MPI_SUM,world);

  // itag = 1st tag that my nodes in group should use

  tagint itag = maxtag_sum - maxtag + 1;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int j = 0; j < npe; j++)
        nodetag[i][j] = itag++;
}

/* ----------------------------------------------------------------------
   add unique tags to any elements with tag = 0
   new tags are grouped by proc and start after max current tag
   called after creating new elements
   error if new tags will exceed MAXTAGINT
   ------------------------------------------------------------------------- */

void Element::tag_extend()
{
  // maxtag_all = max tag for all elements

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag,tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);

  // notag = # of elements I own with no tag (tag = 0)
  // notag_sum = # of total elements on procs <= me with no tag

  bigint notag = 0;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) notag++;

  bigint notag_total;
  MPI_Allreduce(&notag,&notag_total,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (notag_total >= MAXTAGINT)
    error->all(FLERR,"New element IDs exceed maximum allowed ID");

  bigint notag_sum;
  MPI_Scan(&notag,&notag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my untagged elements should use

  tagint itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) tag[i] = itag++;

}

/* ----------------------------------------------------------------------
   check if node coords match up with center coords
   called during read_data
   ------------------------------------------------------------------------- */

void Element::check_node_coords() 
{
  double xtmp,ytmp,ztmp;
  for (int i = 0; i < nlocal; i++) {
    xtmp = ytmp = ztmp = 0.0;
    for (int j = 0; j < npe; j++) {
      xtmp += nodex[i][j][0];
      ytmp += nodex[i][j][1];
      ztmp += nodex[i][j][2];
    }
    xtmp /= npe;
    ytmp /= npe;
    ztmp /= npe;
    if (fabs(xtmp-x[i][0]) > EPSILON ||
        fabs(ytmp-x[i][1]) > EPSILON ||
        fabs(ztmp-x[i][2]) > EPSILON)
      error->one(FLERR,"Node coords does not match with element center coords");
  }
}

/* ----------------------------------------------------------------------
   check that element IDs are valid
   error if any element ID < 0 or element ID = MAXTAGINT
   if any element ID > 0, error if any element ID == 0
   if any element ID > 0, error if tag_enable = 0
   if all element IDs = 0, tag_enable must be 0
   if max element IDs < nelements, must be duplicates
   OK if max element IDs > nelements
   NOTE: not fully checking that element IDs are unique
------------------------------------------------------------------------- */

void Element::tag_check()
{
  tagint min = MAXTAGINT;
  tagint max = 0;

  for (int i = 0; i < nlocal; i++) {
    min = MIN(min,tag[i]);
    max = MAX(max,tag[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_CAC_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_CAC_TAGINT,MPI_MAX,world);

  if (minall < 0) error->all(FLERR,"One or more element IDs is negative");
  if (maxall >= MAXTAGINT) error->all(FLERR,"One or more element IDs is too big");
  if (maxall > 0 && minall == 0)
    error->all(FLERR,"One or more element IDs is zero");
  if (maxall > 0 && tag_enable == 0)
    error->all(FLERR,"Non-zero element IDs with element_modify id = no");
  if (maxall == 0 && nelements && tag_enable)
    error->all(FLERR,"All element IDs = 0 but element_modify id = yes");
  if (tag_enable && maxall < nelements)
    error->all(FLERR,"Duplicate element IDs exist");
}

/* ----------------------------------------------------------------------
   count number of interpolated atoms
   ------------------------------------------------------------------------- */

void Element::count_intpl() 
{
  bigint n = 0;
  for (int i = 0; i < nlocal; i++)
    n += nintpl[etype[i]];
  MPI_Allreduce(&n,&nintpls,1,MPI_CAC_BIGINT,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   init per-element fix/compute/variable values for newly created elements
   called from create_elements, read_data, read_dump,
   fixes, computes, variables may or may not exist when called
------------------------------------------------------------------------- */

void Element::data_fix_compute_variable(int nprev, int nnew)
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

  for (int i = nprev; i < nnew; i++)
    input->variable->set_elem_arrays(i);
}


