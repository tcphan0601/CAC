#include <stdlib.h>
#include <string.h>
#include "element_vec_cac.h"
#include "element.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "group.h"
#include "error.h"
#include "math_extra.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathExtra;

// same as in element.cpp

enum{QUADRILATERAL,TRIANGLE,HEXAHEDRON,TETRAHEDRON,OCTAHEDRON,PYRAMID,WEDGE};

#define EPSILON 1e-6
#define BIG     1e30
#define MAXSUB  5

/*----------------------------------------------------------------*/

ElementVecCAC::ElementVecCAC(CAC *cac) : ElementVec(cac)
{

  ncells = nintgs = NULL;
  intg = NULL;

  // size of data line in write_data_elem command

  size_data_elem = 5;
  size_data_elem_vel = 5;

  // size of data line in write_data command

  size_data_element = 6;
  size_data_node = 6;
  size_data_vel = 5;

  // position of first x column in read_data

  xcol_data = 4;

  // size of data line in write_tecplot command

  size_tecplot_node = 3;

  nrequested_etype = 0;
  maxrequested_etype = 4;
  memory->create(requested_etype,4*8,"evec:requested_etype");

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
  node_set_3D[1][1][3] = 8;

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

}

/*----------------------------------------------------------------*/

ElementVecCAC::~ElementVecCAC()
{
  memory->destroy(ncells);
  memory->destroy(nintgs);
  memory->destroy(intg);
  memory->destroy(requested_etype);
}

/* ----------------------------------------------------------------------
   set communication sizes that are dependent on max_npw
   ------------------------------------------------------------------------- */

void ElementVecCAC::set_data_size() 
{
  int max_npe = element->max_npe;
  int max_apc = element->max_apc;

  size_data_cluster = max_apc;

  // size of data line in write_tecplot command

  if (domain->dimension == 3) 
    size_tecplot_node_connect = 8;
  else
    size_tecplot_node_connect = 4;

  // size of communication

  size_velocity = 3*max_npe;
  size_forward = 10 + 3*max_npe;
  size_reverse = 3*max_npe;
  size_border = 14 + 5*max_npe;
  if (max_apc > 1) size_border += max_apc;
}

/* ----------------------------------------------------------------------
   copy element I info to element J
   ------------------------------------------------------------------------- */

void ElementVecCAC::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  ctype[j] = ctype[i];
  etype[j] = etype[i];
  mask[j] = mask[i];
  image[j] = image[i];
  initial_box_size[j][0] = initial_box_size[i][0];
  initial_box_size[j][1] = initial_box_size[i][1];
  initial_box_size[j][2] = initial_box_size[i][2];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  for (int k = 0; k < element->npe[etype[i]]; k++) {
    nodemask[j][k] = nodemask[i][k];
    nodetag[j][k] = nodetag[i][k];
    nodex[j][k][0] = nodex[i][k][0];
    nodex[j][k][1] = nodex[i][k][1];
    nodex[j][k][2] = nodex[i][k][2];
    nodev[j][k][0] = nodev[i][k][0];
    nodev[j][k][1] = nodev[i][k][1];
    nodev[j][k][2] = nodev[i][k][2];
  }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      modify->fix[element->extra_grow[iextra]]->copy_elem_arrays(i,j,delflag);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
   ------------------------------------------------------------------------- */

void ElementVecCAC::grow(int n)
{
  int max_npe = element->max_npe;
  if (n == 0) grow_nmax();
  else {
    nmax = n;
    element->nmaxintg = nmax * element->max_nintg;
    element->nmaxintpl = nmax * element->max_nintpl;
  }
  element->nmax = nmax;

  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(element->tag,nmax,"element:tag");
  ctype = memory->grow(element->ctype,nmax,"element:ctype");
  etype = memory->grow(element->etype,nmax,"element:etype");
  mask = memory->grow(element->mask,nmax,"element:mask");
  nodemask = memory->grow(element->nodemask,nmax,max_npe,"element:nodemask");
  image = memory->grow(element->image,nmax,"element:image");
  x = memory->grow(element->x,nmax,3,"element:x");
  slip_plane = memory->grow(element->slip_plane,nmax,3,3,"element:slip_plane");
  cell_size = memory->grow(element->cell_size,nmax,3,"element:cell_size");
  nodex = memory->grow(element->nodex,nmax,max_npe,3,"element:nodex");
  nodev = memory->grow(element->nodev,nmax,max_npe,3,"element:nodev");
  nodef = memory->grow(element->nodef,nmax,max_npe,3,"element:nodef");
  nodetag = memory->grow(element->nodetag,nmax,max_npe,"element:nodetag");
  subelem_size = memory->grow(element->subelem_size,nmax,"element:subelem_size");
  initial_box_size = memory->grow(element->initial_box_size,nmax,3,"element:initial_box_size");
  element_bound_box = memory->grow(element->element_bound_box,nmax,6,"element:element_bound_box");
  if (element->max_apc > 1) element_clusters = memory->grow(element->element_clusters,nmax,element->max_apc,"element:element_clusters");

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      modify->fix[element->extra_grow[iextra]]->grow_elem_arrays(nmax);

}

/* ---------------------------------------------------------
   grow etype arrays to n types
 * ------------------------------------------------------*/

void ElementVecCAC::grow_etype_arrays(int n)
{
  int max_npe = element->max_npe;
  int old_netypes = element->netypes;
  if (n <= old_netypes) return;
  else element->netypes = n;

  // grow style specific arrays

  memory->grow(element_type_setflag,n+1,"evec:element_type_setflag");

  // grow integration point arrays

  memory->grow(nintgs,n+1,3,"evec:nintgs");

  // interpolated atom arrays

  memory->grow(ncells,n+1,3,"element:ncells");

  // grow arrays common to all element styles

  npe = memory->grow(element->npe,n+1,"element:npe");
  apc = memory->grow(element->apc,n+1,"element:apc");
  element_shape_ids = memory->grow(element->element_shape_ids,n+1,"element:element_shape_ids");
  element_shape_names = memory->grow(element->element_shape_names,n+1,256,"element:element_shape_names");

  // grow integration point arrays

  nintg = memory->grow(element->nintg,n+1,"element:nintg");
  weighted_shape_array = memory->grow(element->weighted_shape_array,
      n+1,element->maxintg,max_npe,"element:weighted_shape_array");
  i2ia = memory->grow(element->i2ia,n+1,element->maxintg,"element:i2ia");
  weight = memory->grow(element->weight,n+1,element->maxintg,"element:weight");
  i2n = memory->grow(element->i2n,n+1,element->maxintg,"element:i2n");

  // grow node arrays

  n2i = memory->grow(element->n2i,n+1,max_npe,"element:n2i");
  n2ia = memory->grow(element->n2ia,n+1,max_npe,"element:n2ia");

  // interpolated atom arrays

  nintpl = memory->grow(element->nintpl,n+1,"element:nintpl");
  surface_intpl = memory->grow(element->surface_intpl,n+1,element->max_surface,element->max_surface_intpl,"element:surface_intpl");
  nsurface_intpl = memory->grow(element->nsurface_intpl,n+1,element->max_surface,"element:nsurface_intpl");
  edge_intpl = memory->grow(element->edge_intpl,n+1,element->max_edge,element->max_edge_intpl,"element:edge_intpl");
  nedge_intpl = memory->grow(element->nedge_intpl,n+1,element->max_edge,"element:nedge_intpl");

  shape_array = memory->grow(element->shape_array,n+1,element->maxintpl,max_npe,"element:shape_array");
  ia2i = memory->grow(element->ia2i,n+1,element->maxintpl,"element:ia2i");

  // sub-element array

  natom_subelem = memory->grow(element->natom_subelem,n+1,element->maxsubelem,"element:natom_subelem");
  ias2ia = memory->grow(element->ias2ia,n+1,element->maxsubelem,element->maxsubintpl,"element:ias2ia");
  shape_array_center_subelem = memory->grow(element->shape_array_center_subelem,
      n+1,element->maxsubelem,max_npe,"element:shape_array_center_subelem");
  nsubelem = memory->grow(element->nsubelem,n+1,"element:nsubelem");
  subsplit = memory->grow(element->subsplit,n+1,"element:subsplit");

  for (int itype = old_netypes + 1; itype <= n; itype++) {
    element_type_setflag[itype] = 0;
    nintg[itype] = 0; 
    nintpl[itype] = 0; 
    for (int i = 0; i < element->maxintpl; i++) 
      ia2i[itype][i] = -1;
    for (int i = 0; i < element->maxintg; i++) {
      i2ia[itype][i] = -1;
      i2n[itype][i] = -1;
    }
    for (int i = 0; i < max_npe; i++)
      n2i[itype][i] = -1;
  }
}

/* ----------------------------------------------------------------------
   unpack one line from Elements section of data file
   initialize other element quantities
   ------------------------------------------------------------------------- */

void ElementVecCAC::data_element(double *coord, imageint imagetmp, char **values)
{
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);  
  while ((nlocal+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
  while ((nlocal+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();

  tag[nlocal] = ATOTAGINT(values[0]);
  etype[nlocal] = atoi(values[1]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid element ID in Elements section of data file");
  if (etype[nlocal] <= 0 || etype[nlocal] > element->netypes)
    error->one(FLERR,"Invalid element type in Elements section of data file"); 

  ctype[nlocal] = atoi(values[2]);
  if (ctype[nlocal] <= 0 || ctype[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Elements section of data file"); 
  image[nlocal] = imagetmp;

  // 1 equal to the all group, 2 equal to the element group

  mask[nlocal] = 1|2;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  initial_box_size[nlocal][0] = -1.0;
  initial_box_size[nlocal][1] = -1.0;
  initial_box_size[nlocal][2] = -1.0;
  element->nintpls += nintpl[etype[nlocal]];
  element->nlocal++;

}

/* ----------------------------------------------------------------------
   unpack one line from Node Velocities section of data file
   ------------------------------------------------------------------------- */

void ElementVecCAC::data_vel(int m, char **values)
{

  int i = atoi(values[0])-1; 

  if (i >= npe[etype[m]]) error->all(FLERR,"Invalid node index value in Node section of data file");
  nodev[m][i][0] = atof(values[1]);
  nodev[m][i][1] = atof(values[2]);
  nodev[m][i][2] = atof(values[3]);
}

/* ----------------------------------------------------------------------
   unpack one line from Node of data file
   ------------------------------------------------------------------------- */

void ElementVecCAC::data_node(int m, char **values, int shiftflag, double *shift)
{
  int i = atoi(values[0])-1; 
  if (i >= npe[etype[m]]) error->all(FLERR,"Invalid node index value in Node section of data file");
  nodex[m][i][0] = atof(values[1]);
  nodex[m][i][1] = atof(values[2]);
  nodex[m][i][2] = atof(values[3]);
  if (shiftflag) {
    nodex[m][i][0] += shift[0];
    nodex[m][i][1] += shift[1];
    nodex[m][i][2] += shift[2];
  }
  nodetag[m][i] = ATOTAGINT(values[4]);

  // 1 equal to the all group, 2 equal to the element group

  nodemask[m][i] = 1|2;
  nodev[m][i][0] = 0.0;
  nodev[m][i][1] = 0.0;
  nodev[m][i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   pack data for element I for sending to another proc
   ------------------------------------------------------------------------- */

int ElementVecCAC::pack_exchange(int i, double *buf)
{
  int m = 1;
  int inpe = npe[etype[i]];
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = initial_box_size[i][0];
  buf[m++] = initial_box_size[i][1];
  buf[m++] = initial_box_size[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(ctype[i]).d;
  buf[m++] = ubuf(etype[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = element_bound_box[i][0];
  buf[m++] = element_bound_box[i][1];
  buf[m++] = element_bound_box[i][2];
  buf[m++] = element_bound_box[i][3];
  buf[m++] = element_bound_box[i][4];
  buf[m++] = element_bound_box[i][5];
  buf[m++] = subelem_size[i];
  for (int j = 0; j < inpe; j++) {
    buf[m++] = ubuf(nodemask[i][j]).d;
    buf[m++] = ubuf(nodetag[i][j]).d;
    buf[m++] = nodex[i][j][0];
    buf[m++] = nodex[i][j][1];
    buf[m++] = nodex[i][j][2];
    buf[m++] = nodev[i][j][0];
    buf[m++] = nodev[i][j][1];
    buf[m++] = nodev[i][j][2];
  }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      m += modify->fix[element->extra_grow[iextra]]->pack_elem_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------
   unpack element from exchanged proc
   -------------------------------------------------------------- */

int ElementVecCAC::unpack_exchange(double *buf)
{

  // grow per-element arrays if needed

  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);
  while ((nlocal+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
  while ((nlocal+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();

  int m = 1;

  // unpack

  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  initial_box_size[nlocal][0] = buf[m++];
  initial_box_size[nlocal][1] = buf[m++];
  initial_box_size[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  ctype[nlocal] = (int) ubuf(buf[m++]).i;
  etype[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  element_bound_box[nlocal][0] = buf[m++];
  element_bound_box[nlocal][1] = buf[m++];
  element_bound_box[nlocal][2] = buf[m++];
  element_bound_box[nlocal][3] = buf[m++];
  element_bound_box[nlocal][4] = buf[m++];
  element_bound_box[nlocal][5] = buf[m++];
  subelem_size[nlocal] = buf[m++];
  int inpe = npe[etype[nlocal]];
  for (int j = 0; j < inpe; j++) {
    nodemask[nlocal][j] = (int) ubuf(buf[m++]).i;
    nodetag[nlocal][j] = (tagint) ubuf(buf[m++]).i;
    nodex[nlocal][j][0] = buf[m++];
    nodex[nlocal][j][1] = buf[m++];
    nodex[nlocal][j][2] = buf[m++];
    nodev[nlocal][j][0] = buf[m++];
    nodev[nlocal][j][1] = buf[m++];
    nodev[nlocal][j][2] = buf[m++];
  }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      m += modify->fix[element->extra_grow[iextra]]->
        unpack_elem_exchange(nlocal,&buf[m]);

  element->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   pack data for elements in list for border exchange
   ------------------------------------------------------------------------- */

int ElementVecCAC::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k,m,inpe;
  double dx,dy,dz;

  m = 0;

  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(ctype[j]).d;
      buf[m++] = ubuf(etype[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = element_bound_box[j][0];
      buf[m++] = element_bound_box[j][1];
      buf[m++] = element_bound_box[j][2];
      buf[m++] = element_bound_box[j][3];
      buf[m++] = element_bound_box[j][4];
      buf[m++] = element_bound_box[j][5];
      buf[m++] = subelem_size[j];
      inpe = npe[etype[j]];
      for (k = 0; k < inpe; k++) {
        buf[m++] = ubuf(nodemask[j][k]).d;
        buf[m++] = ubuf(nodetag[j][k]).d;
        buf[m++] = nodex[j][k][0];
        buf[m++] = nodex[j][k][1];
        buf[m++] = nodex[j][k][2];
      }
      if (element->element_cluster_flag) 
        for (k = 0; k < element->max_apc; k++)
          buf[m++] = ubuf(element_clusters[j][k]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(ctype[j]).d;
      buf[m++] = ubuf(etype[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      inpe = npe[etype[j]];
      buf[m++] = element_bound_box[j][0] + dx;
      buf[m++] = element_bound_box[j][1] + dy;
      buf[m++] = element_bound_box[j][2] + dz;
      buf[m++] = element_bound_box[j][3] + dx;
      buf[m++] = element_bound_box[j][4] + dy;
      buf[m++] = element_bound_box[j][5] + dz;
      buf[m++] = subelem_size[j];
      for (k = 0; k < inpe; k++) {
        buf[m++] = ubuf(nodemask[j][k]).d;
        buf[m++] = ubuf(nodetag[j][k]).d;
        buf[m++] = nodex[j][k][0] + dx;
        buf[m++] = nodex[j][k][1] + dy;
        buf[m++] = nodex[j][k][2] + dz;
      }
      if (element->element_cluster_flag) 
        for (k = 0; k < element->max_apc; k++)
          buf[m++] = ubuf(element_clusters[j][k]).d;
    }
  }

  if (element->nextra_border)
    for (int iextra = 0; iextra < element->nextra_border; iextra++)
      m += modify->fix[element->extra_border[iextra]]->pack_elem_border(n,list,&buf[m]);

  return m;
}

/* ----------------------------------------------------------------------
   unpack data for elements received from border exchange
   ------------------------------------------------------------------------- */

void ElementVecCAC::unpack_border(int n, int first, double *buf)
{
  int m,last,inpe;
  m = 0;
  last = first + n;
  for (int i = first; i < last; i++) {
    while (i >= nmax) grow(0);
    while ((i+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
    while ((i+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    ctype[i] = (int) ubuf(buf[m++]).i;
    etype[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    element_bound_box[i][0] = buf[m++];
    element_bound_box[i][1] = buf[m++];
    element_bound_box[i][2] = buf[m++];
    element_bound_box[i][3] = buf[m++];
    element_bound_box[i][4] = buf[m++];
    element_bound_box[i][5] = buf[m++];
    subelem_size[i] = buf[m++];
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      nodemask[i][j] = (int) ubuf(buf[m++]).i;
      nodetag[i][j] = (tagint) ubuf(buf[m++]).i;
      nodex[i][j][0] = buf[m++];
      nodex[i][j][1] = buf[m++];
      nodex[i][j][2] = buf[m++];
    }
    if (element->element_cluster_flag) 
      for (int j = 0; j < element->max_apc; j++) 
        element_clusters[i][j] = (tagint) ubuf(buf[m++]).i;

  }

  if (element->nextra_border)
    for (int iextra = 0; iextra < element->nextra_border; iextra++)
      m += modify->fix[element->extra_border[iextra]]->
        unpack_elem_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int ElementVecCAC::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k,inpe;
  double dx,dy,dz;

  int m = 0;

  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      inpe = npe[etype[j]];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = element_bound_box[j][0];
      buf[m++] = element_bound_box[j][1];
      buf[m++] = element_bound_box[j][2];
      buf[m++] = element_bound_box[j][3];
      buf[m++] = element_bound_box[j][4];
      buf[m++] = element_bound_box[j][5];
      buf[m++] = subelem_size[j];
      for (k = 0; k < inpe; k++) {
        buf[m++] = nodex[j][k][0];
        buf[m++] = nodex[j][k][1];
        buf[m++] = nodex[j][k][2];
      }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      inpe = npe[etype[j]];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = element_bound_box[j][0] + dx;
      buf[m++] = element_bound_box[j][1] + dy;
      buf[m++] = element_bound_box[j][2] + dz;
      buf[m++] = element_bound_box[j][3] + dx;
      buf[m++] = element_bound_box[j][4] + dy;
      buf[m++] = element_bound_box[j][5] + dz;
      buf[m++] = subelem_size[j];
      for (k = 0; k < inpe; k++) {
        buf[m++] = nodex[j][k][0] + dx;
        buf[m++] = nodex[j][k][1] + dy;
        buf[m++] = nodex[j][k][2] + dz;
      }
    }
  }
  return m;
}


/* ---------------------------------------------------------------------- */

void ElementVecCAC::unpack_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    element_bound_box[i][0] = buf[m++];
    element_bound_box[i][1] = buf[m++];
    element_bound_box[i][2] = buf[m++];
    element_bound_box[i][3] = buf[m++];
    element_bound_box[i][4] = buf[m++];
    element_bound_box[i][5] = buf[m++];
    subelem_size[i] = buf[m++];
    for (j = 0; j < npe[etype[i]]; j++) {
      nodex[i][j][0] = buf[m++];
      nodex[i][j][1] = buf[m++];
      nodex[i][j][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int ElementVecCAC::pack_reverse(int n, int first, double *buf)
{
  int i,j,m,last,inpe;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
    inpe = npe[etype[i]];
  for (j = 0; j < inpe; j++) {
    buf[m++] = nodef[i][j][0];
    buf[m++] = nodef[i][j][1];
    buf[m++] = nodef[i][j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ElementVecCAC::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,k,m,inpe;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    inpe = npe[etype[j]];
    for (k = 0; k < inpe; k++) {
      nodef[j][k][0] += buf[m++];
      nodef[j][k][1] += buf[m++];
      nodef[j][k][2] += buf[m++];
    }
  }
}


/* ----------------------------------------------------------------------
   add a new element
   ------------------------------------------------------------------------- */

void ElementVecCAC::create_element(double *coord, double **nodecoord, int ietype, int ictype, int itag)
{
  int inpe = npe[ietype];
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);
  while ((nlocal+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
  while ((nlocal+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();

  double bound_box[6]; 
  bound_box[0] = BIG;
  bound_box[1] = BIG;
  bound_box[2] = BIG;
  bound_box[3] = -BIG;
  bound_box[4] = -BIG;
  bound_box[5] = -BIG;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  for (int i = 0; i < inpe; i++) {
    nodex[nlocal][i][0] = nodecoord[i][0];
    nodex[nlocal][i][1] = nodecoord[i][1];
    nodex[nlocal][i][2] = nodecoord[i][2];
    bound_box[0] = MIN(bound_box[0],nodecoord[i][0]);
    bound_box[1] = MIN(bound_box[1],nodecoord[i][1]);
    bound_box[2] = MIN(bound_box[2],nodecoord[i][2]);
    bound_box[3] = MAX(bound_box[3],nodecoord[i][0]);
    bound_box[4] = MAX(bound_box[4],nodecoord[i][1]);
    bound_box[5] = MAX(bound_box[5],nodecoord[i][2]);
    nodev[nlocal][i][0] = 0.0;
    nodev[nlocal][i][1] = 0.0;
    nodev[nlocal][i][2] = 0.0;
    nodemask[nlocal][i] = 1|2;
  }

  etype[nlocal] = ietype;
  ctype[nlocal] = ictype;
  initial_box_size[nlocal][0] = -1.0;
  initial_box_size[nlocal][1] = -1.0;
  initial_box_size[nlocal][2] = -1.0;
  element_bound_box[nlocal][0] = bound_box[0];
  element_bound_box[nlocal][1] = bound_box[1];
  element_bound_box[nlocal][2] = bound_box[2];
  element_bound_box[nlocal][3] = bound_box[3];
  element_bound_box[nlocal][4] = bound_box[4];
  element_bound_box[nlocal][5] = bound_box[5];
  tag[nlocal] = itag;
  mask[nlocal] = 1|2; 
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  element->nlocal++;
}

/* ----------------------------------------------------------------------
   add a new element, node coords are calculated from lattice
   lattice must be defined to call this function
   ------------------------------------------------------------------------- */

void ElementVecCAC::create_element(double *coord, int ietype, int ictype, int itag)
{
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);
  while ((nlocal+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
  while ((nlocal+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  double xtmp = coord[0];
  double ytmp = coord[1];
  double ztmp = coord[2];
  double bound_box[6]; 
  bound_box[0] = BIG;
  bound_box[1] = BIG;
  bound_box[2] = BIG;
  bound_box[3] = -BIG;
  bound_box[4] = -BIG;
  bound_box[5] = -BIG;

  if (domain->lattice)
    domain->lattice->box2lattice(xtmp,ytmp,ztmp);
  else 
    error->one(FLERR,"Lattice must be defined to create new elements");

  double *a1 = domain->lattice->a1;
  double *a2 = domain->lattice->a2;
  double *a3 = domain->lattice->a3;


  double xnode,ynode,znode;

  if (element_shape_ids[ietype] == HEXAHEDRON) {
    int c1,c2,c3;
    c1 = ncells[ietype][0]-1; 
    c2 = ncells[ietype][1]-1; 
    c3 = ncells[ietype][2]-1;

    double d[8][3] = {
      {-1,-1,-1},
      { 1,-1,-1},
      { 1, 1,-1},
      {-1, 1,-1},
      {-1,-1, 1},
      { 1,-1, 1},
      { 1, 1, 1},
      {-1, 1, 1}
    };

    for (int i = 0; i < 8; i++) {
      xnode = xtmp + c1*d[i][0]/2.0;
      ynode = ytmp + c2*d[i][1]/2.0;
      znode = ztmp + c3*d[i][2]/2.0;

      domain->lattice->lattice2box(xnode,ynode,znode);

      nodex[nlocal][i][0] = xnode;
      nodex[nlocal][i][1] = ynode;
      nodex[nlocal][i][2] = znode;
      bound_box[0] = MIN(bound_box[0],xnode);
      bound_box[1] = MIN(bound_box[1],ynode);
      bound_box[2] = MIN(bound_box[2],znode);
      bound_box[3] = MAX(bound_box[3],xnode);
      bound_box[4] = MAX(bound_box[4],ynode);
      bound_box[5] = MAX(bound_box[5],znode);
      nodev[nlocal][i][0] = 0.0;
      nodev[nlocal][i][1] = 0.0;
      nodev[nlocal][i][2] = 0.0;
      nodetag[nlocal][i] = 0;
      nodemask[nlocal][i] = 1|2;
    }
  } else if (element_shape_ids[ietype] == QUADRILATERAL) {
    int c1,c2;
    c1 = ncells[ietype][0]-1; 
    c2 = ncells[ietype][1]-1; 
    double d[4][2] = {
      {-1,-1},
      { 1,-1},
      { 1, 1},
      {-1, 1}
    };

    for (int i = 0; i < 4; i++) {
      xnode = xtmp + c1*d[i][0]/2.0;
      ynode = ytmp + c2*d[i][1]/2.0;
      znode = ztmp;

      domain->lattice->lattice2box(xnode,ynode,znode);

      nodex[nlocal][i][0] = xnode;
      nodex[nlocal][i][1] = ynode;
      nodex[nlocal][i][2] = znode;
      bound_box[0] = MIN(bound_box[0],xnode);
      bound_box[1] = MIN(bound_box[1],ynode);
      bound_box[2] = MIN(bound_box[2],znode);
      bound_box[3] = MAX(bound_box[3],xnode);
      bound_box[4] = MAX(bound_box[4],ynode);
      bound_box[5] = MAX(bound_box[5],znode);
      nodev[nlocal][i][0] = 0.0;
      nodev[nlocal][i][1] = 0.0;
      nodev[nlocal][i][2] = 0.0;
      nodetag[nlocal][i] = 0;
      nodemask[nlocal][i] = 1|2;
    }
  }
  etype[nlocal] = ietype;
  ctype[nlocal] = ictype;
  initial_box_size[nlocal][0] = -1.0;
  initial_box_size[nlocal][1] = -1.0;
  initial_box_size[nlocal][2] = -1.0;
  element_bound_box[nlocal][0] = bound_box[0];
  element_bound_box[nlocal][1] = bound_box[1];
  element_bound_box[nlocal][2] = bound_box[2];
  element_bound_box[nlocal][3] = bound_box[3];
  element_bound_box[nlocal][4] = bound_box[4];
  element_bound_box[nlocal][5] = bound_box[5];
  tag[nlocal] = itag;
  mask[nlocal] = 1|2; 
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  element->nlocal++;
}

/*----------------------------------------------------
  pack element info for Atoms section in data file from write_data_elem command
  ------------------------------------------*/

void ElementVecCAC::pack_element_data_elem(double **buf, tagint tag_offset, int ghostflag)
{
  if (!ghostflag) {
    for (int i = 0; i < element->nlocal; i++) {
      buf[i][0] = ubuf(tag[i] + tag_offset).d;
      buf[i][1] = ubuf(ctype[i]).d;
      buf[i][2] = x[i][0];
      buf[i][3] = x[i][1];
      buf[i][4] = x[i][2];
    }
  } else {
    for (int i = 0; i < element->nlocal + element->nghost; i++) {
      buf[i][0] = ubuf(i + 1 + tag_offset).d;
      buf[i][1] = ubuf(ctype[i]).d;
      buf[i][2] = x[i][0];
      buf[i][3] = x[i][1];
      buf[i][4] = x[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file from write_data_elem command
   ------------------------------------------------------------------------- */

void ElementVecCAC::pack_element_vel_data_elem(double **buf, tagint tag_offset)
{
  int k = 0;
  double vel[3];
  tagint mytag;
  int inpe;
  for (int i = 0; i < element->nlocal; i++) {
    buf[k][0] = ubuf(tag[i] + tag_offset).d;
    buf[k][1] = ubuf(ctype[i]).d;
    vel[0] = vel[1] = vel[2] = 0.0;
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      vel[0] += nodev[i][j][0];
      vel[1] += nodev[i][j][1];
      vel[2] += nodev[i][j][2];
    }
    buf[k][2] = vel[0]/inpe;
    buf[k][3] = vel[1]/inpe;
    buf[k][4] = vel[2]/inpe;
    k++;
  }
}

/* ----------------------------------------------------
   pack node info for Atoms section in data file from write_data_elem command
   ------------------------------------------*/

void ElementVecCAC::pack_node_data_elem(double **buf, tagint tag_offset, int ghostflag)
{
  int k = 0;
  int inpe;
  if (!ghostflag) {
    for (int i = 0; i < element->nlocal; i++) {
      inpe = npe[etype[i]];
      for (int j = 0; j < inpe; j++) {
        buf[k][0] = ubuf(nodetag[i][j]+tag_offset).d;
        buf[k][1] = ubuf(atom->ntypes + 1).d;
        buf[k][2] = nodex[i][j][0];
        buf[k][3] = nodex[i][j][1];
        buf[k][4] = nodex[i][j][2];
        k++;
      }
    }
  } else {
    for (int i = 0; i < element->nlocal + element->nghost; i++) {
      inpe = npe[etype[i]];
      for (int j = 0; j < inpe; j++) {
        buf[k][0] = ubuf(k + 1 + tag_offset).d;
        buf[k][1] = ubuf(atom->ntypes + 1).d;
        buf[k][2] = nodex[i][j][0];
        buf[k][3] = nodex[i][j][1];
        buf[k][4] = nodex[i][j][2];
        k++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack node velocity info for data file from write_data_elem command
   ------------------------------------------------------------------------- */

void ElementVecCAC::pack_node_vel_data_elem(double **buf, tagint tag_offset)
{
  int k = 0;
  int inpe;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      buf[k][0] = ubuf(nodetag[i][j]+tag_offset).d;
      buf[k][1] = ubuf(atom->ntypes+1).d;
      buf[k][2] = nodev[i][j][0];
      buf[k][3] = nodev[i][j][1];
      buf[k][4] = nodev[i][j][2];
      k++;
    }
  }
}

/*----------------------------------------------------
  write node info into Atoms section in data file from write_data_elem command
  ------------------------------------------*/

void ElementVecCAC::write_data_elem(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e \n",
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4]);
}

/* ----------------------------------------------------------------------
   write node velocity info to data file from write_data_elem command
   ------------------------------------------------------------------------- */

void ElementVecCAC::write_vel_data_elem(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n",
        (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4]);
}

/* ----------------------------------------------------
   write element info into CAC Elements section in data file from write_data_adrian command
   ------------------------------------------*/

void ElementVecCAC::write_element_data_adrian(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT " Eight_Node %d %d %d %d\n",
        (tagint) ubuf(buf[i][0]).i,apc[1],
        (int) ubuf(buf[i][2]).i,
        (int) ubuf(buf[i][3]).i,
        (int) ubuf(buf[i][4]).i);

    for (int k = 0; k < apc[1]; k++) {
      for (int j = 0; j < 8; j++) {
        fprintf(fp,"%d %d %d %-1.16e %-1.16e %-1.16e \n",
            j+1,k+1,(int) ubuf(buf[i][1]).i,
            buf[i][5+3*j+24*k],buf[i][6+3*j+24*k],buf[i][7+3*j+24*k]);
      }
    }
  }
}


/* ----------------------------------------------------
   pack element info for CAC Elements section in data file from write_data_adrian command
   ------------------------------------------*/

void ElementVecCAC::pack_element_data_adrian(double **buf)
{
  int ietype;
  double a1[3],a2[3],a3[3],
         a4[3],a5[3],a6[3],
         a7[3],a8[3],a9[3],
         a10[3],a11[3],a12[3];

  int m = -1;
  for (int i = 0; i < element->nlocal; i += element->apc[1] ) {
    m++;
    buf[m][0] = ubuf((tag[i]-1)/apc[1]+1).d;
    buf[m][1] = ubuf(ctype[i]).d;
    ietype = etype[i];
    buf[m][2] = ubuf(ncells[ietype][0]).d;
    buf[m][3] = ubuf(ncells[ietype][1]).d;
    buf[m][4] = ubuf(ncells[ietype][2]).d;

    for (int k = 0; k < apc[1]; k++) {
      a1[0] = nodex[i+k][1][0] - nodex[i+k][0][0]; 
      a1[1] = nodex[i+k][1][1] - nodex[i+k][0][1]; 
      a1[2] = nodex[i+k][1][2] - nodex[i+k][0][2]; 

      a2[0] = nodex[i+k][3][0] - nodex[i+k][0][0]; 
      a2[1] = nodex[i+k][3][1] - nodex[i+k][0][1]; 
      a2[2] = nodex[i+k][3][2] - nodex[i+k][0][2]; 

      a3[0] = nodex[i+k][2][0] - nodex[i+k][3][0]; 
      a3[1] = nodex[i+k][2][1] - nodex[i+k][3][1]; 
      a3[2] = nodex[i+k][2][2] - nodex[i+k][3][2]; 

      a4[0] = nodex[i+k][2][0] - nodex[i+k][1][0]; 
      a4[1] = nodex[i+k][2][1] - nodex[i+k][1][1]; 
      a4[2] = nodex[i+k][2][2] - nodex[i+k][1][2]; 

      a5[0] = nodex[i+k][5][0] - nodex[i+k][1][0]; 
      a5[1] = nodex[i+k][5][1] - nodex[i+k][1][1]; 
      a5[2] = nodex[i+k][5][2] - nodex[i+k][1][2]; 

      a6[0] = nodex[i+k][4][0] - nodex[i+k][0][0]; 
      a6[1] = nodex[i+k][4][1] - nodex[i+k][0][1]; 
      a6[2] = nodex[i+k][4][2] - nodex[i+k][0][2]; 

      a7[0] = nodex[i+k][7][0] - nodex[i+k][3][0]; 
      a7[1] = nodex[i+k][7][1] - nodex[i+k][3][1]; 
      a7[2] = nodex[i+k][7][2] - nodex[i+k][3][2]; 

      a8[0] = nodex[i+k][6][0] - nodex[i+k][2][0]; 
      a8[1] = nodex[i+k][6][1] - nodex[i+k][2][1]; 
      a8[2] = nodex[i+k][6][2] - nodex[i+k][2][2]; 

      a9[0] = nodex[i+k][5][0] - nodex[i+k][4][0]; 
      a9[1] = nodex[i+k][5][1] - nodex[i+k][4][1]; 
      a9[2] = nodex[i+k][5][2] - nodex[i+k][4][2]; 

      a10[0] = nodex[i+k][7][0] - nodex[i+k][4][0]; 
      a10[1] = nodex[i+k][7][1] - nodex[i+k][4][1]; 
      a10[2] = nodex[i+k][7][2] - nodex[i+k][4][2]; 

      a11[0] = nodex[i+k][6][0] - nodex[i+k][7][0]; 
      a11[1] = nodex[i+k][6][1] - nodex[i+k][7][1]; 
      a11[2] = nodex[i+k][6][2] - nodex[i+k][7][2]; 

      a12[0] = nodex[i+k][6][0] - nodex[i+k][5][0]; 
      a12[1] = nodex[i+k][6][1] - nodex[i+k][5][1]; 
      a12[2] = nodex[i+k][6][2] - nodex[i+k][5][2]; 


      int j,ncellx,ncelly,ncellz;
      ncellx = ncells[ietype][0]-1;
      ncelly = ncells[ietype][1]-1;
      ncellz = ncells[ietype][2]-1;
      j = 0;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (- a1[0] - a2[0] - a6[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (- a1[1] - a2[1] - a6[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (- a1[2] - a2[2] - a6[2])/2/ncellz;
      j = 1;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a1[0] - a4[0] - a5[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a1[1] - a4[1] - a5[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a1[2] - a4[2] - a5[2])/2/ncellz;
      j = 2;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a3[0] + a4[0] - a8[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a3[1] + a4[1] - a8[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a3[2] + a4[2] - a8[2])/2/ncellz;
      j = 3;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a2[0] - a3[0] - a7[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a2[1] - a3[1] - a7[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a2[2] - a3[2] - a7[2])/2/ncellz;
      j = 4;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a6[0] - a9[0] - a10[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a6[1] - a9[1] - a10[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a6[2] - a9[2] - a10[2])/2/ncellz;
      j = 5;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a9[0] - a12[0] + a5[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a9[1] - a12[1] + a5[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a9[2] - a12[2] + a5[2])/2/ncellz;
      j = 6;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a12[0] + a11[0] + a8[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a12[1] + a11[1] + a8[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a12[2] + a11[2] + a8[2])/2/ncellz;
      j = 7;
      buf[m][5+j*3+k*3*8] = nodex[i+k][j][0] + (a10[0] - a11[0] + a7[0])/2/ncellx;
      buf[m][6+j*3+k*3*8] = nodex[i+k][j][1] + (a10[1] - a11[1] + a7[1])/2/ncelly;
      buf[m][7+j*3+k*3*8] = nodex[i+k][j][2] + (a10[2] - a11[2] + a7[2])/2/ncellz;
    }
  }
}

/* ----------------------------------------------------
   write element info into Elements section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::write_element_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %d %-1.16e %-1.16e %-1.16e \n",
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, (int) ubuf(buf[i][2]).i,buf[i][3],buf[i][4],buf[i][5]);
}

/* ----------------------------------------------------
   pack element info for Elements section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::pack_element_data(double **buf)
{
  for (int i = 0; i < element->nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(etype[i]).d;
    buf[i][2] = ubuf(ctype[i]).d;
    buf[i][3] = x[i][0];
    buf[i][4] = x[i][1];
    buf[i][5] = x[i][2];
  }
}

/* ----------------------------------------------------
   write element info into Element Clusters  section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::write_element_cluster_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < element->max_apc; j++)
      fprintf(fp,TAGINT_FORMAT " ",(tagint) ubuf(buf[i][j]).i);
    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------
   pack element info for Element Clusters section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::pack_element_cluster_data(double **buf)
{
  int m = 0;
  for (int i = 0; i < element->nlocal; i++) {
    if (tag[i] == element_clusters[i][0]) {
      for (int j = 0; j < element->max_apc; j++)
        buf[m][j] = ubuf(element_clusters[i][j]).d;
      m++;
    }
  }
}

/* ----------------------------------------------------
   write node info into Nodes section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::write_node_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e " TAGINT_FORMAT " \n",
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4],(tagint) ubuf(buf[i][5]).i);
}

/* ----------------------------------------------------
   pack node info for Nodes section in data file from write_data command
   ------------------------------------------*/

void ElementVecCAC::pack_node_data(double **buf)
{
  int k = 0;
  int inpe;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      buf[k][0] = ubuf(tag[i]).d;
      buf[k][1] = ubuf(j+1).d;
      buf[k][2] = nodex[i][j][0];
      buf[k][3] = nodex[i][j][1];
      buf[k][4] = nodex[i][j][2];
      buf[k][5] = ubuf(nodetag[i][j]).d; 
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file from write_data command
   ------------------------------------------------------------------------- */

void ElementVecCAC::pack_vel_data(double **buf)
{
  int k = 0;
  int inpe;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      buf[k][0] = ubuf(tag[i]).d;
      buf[k][1] = ubuf(j+1).d;
      buf[k][2] = nodev[i][j][0];
      buf[k][3] = nodev[i][j][1];
      buf[k][4] = nodev[i][j][2];
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   write velocity info to data file from write_data command
   ------------------------------------------------------------------------- */

void ElementVecCAC::write_vel_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n",
        (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4]);
}

/* ----------------------------------------------------------------------
   write element types information to Element Types section in data file from write_data command
   ------------------------------------------------------------------------- */

void ElementVecCAC::write_element_types(FILE *fp)
{
  for (int i = 1; i <= element->netypes; i++) 
    fprintf(fp,"%d %s %d %d %d %d %d %d %d\n",
        i,element_shape_names[i],apc[i],
        ncells[i][0],ncells[i][1],ncells[i][2],
        nintgs[i][0],nintgs[i][1],nintgs[i][2]);
}

/*----------------------------------------------------
  write node info into Nodes section in dat file from write_tecplot command
  ------------------------------------------*/

void ElementVecCAC::write_node_tecplot(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",
        buf[i][0],buf[i][1],buf[i][2]);
}

/*----------------------------------------------------
  write node connectivity info into Nodes section in .dat file from write_tecplot_command
  ------------------------------------------*/

void ElementVecCAC::write_node_connect_tecplot(FILE *fp, int n, int *buf)
{
  int inpe;
  int k = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < size_tecplot_node_connect; j++)
      fprintf(fp,TAGINT_FORMAT" ",(tagint) buf[k++]);
    fprintf(fp,"\n");
  }
}

/*----------------------------------------------------
  pack node info for Nodes section in .dat file from write_tecplot command
  ------------------------------------------*/

void ElementVecCAC::pack_node_tecplot(double **buf)
{
  int inpe;
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      buf[k][0] = nodex[i][j][0];
      buf[k][1] = nodex[i][j][1];
      buf[k][2] = nodex[i][j][2];
      k++;
    }
  }
}

/*----------------------------------------------------
  pack node info for Nodes section in plt or szplt file from write_tecplot command
  ------------------------------------------*/

void ElementVecCAC::pack_node_tecplot_binary(double *buf, int ival)
{
  int inpe;
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      buf[k] = nodex[i][j][ival];
      k++;
    }
  }
}

/*----------------------------------------------------
  pack node connectivity info for Nodes section in dat file
  specifically set for each element shape
  ------------------------------------------*/

void ElementVecCAC::pack_node_connect_tecplot(int *buf)
{
  int itype;
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    itype = etype[i];

    if (element_shape_ids[itype] == QUADRILATERAL) {
      for (int j = 0; j < 4; j++) 
        buf[k++] = nodetag[i][j];

    } else if (element_shape_ids[itype] == TRIANGLE) {

      // node connectivity scheme: 1 2 3 3

      for (int j = 0; j < 3; j++) 
        buf[k++] = nodetag[i][j];
      buf[k++] = nodetag[i][2];

    } else if (element_shape_ids[itype] == HEXAHEDRON) {

      for (int j = 0; j < 8; j++) 
        buf[k++] = nodetag[i][j];

    } else if (element_shape_ids[itype] == PYRAMID) {

      // node connectivity scheme: 1 2 3 4 5 5 5 5

      for (int j = 0; j < 5; j++) 
        buf[k++] = nodetag[i][j];
      for (int j = 5; j < 8; j++)
        buf[k++] = nodetag[i][4];

    } else if (element_shape_ids[itype] == TETRAHEDRON) {

      // node connectivity scheme: 1 2 3 3 4 4 4 4

      for (int j = 0; j < 3; j++) 
        buf[k++] = nodetag[i][j];
      buf[k++] = nodetag[i][2];
      for (int j = 4; j < 8; j++) 
        buf[k++] = nodetag[i][3];

    } else if (element_shape_ids[itype] == WEDGE) {

      // node connectivity scheme: 1 2 3 3 4 5 6 6

      for (int j = 0; j < 3; j++) 
        buf[k++] = nodetag[i][j];
      for (int j = 3; j < 7; j++) 
        buf[k++] = nodetag[i][j-1];
      buf[k++] = nodetag[i][5];
    }
  }
}

/* ---------------------------------------------
   Setup element types info (shape, interpolation, integration)
   -------------------------------*/

void ElementVecCAC::set_element_types(const char *str, int type_offset)
{
  int itype,iapc;
  int element_shape_id = -1;
  int ncellx,ncelly,ncellz,nintgx,nintgy,nintgz;
  char element_shape_name[256];

  // remove user input for intergration points later and have a fixed scheme

  int n = sscanf(str,"%d %s %d %d %d %d %d %d %d",&itype,element_shape_name,&iapc,&ncellx,&ncelly,&ncellz,&nintgx,&nintgy,&nintgz);
  if (n != 9) error->all(FLERR, "Invalid Element Type line in data file");
  itype += type_offset;

  if (itype < 1 || itype > element->netypes)
    error->all(FLERR,"Invalid element type for interpolate set");
  if (iapc < 1) error->all(FLERR,"Invalid apc value for this element type");
  if (iapc > element->max_apc) error->all(FLERR,"apc for this type was greater than maxapc defined, increase maxapc through element_modify command");
  // determine element shape id

  for (int i = 0; i < element->num_element_shapes; i++) {
    if (strcmp(element_shape_name,element->element_shape_list[i]) == 0)
      element_shape_id = i;
  }

  if (element_shape_id == -1) error->all(FLERR,"Invalid element shape name");

  if (element_type_setflag[itype]) {

    // if identical element type is defined, skip the setup
    // else throw error;

    if (ncellx == ncells[itype][0] &&
        ncelly == ncells[itype][1] &&
        ncellz == ncells[itype][2] && 
        nintgx == nintgs[itype][0] &&
        nintgy == nintgs[itype][1] &&
        nintgz == nintgs[itype][2] &&
        apc[itype] == iapc &&
        element_shape_id == element_shape_ids[itype]) 
      return;

    error->all(FLERR,"This element type has been set");
  }

  apc[itype] = iapc;
  if (iapc > 1 && element->element_cluster_flag == 0) 
    element->element_cluster_flag = 1;
  element_shape_ids[itype] = element_shape_id;
  strcpy(element_shape_names[itype],element_shape_name);
  element_type_setflag[itype] = 1;

  ncells[itype][0] = ncellx;
  ncells[itype][1] = ncelly;
  ncells[itype][2] = ncellz;

  nintgs[itype][0] = nintgx;
  nintgs[itype][1] = nintgy;
  nintgs[itype][2] = nintgz;

  // check input error for each element shapes, specify constants for this type
  // add here to include more element shapes

  if (element_shape_id == QUADRILATERAL) {
    if (domain->dimension == 3) error->all(FLERR,"2D element requires 2D simulation");
    //if (ncellz != 1) error->all(FLERR,"Number of cells along third direction must be one for 2D simulation"); 
    if (ncellx < 2 || ncelly < 2) 
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl || 
        ncelly-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if (nintgx < 2 || nintgy < 2 || nintgx > ncellx || nintgy > ncelly) 
      error->all(FLERR,"Invalid integration values");

    if ((ncellx-2)*(ncelly-2) > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceed max_surface_intpl");
    npe[itype] = 4;
    n = ncellx*ncelly;
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    n = nintgx * nintgy;
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == TRIANGLE) {
    if (domain->dimension == 3) error->all(FLERR,"2D element requires 2D simulation");
    //if (ncellz != 1) error->all(FLERR,"Number of cells along third direction must be one for 2D simulation"); 
    if (ncellx < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if (nintgx < 2 || nintgx > ncellx)
      error->all(FLERR,"Invalid integration values");

    if ((ncellx - 2)*(ncellx - 1)/2 > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceed max_surface_intpl");
    npe[itype] = 3;
    n = ncellx*(ncellx + 1)/2;
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    n = nintgx*(nintgx + 1)/2;
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == HEXAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR,"3D element requires 3D simulation");
    if (ncellx < 2 || ncelly < 2 || ncellz < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl ||
        ncelly-2 > element->max_edge_intpl ||
        ncellz-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if ((ncellx-2)*(ncelly-2) > element->max_surface_intpl ||
        (ncellz-2)*(ncelly-2) > element->max_surface_intpl ||
        (ncellx-2)*(ncellz-2) > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceeds max_surface_intpl");
    if (nintgx < 2 || nintgy < 2 || nintgz < 2 ||
        nintgx > ncellx || nintgy > ncelly || nintgz > ncellz) 
      error->all(FLERR,"Invalid integration values");

    npe[itype] = 8;
    n = ncellx*ncelly*ncellz;
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    n = nintgx * nintgy * nintgz;
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == TETRAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR,"3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl)
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if ((ncellx-2)*(ncellx-1)/2 > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceeds max_surface_intpl");
    if (nintgx < 2 || nintgx > ncellx)
      error->all(FLERR,"Invalid integration values");
    npe[itype] = 4;
    n = sum_tetrahedron(ncellx);
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;
    n = sum_tetrahedron(nintgx);
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == OCTAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR,"3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if ((ncellx-2)*(ncellx-1)/2 > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceeds max_surface_intpl");
    if (nintgx < 2 || nintgx > 4 ||
        (nintgx < 4 && nintgx != ncellx))
      error->all(FLERR,"Invalid integration values");

    npe[itype] = 6;
    n = sum_octahedron(ncellx);
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    n = sum_octahedron(nintgx);
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == PYRAMID) {
    if (domain->dimension == 2) error->all(FLERR,"3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if ((ncellx-2)*(ncellx-2) > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceeds max_surface_intpl");
    if (nintgx < 2 || nintgx > 4 ||
        (nintgx < 4 && nintgx != ncellx))
      error->all(FLERR,"Invalid integration values");

    npe[itype] = 5;
    n = sum_pyramid(ncellx);
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    if (ncellx > 4)
      n = 42;
    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  } else if (element_shape_id == WEDGE) {
    if (domain->dimension == 2) error->all(FLERR,"3D element requires 3D simulation");
    if (ncellx < 2 || ncelly < 2)
      error->all(FLERR,"Invalid interpolate values");
    if (ncellx-2 > element->max_edge_intpl || ncelly-2 > element->max_edge_intpl) 
      error->all(FLERR,"Number of interpolated atoms on an edge exceeds max_edge_intpl");
    if ((ncellx-2)*(ncelly-2) > element->max_surface_intpl ||
        (ncelly-1)*(ncelly-2)/2 > element->max_surface_intpl)
      error->all(FLERR,"Number of interpolated atoms on a surface exceeds max_surface_intpl");
    if (nintgy < 2 || nintgy > 4 ||
        (nintgy < 4 && nintgy != ncelly) ||
        nintgx < 2 || nintgx > ncellx)
      error->all(FLERR,"Invalid integration values");

    npe[itype] = 6;
    n = ncellx*(ncelly + 1)*ncelly/2;
    if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
    element->max_nintpl = MAX(element->max_nintpl,n);
    nintpl[itype] = n;

    if (ncelly > 4)
      n = 12*nintgx;
    else 
      n = (nintgy + 1)*nintgy*nintgx/2;

    if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
    element->max_nintg = MAX(element->max_nintg,n);
    nintg[itype] = n;

  }

  // setup interpolation, sub-element division, and integration scheme for this element type

  setup_interpolate_element(itype);
  setup_integration_point(itype);
  setup_sub_element(itype);
}

/*-----------------------------------------------------------------------------------
  Calculate interpolated atom related arrays for an element type
  - Shape function array to interpolate atoms
  - List of element nodes in natural coordinate system
  + Quadrilateral
  - 1: (-1,-1)
  - 2: (1,-1)
  - 3: (1,1)
  - 4: (-1,1)
  + Triangle
  - 1: (0,0)
  - 2: (1,0)
  - 3: (0,1)
  + Hexahedron
  - 1: (-1,-1,-1)
  - 2: (1,-1,-1)
  - 3: (1,1,-1)
  - 4: (-1,1,-1)
  - 5: (-1,-1,1)
  - 6: (1,-1,1)
  - 7: (1,1,1)
  - 8: (-1,1,1)
  + Tetrahedron
  - 1: (0,0,0)
  - 2: (1,0,0)
  - 3: (0,1,0)
  - 4: (0,0,1)
  + Pyramid
  - 1: (-1,-1,0)
  - 2: (1,-1,0)
  - 3: (1,1,0)
  - 4: (-1,1,0)
  - 5: (0,0,1)
  + Wedge
  - 1: (-1,0,0)
  - 2: (-1,1,0)
  - 3: (-1,0,1)
  - 4: (1,0,0)
  - 5: (1,1,0)
  - 6: (1,0,1)

  - Lists of indices of interpolated atoms on edges and surfaces
  ------------------------------------------------------------------------------------*/

void ElementVecCAC::setup_interpolate_element(int itype)
{
  int ncellx = ncells[itype][0];
  int ncelly = ncells[itype][1];
  int ncellz = ncells[itype][2];

  // interpolated atom index counter;

  int n = 0;

  // calculate shape function arrays based on its shape and size
  // add here to include other element shapes

  if (element_shape_ids[itype] == QUADRILATERAL) {

    int nedge[4];
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = 0;

    double px,py;
    double dx_pc,dy_pc;
    dx_pc = 2.0/(ncellx - 1);
    dy_pc = 2.0/(ncelly - 1);
    for (int i = 0; i < ncellx; i++) {
      px = dx_pc*i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dy_pc*j - 1.0;
        shape_array[itype][n][0] = 0.25*(1 - px)*(1 - py);
        shape_array[itype][n][1] = 0.25*(1 + px)*(1 - py);
        shape_array[itype][n][2] = 0.25*(1 + px)*(1 + py);
        shape_array[itype][n][3] = 0.25*(1 - px)*(1 + py);

        // store indices of interpolated atoms on edges for this element type

        if (i == 0 && j != 0 && j != ncelly - 1)
          edge_intpl[itype][0][nedge[0]++] = n; 
        else if (i == ncellx - 1 && j != 0 && j != ncelly - 1)
          edge_intpl[itype][1][nedge[1]++] = n; 
        else if (j == 0 && i != 0 && i != ncellx - 1)
          edge_intpl[itype][2][nedge[2]++] = n; 
        else if (j == ncelly - 1 && i != 0 && i != ncellx - 1)
          edge_intpl[itype][3][nedge[3]++] = n; 

        n++;
      }
    }
    for (int i = 0; i < 4; i++)
      nedge_intpl[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == TRIANGLE) {
    int nedge[3];
    double px,py;
    nedge[0] = nedge[1] = nedge[2] = 0;
    double dxy_pc = 1.0/(ncellx - 1);

    for (int i = 0; i < ncellx; i++) {
      px = dxy_pc*i;
      for (int j = 0; j < ncellx - i; j++) {
        py = dxy_pc*j;
        shape_array[itype][n][0] = px;
        shape_array[itype][n][1] = py;
        shape_array[itype][n][2] = 1- px - py;

        // store indices of interpolated atoms on edges for this element type

        if (i == 0 && j != 0 && j != ncelly - 1)
          edge_intpl[itype][0][nedge[0]++] = n; 
        if (j == 0 && i != 0 && i != ncellx - 1)
          edge_intpl[itype][1][nedge[1]++] = n; 
        if (i+j == ncelly - 1 && i != 0 && j != 0)
          edge_intpl[itype][2][nedge[2]++] = n; 

        n++;
      }
    }
    for (int i = 0; i < 3; i++)
      nedge_intpl[itype][i] = nedge[i];


  } else if (element_shape_ids[itype] == HEXAHEDRON) {

    int nsur[6],nedge[12];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = nsur[5] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = nedge[9] = nedge[10] = nedge[11] = 0;
    double px,py,pz;
    double dx_pc,dy_pc,dz_pc;
    dx_pc = 2.0/(ncellx - 1);
    dy_pc = 2.0/(ncelly - 1);
    dz_pc = 2.0/(ncellz - 1);

    for (int i = 0; i < ncellx; i++) {
      px = dx_pc*i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dy_pc*j - 1.0;
        for (int k = 0; k < ncellz; k++) {
          pz = dz_pc*k - 1.0;
          shape_array[itype][n][0] = 0.125*(1 - px)*(1 - py)*(1 - pz);
          shape_array[itype][n][1] = 0.125*(1 + px)*(1 - py)*(1 - pz);
          shape_array[itype][n][2] = 0.125*(1 + px)*(1 + py)*(1 - pz);
          shape_array[itype][n][3] = 0.125*(1 - px)*(1 + py)*(1 - pz);
          shape_array[itype][n][4] = 0.125*(1 - px)*(1 - py)*(1 + pz);
          shape_array[itype][n][5] = 0.125*(1 + px)*(1 - py)*(1 + pz);
          shape_array[itype][n][6] = 0.125*(1 + px)*(1 + py)*(1 + pz);
          shape_array[itype][n][7] = 0.125*(1 - px)*(1 + py)*(1 + pz);

          // store indices of interpolated atoms on edges and surfaces for this element type       

          if (k != 0 && k != ncellz-1) {
            if (i == 0 && j == 0) {
              edge_intpl[itype][0][nedge[0]++] = n; 
            } else if (i == 0 && j == ncelly-1) {
              edge_intpl[itype][1][nedge[1]++] = n; 
            } else if (i == ncellx-1 && j == 0) {
              edge_intpl[itype][2][nedge[2]++] = n; 
            } else if (i == ncellx-1 && j == ncelly-1) {
              edge_intpl[itype][3][nedge[3]++] = n; 
            } else if (i == 0) {
              surface_intpl[itype][0][nsur[0]++] = n;
            } else if (i == ncellx-1) {
              surface_intpl[itype][1][nsur[1]++] = n;
            } else if (j == 0) {
              surface_intpl[itype][2][nsur[2]++] = n;
            } else if (j == ncelly-1) {
              surface_intpl[itype][3][nsur[3]++] = n;
            }
          } else if (k == 0) {
            if (j != 0 && j != ncelly-1) {
              if (i == 0) {
                edge_intpl[itype][4][nedge[4]++] = n; 
              } else if (i == ncellx-1) {
                edge_intpl[itype][6][nedge[6]++] = n; 
              } else {
                surface_intpl[itype][4][nsur[4]++] = n;
              }
            } else if (i != 0 && i != ncellx-1) {
              if (j == 0) {
                edge_intpl[itype][8][nedge[8]++] = n; 
              } else {
                edge_intpl[itype][10][nedge[10]++] = n; 
              }
            }
          } else {
            if (j != 0 && j != ncelly-1) {
              if (i == 0) {
                edge_intpl[itype][5][nedge[5]++] = n; 
              } else if (i == ncellx-1) {
                edge_intpl[itype][7][nedge[7]++] = n; 
              } else {
                surface_intpl[itype][5][nsur[5]++] = n;
              }
            } else if (i != 0 && i != ncellx-1) {
              if (j == 0) {
                edge_intpl[itype][9][nedge[9]++] = n; 
              } else {
                edge_intpl[itype][11][nedge[11]++] = n; 
              }
            }
          }
          n++;
        }
      }
    }
    for (int i = 0; i < 6; i++)
      nsurface_intpl[itype][i] = nsur[i];
    for (int i = 0; i < 12; i++)
      nedge_intpl[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == TETRAHEDRON) {

    int nsur[4],nedge[6];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = 0;
    nedge[0] = nedge[1] = nedge[2] = 0;
    nedge[3] = nedge[4] = nedge[5] = 0;
    double px,py,pz;
    double dxyz_pc;
    int xnum,ynum;
    dxyz_pc = 1.0/(ncellx-1);

    for (int k = 0; k < ncellx; k++) {
      pz = dxyz_pc*k;
      xnum = ncellx - k;
      for (int i = 0; i < xnum; i++) {
        px = dxyz_pc*i;
        ynum = xnum - i;
        for (int j = 0; j < ynum; j++) {
          py = dxyz_pc*j;
          shape_array[itype][n][0] = 1 - px - py - pz;
          shape_array[itype][n][1] = px;
          shape_array[itype][n][2] = py;
          shape_array[itype][n][3] = pz;

          // store indices of interpolated atoms on edges and surfaces for this element type       

          if (i != 0 && j != 0 && k == 0 && i+j+k != ncellx-1)          // xy plane
            surface_intpl[itype][0][nsur[0]++] = n;
          if (i == 0 && j != 0 && k != 0 && i+j+k != ncellx-1)          // yz plane
            surface_intpl[itype][1][nsur[1]++] = n;
          if (i != 0 && j == 0 && k != 0 && i+j+k != ncellx-1)          // xz plane
            surface_intpl[itype][2][nsur[2]++] = n;
          if (i != 0 && j != 0 && k != 0 && i+j+k == ncellx-1)          // x+y+z = 1 plane
            surface_intpl[itype][3][nsur[3]++] = n;

          if (i == 0 && j == 0 && k != 0 && k != ncellx-1)              // z axis
            edge_intpl[itype][0][nedge[0]++] = n;
          if (i == 0 && k == 0 && j != 0 && j != ncellx-1)              // y axis
            edge_intpl[itype][1][nedge[1]++] = n;
          if (j == 0 && k == 0 && i != 0 && i != ncellx-1)              // x axis
            edge_intpl[itype][2][nedge[2]++] = n;
          if (i == 0 && j+k == ncellx-1 && j != 0 && k != 0)            // y+z = 1 line
            edge_intpl[itype][3][nedge[3]++] = n;
          if (j == 0 && i+k == ncellx-1 && i != 0 && k != 0)            // x+z = 1 line
            edge_intpl[itype][4][nedge[4]++] = n;
          if (k == 0 && i+j == ncellx-1 && i != 0 && j != 0)            // x+y = 1 line
            edge_intpl[itype][5][nedge[5]++] = n;

          n++;
        }
      }
    }

    for (int i = 0; i < 4; i++) 
      nsurface_intpl[itype][i] = nsur[i];
    for (int i = 0; i < 6; i++)
      nedge_intpl[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == OCTAHEDRON) {

    error->all(FLERR,"Shape function for octahedron element is not found yet");
    int nsur[8],nedge[12];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = nsur[5] = nsur[6] = nsur[7] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = nedge[9] = nedge[10] = nedge[11] = 0;
    double px,py,pz;
    double dxy_pc,dz_pc;
    double xycorner;
    int xynum;
    dxy_pc = 2.0/(ncellx - 1);
    dz_pc = 1.0/(ncellx - 1);

    for (int k = 1 - ncellx; k <= ncellx - 1; k++) {
      pz = dz_pc*k;
      xynum = ncellx - abs(k);            // size of current xy plane
      xycorner = -(xynum - 1)*dz_pc;
      for (int i = 0; i < xynum; i++) {
        px = dxy_pc*i + xycorner;
        for (int j = 0; j < xynum; j++) {
          py = dxy_pc*j + xycorner;

          // THIS IS INCORRECT, NEED TO FIND NEW SHAPE FUNCTION
          shape_array[itype][n][0] = 0.5*pz*(pz-1);
          shape_array[itype][n][1] = 0.25*(1-pz*pz)*(1-px)*(1-py);
          shape_array[itype][n][2] = 0.25*(1-pz*pz)*(1+px)*(1-py);
          shape_array[itype][n][3] = 0.25*(1-pz*pz)*(1+px)*(1+py);
          shape_array[itype][n][4] = 0.25*(1-pz*pz)*(1-px)*(1+py);
          shape_array[itype][n][5] = 0.5*pz*(pz+1);

          // store indices of interpolated atoms on edges and surfaces for this element type       

          if (k > 1 - ncellx && k < 0) {          // lower half
            if (i == 0 && j == 0)                       // on edge
              edge_intpl[itype][0][nedge[0]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_intpl[itype][1][nedge[1]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_intpl[itype][2][nedge[2]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_intpl[itype][3][nedge[3]++] = n;
            else if (i == 0)                            // on surface
              surface_intpl[itype][0][nsur[0]++] = n;
            else if (j == 0)                            // on surface
              surface_intpl[itype][1][nsur[1]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_intpl[itype][2][nsur[2]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_intpl[itype][3][nsur[3]++] = n;
          } else if (k == 0) {                    // middle plane
            if (!(i == 0 && j == 0) &&       
                !(i == xynum - 1 && j == 0) &&
                !(i == 0 && j == xynum - 1) &&
                !(i == xynum - 1 && j == xynum - 1)) {  // exclude 4 corners
              if (i == 0)                               // on edge
                edge_intpl[itype][4][nedge[4]++] = n;
              else if (j == 0)                          // on edge
                edge_intpl[itype][5][nedge[5]++] = n;
              else if (i == xynum - 1)                  // on edge
                edge_intpl[itype][6][nedge[6]++] = n;
              else if (j == xynum - 1)                  // on edge
                edge_intpl[itype][7][nedge[7]++] = n;
            }
          } else if (k > 0 && k < ncellx - 1) {   // upper half
            if (i == 0 && j == 0)                       // on edge
              edge_intpl[itype][8][nedge[8]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_intpl[itype][9][nedge[9]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_intpl[itype][10][nedge[10]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_intpl[itype][11][nedge[11]++] = n;
            else if (i == 0)                            // on surface
              surface_intpl[itype][4][nsur[4]++] = n;
            else if (j == 0)                            // on surface
              surface_intpl[itype][5][nsur[5]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_intpl[itype][6][nsur[6]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_intpl[itype][7][nsur[7]++] = n;
          }
          n++;
        }
      }
    }
    for (int i = 0; i < 8; i++)
      nsurface_intpl[itype][i] = nsur[i];
    for (int i = 0; i < 12; i++)
      nedge_intpl[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == PYRAMID) {

    int nsur[5],nedge[8];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = 0;
    nedge[4] = nedge[5] = nedge[6] = nedge[7] = 0;
    double px,py,pz;
    double dxy_pc,dz_pc;
    double xycorner;
    int xynum;
    dxy_pc = 2.0/(ncellx-1);
    dz_pc = 1.0/(ncellx-1);

    for (int k = 0; k < ncellx-1; k++) {
      pz = dz_pc*k;
      xynum = ncellx - k;                // size of current xy plane
      xycorner = -(xynum - 1)*dz_pc;
      for (int i = 0; i < xynum; i++) {
        px = dxy_pc*i + xycorner;
        for (int j = 0; j < xynum; j++) {
          py = dxy_pc*j + xycorner;
          shape_array[itype][n][0] = 0.25*(1.0 - px - pz)*(1.0 - py - pz)/(1 - pz);
          shape_array[itype][n][1] = 0.25*(1.0 + px - pz)*(1.0 - py - pz)/(1 - pz);
          shape_array[itype][n][2] = 0.25*(1.0 + px - pz)*(1.0 + py - pz)/(1 - pz);
          shape_array[itype][n][3] = 0.25*(1.0 - px - pz)*(1.0 + py - pz)/(1 - pz);
          shape_array[itype][n][4] = pz;

          // store indices of interpolated atoms on edges and surfaces for this element type       

          if (k == 0) {                                 // base plane
            if (!(i == 0 && j == 0) &&       
                !(i == xynum - 1 && j == 0) &&
                !(i == 0 && j == xynum - 1) &&
                !(i == xynum - 1 && j == xynum - 1)) {  // exclude 4 corners
              if (i == 0)                               // on edge
                edge_intpl[itype][0][nedge[0]++] = n;
              else if (j == 0)                          // on edge
                edge_intpl[itype][1][nedge[1]++] = n;
              else if (i == xynum - 1)                  // on edge
                edge_intpl[itype][2][nedge[2]++] = n;
              else if (j == xynum - 1)                  // on edge
                edge_intpl[itype][3][nedge[3]++] = n;
              else                                      // on surface
                surface_intpl[itype][0][nsur[0]++] = n;
            }
          } else if (k < ncellx - 1) {                  
            if (i == 0 && j == 0)                       // on edge
              edge_intpl[itype][4][nedge[4]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_intpl[itype][5][nedge[5]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_intpl[itype][6][nedge[6]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_intpl[itype][7][nedge[7]++] = n;
            else if (i == 0)                            // on surface
              surface_intpl[itype][1][nsur[1]++] = n;
            else if (j == 0)                            // on surface
              surface_intpl[itype][2][nsur[2]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_intpl[itype][3][nsur[3]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_intpl[itype][4][nsur[4]++] = n;
          }
          n++;
        }
      }
    }

    // do peak point separately to avoid singularity (pz = 1)

    shape_array[itype][n][0] = 0.0;
    shape_array[itype][n][1] = 0.0;
    shape_array[itype][n][2] = 0.0;
    shape_array[itype][n][3] = 0.0;
    shape_array[itype][n][4] = 1.0;

    n++;

    for (int i = 0; i < 5; i++) 
      nsurface_intpl[itype][i] = nsur[i];
    for (int i = 0; i < 8; i++)
      nedge_intpl[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == WEDGE) {

    int nsur[5],nedge[9];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = 0;
    nedge[0] = nedge[1] = nedge[2] = 0;
    nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = 0;
    double px,py,pz;
    int znum;
    double dx_pc = 2.0/(ncellx-1);
    double dyz_pc = 1.0/(ncelly-1);

    for (int i = 0; i < ncellx; i++) {
      px = dx_pc*i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dyz_pc*j;
        znum = ncelly - j;
        for (int k = 0; k < znum; k++) {
          pz = dyz_pc*k;

          shape_array[itype][n][0] = 0.5*(1.0 - px)*(1.0 - py - pz);
          shape_array[itype][n][1] = 0.5*(1.0 - px)*py;
          shape_array[itype][n][2] = 0.5*(1.0 - px)*pz;
          shape_array[itype][n][3] = 0.5*(1.0 + px)*(1.0 - py - pz);
          shape_array[itype][n][4] = 0.5*(1.0 + px)*py;
          shape_array[itype][n][5] = 0.5*(1.0 + px)*pz;

          // store indices of interpolated atoms on edges and surfaces for this element type       

          if (i == ncellx-1 && j != 0 && k != 0 && j+k != ncelly-1)             // +yz plane
            surface_intpl[itype][0][nsur[0]++] = n;
          if (i == 0 && j != 0 && k != 0 && j+k != ncelly-1)                    // -yz plane
            surface_intpl[itype][1][nsur[1]++] = n;
          if (i != 0 && i != ncellx-1 && j != 0 && k == 0 && j+k != ncelly-1)   // xy plane
            surface_intpl[itype][2][nsur[2]++] = n;
          if (i != 0 && i != ncellx-1 && j == 0 && k != 0 && j+k != ncelly-1)   // xz plane
            surface_intpl[itype][3][nsur[3]++] = n;
          if (i != 0 && i != ncellx-1 && j != 0 && k != 0 && j+k == ncelly-1)   // y + z = 1 plane
            surface_intpl[itype][4][nsur[4]++] = n;

          if (i == ncellx-1 && j == 0 && k != 0 && k != ncelly-1)               // x = 1, y = 0 line
            edge_intpl[itype][0][nedge[0]++] = n;
          if (i == ncellx-1 && k == 0 && j != 0 && j != ncelly-1)               // x = 1, z = 0 line
            edge_intpl[itype][1][nedge[1]++] = n;
          if (i == ncellx-1 && j+k == ncelly-1 && j != 0 && k != 0)             // x = 1, y+z = 1 line
            edge_intpl[itype][2][nedge[2]++] = n;
          if (i == 0 && j == 0 && k != 0 && k != ncelly-1)                      // x = -1, y = 0 line
            edge_intpl[itype][3][nedge[3]++] = n;
          if (i == 0 && k == 0 && j != 0 && j != ncelly-1)                      // x = -1, z = 0 line
            edge_intpl[itype][4][nedge[4]++] = n;
          if (i == 0 && j+k == ncelly-1 && j != 0 && k != 0)                    // x = -1, y+z = 1 line
            edge_intpl[itype][5][nedge[5]++] = n;
          if (i != 0 && i != ncellx-1 && j == 0 && k == 0)                      // y = 0, z = 0 line
            edge_intpl[itype][6][nedge[6]++] = n;
          if (i != 0 && i != ncellx-1 && j == ncelly-1)                         // y = 0 line
            edge_intpl[itype][7][nedge[7]++] = n;
          if (i != 0 && i != ncellx-1 && k == ncelly-1)                         // z = 0 line
            edge_intpl[itype][8][nedge[8]++] = n;

          n++;
        }
      }
    }

    for (int i = 0; i < 4; i++) 
      nsurface_intpl[itype][i] = nsur[i];
    for (int i = 0; i < 6; i++)
      nedge_intpl[itype][i] = nedge[i];


  } else if (0) {
    // add new interpolation scheme for new element shape here
  }
  if (n != nintpl[itype]) error->all(FLERR,"TEST"); 
}

/* -----------------------------------------------------------------------------------
   Calculate integration point related arrays:
   - Weight of integration points
   - i2ia to convert integration point index 
   to interpolated atom index in an element,
   used for interpolating integration points
   - n2i to convert node index to integration point index
   Must be done after setting up interpolate element
   ------------------------------------------------------------------------------------*/

void ElementVecCAC::setup_integration_point(int itype)
{
  int nintgx,nintgy,nintgz;
  int ncellx,ncelly,ncellz;
  int iintpl;

  int iintg = 0;                      // integration point index counter

  double **xgauss = element->xgauss;
  double **wgauss = element->wgauss;

  // number of interpolated atoms along each edge

  ncellx = ncells[itype][0];
  ncelly = ncells[itype][1];
  ncellz = ncells[itype][2];

  // number of integration points along each edge

  nintgx = nintgs[itype][0];
  nintgy = nintgs[itype][1];
  nintgz = nintgs[itype][2];

  // determine the indices of integration points and their weights for this element type

  if (element_shape_ids[itype] == QUADRILATERAL) {

    // calculate weight fraction for integration points

    int a,b,w;
    int px,py;
    double w_edge_x;
    double w_edge_y;
    double w_inner;
    double scale = ncellx * ncelly;
    w_edge_x = static_cast<double> (ncellx - 2)/2.0;
    w_edge_y = static_cast<double> (ncelly - 2)/2.0;
    w_inner = w_edge_x * w_edge_y;
    double dx = (ncellx-1) / static_cast<double> (nintgx - 1);
    double dy = (ncelly-1) / static_cast<double> (nintgy - 1);
    for (int i = 0; i < nintgx; i++) {
      if (i == 0) px = 0;
      else if (i == nintgx-1) px = ncellx-1;
      else px = (int) ((1.0 + xgauss[nintgx-2][i])/2.0*(ncellx-2)+1.0);
      for (int j = 0; j < nintgy; j++) {
        if (j == 0) py = 0;
        else if (j == nintgy-1) py = ncelly-1;
        else py = (int) ((1.0 + xgauss[nintgy-2][j])/2.0*(ncelly-2)+1.0);

        a = (i == 0 || i == (nintgx - 1));
        b = (j == 0 || j == (nintgy - 1));
        w = a + b;

        // integration point is a node

        if (w == 2) weight[itype][iintg] = 1;

        // integration point is on an edge

        else if (w == 1) {

          // integration point is on edge in x direction

          if (!a) weight[itype][iintg] = w_edge_x;

          // integration point is on edge in y direction

          else if (!b) weight[itype][iintg] = w_edge_y;

          // integration point is inside element

        } else weight[itype][iintg] = w_inner;

        // list to convert integration point index to 
        // interpolated atom index

        iintpl = px*ncelly + py;
        i2ia[itype][iintg] = iintpl;
        ia2i[itype][iintpl] = iintg;

        for (int inode = 0; inode < 4; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = 
            shape_array[itype][iintpl][inode] * weight[itype][iintg] * 4.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
        iintg++;
      }
    }

  } else if (element_shape_ids[itype] == HEXAHEDRON) {

    // calculate weight fraction for integration points

    int a,b,c,w;
    int px,py,pz;
    double w_edge_x,w_face_xy;
    double w_edge_y,w_face_yz;
    double w_edge_z,w_face_xz;
    double w_inner;
    w_edge_x = static_cast<double> (ncellx - 2)/2.0;
    w_edge_y = static_cast<double> (ncelly - 2)/2.0;
    w_edge_z = static_cast<double> (ncellz - 2)/2.0;
    w_face_xy = w_edge_x * w_edge_y;
    w_face_yz = w_edge_y * w_edge_z;
    w_face_xz = w_edge_x * w_edge_z;
    w_inner = w_edge_x * w_edge_y * w_edge_z;
    double dx = (ncellx-1) / static_cast<double> (nintgx - 1);
    double dy = (ncelly-1) / static_cast<double> (nintgy - 1);
    double dz = (ncellz-1) / static_cast<double> (nintgz - 1);
    for (int i = 0; i < nintgx; i++) {
      if (i == 0) px = 0;
      else if (i == nintgx-1) px = ncellx-1;
      else px = (int) ((1.0 + xgauss[nintgx-2][i])/2.0*(ncellx-2)+1.0);
      for (int j = 0; j < nintgy; j++) {
        if (j == 0) py = 0;
        else if (j == nintgy-1) py = ncelly-1;
        else py = (int) ((1.0 + xgauss[nintgy-2][j])/2.0*(ncelly-2)+1.0);
        for (int k = 0; k < nintgz; k++) {
          if (k == 0) pz = 0;
          else if (k == nintgz-1) pz = ncellz-1;
          else pz = (int) ((1.0 + xgauss[nintgz-2][k])/2.0*(ncellz-2)+1.0);

          // calculate weight of integration points

          a = (i == 0 || i == (nintgx - 1));
          b = (j == 0 || j == (nintgy - 1));
          c = (k == 0 || k == (nintgz - 1));
          w = a + b + c;

          // integration point is a node

          if (w == 3) weight[itype][iintg] = 1.0;

          // integration point is on an edge

          else if (w == 2) {

            // integration point is on edge in x direction

            if (!a) weight[itype][iintg] = w_edge_x*wgauss[nintgx-2][i];

            // integration point is on edge in y direction

            else if (!b) weight[itype][iintg] = w_edge_y*wgauss[nintgy-2][j];

            // integration point is on edge in z direction

            else if (!c) weight[itype][iintg] = w_edge_z*wgauss[nintgy-2][k];

            // integration point is on a face

          } else if (w == 1) {

            // integration point is on yz face 

            if (a) weight[itype][iintg] = w_face_yz*wgauss[nintgy-2][j]*wgauss[nintgz-2][k];

            // integration point is on xz face 

            else if (b) weight[itype][iintg] = w_face_xz*wgauss[nintgx-2][i]*wgauss[nintgz-2][k];

            // integration point is on xy face 

            else if (c) weight[itype][iintg] = w_face_xy*wgauss[nintgx-2][i]*wgauss[nintgy-2][j];

            // integration point is inside element

          } else weight[itype][iintg] = w_inner*wgauss[nintgx-2][i]*wgauss[nintgy-2][j]*wgauss[nintgz-2][k];

          // list to convert integration point index to 
          // interpolated atom index

          iintpl = px*ncelly*ncellz + py*ncellz + pz;
          i2ia[itype][iintg] = iintpl;
          ia2i[itype][iintpl] = iintg;

          for (int inode = 0; inode < 8; inode++) {

            // weighted shape array to distribute force from integration points to nodes

            weighted_shape_array[itype][iintg][inode] = 
              shape_array[itype][iintpl][inode] * weight[itype][iintg] * 8.0 / nintpl[itype];

            // list to convert node index to integration point index

            if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
              if (n2i[itype][inode] < 0) {
                n2i[itype][inode] = iintg;
                i2n[itype][iintg] = inode;
              } else error->all(FLERR,"Node to integration point mapping is not correct");
            }
          }
          iintg++;
        }
      }
    }

  } else if (element_shape_ids[itype] == OCTAHEDRON) {

    // calculate weight fraction for integration points
    // for element size <= 4 unit cells, all atoms are integration points

    if (ncellx <= 4) {
      for (iintg = 0; iintg < nintpl[itype]; iintg++) {

        // interpolated atom index
        iintpl = iintg; 
        i2ia[itype][iintg] = iintg;
        ia2i[itype][iintg] = iintg;
        weight[itype][iintg] = 1.0;
        for (int inode = 0; inode < 6; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = shape_array[itype][iintpl][inode] * 6.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    } else {
      int xynum;
      double dxy_pc = 2.0/(ncellx-1);
      double dz_pc = 1.0/(ncellx-1);
      double w_inner = static_cast<double> (sum_octahedron(ncellx - 1))/6.0;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iintpl = 0;

      int nsur = 0;
      int nedge = 0;
      int ninner = 0;
      int nnode = 0;

      // position of integration points on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0*(ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0*(ncellx - 2) + 1.0);

      // position of integration points on surface

      int isur = (int) (1.0/6.0*(ncellx - 4.0) + 2);

      for (int k = 1 - ncellx; k <= ncellx - 1; k++) {
        xynum = ncellx - abs(k);            // size of current xy plane
        for (int i = 0; i < xynum; i++) {
          for (int j = 0; j < xynum; j++) {

            // nodal integration point

            if (k == 1 - ncellx || k == ncellx - 1 || 
                (k == 0 && (i == 0 || i == xynum - 1) && (j == 0 || j == xynum -1))) {
              weight[itype][iintg] = 1.0;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
              nnode++;
            }

            // edge integration point 

            if (k != 0 && k != 1 - ncellx && k != ncellx - 1) {     // upper and lower half edges
              if ((abs(k) == i1 || abs(k) == i2) &&
                  (i == 0 || i == xynum - 1) &&
                  (j == 0 || j == xynum - 1)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
                nedge++;
              }
            } else if (k == 0) {                                    // middle plane edges
              if ((i == 0 || i == xynum - 1) && (j == i1 || j == i2) ||
                  (j == 0 || j == xynum - 1) && (i == i1 || i == i2)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
                //nedge++;
              }
            }

            // surface integration point

            if (abs(k) == ncellx - 1 - isur*2) {
              if (((i == 0 || i == xynum - 1) && j == (xynum-1)/2) ||
                  ((j == 0 || j == xynum - 1) && i == (xynum-1)/2)) {
                weight[itype][iintg] = w_face;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
                nsur++;
              }
            } else if (abs(k) == isur) {
              if (((i == 0 || i == xynum - 1) && (j == isur || j == xynum - 1 - isur)) ||
                  ((j == 0 || j == xynum - 1) && (i == isur || i == xynum - 1 - isur))) {
                weight[itype][iintg] = w_face;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
                nsur++;
              }
            }

            // inner integration points

            if ((k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) ||                 // on either inner half
                (abs(k) == ncellx - 1 - i1*2 && i == (xynum-1)/2 && j == (xynum-1)/2)) {    // on middle plane
              weight[itype][iintg] = w_inner;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
              ninner++;
            }


            iintpl++;
          }
        }
      }
    }
    if (iintg != nintg[itype]) error->all(FLERR,"Wrong integration setup");

    for (iintg = 0; iintg < nintg[itype]; iintg++) {
      iintpl == i2ia[itype][iintg];
      for (int inode = 0; inode < 6; inode++) {

        // weighted shape array to distribute force from integration points to nodes

        weighted_shape_array[itype][iintg][inode] = 
          shape_array[itype][iintpl][inode] * weight[itype][iintg] * 6.0 / nintpl[itype];

        // list to convert node index to integration point index

        if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
          if (n2i[itype][inode] < 0) {
            n2i[itype][inode] = iintg;
            i2n[itype][iintg] = inode;
          } else error->all(FLERR,"Node to integration point mapping is not correct");
        }
      }
    }

  } else if (element_shape_ids[itype] == PYRAMID) {

    // for element size <= 4 unit cells, all atoms are integration points

    if (ncellx <= 4) {
      for (iintg = 0; iintg < nintpl[itype]; iintg++) {

        // interpolated atom index
        iintpl = iintg; 
        i2ia[itype][iintg] = iintg;
        ia2i[itype][iintg] = iintg;
        weight[itype][iintg] = 1.0;
        for (int inode = 0; inode < 5; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = shape_array[itype][iintpl][inode] * 5.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    } else {

      double px,py,pz;
      double xycorner;
      int xynum;
      double dxy_pc = 2.0/(ncellx-1);
      double dz_pc = 1.0/(ncellx-1);
      double w_inner_total = static_cast<double> (sum_pyramid(ncellx - 3));
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face_tri = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      double w_face_quad = w_edge*w_edge;
      iintpl = 0;

      // position of integration points on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0*(ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0*(ncellx - 2) + 1.0);

      // position of integration points on surface

      int isur = (int) (1.0/6.0*(ncellx - 4.0) + 2.0);


      // for inner integration points, locate the closest atoms
      // parameters from Liping Liu - Numerical Integration over Pyramids
      // map the coordinates to the inner pyramid

      // position parameters of lower 4 inner integration points

      double z1 = (35.0 - sqrt(140.0))/140.0;
      z1 = (z1*(ncellx - 4.0) + 1.0)/(ncellx - 1.0);
      double a = sqrt(5.0/21.0);
      a = a*(ncellx - 4.0)/(ncellx - 1.0);
      double w1 = 0.21;

      // position parameter of upper inner integration points

      double z0 = (70.0 + 21.0*sqrt(35.0))/280.0;
      z0 = (z0*(ncellx - 4.0) + 1.0)/(ncellx - 1.0);
      double w0 = 0.16;

      // list of minimum distance and indexes of closest atoms for 5 inner integration points

      double mindistsq[5];
      int minintpl[5];
      mindistsq[0] = mindistsq[1] = mindistsq[2] = mindistsq[3] = mindistsq[4] = BIG; 
      double rsq;

      int nnode = 0;
      int nedge = 0;
      int nsur = 0;

      for (int k = 0; k < ncellx; k++) {
        xynum = ncellx - k;            // size of current xy plane
        pz = k*dz_pc;
        xycorner = -(xynum - 1)*dz_pc;
        for (int i = 0; i < xynum; i++) {
          px = dxy_pc*i + xycorner;
          for (int j = 0; j < xynum; j++) {
            py = dxy_pc*j + xycorner;

            // nodal integration point

            if (k == ncellx - 1 || 
                (k == 0 && (i == 0 || i == xynum - 1) && (j == 0 || j == xynum -1))) {
              weight[itype][iintg] = 1.0;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              //              fprintf(screen,"px py pz = %g %g %g iintg = %d iintpl = %d\n",px,py,pz,iintg,iintpl);
              iintg++;
            }

            // edge integration point 

            if (k != 0 && k != ncellx - 1) {     // upper edges
              if ((k == i1 || k == i2) &&
                  (i == 0 || i == xynum - 1) &&
                  (j == 0 || j == xynum - 1)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == 0) {                 // base plane edges
              if ((i == 0 || i == xynum - 1) && (j == i1 || j == i2) ||
                  (j == 0 || j == xynum - 1) && (i == i1 || i == i2)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            }

            // surface integration point

            if (k == ncellx - 1 - isur*2) {
              if (((i == 0 || i == xynum - 1) && j == (xynum-1)/2) ||
                  ((j == 0 || j == xynum - 1) && i == (xynum-1)/2)) {
                weight[itype][iintg] = w_face_tri;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == isur) {
              if (((i == 0 || i == xynum - 1) && (j == isur || j == xynum - 1 - isur)) ||
                  ((j == 0 || j == xynum - 1) && (i == isur || i == xynum - 1 - isur))) {
                weight[itype][iintg] = w_face_tri;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) {
              weight[itype][iintg] = w_face_quad;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            }


            // inner integration points

            // upper inner integration point

            rsq = px*px + py*py + (pz - z0)*(pz - z0);
            if (rsq < mindistsq[0]) {
              mindistsq[0] = rsq;
              minintpl[0] = iintpl; 
            }

            // lower 4 inner integration points

            rsq = (px + a)*(px + a) + (py + a)*(py + a) + (pz - z1)*(pz - z1);    // (-a,-a,z1)
            if (rsq < mindistsq[1]) {
              mindistsq[1] = rsq;
              minintpl[1] = iintpl; 
            }
            rsq = (px - a)*(px - a) + (py + a)*(py + a) + (pz - z1)*(pz - z1);    // (+a,-a,z1)
            if (rsq < mindistsq[2]) {
              mindistsq[2] = rsq;
              minintpl[2] = iintpl; 
            }
            rsq = (px - a)*(px - a) + (py - a)*(py - a) + (pz - z1)*(pz - z1);    // (+a,+a,z1)
            if (rsq < mindistsq[3]) {
              mindistsq[3] = rsq;
              minintpl[3] = iintpl; 
            }
            rsq = (px + a)*(px + a) + (py - a)*(py - a) + (pz - z1)*(pz - z1);    // (-a,+a,z1)
            if (rsq < mindistsq[4]) {
              mindistsq[4] = rsq;
              minintpl[4] = iintpl; 
            }

            iintpl++;
          }
        }
      }

      // save indices and weights of 5 inner integration points

      iintpl = minintpl[0];
      weight[itype][iintg] = w_inner_total*w0;
      i2ia[itype][iintg] = iintpl;
      ia2i[itype][iintpl] = iintg;
      iintg++;

      for (int i = 1; i < 5; i++) {
        iintpl = minintpl[i];
        weight[itype][iintg] = w_inner_total*w1;
        i2ia[itype][iintg] = iintpl;
        ia2i[itype][iintpl] = iintg;
        iintg++;
      }

      if (iintg != nintg[itype]) error->all(FLERR,"Wrong integration setup");

      // setup weighted shape array and node to integration point mapping

      for (iintg = 0; iintg < nintg[itype]; iintg++) {
        iintpl = i2ia[itype][iintg];
        for (int inode = 0; inode < 5; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = 
            shape_array[itype][iintpl][inode] * weight[itype][iintg] * 5.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            //            fprintf(screen,"n2i = %d inode = %d iintg = %d shape_array = %g\n",n2i[itype][inode],inode,iintg,shape_array[itype][iintpl][inode]);
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    }

  } else if (element_shape_ids[itype] == TRIANGLE) {

    // for element size <= 4 unit cells, all atoms are integration points

    if (ncellx <= 4) {
      for (iintg = 0; iintg < nintpl[itype]; iintg++) {

        // interpolated atom index
        iintpl = iintg; 
        i2ia[itype][iintg] = iintg;
        ia2i[itype][iintg] = iintg;
        weight[itype][iintg] = 1.0;
        for (int inode = 0; inode < 3; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = shape_array[itype][iintpl][inode] * 3.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    } else {

      int ynum;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iintpl = 0;

      // position of integration points on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0*(ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0*(ncellx - 2) + 1.0);

      // position of integration points on surface

      int isur = (int) (1.0/6.0*(ncellx - 4.0) + 1.0);

      int nnode = 0;
      int nedge = 0;
      int nsur = 0;

      for (int i = 0; i < ncellx; i++) {
        ynum = ncellx - i;
        for (int j = 0; j < ynum; j++) {

          // nodal integration points


          if (i == ncellx - 1 ||
              j == ncellx - 1 ||
              (i == 0 && j == 0)) {
            weight[itype][iintg] = 1.0;
            i2ia[itype][iintg] = iintpl;
            ia2i[itype][iintpl] = iintg;
            iintg++;
          }

          // edge integration points

          if (i == 0 && (j == i1 || j == i2) ||
              (j == 0 || j == ynum - 1) && (i == i1 || i == i2)) {
            weight[itype][iintg] = w_edge;
            i2ia[itype][iintg] = iintpl;
            ia2i[itype][iintpl] = iintg;
            iintg++;
          }

          // surface integration point

          if ((i == isur || i == i2) && (j == i1 || j == i2)) {
            weight[itype][iintg] = w_face;
            i2ia[itype][iintg] = iintpl;
            ia2i[itype][iintpl] = iintg;
            iintg++;
          }

          iintpl++;
        }
      }
    }

    if (iintg != nintg[itype]) error->all(FLERR,"Wrong integration setup");

    // setup weighted shape array and node to integration point mapping

    for (iintg = 0; iintg < nintg[itype]; iintg++) {
      iintpl = i2ia[itype][iintg];
      for (int inode = 0; inode < 3; inode++) {

        // weighted shape array to distribute force from integration points to nodes

        weighted_shape_array[itype][iintg][inode] = 
          shape_array[itype][iintpl][inode] * weight[itype][iintg] * 3.0 / nintpl[itype];

        // list to convert node index to integration point index

        if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
          //            fprintf(screen,"n2i = %d inode = %d iintg = %d shape_array = %g\n",n2i[itype][inode],inode,iintg,shape_array[itype][iintpl][inode]);
          if (n2i[itype][inode] < 0) {
            n2i[itype][inode] = iintg;
            i2n[itype][iintg] = inode;
          } else error->all(FLERR,"Node to integration point mapping is not correct");
        }

      }
    }

  } else if (element_shape_ids[itype] == TETRAHEDRON) {

    // for element size <= 4 unit cells, all atoms are integration points

    if (ncellx <= 4) {
      for (iintg = 0; iintg < nintpl[itype]; iintg++) {

        // interpolated atom index
        iintpl = iintg; 
        i2ia[itype][iintg] = iintg;
        ia2i[itype][iintg] = iintg;
        weight[itype][iintg] = 1.0;
        for (int inode = 0; inode < 4; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = shape_array[itype][iintpl][inode] * 4.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    } else {

      int xnum,ynum;
      double px,py,pz;
      double dxyz_pc = 1.0/(ncellx-1);
      double w_inner = static_cast<double> (sum_tetrahedron(ncellx - 3))/4.0;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iintpl = 0;

      // position of integration points on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0*(ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0*(ncellx - 2) + 1.0);

      // position of integration points on surface

      int isur = (int) (1.0/6.0*(ncellx - 4.0) + 1.0);

      int nnode = 0;
      int nedge = 0;
      int nsur = 0;

      // list of minimum distance and indexes of closest atoms for 4 inner integration points

      double mindistsq[4];
      int minintpl[4];
      mindistsq[0] = mindistsq[1] = mindistsq[2] = mindistsq[3] = BIG; 
      double rsq;

      // coordinates of inner integration points
      // map coordinates to the inner tetrahedron domain

      double a = 0.1381966011250110;
      a = (a*(ncellx - 5.0) + 1.0)/(ncellx - 1.0);
      double b = 0.5854101966249680;
      b = (b*(ncellx - 5.0) + 1.0)/(ncellx - 1.0);

      for (int k = 0; k < ncellx; k++) {
        xnum = ncellx - k;            // size of current xy plane
        pz = k*dxyz_pc;
        for (int i = 0; i < xnum; i++) {
          ynum = xnum - i;
          px = dxyz_pc*i;
          for (int j = 0; j < ynum; j++) {
            py = dxyz_pc*j;

            // nodal integration points

            if (k == ncellx - 1 || 
                i == ncellx - 1 ||
                j == ncellx - 1 ||
                (k == 0 && i == 0 && j == 0)) {
              weight[itype][iintg] = 1.0;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            }

            // edge integration points

            if (k != 0 && k != ncellx - 1) {     // upper edges
              if ((k == i1 || k == i2) &&
                  (i == 0 || i == xnum - 1) &&
                  (j == 0 || j == ynum - 1)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == 0) {                 // base plane edges
              if (i == 0 && (j == i1 || j == i2) ||
                  (j == 0 || j == ynum - 1) && (i == i1 || i == i2)) {
                weight[itype][iintg] = w_edge;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            }

            // surface integration point

            if (k == ncellx - 1 - isur*2) {
              if ((i == 0 && j == (xnum-1)/2) ||
                  (j == 0 && i == (ynum-1)/2) ||
                  (i == ynum-1)) {
                weight[itype][iintg] = w_face;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == isur) {
              if ((i == 0 && (j == isur || j == ynum - 1 - isur)) ||
                  ((j == 0 || j == ynum - 1) && (i == isur || i == xnum - 1 - isur))) {
                weight[itype][iintg] = w_face;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) {
              weight[itype][iintg] = w_face;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            }


            // inner integration points

            rsq = (px - a)*(px - a) + (py - a)*(py - a) + (pz - a)*(pz - a);    // (a,a,a)
            if (rsq < mindistsq[0]) {
              mindistsq[0] = rsq;
              minintpl[0] = iintpl; 
            }
            rsq = (px - b)*(px - b) + (py - a)*(py - a) + (pz - a)*(pz - a);    // (b,a,a)
            if (rsq < mindistsq[1]) {
              mindistsq[1] = rsq;
              minintpl[1] = iintpl; 
            }
            rsq = (px - a)*(px - a) + (py - b)*(py - b) + (pz - a)*(pz - a);    // (a,b,a)
            if (rsq < mindistsq[2]) {
              mindistsq[2] = rsq;
              minintpl[2] = iintpl; 
            }
            rsq = (px - a)*(px - a) + (py - a)*(py - a) + (pz - b)*(pz - b);    // (a,a,b)
            if (rsq < mindistsq[3]) {
              mindistsq[3] = rsq;
              minintpl[3] = iintpl; 
            }

            iintpl++;
          }
        }
      }

      // save indices and weights of 4 inner integration points

      for (int i = 0; i < 4; i++) {
        iintpl = minintpl[i];
        weight[itype][iintg] = w_inner;
        i2ia[itype][iintg] = iintpl;
        ia2i[itype][iintpl] = iintg;
        iintg++;
      }

      if (iintg != nintg[itype]) error->all(FLERR,"Wrong integration setup");

      // setup weighted shape array and node to integration point mapping

      for (iintg = 0; iintg < nintg[itype]; iintg++) {
        iintpl = i2ia[itype][iintg];
        for (int inode = 0; inode < 4; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = 
            shape_array[itype][iintpl][inode] * weight[itype][iintg] * 4.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            //            fprintf(screen,"n2i = %d inode = %d iintg = %d shape_array = %g\n",n2i[itype][inode],inode,iintg,shape_array[itype][iintpl][inode]);
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    }

  } else if (element_shape_ids[itype] == WEDGE) {

    // position of integration points on edge

    int i1x = (int) ((1.0 + xgauss[2][1])/2.0*(ncellx - 2) + 1.0);
    int i2x = (int) ((1.0 + xgauss[2][2])/2.0*(ncellx - 2) + 1.0);
    int i1y = (int) ((1.0 + xgauss[2][1])/2.0*(ncelly - 2) + 1.0);
    int i2y = (int) ((1.0 + xgauss[2][2])/2.0*(ncelly - 2) + 1.0);

    // for element size in yz <= 4 unit cells, all atoms on yz plane are integration points

    if (ncellx <= 4) {
      int znum;
      double dx_pc = 2.0/(ncellx-1);
      double dyz_pc = 1.0/(ncelly-1);
      double w_edge_x = (ncellx-2.0)/(nintgx-2.0);

      iintpl = 0;
      for (int i = 0; i < ncellx; i++) {
        for (int j = 0; j < ncelly; j++) {
          znum = ncelly - j;
          for (int k = 0; k < znum; k++) {
            if (i == 0 || i == ncellx-1) {
              weight[itype][iintg] = 1;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            } else if (i == i1x || i == i2x) {
              weight[itype][iintg] = w_edge_x;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            }
            iintpl++;
          }
        }
      }

      // setup weighted shape array and node to integration point mapping

      for (iintg = 0; iintg < nintg[itype]; iintg++) {
        iintpl = i2ia[itype][iintg];
        for (int inode = 0; inode < 6; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = 
            shape_array[itype][iintpl][inode] * weight[itype][iintg] * 6.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    } else {

      int znum;
      double px,py,pz;
      double dxyz_pc = 1.0/(ncellx-1);
      double w_inner = ((ncellx - 2.0)*(ncelly - 3.0)*(ncelly - 2.0))/12.0;
      double w_edge_x = (ncellx - 2.0)/2.0;
      double w_edge_yz = (ncelly - 2.0)/2.0;
      double w_face_tri = static_cast<double> (sum_triangle(ncelly - 3))/3.0;
      double w_face_quad = w_edge_x*w_edge_yz;

      iintpl = 0;

      // position of integration points on surface

      int isur = (int) (1.0/6.0*(ncelly - 4.0) + 1.0);

      for (int i = 0; i < ncellx; i++) {
        for (int j = 0; j < ncelly; j++) {
          znum = ncelly - j;
          for (int k = 0; k < znum; k++) {

            // nodal integration points

            if ((i == ncellx - 1 || i == 0) &&
                (j == ncelly - 1 || j == 0) &&
                (k == ncelly - 1 || k == 0)) {
              weight[itype][iintg] = 1.0;
              i2ia[itype][iintg] = iintpl;
              ia2i[itype][iintpl] = iintg;
              iintg++;
            }

            // edge integration points

            if (i == i1x || i == i2x) {                   // x edge
              if ((j == 0 && k == 0) ||
                  j == ncelly-1 || k == ncelly-1) {
                weight[itype][iintg] = w_edge_x;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;

              }
            } else if (i == 0 || i == ncellx-1) {         // y,z edges
              if ((j == 0 && (k == i1y || k == i2y)) ||
                  ((j == i1y || j == i2y) && (k == 0 || k == znum-1))) {
                weight[itype][iintg] = w_edge_yz;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            }

            // surface integration point

            if (i == ncellx - 1 || i == 0) {             // x = +-1 planes
              if ((j == isur && k == isur) ||
                  (j == isur && k == znum-1-isur) ||
                  (j == ncelly-1-2*isur && k == znum/2)) {
                weight[itype][iintg] = w_face_tri;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (k == 0 || k == znum-1) {          // xy, y + z = 1 planes
              if ((i == i1x || i == i2x) && (j == i1y || j == i2y)) { 
                weight[itype][iintg] = w_face_quad;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            } else if (j == 0) {          // xz plane
              if ((i == i1x || i == i2x) && (k == i1y || k == i2y)) { 
                weight[itype][iintg] = w_face_quad;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            }

            // inner integration points

            if (i == i1x || i == i2x) {      
              if ((j == isur && k == isur) ||
                  (j == isur && k == znum-1-isur) ||
                  (j == ncelly-1-2*isur && k == znum/2)) {
                weight[itype][iintg] = w_inner;
                i2ia[itype][iintg] = iintpl;
                ia2i[itype][iintpl] = iintg;
                iintg++;
              }
            }

            iintpl++;
          }
        }
      }

      if (iintg != nintg[itype]) error->all(FLERR,"Wrong integration setup");

      // setup weighted shape array and node to integration point mapping

      for (iintg = 0; iintg < nintg[itype]; iintg++) {
        iintpl = i2ia[itype][iintg];
        for (int inode = 0; inode < 6; inode++) {

          // weighted shape array to distribute force from integration points to nodes

          weighted_shape_array[itype][iintg][inode] = 
            shape_array[itype][iintpl][inode] * weight[itype][iintg] * 6.0 / nintpl[itype];

          // list to convert node index to integration point index

          if (fabs(shape_array[itype][iintpl][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        }
      }
    }
  } else {

    // add new integration scheme for new element shape here
  }


  // check if all node to integration point mapping is set

  for (int inode = 0; inode < npe[itype]; inode++) {
    if (n2i[itype][inode] < 0) 
      error->all(FLERR,"Node to integration point mapping is not complete");
    if (i2n[itype][n2i[itype][inode]] < 0)
      error->all(FLERR,"Integration point to node mapping is not complete");
    n2ia[itype][inode] = i2ia[itype][n2i[itype][inode]];
  }

  // check if weight is assigned correctly

  double total_weight = 0.0;
  for (int l = 0; l < nintg[itype]; l++) total_weight += weight[itype][l];
  //fprintf(screen,"total_weight = %g nintpl = %d\n",total_weight,nintpl[itype]);
  if (fabs(total_weight - nintpl[itype]) > EPSILON) error->all(FLERR,"Integraion weight is not assigned correctly");
}

/*-----------------------------------------------------------------------------------
  Calculate sub-element related arrays for an element type
  - ias2ia to map atom index inside sub-element to atom index inside parent element
  - natom_subelem: number of atoms inside each sub element
  ------------------------------------------------------------------------------------*/

void ElementVecCAC::setup_sub_element(int itype)
{
  int m,n,ntotal;
  int isub = 0;
  int ncellx,ncelly,ncellz;
  int min_split_size = element->min_split_size;

  ncellx = ncells[itype][0];
  ncelly = ncells[itype][1];
  ncellz = ncells[itype][2];

  nsubelem[itype] = 0;

  // total interpolated atoms counter

  ntotal = 0;

  // makes sure atoms on element boundaries are included
  // atoms might sit right on sub-element boundaries
  // so subtract a small number from the boundaries to make sure
  // all atoms are assigned to only one sub-element

  if (element_shape_ids[itype] == QUADRILATERAL) {

    double px,py;
    double xlo,xhi,ylo,yhi;
    int nsplitx = ncellx/min_split_size;
    int nsplity = ncelly/min_split_size;
    int esplit = MAX(nsplitx,nsplity);
    subsplit[itype] = esplit;
    double ds = 2.0/esplit;
    double dx_pc = 2.0/(ncellx-1);
    double dy_pc = 2.0/(ncelly-1);

    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix*ds - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds;

      for (int iy = 0; iy < esplit; iy++) {

        ylo = -1.0 + iy*ds - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds;

        // calculate sub-element center shape function array

        px = (ix + 0.5) * ds - 1.0;
        py = (iy + 0.5) * ds - 1.0;

        shape_array_center_subelem[itype][isub][0] = 0.25*(1 - px)*(1 - py);
        shape_array_center_subelem[itype][isub][1] = 0.25*(1 + px)*(1 - py);
        shape_array_center_subelem[itype][isub][2] = 0.25*(1 + px)*(1 + py);
        shape_array_center_subelem[itype][isub][3] = 0.25*(1 - px)*(1 + py);

        // place interpolated atoms into sub element isub

        n = 0; // number of atoms inside a sub element counter
        m = 0; // index of interpolated atom in element

        for (int i = 0; i < ncellx; i++) {
          px = dx_pc*i - 1.0;
          for (int j = 0; j < ncelly; j++) {
            py = dy_pc*j - 1.0;
            if (px >= xlo && px < xhi &&
                py >= ylo && py < yhi) {
              ias2ia[itype][isub][n++] = m;
              if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
            }
            m++;
          }
        }

        // don't store info if sub-element contains no atom

        if (n == 0) continue;

        nsubelem[itype]++;
        if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");
        natom_subelem[itype][isub++] = n;
        ntotal += n;
      }
    }

  } else if (element_shape_ids[itype] == HEXAHEDRON) {

    double px,py,pz;
    double dx_pc,dy_pc,dz_pc;
    double xlo,xhi,ylo,yhi,zlo,zhi;
    int nsplitx = ncellx/min_split_size;
    int nsplity = ncelly/min_split_size;
    int nsplitz = ncellz/min_split_size;
    int esplit = MAX(nsplitx,nsplity);
    esplit = MAX(esplit,nsplitz);
    subsplit[itype] = esplit;
    if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");
    double ds = 2.0/esplit;
    dx_pc = 2.0/(ncellx-1);
    dy_pc = 2.0/(ncelly-1);
    dz_pc = 2.0/(ncellz-1);
    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix*ds - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = -1.0 + iy*ds - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = -1.0 + iz*ds - EPSILON;
          if (iz == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds;

          // calculate sub-element center shape function array

          px = (ix + 0.5) * ds - 1.0;
          py = (iy + 0.5) * ds - 1.0;
          pz = (iz + 0.5) * ds - 1.0;

          shape_array_center_subelem[itype][isub][0] = 0.125*(1 - px)*(1 - py)*(1 - pz);
          shape_array_center_subelem[itype][isub][1] = 0.125*(1 + px)*(1 - py)*(1 - pz);
          shape_array_center_subelem[itype][isub][2] = 0.125*(1 + px)*(1 + py)*(1 - pz);
          shape_array_center_subelem[itype][isub][3] = 0.125*(1 - px)*(1 + py)*(1 - pz);
          shape_array_center_subelem[itype][isub][4] = 0.125*(1 - px)*(1 - py)*(1 + pz);
          shape_array_center_subelem[itype][isub][5] = 0.125*(1 + px)*(1 - py)*(1 + pz);
          shape_array_center_subelem[itype][isub][6] = 0.125*(1 + px)*(1 + py)*(1 + pz);
          shape_array_center_subelem[itype][isub][7] = 0.125*(1 - px)*(1 + py)*(1 + pz);

          // place interpolated atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of interpolated atom in element

          for (int i = 0; i < ncellx; i++) {
            px = dx_pc*i - 1.0;
            for (int j = 0; j < ncelly; j++) {
              py = dy_pc*j - 1.0;
              for (int k = 0; k < ncellz; k++) {
                pz = dz_pc*k - 1.0;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  ias2ia[itype][isub][n++] = m;
                  if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          nsubelem[itype]++;
          if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");

          natom_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }

  } else if (element_shape_ids[itype] == PYRAMID) {
    double xlo,xhi,ylo,yhi,zlo,zhi;
    double px,py,pz;
    double xycorner;
    int xynum;
    int esplit = ncellx/min_split_size;
    subsplit[itype] = esplit;
    double ds_xy = 2.0/esplit;
    double ds_z = 1.0/esplit;
    double dxy_pc = 2.0/(ncellx-1);
    double dz_pc = 1.0/(ncellx-1);

    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix*ds_xy - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds_xy;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = -1.0 + iy*ds_xy - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds_xy;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = iz*ds_z - EPSILON;
          if (iz == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds_z;

          // calculate sub-element center shape function array

          px = (ix + 0.5) * ds_xy - 1.0;
          py = (iy + 0.5) * ds_xy - 1.0;
          pz = (iz + 0.5) * ds_z;

          shape_array_center_subelem[itype][isub][0] = 0.25*(1.0 - px - pz)*(1.0 - py - pz)/(1 - pz);
          shape_array_center_subelem[itype][isub][1] = 0.25*(1.0 + px - pz)*(1.0 - py - pz)/(1 - pz);
          shape_array_center_subelem[itype][isub][2] = 0.25*(1.0 + px - pz)*(1.0 + py - pz)/(1 - pz);
          shape_array_center_subelem[itype][isub][3] = 0.25*(1.0 - px - pz)*(1.0 + py - pz)/(1 - pz);
          shape_array_center_subelem[itype][isub][4] = pz;

          // place interpolated atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of interpolated atom in element

          for (int k = 0; k <= ncellx - 1; k++) {
            pz = dz_pc*k;
            xynum = ncellx - k;            // size of current xy plane
            xycorner = -(xynum - 1)*dz_pc;
            for (int i = 0; i < xynum; i++) {
              px = dxy_pc*i + xycorner;
              for (int j = 0; j < xynum; j++) {
                py = dxy_pc*j + xycorner;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  ias2ia[itype][isub][n++] = m;
                  if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom


          if (n == 0) continue;

          nsubelem[itype]++;
          if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");

          natom_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }
  } else if (element_shape_ids[itype] == TETRAHEDRON) {
    double xlo,xhi,ylo,yhi,zlo,zhi;
    double px,py,pz;
    int xnum,ynum;
    int esplit = ncellx/min_split_size;
    subsplit[itype] = esplit;
    double ds_xyz = 1.0/esplit;
    double dxyz_pc = 1.0/(ncellx-1);

    for (int ix = 0; ix < esplit; ix++) {
      xlo = ix*ds_xyz - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds_xyz;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = iy*ds_xyz - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds_xyz;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = iz*ds_xyz - EPSILON;
          if (iz == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds_xyz;

          // calculate sub-element center shape function array

          px = (ix + 0.5) * ds_xyz;
          py = (iy + 0.5) * ds_xyz;
          pz = (iz + 0.5) * ds_xyz;

          shape_array_center_subelem[itype][isub][0] = px;
          shape_array_center_subelem[itype][isub][1] = py;
          shape_array_center_subelem[itype][isub][2] = pz;
          shape_array_center_subelem[itype][isub][3] = 1.0 - px - py - pz;

          // place interpolated atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of interpolated atom in element

          for (int k = 0; k < ncellx; k++) {
            pz = dxyz_pc*k;
            xnum = ncellx - k;
            for (int i = 0; i < xnum; i++) {
              px = dxyz_pc*i;
              ynum = xnum - i;
              for (int j = 0; j < ynum; j++) {
                py = dxyz_pc*j;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  ias2ia[itype][isub][n++] = m;
                  if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom


          if (n == 0) continue;

          nsubelem[itype]++;
          if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");

          natom_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }

  } else if (element_shape_ids[itype] == WEDGE) {

    double px,py,pz;
    int znum;
    double xlo,xhi,ylo,yhi,zlo,zhi;
    int nsplitx = ncellx/min_split_size;
    int nsplity = ncelly/min_split_size;
    int esplit = MAX(nsplitx,nsplity);
    subsplit[itype] = esplit;
    if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");
    double ds_x = 2.0/esplit;
    double ds_yz = 1.0/esplit;
    double dx_pc = 2.0/(ncellx-1);
    double dyz_pc = 1.0/(ncelly-1);

    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix*ds_x - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds_x;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = iy*ds_yz - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds_yz;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = iz*ds_yz - EPSILON;
          if (iz == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds_yz;

          // calculate sub-element center shape function array

          px = (ix + 0.5) * ds_x - 1.0;
          py = (iy + 0.5) * ds_yz;
          pz = (iz + 0.5) * ds_yz;

          shape_array_center_subelem[itype][isub][0] = 0.5*(1 - px)*(1 - py - pz);
          shape_array_center_subelem[itype][isub][1] = 0.5*(1 - px)*py;
          shape_array_center_subelem[itype][isub][2] = 0.5*(1 - px)*pz;
          shape_array_center_subelem[itype][isub][3] = 0.5*(1 + px)*(1 - py - pz);
          shape_array_center_subelem[itype][isub][4] = 0.5*(1 + px)*py;
          shape_array_center_subelem[itype][isub][5] = 0.5*(1 + px)*pz;

          // place interpolated atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of interpolated atom in element

          for (int i = 0; i < ncellx; i++) {
            px = dx_pc*i - 1.0;
            for (int j = 0; j < ncelly; j++) {
              py = dyz_pc*j;
              znum = ncelly - j;
              for (int k = 0; k < znum; k++) {
                pz = dyz_pc*k;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  ias2ia[itype][isub][n++] = m;
                  if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          nsubelem[itype]++;
          if (nsubelem[itype] > element->maxsubelem) error->all(FLERR,"Too many sub-elements, boost maxsubelem or min_split_size through element_modify command");

          natom_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }

  } else if (0) {
    // add new split scheme for new element shape here

  } else {

    // for element shape not yet having a split scheme

    subsplit[itype] = 1;
    nsubelem[itype] = 1;
    natom_subelem[itype][0] = nintpl[itype];
    if (nintpl[itype] > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
    for (int i = 0; i < nintpl[itype]; i++)
      ias2ia[itype][0][i] = i;

  }

  if (ntotal != nintpl[itype])
    error->all(FLERR,"Total atom count in sub-elements do not match");
}

/* -----------------------------------
   update all local element center coordinates from node coordinates
   to be called by Fix after integrating node positions
   ------------------------------------*/

void ElementVecCAC::update_center_coord()
{

  // zero center coords of elements

  size_t nbytes = sizeof(double) * element->nlocal;
  if (nbytes) {
    memset(&x[0][0],0,3*nbytes);
  }

  // update center coords of elements

  for (int i = 0; i < element->nlocal; i++) {
    int inpe = npe[etype[i]];
    for (int j = 0; j < inpe; j++) {
      x[i][0] += nodex[i][j][0];
      x[i][1] += nodex[i][j][1];
      x[i][2] += nodex[i][j][2];
    }
    x[i][0] = x[i][0]/inpe;
    x[i][1] = x[i][1]/inpe;
    x[i][2] = x[i][2]/inpe;
  }
}

/* -----------------------------------
   update local element I node idim coordinate 
   from center idim coordinate
   ------------------------------------*/

void ElementVecCAC::update_node_coord(int i, int idim)
{
  double xtmp = 0.0;
  int inpe = npe[etype[i]];
  for (int j = 0; j < inpe; j++)
    xtmp += nodex[i][j][idim];
  xtmp = x[i][idim] - xtmp/inpe;
  for (int j = 0; j < inpe; j++) 
    nodex[i][j][idim] += xtmp;
}

/* -----------------------------------
   update all local element node coordinates 
   from center coordinates
   ------------------------------------*/

void ElementVecCAC::update_node_coord()
{
  double xtmp,ytmp,ztmp;
  for (int i = 0; i < element->nlocal; i++) {
    int inpe = npe[etype[i]];
    xtmp = ytmp = ztmp = 0.0;
    for (int j = 0; j < inpe; j++) {
      xtmp += nodex[i][j][0];
      ytmp += nodex[i][j][1];
      ztmp += nodex[i][j][2];
    }
    xtmp = x[i][0] - xtmp/inpe;
    ytmp = x[i][1] - ytmp/inpe;
    ztmp = x[i][2] - ztmp/inpe;
    for (int j = 0; j < inpe; j++) {
      nodex[i][j][0] += xtmp;
      nodex[i][j][1] += ytmp;
      nodex[i][j][2] += ztmp;
    }
  }
}

/* --------------------------------------------------------------
   discritize element I into discrete atoms
   to be called by disc_element command
   return number of atoms created
   ---------------------------------------------------------------*/

int ElementVecCAC::element2atom(int i) 
{

  // get the coordinates of interpolated atom k
  // move last element to element i

  int ietype = etype[i];
  int ictype = ctype[i];
  double coord[3];
  for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    for (int inode = 0; inode < npe[ietype]; inode++) {
      coord[0] += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
      coord[1] += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
      coord[2] += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
    }
    atom->avec->create_atom(coord,ictype,0); 
  }

  copy(--element->nlocal,i,1);
  return nintpl[ietype];
}

/* --------------------------------------------------------------
   discritize element I into smaller elements with new etype and same ctype
   to be called by disc_element command
   return number of elements added
 ****NEED TO BE FIXED FOR NEW DATA STRUCTURE****
 ---------------------------------------------------------------*/

int ElementVecCAC::element2element(int i, int ietype_new, tagint id_offset) 
{

  // get the coordinates of interpolated atom k
  // move last element to element i

  int ietype_old = etype[i];
  int ictype = ctype[i];
  if (apc[ietype_old] != apc[ietype_new]) 
    error->all(FLERR,"New and old etype must have the same apc");

  // number of new elements along each direction

  int ne[3];
  tagint itag = tag[i];
  tagint mytag;
  double *coord = new double[3];            // coords of center
  int iintpl;

  int n = 0;

  if (element_shape_ids[ietype_old] == HEXAHEDRON) {

    double **nodecoord;
    memory->create(nodecoord,8,3,"evec:nodecoord");
    int d[8][3] = {
      {0,0,0},
      {1,0,0},
      {1,1,0},
      {0,1,0},
      {0,0,1},
      {1,0,1},
      {1,1,1},
      {0,1,1}};

    if (element_shape_ids[ietype_new] != HEXAHEDRON) 
      error->one(FLERR,"Rhombohedron elements can only be divided into rhombohedron elements");

    // check if element i can be divided evenly

    if (ncells[ietype_old][0] % ncells[ietype_new][0] ||
        ncells[ietype_old][1] % ncells[ietype_new][1] ||
        ncells[ietype_old][2] % ncells[ietype_new][2]) 
      error->one(FLERR,"Cannot divide elements evenly");

    for (int ix = 0; ix < ncells[ietype_old][0]; ix += ncells[ietype_new][0])
      for (int iy = 0; iy < ncells[ietype_old][1]; iy += ncells[ietype_new][1])
        for (int iz = 0; iz < ncells[ietype_old][2]; iz += ncells[ietype_new][2]) {

          // interpolate new node coords from element I

          coord[0] = coord[1] = coord[2] = 0.0;

          for (int new_node = 0; new_node < 8; new_node++) { 

            iintpl = (ix + d[new_node][0]*(ncells[ietype_new][0]-1))*ncells[ietype_old][1]*ncells[ietype_old][2] 
              + (iy + d[new_node][1]*(ncells[ietype_new][1]-1))*ncells[ietype_old][2] 
              + (iz + d[new_node][2]*(ncells[ietype_new][2]-1));

            nodecoord[new_node][0] = 0.0;
            nodecoord[new_node][1] = 0.0;
            nodecoord[new_node][2] = 0.0;
            for (int inode = 0; inode < 8; inode++) {
              nodecoord[new_node][0] += nodex[i][inode][0]*shape_array[ietype_old][iintpl][inode];
              nodecoord[new_node][1] += nodex[i][inode][1]*shape_array[ietype_old][iintpl][inode];
              nodecoord[new_node][2] += nodex[i][inode][2]*shape_array[ietype_old][iintpl][inode];
            }
            coord[0] += nodecoord[new_node][0];
            coord[1] += nodecoord[new_node][1];
            coord[2] += nodecoord[new_node][2];

          } 

          coord[0] /= 8.0; 
          coord[1] /= 8.0; 
          coord[2] /= 8.0; 

          // first element will have the same tag as old element
          // other elements will start from id_offset and increment by apc
          // this will ensure cluster ordering for multiple atoms per node cluster case
          // if id_offset = 0, set mytag = 0;

          if (id_offset) {
            if (n) mytag = id_offset + (n-1)*apc[ietype_old];
            else mytag = itag;
          } else mytag = 0;

          create_element(coord,nodecoord,ietype_new,ictype,mytag);
          n++;
        }

    memory->destroy(nodecoord);
  } else if (element_shape_ids[ietype_old] == QUADRILATERAL) {

    double **nodecoord;
    memory->create(nodecoord,4,3,"evec:nodecoord");
    int d[4][2] = {
      {0,0},
      {1,0},
      {1,1},
      {0,1}};

    if (element_shape_ids[ietype_new] != QUADRILATERAL)
      error->one(FLERR,"Rhombus2D elements can only be divided into rhombus2D elements");

    // check if element i can be divided evenly

    if (ncells[ietype_old][0] % ncells[ietype_new][0] ||
        ncells[ietype_old][1] % ncells[ietype_new][1])
      error->one(FLERR,"Cannot divide elements evenly");

    // interpolate new node coords from element I
    //
    for (int ix = 0; ix < ncells[ietype_old][0]; ix += ncells[ietype_new][0])
      for (int iy = 0; iy < ncells[ietype_old][1]; iy += ncells[ietype_new][1]) {

        coord[0] = coord[1] = 0.0;

        for (int new_node = 0; new_node < 4; new_node++) { 
          nodecoord[new_node][0] = 0.0;
          nodecoord[new_node][1] = 0.0;
          nodecoord[new_node][2] = 0.0;

          // interpolated index of this new node in old element

          iintpl = (ix + d[new_node][0]*(ncells[ietype_new][0]-1))*ncells[ietype_old][1]
            + (iy + d[new_node][1]*(ncells[ietype_new][1]-1));

          for (int inode = 0; inode < 4; inode++) {
            nodecoord[new_node][0] += nodex[i][inode][0]*shape_array[ietype_old][iintpl][inode];
            nodecoord[new_node][1] += nodex[i][inode][1]*shape_array[ietype_old][iintpl][inode];
          }
          coord[0] += nodecoord[new_node][0];
          coord[0] += nodecoord[new_node][1];

          nodecoord[new_node][3] = 0.0;
        }

        coord[0] /= 4.0; 
        coord[1] /= 4.0; 

        // first element will have the same tag as old element
        // other elements will start from id_offset and increment by apc
        // this will ensure cluster ordering for multiple atoms per node cluster case
        // if id_offset = 0, set mytag = 0;

        if (id_offset) {
          if (n) mytag = id_offset + (n-1)*apc[ietype_old];
          else mytag = itag;
        } else mytag = 0;

        create_element(coord,nodecoord,ietype_new,ictype,mytag);
        n++;
      }
    memory->destroy(nodecoord);

  } else {
    // add new element shape here
    error->one(FLERR,"Element split for this type has not been set yet");
  }

  // init per-element fix/compute/variable values for created elements

  element->data_fix_compute_variable(element->nlocal-n,element->nlocal);

  // remove element I by replacing it with last element in list

  copy(--element->nlocal,i,1);

  // clean up and return number of additional elements

  delete [] coord;
  return n;
}


/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint ElementVecCAC::memory_usage()
{
  bigint bytes = 0;

  int max_npe = element->max_npe;
  if (element->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (element->memcheck("ctype")) bytes += memory->usage(ctype,nmax);
  if (element->memcheck("etype")) bytes += memory->usage(etype,nmax);
  if (element->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (element->memcheck("image")) bytes += memory->usage(image,nmax);
  if (element->memcheck("nodemask")) bytes += memory->usage(nodemask,nmax,max_npe);
  if (element->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (element->memcheck("slip_plane")) bytes += memory->usage(slip_plane,nmax,3,3);
  if (element->memcheck("cell_size")) bytes += memory->usage(cell_size,nmax,3);
  if (element->memcheck("nodex")) bytes += memory->usage(nodex,nmax,max_npe,3);
  if (element->memcheck("nodev")) bytes += memory->usage(nodev,nmax,max_npe,3);
  if (element->memcheck("nodef")) bytes += memory->usage(nodef,nmax,max_npe,3);
  if (element->memcheck("nodetag")) bytes += memory->usage(nodetag,nmax,max_npe);
  if (element->memcheck("subelem_size")) bytes += memory->usage(subelem_size,nmax);
  if (element->memcheck("initial_box_size")) bytes += memory->usage(initial_box_size,nmax,3);
  if (element->memcheck("element_bound_box")) bytes += memory->usage(element_bound_box,nmax,6);
  if (element->max_apc > 1)
    if (element->memcheck("element_clusters")) bytes += memory->usage(element_clusters,nmax,element->max_apc);
  return bytes;
}

/* ----------------------------------------------------------------------
   update internal slip_plane (axes) of all elements
   slip_plane is defined as unit normal vector of the plane
   ------------------------------------------------------------------------- */

void ElementVecCAC::update_slip_plane()
{
  double xnorm,ynorm,znorm;
  double xdir[3],ydir[3],zdir[3];
  for (int i = 0; i < element->nlocal; i++) {
    int itype = etype[i];
    if (element_shape_ids[itype] == QUADRILATERAL) {
      for (int j = 0; j < 2; j++) {
        xdir[j] = 0.5*
          (nodex[i][2][j] - nodex[i][1][j] + 
           nodex[i][3][j] - nodex[i][4][j]);
        ydir[j] = 0.5*
          (nodex[i][4][j] - nodex[i][1][j] +
           nodex[i][3][j] - nodex[i][2][j]);
        zdir[j] = 0.0;
      }
      xdir[2] = 0.0;
      ydir[2] = 0.0;
      zdir[2] = 1.0;
    } else if (element_shape_ids[itype] == HEXAHEDRON) {
      for (int j = 0; j < 3; j++) {
        xdir[j] = 0.25*
          (nodex[i][1][j] - nodex[i][0][j] + 
           nodex[i][2][j] - nodex[i][3][j] +
           nodex[i][5][j] - nodex[i][4][j] +
           nodex[i][6][j] - nodex[i][7][j]);
        ydir[j] = 0.25*
          (nodex[i][3][j] - nodex[i][0][j] +
           nodex[i][2][j] - nodex[i][1][j] +
           nodex[i][7][j] - nodex[i][4][j] +
           nodex[i][6][j] - nodex[i][5][j]);
        zdir[j] = 0.25*
          (nodex[i][4][j] - nodex[i][0][j] +
           nodex[i][7][j] - nodex[i][3][j] +
           nodex[i][5][j] - nodex[i][1][j] +
           nodex[i][6][j] - nodex[i][2][j]);
      }
    } else {
      // add new element shape here
      error->one(FLERR,"Slip plane for this element shape has not been set");
    }

    // calculate cell sizes and
    // normalize vectors

    xnorm = norm3(xdir);
    ynorm = norm3(ydir);
    znorm = norm3(zdir);
    cross3(ydir,zdir,slip_plane[i][0]);
    cross3(zdir,xdir,slip_plane[i][1]);
    cross3(xdir,ydir,slip_plane[i][2]);
    norm3(slip_plane[i][0]);
    norm3(slip_plane[i][1]);
    norm3(slip_plane[i][2]);
    cell_size[i][0] = xnorm*fabs(dot3(xdir,slip_plane[i][0]));
    cell_size[i][1] = ynorm*fabs(dot3(ydir,slip_plane[i][1]));
    cell_size[i][2] = znorm*fabs(dot3(zdir,slip_plane[i][2]));
  }
}

/* ----------------------------------------------------------------------
   interpolate values from nodevalues in element I at interpolated atom J
   n = number of values
   ------------------------------------------------------------------------- */

void ElementVecCAC::interpolate(double *value, double ***nodevalue, int i, int j, int n) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  for (int k = 0; k < n; k++) {
    value[k] = 0.0;
    for (int l = 0; l < inpe; l++) 
      value[k] += shape_array[ietype][j][l]*nodevalue[i][l][k];
  }
}

/* ----------------------------------------------------------------------
   add atom plane next to elements
   Called from input add_atoms command
   ------------------------------------------------------------------------- */

void ElementVecCAC::add_atoms(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR,"Illegal add_atoms command");

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find add_atoms group ID");
  int groupbit = group->bitmask[igroup];

  // check which type to modify

  int n1list[4],n2list[4],dim0,dim1,dim2;
  if (strcmp(arg[1],"+x") == 0) {
    n1list[0] = 1; n2list[0] = 0;
    n1list[1] = 2; n2list[1] = 3;
    n1list[2] = 6; n2list[2] = 7;
    n1list[3] = 5; n2list[3] = 4;
    dim0 = 0; dim1 = 1; dim2 = 2;
  } else if (strcmp(arg[1],"-x") == 0) { 
    n1list[0] = 0; n2list[0] = 1;
    n1list[1] = 3; n2list[1] = 2;
    n1list[2] = 7; n2list[2] = 6;
    n1list[3] = 4; n2list[3] = 5;
    dim0 = 0; dim1 = 1; dim2 = 2;
  } else if (strcmp(arg[1],"+y") == 0) { 
    n1list[0] = 2; n2list[0] = 1;
    n1list[1] = 3; n2list[1] = 0;
    n1list[2] = 7; n2list[2] = 4;
    n1list[3] = 6; n2list[3] = 5;
    dim0 = 1; dim1 = 0; dim2 = 2;
  } else if (strcmp(arg[1],"-y") == 0) { 
    n1list[0] = 1; n2list[0] = 2;
    n1list[1] = 0; n2list[1] = 3;
    n1list[2] = 4; n2list[2] = 7;
    n1list[3] = 5; n2list[3] = 6;
    dim0 = 1; dim1 = 0; dim2 = 2;
  } else if (strcmp(arg[1],"+z") == 0) { 
    n1list[0] = 4; n2list[0] = 0;
    n1list[1] = 5; n2list[1] = 1;
    n1list[2] = 6; n2list[2] = 2;
    n1list[3] = 7; n2list[3] = 3;
    dim0 = 2; dim1 = 0; dim2 = 1;
  } else if (strcmp(arg[1],"-z") == 0) { 
    n1list[0] = 0; n2list[0] = 4;
    n1list[1] = 1; n2list[1] = 5;
    n1list[2] = 2; n2list[2] = 6;
    n1list[3] = 3; n2list[3] = 7;
    dim0 = 2; dim1 = 0; dim2 = 1;
  } else error->all(FLERR,"Illegal modify_elements command");

  int nlayers = universe->inumeric(FLERR,arg[2]);

  // modify elements and add atoms

  double coord[3];
  int iintpl,itype,ncellx,ncelly,ncellz;
  int natoms_added = 0;
  int nlocal_previous = atom->nlocal;
  for (int i = 0; i < element->nlocal; i++) 
    if (mask[i] & groupbit) {
      itype = etype[i];
      if (element_shape_ids[itype] != HEXAHEDRON)
        error->all(FLERR,"add_atoms command only works for rhombohedron element for now");

      ncellx = ncells[etype[i]][0];
      ncelly = ncells[etype[i]][1];
      ncellz = ncells[etype[i]][2];
      for (int n = 0; n < nlayers; n++) {
        for (int j = 0; j < ncells[itype][dim1]; j++) 
          for (int k = 0; k < ncells[itype][dim2]; k++) {
            iintpl = (j*(dim1 == 2) + k*(dim2 == 2) + (ncellz-1)*((dim1!=2)&&(dim2!=2)&&(n1list[0]==4)))
              + ncellz * (j*(dim1 == 1) + k*(dim2 == 1) + (ncelly-1)*((dim1!=1)&&(dim2!=1)&&(n1list[0]==2)))
              + ncelly * ncellz * (j*(dim1 == 0) + k*(dim2 == 0) + (ncellx-1)*((dim1!=0)&&(dim2!=0)&&(n1list[0]==1)));
            coord[0] = coord[1] = coord[2] = 0.0;
            for (int l = 0; l < npe[itype]; l++) {
              coord[0] += shape_array[itype][iintpl][l]*nodex[i][l][0];
              coord[1] += shape_array[itype][iintpl][l]*nodex[i][l][1];
              coord[2] += shape_array[itype][iintpl][l]*nodex[i][l][2];
            }
            atom->avec->create_atom(coord,ctype[i],0); 
            natoms_added++;
          }

        for (int j = 0; j < 4; j++) 
          for (int k = 0; k < 3; k++) 
            nodex[i][n1list[j]][k] += (nodex[i][n2list[j]][k] - nodex[i][n1list[j]][k])/ncells[etype[i]][dim0];
      }
    }

  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_variable(nlocal_previous,atom->nlocal);

  // add new tags to atoms

  atom->tag_extend();
  int total_atoms_added;
  MPI_Allreduce(&natoms_added,&total_atoms_added,1,MPI_INT,MPI_SUM,world);
  atom->natoms += total_atoms_added;
  if (comm->me == 0) fprintf(screen," %d atoms added by add_atoms command.\n",total_atoms_added);
}

/* ----------------------------------------------------------------------
   create a list of template elements from lattice within given bounds
   return total number of template elements
   called by create_elements command
   ------------------------------------------------------------------------- */

bigint ElementVecCAC::create_element_template(
    int klo, int khi, int jlo, int jhi, int ilo, int ihi, 
    int ietype, int *basisctype, double **xlist, int *clist)
{
  if (element_shape_ids[ietype] != HEXAHEDRON)
    error->all(FLERR,"create_elements command can only create rhombohedron element for now");
  bigint n = 0;

  int nbasis = domain->lattice->nbasis;
  double **basis = domain->lattice->basis;
  int i,j,k,m;
  klo = (klo/ncells[ietype][2]-1)*ncells[ietype][2];
  jlo = (jlo/ncells[ietype][1]-1)*ncells[ietype][1];
  ilo = (ilo/ncells[ietype][0]-1)*ncells[ietype][0];
  for (k = klo; k <= khi; k += ncells[ietype][2]) 
    for (j = jlo; j <= jhi; j += ncells[ietype][1]) 
      for (i = ilo; i <= ihi; i += ncells[ietype][0]) 
        for (m = 0; m < nbasis; m++) {
          xlist[n][0] = i + basis[m][0];
          xlist[n][1] = j + basis[m][1];
          xlist[n][2] = k + basis[m][2];
          clist[n] = basisctype[m];

          // convert from lattice coords to box coords

          domain->lattice->lattice2box(xlist[n][0],xlist[n][1],xlist[n][2]);
          n++;
        }

  return n;
}

/* ----------------------------------------------------------------------
   called from fix_adaptive to split element I along 'idim' direction 
   flag = 0: first pass, check if element need to be splitted and request new etype. Element is not splitted yet
   flag = 1: second pass, new etype has been created and element can now be splitted
   return 0 if element does not need splitted or discretized to all atoms
   1 if 1 new element created
   2 if 2 new elements created
   ------------------------------------------------------------------------- */

int ElementVecCAC::split_element(int i, int idim, double iplane, int split, int flag, double *plane_id)
{
  int split_flag = 0;
  int ietype = etype[i];
  int ictype = ctype[i];
  tagint itag = tag[i];
  double coord[3];
  int ncellx = ncells[ietype][0]; 
  int ncelly = ncells[ietype][1]; 
  int ncellz = ncells[ietype][2]; 
  int index_offset = (idim==2) + ncellz*(idim==1) + ncelly*ncellz*(idim==0);
  int type_info[8];
  type_info[0] = ncellx;
  type_info[1] = ncelly;
  type_info[2] = ncellz;
  type_info[3] = nintgs[ietype][0]; 
  type_info[4] = nintgs[ietype][1]; 
  type_info[5] = nintgs[ietype][2]; 
  type_info[6] = element_shape_ids[ietype];
  type_info[7] = apc[ietype];
  int clo = static_cast<int> (iplane - (ncells[ietype][idim] + 1)/2);
  int chi = static_cast<int> (iplane + (ncells[ietype][idim] - 1)/2);

  // check if element doesn't need splitting 

  if (split < clo || split >= chi) 
    return split_flag;

  // if flag = 0, request new etype if needed,
  // not actually splitting yet

  if (flag == 0) {

    split_flag = 1;

    // if has only one plane, discretize into atoms
    // 2D and 3D elements can't be in the same simulation yet

    // lower element   

    if (split > clo) {

      type_info[idim] = split - clo + 1;

      // limit the number of integration points to the number of unit cells

      type_info[3+idim] = MIN(type_info[idim],nintgs[ietype][idim]);

      // request new etype

      request_new_etype(type_info);
    }

    // upper element

    if (split+1 < chi) {

      type_info[idim] = chi - split;

      // limit the number of integration points to the number of unit cells

      type_info[3+idim] = MIN(type_info[idim],nintgs[ietype][idim]);

      // request new etype

      request_new_etype(type_info);
    }

    // if flag != 0, split elements

  } else {
    double **nodecoord;
    memory->create(nodecoord,element->max_npe,3,"evec:nodecoord");
    int *iintpl_list = new int[element->max_npe];
    int *n2ia = element->n2ia[ietype];

    // lower element
    // if has only one plane, discretize into atoms

    if (split == clo) {
      int hi[3],iintpl;
      hi[0] = ncellx;
      hi[1] = ncelly;
      hi[2] = ncellz;
      hi[idim] = 1;
      for (int ii = 0; ii < hi[0]; ii ++)
        for (int jj = 0; jj < hi[1]; jj ++)
          for (int kk = 0; kk < hi[2]; kk ++) {
            iintpl = kk + jj*ncellz + ii*ncellz*ncelly;
            interpolate(coord,nodex,i,iintpl,3);
            atom->avec->create_atom(coord,ictype,0);

            // pass on stored data for new atom J from fix/compute/variable and velocity

            int j = atom->nlocal-1;
            atom->data_pass_on_fix_compute_variable(j,i,iintpl);
            interpolate(atom->v[j],nodev,i,iintpl,3);
          } 
    } else {
      split_flag++;
      type_info[idim] = split - clo + 1;

      // limit the number of integration points to the number of unit cells

      type_info[3+idim] = MIN(type_info[idim],nintgs[ietype][idim]);

      int myetype = find_etype(type_info);
      if (myetype == 0) error->one(FLERR,"Etype not defined yet");

      coord[0] = coord[1] = coord[2] = 0;

      if (domain->dimension == 3) {
        for (int ii = 0; ii < 4; ii++) {

          // lower set of nodes stays are the same as the origional element
          // calculate upper set of nodes by scaling along edges

          int lnode = node_set_3D[idim][0][ii];
          int unode = node_set_3D[idim][1][ii]; 

          iintpl_list[lnode] = n2ia[lnode];
          iintpl_list[unode] = n2ia[lnode] + (type_info[idim]-1)*index_offset;

          for (int jj = 0; jj < 3; jj++) {
            nodecoord[lnode][jj] = nodex[i][lnode][jj];
            nodecoord[unode][jj] = nodecoord[lnode][jj] + 
              (nodex[i][unode][jj]-nodex[i][lnode][jj])
              *(type_info[idim]-1)/(ncells[ietype][idim]-1);
            coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
          }
        }
        coord[0] /= 8;
        coord[1] /= 8;
        coord[2] /= 8;
        create_element(coord,nodecoord,myetype,ictype,0);

        // pass on stored data for new element in fix/compute/variable and velocity

        int j = element->nlocal-1;
        element->data_pass_on_fix_compute_variable(j,i,iintpl_list);
        for (int ii = 0; ii < 8; ii++)
          interpolate(nodev[j][ii],nodev,i,iintpl_list[ii],3);

      } else {
        for (int ii = 0; ii < 2; ii++) {

          // lower set of nodes stays are the same as the origional element
          // calculate upper set of nodes by scaling along edges

          int lnode = node_set_2D[idim][0][ii];
          int unode = node_set_2D[idim][1][ii]; 

          iintpl_list[lnode] = n2ia[lnode];
          iintpl_list[unode] = n2ia[lnode] + (type_info[idim]-1)*index_offset;

          for (int jj = 0; jj < 2; jj++) {
            nodecoord[lnode][jj] = nodex[i][lnode][jj];
            nodecoord[unode][jj] = nodecoord[lnode][jj] + 
              (nodex[i][unode][jj]-nodex[i][lnode][jj])
              *(type_info[idim]-1)/(ncells[ietype][idim]-1);
            coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
          }
        }
        coord[0] /= 4;
        coord[1] /= 4;

        // z direction coord stay the same

        coord[2] = x[i][2];       
        for (int ii = 0; ii < 4; ii++) 
          nodecoord[ii][2] = nodex[i][ii][2];
        create_element(coord,nodecoord,ietype,ictype,0);
        int j = element->nlocal-1;
        element->data_pass_on_fix_compute_variable(j,i,iintpl_list);
        for (int ii = 0; ii < 4; ii++)
          interpolate(nodev[j][ii],nodev,i,iintpl_list[ii],3);

      }

      // assign new plane_id for new element

      plane_id[element->nlocal-1] = static_cast<double> (split+clo);
      plane_id[element->nlocal-1] /= 2;
    }

    // upper element
    // if has only one plane, discretize into atoms

    if (split+1 == chi) {
      int lo[3],iintpl;
      lo[0] = lo[1] = lo[2] = 0;
      lo[idim] = ncells[ietype][idim]-1;
      for (int ii = lo[0]; ii < ncellx; ii ++)
        for (int jj = lo[1]; jj < ncelly; jj ++)
          for (int kk = lo[2]; kk < ncellz; kk ++) {
            iintpl = kk + jj*ncellz + ii*ncellz*ncelly;
            interpolate(coord,nodex,i,iintpl,3);
            atom->avec->create_atom(coord,ictype,0);
            int j = atom->nlocal-1;
            atom->data_pass_on_fix_compute_variable(j,i,iintpl);
            interpolate(atom->v[j],nodev,i,iintpl,3);
          } 

    } else {

      split_flag++;
      type_info[idim] = chi - split;

      // limit the number of integration points to the number of unit cells

      type_info[3+idim] = MIN(type_info[idim],nintgs[ietype][idim]);

      int myetype = find_etype(type_info);
      if (myetype == 0) error->one(FLERR,"Etype not defined yet");

      coord[0] = coord[1] = coord[2] = 0;
      if (domain->dimension == 3) {
        for (int ii = 0; ii < 4; ii++) {

          // upper set of nodes stays are the same as the origional element
          // calculate lower set of nodes by scaling along edges

          int lnode = node_set_3D[idim][0][ii];
          int unode = node_set_3D[idim][1][ii]; 

          iintpl_list[unode] = n2ia[unode];
          iintpl_list[lnode] = n2ia[unode] - (type_info[idim]-1)*index_offset;

          for (int jj = 0; jj < 3; jj++) {
            nodecoord[unode][jj] = nodex[i][unode][jj];
            nodecoord[lnode][jj] = nodecoord[unode][jj] + 
              (nodex[i][lnode][jj]-nodex[i][unode][jj])
              *(type_info[idim]-1)/(ncells[ietype][idim]-1);
            coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
          }
        }
        coord[0] /= 8;
        coord[1] /= 8;
        coord[2] /= 8;
        create_element(coord,nodecoord,myetype,ictype,0);
        int j = element->nlocal-1;
        element->data_pass_on_fix_compute_variable(j,i,iintpl_list);
        for (int ii = 0; ii < 8; ii++)
          interpolate(nodev[j][ii],nodev,i,iintpl_list[ii],3);

      } else {
        for (int ii = 0; ii < 2; ii++) {

          // lower set of nodes stays are the same as the origional element
          // calculate upper set of nodes by scaling along edges

          int lnode = node_set_2D[idim][0][ii];
          int unode = node_set_2D[idim][1][ii]; 

          iintpl_list[unode] = n2ia[unode];
          iintpl_list[lnode] = n2ia[unode] - (type_info[idim]-1)*index_offset;

          for (int jj = 0; jj < 2; jj++) {
            nodecoord[unode][jj] = nodex[i][unode][jj];
            nodecoord[lnode][jj] = nodecoord[unode][jj] + 
              (nodex[i][lnode][jj]-nodex[i][unode][jj])
              *(type_info[idim]-1)/(ncells[ietype][idim]-1);
            coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
          }
        }
        coord[0] /= 4;
        coord[1] /= 4;

        // z direction coord stay the same

        coord[2] = x[i][2];       
        for (int ii = 0; ii < 4; ii++) 
          nodecoord[ii][2] = nodex[i][ii][2];

        create_element(coord,nodecoord,ietype,ictype,0);
        int j = element->nlocal-1;
        element->data_pass_on_fix_compute_variable(j,i,iintpl_list);
        for (int ii = 0; ii < 4; ii++)
          interpolate(nodev[j][ii],nodev,i,iintpl_list[ii],3);

      }

      // assign new plane_id for new element

      plane_id[element->nlocal-1] = static_cast<double> (split+chi+1);
      plane_id[element->nlocal-1] /= 2;

    }

    // delete element I

    copy(element->nlocal-1,i,1);
    element->nlocal--;
    delete [] iintpl_list;
    memory->destroy(nodecoord);
  }
  return split_flag;
}

/* ----------------------------------------------------------------------
   request new etype
   return if etype already exist or requested
   otherwise add to list of requested etype
   ------------------------------------------------------------------------- */

void ElementVecCAC::request_new_etype(int *type_info)
{
  if (find_etype(type_info)) return;

  int n = nrequested_etype;

  // check if requested etype exists

  for (int i = 0; i < n; i++) 
    if (requested_etype[i*8] == type_info[0] &&
        requested_etype[i*8+1] == type_info[1] &&
        requested_etype[i*8+2] == type_info[2] &&
        requested_etype[i*8+3] == type_info[3] &&
        requested_etype[i*8+4] == type_info[4] &&
        requested_etype[i*8+5] == type_info[5] &&
        requested_etype[i*8+6] == type_info[6] &&
        requested_etype[i*8+7] == type_info[7])
      return;

  // grow requested_etype array if needed

  if (n == maxrequested_etype) {
    maxrequested_etype *= 2;
    memory->grow(requested_etype,maxrequested_etype*8,"evec:requested_etype");
  }

  // add requested etype to the list

  for (int i = 0; i < 8; i++)
    requested_etype[n*8+i] = type_info[i];

  nrequested_etype++;

}

/* ----------------------------------------------------------------------
   find existing etype
   return the index of etype if match
   return 0 if no etype found
   ------------------------------------------------------------------------- */

int ElementVecCAC::find_etype(int *type_info)
{
  for (int i = 1; i <= element->netypes; i++)
    if (type_info[0] == ncells[i][0] && 
        type_info[1] == ncells[i][1] && 
        type_info[2] == ncells[i][2] && 
        type_info[3] == nintgs[i][0] && 
        type_info[4] == nintgs[i][1] && 
        type_info[5] == nintgs[i][2] &&
        type_info[6] == element_shape_ids[i] &&
        type_info[7] == apc[i])

      return i;
  return 0;
}

/* ----------------------------------------------------------------------
   union all requested etype from all procs
   add the requested etype
   ------------------------------------------------------------------------- */

void ElementVecCAC::add_requested_etype()
{
  int nprocs = comm->nprocs;
  int *recvcount = new int[nprocs];
  int *recvbuf = new int[nprocs];
  int *displs = new int[nprocs];
  for (int i = 0; i < nprocs; i++) {
    recvcount[i] = 1;
    displs[i] = i;
  }

  MPI_Allgatherv(&nrequested_etype,1,MPI_INT,recvbuf,recvcount,displs,MPI_INT,world);

  int sum = 0;
  for (int i = 0; i < nprocs; i++) {
    recvcount[i] = recvbuf[i]*8;
    displs[i] = sum;
    sum += recvcount[i];
  }
  int nrequests_all = sum/8;
  int *requests_all = new int[sum];

  MPI_Allgatherv(requested_etype,nrequested_etype*8,MPI_INT,
      requests_all,recvcount,displs,MPI_INT,world);

  char **newarg = new char*[9];
  for (int i = 0; i < nrequests_all; i++) {
    int *type_info = requests_all+i*8;
    if (!find_etype(type_info)) {
      for (int j = 0; j < 9; j++)
        newarg[j] = new char[256];
      sprintf(newarg[0],"%d",element->netypes+1);
      sprintf(newarg[1],"%s",element->element_shape_list[type_info[6]]);
      sprintf(newarg[2],"%d",type_info[7]);
      sprintf(newarg[3],"%d",type_info[0]);
      sprintf(newarg[4],"%d",type_info[1]);
      sprintf(newarg[5],"%d",type_info[2]);
      sprintf(newarg[6],"%d",type_info[3]);
      sprintf(newarg[7],"%d",type_info[4]);
      sprintf(newarg[8],"%d",type_info[5]);
      element->add_etype(9,newarg);
      for (int j = 0; j < 9; j++)
        delete [] newarg[j];
    }
  }

  // reset number of requested etype

  nrequested_etype = 0;

  // update elemscale in pair

  force->pair->update_elemscale();

  // clean up

  delete [] newarg; 
  delete [] requests_all;
  delete [] displs;
  delete [] recvbuf;
  delete [] recvcount;
}

/* ----------------------------------------------------------------------
   check if element J surface K cut through element I
   ------------------------------------------------------------------------- */

int ElementVecCAC::check_split_elem(int i, int j, int surface, int &cut, int &cutdim)
{
  // look for node closed to element I
  // project node coord onto each slip plane normal line
  // if projection falls inside element --> need splitting

  int jnode_closest,jnode,jj;
  double delx,dely,delz,rsq;
  double rsqmin = BIG;
  double ix = x[i][0];
  double iy = x[i][1];
  double iz = x[i][2];
  int ietype = etype[i];

  cutdim = surface/2;
  int dir = surface%2;
  if (domain->dimension == 3) {
    for (jj = 0; jj < 4; jj++) {
      jnode = node_set_3D[cutdim][dir][jj];
      delx = ix - nodex[j][jnode][0];
      dely = iy - nodex[j][jnode][1];
      delz = iz - nodex[j][jnode][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) { 
        rsqmin = rsq;
        jnode_closest = jnode;
      }
    }
  } else {
    for (jj = 0; jj < 2; jj++) {
      jnode = node_set_2D[cutdim][dir][jj];
      delx = ix - nodex[j][jnode][0];
      dely = iy - nodex[j][jnode][1];
      rsq = delx*delx + dely*dely;
      if (rsq < rsqmin) { 
        rsqmin = rsq;
        jnode_closest = jnode;
      }
    }
  }
  double dist = dot3(slip_plane[i][cutdim],nodex[j][jnode_closest]) 
    - dot3(slip_plane[i][cutdim],x[i]);
  double a = cell_size[i][cutdim]/(ncells[ietype][cutdim]-1);
  double hi = cell_size[i][cutdim]/2 - a/2;
  double lo = -hi;

  if (dir) lo -= a;
  else hi += a;
  if (dist < lo || dist > hi) return 0;

  cut = static_cast<int> ((dist-lo)/a);
  cut++;
  return 1;
}

/* ----------------------------------------------------------------------
   called from fix_adaptive to split elements in list
   return number of elements being splitted
   ------------------------------------------------------------------------- */

int ElementVecCAC::split_element(int **splitlist, int nsplit)
{
  int ietype,ictype;
  int type_info[8];
  //***NEED FIX
  type_info[6]=type_info[7]=1;
  int *nelem,**cellsize;
  int ii = 0;
  int i,icut,icutdim;
  int *ilist,**nelem_array,***cellsize_array;

  // counter for number of elements being splitted

  int n = 0;

  if (nsplit > 0) {
    memory->create(ilist,nsplit,"evec:ilist");
    memory->create(nelem_array,nsplit,3,"evec:nelem");
    memory->create(cellsize_array,nsplit,3,MAXSUB,"evec:cellsize");

    // first pass: check for new etypes
    // loop through all splits in list
    // each element might have several splits and
    // they are all next to each other in splitlist

    while (ii < nsplit) {
      nelem = nelem_array[n];
      cellsize = cellsize_array[n];
      nelem[0] = nelem[1] = nelem[2] = 1;

      i = splitlist[ii][0];

      // apply first split to element I

      icut = splitlist[ii][1];
      icutdim = splitlist[ii][2];
      ilist[n] = i;
      ietype = etype[i];
      cellsize[0][0] = ncells[ietype][0];
      cellsize[1][0] = ncells[ietype][1];
      cellsize[2][0] = ncells[ietype][2];
      cellsize[icutdim][0] = icut;
      cellsize[icutdim][1] = ncells[ietype][icutdim]-icut;
      nelem[icutdim]++;

      // apply each additional split to element cellsize list if next split 
      // still cut element I

      if (ii < nsplit-1) {
        while (i == splitlist[ii+1][0]) {

          icut = splitlist[ii+1][1];
          icutdim = splitlist[ii+1][2];

          if (nelem[icutdim] == MAXSUB) 
            error->one(FLERR,"Too many split in one direction");
          int sum = 0;

          for (int l = 0; l < nelem[icutdim]; l++) {
            if (icut > sum + cellsize[icutdim][l])
              sum += cellsize[icutdim][l];
            else if (icut < sum + cellsize[icutdim][l]) {
              for (int m = l+1; m < nelem[icutdim]; m++) 
                cellsize[icutdim][m+1] = cellsize[icutdim][m];
              cellsize[icutdim][l] = icut-sum;
              cellsize[icutdim][l+1] -= icut-sum;
              nelem[icutdim]++;
              break;
            } else {
              break;
            }
          }
          ii++;
          if (ii == nsplit-1) break;
        }
      }
      for (int ix = 0; ix < nelem[0]; ix++) {
        type_info[0] = cellsize[0][ix];
        if (type_info[0] == 1) continue;
        type_info[3] = MIN(type_info[0],nintgs[ietype][0]);
        for (int iy = 0; iy < nelem[1]; iy++) {
          type_info[1] = cellsize[1][iy];
          if (type_info[1] == 1) continue;
          type_info[4] = MIN(type_info[1],nintgs[ietype][1]);
          for (int iz = 0; iz < nelem[2]; iz++) {
            type_info[2] = cellsize[2][iz];
            type_info[5] = MIN(type_info[2],nintgs[ietype][2]);
            if (type_info[2] == 1 && domain->dimension == 3) 
              continue;
            request_new_etype(type_info);
          }
        }
      }
      n++;
      ii++;
    }
  }

  add_requested_etype();

  // second pass: split elements

  if (nsplit > 0) {
    int index_offset[3],index_node0,iintpl;
    int myetype,ncellx,ncelly,ncellz;
    double coord[3],**nodecoord;
    memory->create(nodecoord,element->max_npe,3,"evec:nodecoord");

    // iintpl_list = list of intpl index in old element for each node of new element

    int *iintpl_list = new int[element->max_npe]; 

    for (ii = 0; ii < n; ii++) {
      i = ilist[ii];
      ietype = etype[i];
      ncellx = ncells[ietype][0];
      ncelly = ncells[ietype][1];
      ncellz = ncells[ietype][2];
      ictype = ctype[i];
      nelem = nelem_array[ii];
      cellsize = cellsize_array[ii];
      index_offset[0] = index_offset[1] = index_offset[2] = 0;
      index_node0 = 0;
      for (int ix = 0; ix < nelem[0]; ix++) {
        type_info[0] = cellsize[0][ix];
        type_info[3] = MIN(type_info[0],nintgs[ietype][0]);
        index_offset[0] = ncelly*ncellz*(cellsize[0][ix]-1);
        for (int iy = 0; iy < nelem[1]; iy++) {
          type_info[1] = cellsize[1][iy];
          type_info[4] = MIN(type_info[1],nintgs[ietype][1]);
          index_offset[1] = ncellz*(cellsize[1][iy]-1);
          if (domain->dimension == 3) {
            for (int iz = 0; iz < nelem[2]; iz++) {
              type_info[2] = cellsize[2][iz];
              type_info[5] = MIN(type_info[2],nintgs[ietype][2]);
              index_offset[2] = (cellsize[2][iz]-1);
              if (type_info[0] == 1) {
                for (int iiy = 0; iiy < type_info[1]; iiy++)
                  for (int iiz = 0; iiz < type_info[2]; iiz++) 
                    create_pass_on_atom(i,index_node0+iiz+iiy*ncellz,ictype,0);
              } else if (type_info[1] == 1) {
                for (int iix = 0; iix < type_info[0]; iix++) 
                  for (int iiz = 0; iiz < type_info[2]; iiz++) 
                    create_pass_on_atom(i,index_node0+iiz+iix*ncelly*ncellz,ictype,0);
              } else if (type_info[2] == 1) {
                for (int iix = 0; iix < type_info[0]; iix++) 
                  for (int iiy = 0; iiy < type_info[1]; iiy++) 
                    create_pass_on_atom(i,index_node0+iiy*ncellz+iix*ncelly*ncellz,ictype,0);
              } else {
                myetype = find_etype(type_info);
                if (myetype == 0) error->one(FLERR,"Etype not defined yet");

                // calculate intpl index in element I for each node of new element

                iintpl_list[0] = index_node0; 
                iintpl_list[1] = index_node0 + index_offset[0];
                iintpl_list[2] = index_node0 + index_offset[0] 
                  + index_offset[1];
                iintpl_list[3] = index_node0 + index_offset[1];
                iintpl_list[4] = iintpl_list[0] + index_offset[2];
                iintpl_list[5] = iintpl_list[1] + index_offset[2];
                iintpl_list[6] = iintpl_list[2] + index_offset[2];
                iintpl_list[7] = iintpl_list[3] + index_offset[2];
                create_pass_on_element(i,iintpl_list,myetype,ictype,0);
              }
              index_node0 += index_offset[2]+1;
            }
            index_node0 += index_offset[1];
          } else {
            type_info[2] = type_info[5] = 1;
            if (type_info[0] == 1) {
              for (int iiy = 0; iiy < type_info[1]; iiy++)
                create_pass_on_atom(i,index_node0+iiy,ictype,0);
            } else if (type_info[1] == 1) {
              error->one(FLERR,"TEST"); 
              //            for (int iix = 0; iix < info[0]; iix++) 
              //              create_pass_on_atom(i,index_node0+iix*ncellyz,ictype,0);
            } else {
              myetype = find_etype(type_info);
              if (myetype == 0) error->one(FLERR,"Etype not defined yet");
              iintpl_list[0] = index_node0; 
              iintpl_list[1] = index_node0 + index_offset[0];
              iintpl_list[2] = index_node0 + index_offset[0] 
                + index_offset[1];
              iintpl_list[3] = index_node0 + index_offset[1];
              create_pass_on_element(i,iintpl_list,myetype,ictype,0);
            }
            index_node0 += index_offset[1]+1;
          }
        }
        index_node0 += index_offset[0];
      }

      // delete element I

      copy(element->nlocal-1,i,1);
      element->nlocal--;
    }
    memory->destroy(nelem_array); 
    memory->destroy(cellsize_array);
    memory->destroy(nodecoord);
    delete [] iintpl_list;
  }

  return n;
}

/* ----------------------------------------------------------------------
   create new element from an old element
   ------------------------------------------------------------------------- */

void ElementVecCAC::create_pass_on_element(int i, int *iintpl_list, int ietype, int ictype, tagint itag)
{
  double **nodecoord,coord[3];
  int inpe = etype[ietype];
  memory->create(nodecoord,inpe,3,"evec:nodecoord");
  coord[0] = coord[1] = coord[2] = 0;
  for (int jj = 0; jj < inpe; jj++) {
    interpolate(nodecoord[jj],nodex,i,iintpl_list[jj],3);
    coord[0] += nodecoord[jj][0];
    coord[1] += nodecoord[jj][1];
    coord[2] += nodecoord[jj][2];
  }
  coord[0] /= inpe;
  coord[1] /= inpe;
  coord[2] /= inpe;
  create_element(coord,nodecoord,ietype,ictype,itag);

  // pass on stored data for new element in fix/compute/variable and velocity
  // new elements will be in the same group as old element

  int j = element->nlocal-1;
  element->data_pass_on_fix_compute_variable(j,i,iintpl_list);
  mask[j] = mask[i];
  for (int jj = 0; jj < inpe; jj++)
    interpolate(nodev[j][jj],nodev,i,iintpl_list[jj],3);
  memory->destroy(nodecoord);
}

/* ----------------------------------------------------------------------
   create new atom from an old element
   ------------------------------------------------------------------------- */

void ElementVecCAC::create_pass_on_atom(int i, int iintpl, int itype, tagint itag)
{
  double coord[3];
  interpolate(coord,nodex,i,iintpl,3);
  atom->avec->create_atom(coord,itype,itag);

  // pass on stored data for new element in fix/compute/variable and velocity
  // new atoms will be in the same group as old element

  int j = atom->nlocal-1;
  atom->data_pass_on_fix_compute_variable(j,i,iintpl);
  atom->mask[j] = mask[i];
  interpolate(atom->v[j],nodev,i,iintpl,3);
}

