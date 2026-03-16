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
#include "update.h"
#include "pair.h"
#include "group.h"
#include "error.h"
#include "math_extra.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathExtra;

#define EPSILON 1e-6
#define BIG     1e30
#define MAXSUB  20

/* ---------------------------------------------------------------- */

ElementVecCAC::ElementVecCAC(CAC *cac) : ElementVec(cac)
{

  ncells = ngcells = NULL;

  // size of data line in write_data_elem command

  size_data_elem = 5;
  size_data_elem_vel = 5;

  // size of data line in write_data command

  size_data_element = 5;
  size_data_node = 7;
  size_data_vel = 6;

  // position of first x column in read_data

  xcol_data = 3;

  // size of data line in write_tecplot command

  size_tecplot_node = 6;

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

}

/* ---------------------------------------------------------------- */

ElementVecCAC::~ElementVecCAC()
{

}

/*  ----------------------------------------------------------------------
   set communication sizes that are dependent on maxnpe and maxapc
   -------------------------------------------------------------------------  */

void ElementVecCAC::set_data_size() 
{
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;

  // size of data line in write_tecplot command

  if (domain->dimension == 3) 
    size_tecplot_node_connect = 8;
  else
    size_tecplot_node_connect = 4;

  // size of communication

  size_velocity = 3 * maxnpe * maxapc;
  size_forward = 9 + 3 * maxnpe * maxapc;
  size_reverse = 3 * element->maxgcell * maxapc;
  if (atom->atom_strain_flag) {
    size_border = 13 + (7 * maxnpe + 1) * maxapc;
  } else {
    size_border = 12 + (4 * maxnpe + 1) * maxapc;
  }

}

/*  ----------------------------------------------------------------------
   copy element I info to element J
   -------------------------------------------------------------------------  */

void ElementVecCAC::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  for (int l = 0; l < element->apc[etype[i]]; l++) 
    ctype[j][l] = ctype[i][l];
  etype[j] = etype[i];
  mask[j] = mask[i];
  image[j] = image[i];
  initial_box_size[j][0] = initial_box_size[i][0];
  initial_box_size[j][1] = initial_box_size[i][1];
  initial_box_size[j][2] = initial_box_size[i][2];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  for (int l = 0; l < element->apc[etype[i]]; l++) 
    for (int k = 0; k < element->npe[etype[i]]; k++)
      nodemask[j][l][k] = nodemask[i][l][k];

  for (int k = 0; k < element->npe[etype[i]]; k++) 
    for (int l = 0; l < element->apc[etype[i]]; l++) {
      nodex[j][l][k][0] = nodex[i][l][k][0];
      nodex[j][l][k][1] = nodex[i][l][k][1];
      nodex[j][l][k][2] = nodex[i][l][k][2];
      if (atom->atom_strain_flag) {
        nodex_current[j][l][k][0] = nodex_current[i][l][k][0];
        nodex_current[j][l][k][1] = nodex_current[i][l][k][1];
        nodex_current[j][l][k][2] = nodex_current[i][l][k][2];
      }
      nodev[j][l][k][0] = nodev[i][l][k][0];
      nodev[j][l][k][1] = nodev[i][l][k][1];
      nodev[j][l][k][2] = nodev[i][l][k][2];
    }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      modify->fix[element->extra_grow[iextra]]->copy_elem_arrays(i, j, delflag);
}

/*  ----------------------------------------------------------------------
   grow element arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
   -------------------------------------------------------------------------  */

void ElementVecCAC::grow(int n)
{
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;
  int maxgcell = element->maxgcell;
  if (n == 0) grow_nmax();
  else {
    nmax = n;
  }
  element->nmax = nmax;

  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR, "Per-processor system is too big");

  tag = memory->grow(element->tag, nmax, "element:tag");
  ctype = memory->grow(element->ctype, nmax, maxapc, "element:ctype");
  etype = memory->grow(element->etype, nmax, "element:etype");
  mask = memory->grow(element->mask, nmax, "element:mask");
  nodemask = memory->grow(element->nodemask, nmax, maxapc, maxnpe, "element:nodemask");
  image = memory->grow(element->image, nmax, "element:image");
  x = memory->grow(element->x, nmax, 3, "element:x");
  surface_plane = memory->grow(element->surface_plane, nmax, 3, 3, "element:surface_plane");
  cell_size = memory->grow(element->cell_size, nmax, 3, "element:cell_size");
  nodex = memory->grow(element->nodex, nmax, maxapc, maxnpe, 3, "element:nodex");
  if (atom->atom_strain_flag)
    nodex_current = memory->grow(element->nodex_current, nmax, maxapc, maxnpe, 3, "element:nodex_current");
  nodev = memory->grow(element->nodev, nmax, maxapc, maxnpe, 3, "element:nodev");
  nodef = memory->grow(element->nodef, nmax, maxapc, maxnpe, 3, "element:nodef");
  gaussf = memory->create(element->gaussf, nmax, maxapc, maxgcell, 3, "element:gaussf");
  subelem_size = memory->grow(element->subelem_size, nmax, "element:subelem_size");
  initial_box_size = memory->grow(element->initial_box_size, nmax, 3, "element:initial_box_size");
  element_bound_box = memory->grow(element->element_bound_box, nmax, 6, "element:element_bound_box");

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      modify->fix[element->extra_grow[iextra]]->grow_elem_arrays(nmax);

}

/*  ---------------------------------------------------------
    grow etype arrays to n types
    if n = 0, destroy old arrays
 * ------------------------------------------------------ */

void ElementVecCAC::grow_etype_arrays(int n)
{
  if (n == 0) {
    element->netypes = 0;
    memory->destroy(element_type_setflag);
    memory->destroy(ngcells);
    memory->destroy(ncells);
    element->destroy_etype_arrays();

  } else {
    int maxnpe = element->maxnpe;
    int old_netypes = element->netypes;
    if (n <= old_netypes) return;
    else element->netypes = n;


    // grow style specific arrays

    memory->grow(element_type_setflag, n+1, "evec:element_type_setflag");

    // grow gaussian cell arrays

    memory->grow(ngcells, n+1, 3, "evec:ngcells");

    // virtual atom arrays

    memory->grow(ncells, n+1, 3, "evec:ncells");

    // grow arrays common to all element styles

    npe = memory->grow(element->npe, n+1, "element:npe");
    apc = memory->grow(element->apc, n+1, "element:apc");
    nodal_weight = memory->grow(element->nodal_weight, n+1, "element:nodal_weight");
    element_shape_ids = memory->grow(element->element_shape_ids, n+1, "element:element_shape_ids");
    element_shape_names = memory->grow(element->element_shape_names, n+1, 256, "element:element_shape_names");

    // grow gaussian cell arrays

    ngcell = memory->grow(element->ngcell, n+1, "element:ngcell");
    weighted_shape_array = memory->grow(element->weighted_shape_array, 
        n+1, element->maxgcell, maxnpe, "element:weighted_shape_array");
    g2u = memory->grow(element->g2u, n+1, element->maxgcell, "element:g2u");
    weight = memory->grow(element->weight, n+1, element->maxgcell, "element:weight");
    g2n = memory->grow(element->g2n, n+1, element->maxgcell, "element:g2n");

    debug_gcell_type = memory->grow(element->debug_gcell_type, n+1, element->maxgcell, "evec:debug_gcell_type");
    debug_gcell_inactive = memory->grow(element->debug_gcell_inactive, n+1, 4, "evec:debug_gcell_inactive");

    // grow node arrays

    n2g = memory->grow(element->n2g, n+1, maxnpe, "element:n2g");
    n2u = memory->grow(element->n2u, n+1, maxnpe, "element:n2u");

    // virtual atom arrays

    nucell = memory->grow(element->nucell, n+1, "element:nucell");
    surface_ucell = memory->grow(element->surface_ucell, n+1, element->max_surface, element->max_surface_ucell, "element:surface_ucell");
    nsurface_ucell = memory->grow(element->nsurface_ucell, n+1, element->max_surface, "element:nsurface_ucell");
    edge_ucell = memory->grow(element->edge_ucell, n+1, element->max_edge, element->max_edge_ucell, "element:edge_ucell");
    nedge_ucell = memory->grow(element->nedge_ucell, n+1, element->max_edge, "element:nedge_ucell");
    is_outer = memory->grow(element->is_outer, n+1, element->maxucell, "element:is_outer");
    shape_array = memory->grow(element->shape_array, n+1, element->maxucell, maxnpe, "element:shape_array");
    u2g = memory->grow(element->u2g, n+1, element->maxucell, "element:u2g");

    // sub-element array

    nucell_subelem = memory->grow(element->nucell_subelem, n+1, element->maxsubelem, "element:nucell_subelem");
    us2u = memory->grow(element->us2u, n+1, element->maxsubelem, element->maxsubucell, "element:us2u");
    shape_array_center_subelem = memory->grow(element->shape_array_center_subelem, 
        n+1, element->maxsubelem, maxnpe, "element:shape_array_center_subelem");
    nsubelem = memory->grow(element->nsubelem, n+1, "element:nsubelem");
    subsplit = memory->grow(element->subsplit, n+1, "element:subsplit");

    for (int itype = old_netypes + 1; itype <= n; itype++) {
      element_type_setflag[itype] = 0;
      ngcell[itype] = 0; 
      nucell[itype] = 0; 
      for (int i = 0; i < element->maxucell; i++) {
        u2g[itype][i] = -1;
        is_outer[itype][i] = 0;
      }
      for (int i = 0; i < element->maxgcell; i++) {
        g2u[itype][i] = -1;
        g2n[itype][i] = -1;
        debug_gcell_type[itype][i] = -1;
      }
      for (int i = 0; i < maxnpe; i++)
        n2g[itype][i] = -1;
      for (int i = 0; i < 4; i++)
        debug_gcell_inactive[itype][i] = 0;
    }
  }
}

/*  ----------------------------------------------------------------------
    unpack one line from Elements section of data file
    initialize other element quantities
    -------------------------------------------------------------------------  */

void ElementVecCAC::data_element(double *coord, imageint imagetmp, char **values)
{
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);  

  tag[nlocal] = ATOTAGINT(values[0]);
  etype[nlocal] = atoi(values[1]);
  if (tag[nlocal] <= 0) {
    error->one(FLERR, "Invalid element ID in Elements section of data file");
  }
  if (etype[nlocal] <= 0 || etype[nlocal] > element->netypes)
    error->one(FLERR, "Invalid element type in Elements section of data file"); 

  image[nlocal] = imagetmp;

  // 1 equal to the all group, 2 equal to the element group

  mask[nlocal] = 1|2;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  initial_box_size[nlocal][0] = -1.0;
  initial_box_size[nlocal][1] = -1.0;
  initial_box_size[nlocal][2] = -1.0;
  element->nucells += nucell[etype[nlocal]];
  element->nvatoms += nucell[etype[nlocal]] * apc[etype[nlocal]];
  element->nlocal++;

}

/*  ----------------------------------------------------------------------
    unpack one line from Node Velocities section of data file
    -------------------------------------------------------------------------  */

void ElementVecCAC::data_vel(int m, char **values)
{
  int inode = atoi(values[0])-1; 
  int ibasis = atoi(values[1])-1; 
  if (inode >= npe[etype[m]] || inode < 0) 
    error->all(FLERR, "Invalid node index value in Node section of data file");
  if (ibasis >= apc[etype[m]] || ibasis < 0) 
    error->all(FLERR, "Invalid basis index value in Node section of data file");
  nodev[m][ibasis][inode][0] = atof(values[1]);
  nodev[m][ibasis][inode][1] = atof(values[2]);
  nodev[m][ibasis][inode][2] = atof(values[3]);
}

/*  ----------------------------------------------------------------------
    unpack one line from Node of data file
    -------------------------------------------------------------------------  */

void ElementVecCAC::data_node(int m, char **values, int ctype_offset, int shiftflag, double *shift)
{
  int inode = atoi(values[0])-1; 
  int ibasis = atoi(values[1])-1; 
  ctype[m][ibasis] = atoi(values[2]);
  if (ctype_offset) {
    ctype[m][ibasis] += ctype_offset;
    if (ctype[m][ibasis] > atom->ntypes)
      error->one(FLERR, "Invalid atom type in Nodes section of data file");
  }

  if (inode >= npe[etype[m]] || inode < 0) 
    error->all(FLERR, "Invalid node index value in Node section of data file");
  if (ibasis >= apc[etype[m]] || ibasis < 0) 
    error->all(FLERR, "Invalid basis index value in Node section of data file");

  double coord[3];

  coord[0] = atof(values[3]);
  coord[1] = atof(values[4]);
  coord[2] = atof(values[5]);
  domain->unmap_reverse(coord, image[m]);
  nodex[m][ibasis][inode][0] = coord[0];
  nodex[m][ibasis][inode][1] = coord[1];
  nodex[m][ibasis][inode][2] = coord[2];
  if (shiftflag) {
    nodex[m][ibasis][inode][0] += shift[0];
    nodex[m][ibasis][inode][1] += shift[1];
    nodex[m][ibasis][inode][2] += shift[2];
  }

  // 1 equal to the all group, 2 equal to the element group

  nodemask[m][ibasis][inode] = 1|2;
  nodev[m][ibasis][inode][0] = 0.0;
  nodev[m][ibasis][inode][1] = 0.0;
  nodev[m][ibasis][inode][2] = 0.0;
}

/*  ----------------------------------------------------------------------
    unpack one line from Node of data file
    -------------------------------------------------------------------------  */

void ElementVecCAC::data_node_strain(int m, char **values)
{
  int inode = atoi(values[0])-1; 
  int ibasis = atoi(values[1])-1; 
  if (inode >= npe[etype[m]] || inode < 0) 
    error->all(FLERR, "Invalid node index value in Node section of data file");
  if (ibasis >= apc[etype[m]] || ibasis < 0) 
    error->all(FLERR, "Invalid basis index value in Node section of data file");
  double coord[3];
  nodex[m][ibasis][inode][0] = atof(values[2]);
  nodex[m][ibasis][inode][1] = atof(values[3]);
  nodex[m][ibasis][inode][2] = atof(values[4]);
  domain->unmap_reverse(nodex[m][ibasis][inode], image[m]);
  nodex_current[m][ibasis][inode][0] = atof(values[5]);
  nodex_current[m][ibasis][inode][1] = atof(values[6]);
  nodex_current[m][ibasis][inode][2] = atof(values[7]);
  domain->unmap_reverse_current(nodex_current[m][ibasis][inode], image[m]);

  // 1 equal to the all group, 2 equal to the element group

  nodemask[m][ibasis][inode] = 1|2;
  nodev[m][ibasis][inode][0] = 0.0;
  nodev[m][ibasis][inode][1] = 0.0;
  nodev[m][ibasis][inode][2] = 0.0;
}

/*  ----------------------------------------------------------------------
    pack data for element I for sending to another proc
    -------------------------------------------------------------------------  */

int ElementVecCAC::pack_exchange(int i, double *buf)
{
  int m = 1;
  int inpe = npe[etype[i]];
  int iapc = apc[etype[i]];
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = initial_box_size[i][0];
  buf[m++] = initial_box_size[i][1];
  buf[m++] = initial_box_size[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(etype[i]).d;
  for (int j = 0; j < iapc; j++) 
    buf[m++] = ubuf(ctype[i][j]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = element_bound_box[i][0];
  buf[m++] = element_bound_box[i][1];
  buf[m++] = element_bound_box[i][2];
  buf[m++] = element_bound_box[i][3];
  buf[m++] = element_bound_box[i][4];
  buf[m++] = element_bound_box[i][5];
  buf[m++] = subelem_size[i];

  for (int j = 0; j < iapc; j++) 
    for (int k = 0; k < inpe; k++) 
      buf[m++] = ubuf(nodemask[i][j][k]).d;

  for (int j = 0; j < iapc; j++) 
    for (int k = 0; k < inpe; k++) {
      buf[m++] = nodex[i][j][k][0];
      buf[m++] = nodex[i][j][k][1];
      buf[m++] = nodex[i][j][k][2];
      if (atom->atom_strain_flag) {
        buf[m++] = nodex_current[i][j][k][0];
        buf[m++] = nodex_current[i][j][k][1];
        buf[m++] = nodex_current[i][j][k][2];
      }
      buf[m++] = nodev[i][j][k][0];
      buf[m++] = nodev[i][j][k][1];
      buf[m++] = nodev[i][j][k][2];
    }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      m += modify->fix[element->extra_grow[iextra]]->pack_elem_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

/*  ---------------------------------------------------------------
    unpack element from exchanged proc
    --------------------------------------------------------------  */

int ElementVecCAC::unpack_exchange(double *buf)
{

  // grow per-element arrays if needed

  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;

  // unpack

  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  initial_box_size[nlocal][0] = buf[m++];
  initial_box_size[nlocal][1] = buf[m++];
  initial_box_size[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  etype[nlocal] = (int) ubuf(buf[m++]).i;
  for (int j = 0; j < apc[etype[nlocal]]; j++) 
    ctype[nlocal][j] = (int) ubuf(buf[m++]).i;
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
  int iapc = apc[etype[nlocal]];
  for (int j = 0; j < iapc; j++) 
    for (int k = 0; k < inpe; k++) 
      nodemask[nlocal][j][k] = (int) ubuf(buf[m++]).i;

  for (int j = 0; j < iapc; j++)
    for (int k = 0; k < inpe; k++) {
      nodex[nlocal][j][k][0] = buf[m++];
      nodex[nlocal][j][k][1] = buf[m++];
      nodex[nlocal][j][k][2] = buf[m++];
      if (atom->atom_strain_flag) {
        nodex_current[nlocal][j][k][0] = buf[m++];
        nodex_current[nlocal][j][k][1] = buf[m++];
        nodex_current[nlocal][j][k][2] = buf[m++];
      }
      nodev[nlocal][j][k][0] = buf[m++];
      nodev[nlocal][j][k][1] = buf[m++];
      nodev[nlocal][j][k][2] = buf[m++];
    }

  if (element->nextra_grow)
    for (int iextra = 0; iextra < element->nextra_grow; iextra++)
      m += modify->fix[element->extra_grow[iextra]]->
        unpack_elem_exchange(nlocal, &buf[m]);

  element->nlocal++;
  return m;
}

/*  ----------------------------------------------------------------------
    pack data for elements in list for border exchange
    -------------------------------------------------------------------------  */

int ElementVecCAC::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, k, l, m, inpe, iapc;
  double dx, dy, dz;
  double dx_current, dy_current, dz_current;
  dx = dy = dz = 0.0;
  dx_current = dy_current = dz_current = 0.0;
  m = 0;

  if (pbc_flag) {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
      dx_current = pbc[0] * domain->xprd_current;
      dy_current = pbc[1] * domain->yprd_current;
      dz_current = pbc[2] * domain->zprd_current;
    } else {
      dx = dx_current = pbc[0];
      dy = dy_current = pbc[1];
      dz = dz_current = pbc[2];
    }
  }
  for (i = 0; i < n; i++) {
    j = list[i];
    inpe = npe[etype[j]];
    iapc = apc[etype[j]];
    buf[m++] = x[j][0] + dx;
    buf[m++] = x[j][1] + dy;
    buf[m++] = x[j][2] + dz;
    buf[m++] = ubuf(tag[j]).d;
    buf[m++] = ubuf(etype[j]).d;
    for (k = 0; k < iapc; k++) 
      buf[m++] = ubuf(ctype[j][k]).d;
    buf[m++] = ubuf(mask[j]).d;
    buf[m++] = element_bound_box[j][0] + dx;
    buf[m++] = element_bound_box[j][1] + dy;
    buf[m++] = element_bound_box[j][2] + dz;
    buf[m++] = element_bound_box[j][3] + dx;
    buf[m++] = element_bound_box[j][4] + dy;
    buf[m++] = element_bound_box[j][5] + dz;
    if (atom->atom_strain_flag)
      buf[m++] = ubuf(image[j]).d;
    for (k = 0; k < iapc; k++) 
      for (l = 0; l < inpe; l++) {
        buf[m++] = ubuf(nodemask[j][k][l]).d;
        buf[m++] = nodex[j][k][l][0] + dx;
        buf[m++] = nodex[j][k][l][1] + dy;
        buf[m++] = nodex[j][k][l][2] + dz;
        if (atom->atom_strain_flag) {
          buf[m++] = nodex_current[j][k][l][0] + dx_current;
          buf[m++] = nodex_current[j][k][l][1] + dy_current;
          buf[m++] = nodex_current[j][k][l][2] + dz_current;
        }
      }
  }

  if (element->nextra_border)
    for (int iextra = 0; iextra < element->nextra_border; iextra++)
      m += modify->fix[element->extra_border[iextra]]->pack_elem_border(n, list, &buf[m]);

  return m;
}

/*  ----------------------------------------------------------------------
    unpack data for elements received from border exchange
    -------------------------------------------------------------------------  */

void ElementVecCAC::unpack_border(int n, int first, double *buf)
{
  int m, last, inpe, iapc;
  m = 0;
  last = first + n;
  for (int i = first; i < last; i++) {
    while (i >= nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    etype[i] = (int) ubuf(buf[m++]).i;
    inpe = npe[etype[i]];
    iapc = apc[etype[i]];
    for (int j = 0; j < iapc; j++) 
      ctype[i][j] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    element_bound_box[i][0] = buf[m++];
    element_bound_box[i][1] = buf[m++];
    element_bound_box[i][2] = buf[m++];
    element_bound_box[i][3] = buf[m++];
    element_bound_box[i][4] = buf[m++];
    element_bound_box[i][5] = buf[m++];
    if (atom->atom_strain_flag) 
      image[i] = (imageint) ubuf(buf[m++]).i;
    for (int j = 0; j < iapc; j++) {
      for (int k = 0; k < inpe; k++) {
        nodemask[i][j][k] = (int) ubuf(buf[m++]).i;
        nodex[i][j][k][0] = buf[m++];
        nodex[i][j][k][1] = buf[m++];
        nodex[i][j][k][2] = buf[m++];
        if (atom->atom_strain_flag) {
          nodex_current[i][j][k][0] = buf[m++];
          nodex_current[i][j][k][1] = buf[m++];
          nodex_current[i][j][k][2] = buf[m++];
        }
      }
    }
  }

  if (element->nextra_border)
    for (int iextra = 0; iextra < element->nextra_border; iextra++)
      m += modify->fix[element->extra_border[iextra]]->
        unpack_elem_border(n, first, &buf[m]);
}

/*  ----------------------------------------------------------------------  */

int ElementVecCAC::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, k, l, inpe, iapc;
  double dx, dy, dz;
  dx = dy = dz = 0.0;

  int m = 0;

  if (pbc_flag) {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
  }
  for (i = 0; i < n; i++) {
    j = list[i];
    inpe = npe[etype[j]];
    iapc = apc[etype[j]];
    buf[m++] = x[j][0] + dx;
    buf[m++] = x[j][1] + dy;
    buf[m++] = x[j][2] + dz;
    buf[m++] = element_bound_box[j][0] + dx;
    buf[m++] = element_bound_box[j][1] + dy;
    buf[m++] = element_bound_box[j][2] + dz;
    buf[m++] = element_bound_box[j][3] + dx;
    buf[m++] = element_bound_box[j][4] + dy;
    buf[m++] = element_bound_box[j][5] + dz;
    for (k = 0; k < inpe; k++) 
      for (l = 0; l < iapc; l++) {
        buf[m++] = nodex[j][k][l][0] + dx;
        buf[m++] = nodex[j][k][l][1] + dy;
        buf[m++] = nodex[j][k][l][2] + dz;
      }
  }

  return m;
}


/*  ----------------------------------------------------------------------  */

void ElementVecCAC::unpack_comm(int n, int first, double *buf)
{
  int i, j, k, m, last;

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
    for (j = 0; j < npe[etype[i]]; j++) 
      for (k = 0; k < apc[etype[i]]; k++) {
        nodex[i][j][k][0] = buf[m++];
        nodex[i][j][k][1] = buf[m++];
        nodex[i][j][k][2] = buf[m++];
      }
  }
}

/*  ----------------------------------------------------------------------  */

int ElementVecCAC::pack_reverse(int n, int first, double *buf)
{
  int i, j, k, m, last, ingcell, iapc;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    ingcell = ngcell[etype[i]];
    iapc = apc[etype[i]];
    for (j = 0; j < iapc; j++) 
      for (k = 0; k < ingcell; k++) {
        buf[m++] = gaussf[i][j][k][0];
        buf[m++] = gaussf[i][j][k][1];
        buf[m++] = gaussf[i][j][k][2];
      }
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

void ElementVecCAC::unpack_reverse(int n, int *list, double *buf)
{
  int i, j, k, l, m, ingcell, iapc;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ingcell = ngcell[etype[i]];
    iapc = apc[etype[i]];
    for (k = 0; k < iapc; k++) 
      for (l = 0; l < ingcell; l++) {
        gaussf[j][k][l][0] += buf[m++];
        gaussf[j][k][l][1] += buf[m++];
        gaussf[j][k][l][2] += buf[m++];
      }
  }
}


/*  ----------------------------------------------------------------------
    add a new element from node coords
    element center coords are calculated from node coords
    -------------------------------------------------------------------------  */

void ElementVecCAC::create_element(double ***nodecoord, int ietype, int *ictype, int itag)
{
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);

  double bound_box[6]; 
  bound_box[0] = BIG;
  bound_box[1] = BIG;
  bound_box[2] = BIG;
  bound_box[3] = -BIG;
  bound_box[4] = -BIG;
  bound_box[5] = -BIG;

  x[nlocal][0] = 0.0;
  x[nlocal][1] = 0.0;
  x[nlocal][2] = 0.0;
  for (int i = 0; i < inpe; i++) {
    for (int j = 0; j < iapc; j++) {
      nodex[nlocal][i][j][0] = nodecoord[i][j][0];
      nodex[nlocal][i][j][1] = nodecoord[i][j][1];
      nodex[nlocal][i][j][2] = nodecoord[i][j][2];
      x[nlocal][0] += nodecoord[i][j][0];
      x[nlocal][1] += nodecoord[i][j][1];
      x[nlocal][2] += nodecoord[i][j][2];
      bound_box[0] = MIN(bound_box[0], nodecoord[i][j][0]);
      bound_box[1] = MIN(bound_box[1], nodecoord[i][j][1]);
      bound_box[2] = MIN(bound_box[2], nodecoord[i][j][2]);
      bound_box[3] = MAX(bound_box[3], nodecoord[i][j][0]);
      bound_box[4] = MAX(bound_box[4], nodecoord[i][j][1]);
      bound_box[5] = MAX(bound_box[5], nodecoord[i][j][2]);
      nodev[nlocal][i][j][0] = 0.0;
      nodev[nlocal][i][j][1] = 0.0;
      nodev[nlocal][i][j][2] = 0.0;
      nodemask[nlocal][i][j] = 1|2;
    }
  }
  x[nlocal][0] /= inpe * iapc;
  x[nlocal][1] /= inpe * iapc;
  x[nlocal][2] /= inpe * iapc;

  etype[nlocal] = ietype;
  for (int i = 0; i < iapc; i++) 
    ctype[nlocal][i] = ictype[i];
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

/*  ----------------------------------------------------------------------
    add a new element from centroid position, node coords are calculated from lattice
    lattice must be defined to call this function
    coord is in lattice coordinate
    -------------------------------------------------------------------------  */

void ElementVecCAC::create_element(double *coord, int ietype, int *ictype, int itag)
{
  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  domain->lattice->lattice2box(x[nlocal][0], x[nlocal][1], x[nlocal][2]);

  double bound_box[6]; 
  bound_box[0] = BIG;
  bound_box[1] = BIG;
  bound_box[2] = BIG;
  bound_box[3] = -BIG;
  bound_box[4] = -BIG;
  bound_box[5] = -BIG;

  if (!domain->lattice)
    error->one(FLERR, "Lattice must be defined to create new elements");

  double xnode[3];
  int nbasis = domain->lattice->nbasis;
  if (apc[ietype] != nbasis) error->one(FLERR, "Number of atoms per unit cell for this element type does not match with lattice's number of basis atoms");
  double **basis = domain->lattice->basis;

  // offset centroid coordinate to corner
  for (int j = 0; j < nbasis; j++) {
    coord[0] -= basis[j][0] / nbasis; 
    coord[1] -= basis[j][1] / nbasis; 
    coord[2] -= basis[j][2] / nbasis; 
  }

  if (element_shape_ids[ietype] == Element::HEXAHEDRON) {
    if (domain->dimension != 3) error->one(FLERR, "Can only create Hex element in 3D");
    int c1, c2, c3;
    c1 = ncells[ietype][0] - 1; 
    c2 = ncells[ietype][1] - 1; 
    c3 = ncells[ietype][2] - 1;

    double d[8][3] = {
      {-1, -1, -1}, 
      { 1, -1, -1}, 
      { 1, 1, -1}, 
      {-1, 1, -1}, 
      {-1, -1, 1}, 
      { 1, -1, 1}, 
      { 1, 1, 1}, 
      {-1, 1, 1}
    };

    for (int k = 0; k < 8; k++) {
      for (int j = 0; j < nbasis; j++) {
        xnode[0] = coord[0] + basis[j][0] + c1 * d[k][0] / 2.0;
        xnode[1] = coord[1] + basis[j][1] + c2 * d[k][1] / 2.0;
        xnode[2] = coord[2] + basis[j][2] + c3 * d[k][2] / 2.0;

        domain->lattice->lattice2box(xnode[0], xnode[1], xnode[2]);

        nodex[nlocal][j][k][0] = xnode[0];
        nodex[nlocal][j][k][1] = xnode[1];
        nodex[nlocal][j][k][2] = xnode[2];
        if (domain->triclinic) 
          domain->x2lamda(xnode, xnode);
        bound_box[0] = MIN(bound_box[0], xnode[0]);
        bound_box[1] = MIN(bound_box[1], xnode[1]);
        bound_box[2] = MIN(bound_box[2], xnode[2]);
        bound_box[3] = MAX(bound_box[3], xnode[0]);
        bound_box[4] = MAX(bound_box[4], xnode[1]);
        bound_box[5] = MAX(bound_box[5], xnode[2]);

        nodev[nlocal][j][k][0] = 0.0;
        nodev[nlocal][j][k][1] = 0.0;
        nodev[nlocal][j][k][2] = 0.0;
        nodemask[nlocal][j][k] = 1 | 2;
      }
    }
  } else if (element_shape_ids[ietype] == Element::QUADRILATERAL) {
    if (domain->dimension != 2) error->one(FLERR, "Can only create Quad element in 2D");
    int c1, c2;
    c1 = ncells[ietype][0]-1; 
    c2 = ncells[ietype][1]-1; 
    double d[4][2] = {
      {-1, -1}, 
      { 1, -1}, 
      { 1, 1}, 
      {-1, 1}
    };

    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < nbasis; j++) {
        xnode[0] = coord[0] + basis[j][0] + c1 * d[k][0] / 2.0;
        xnode[1] = coord[1] + basis[j][1] + c2 * d[k][1] / 2.0;
        xnode[2] = coord[2];

        domain->lattice->lattice2box(xnode[0], xnode[1], xnode[2]);

        nodex[nlocal][j][k][0] = xnode[0];
        nodex[nlocal][j][k][1] = xnode[1];
        nodex[nlocal][j][k][2] = xnode[2];
        if (domain->triclinic) 
          domain->x2lamda(xnode, xnode);
        bound_box[0] = MIN(bound_box[0], xnode[0]);
        bound_box[1] = MIN(bound_box[1], xnode[1]);
        bound_box[2] = MIN(bound_box[2], xnode[2]);
        bound_box[3] = MAX(bound_box[3], xnode[0]);
        bound_box[4] = MAX(bound_box[4], xnode[1]);
        bound_box[5] = MAX(bound_box[5], xnode[2]);
        nodev[nlocal][j][k][0] = 0.0;
        nodev[nlocal][j][k][1] = 0.0;
        nodev[nlocal][j][k][2] = 0.0;
        nodemask[nlocal][j][k] = 1 | 2;
      }
    }
  }
  etype[nlocal] = ietype;
  for (int j = 0; j < nbasis; j++) 
    ctype[nlocal][j] = ictype[j];
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

/* ----------------------------------------------------
   pack element info for Atoms section in data file from write_data_elem command
   ------------------------------------------ */

void ElementVecCAC::pack_element_data_elem(double **buf, tagint tag_offset, int ghostflag)
{
  /*  NEED TO FIX
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
      */
}

/*  ----------------------------------------------------------------------
    pack velocity info for data file from write_data_elem command
    -------------------------------------------------------------------------  */

void ElementVecCAC::pack_element_vel_data_elem(double **buf, tagint tag_offset)
{
  /*  NEED TO FIX
      int m = 0;
      double vel[3];
      tagint mytag;
      int inpe, iapc;
      for (int i = 0; i < element->nlocal; i++) {
      buf[m][0] = ubuf(tag[i] + tag_offset).d;
      buf[m][1] = ubuf(ctype[i]).d;
      vel[0] = vel[1] = vel[2] = 0.0;
      inpe = npe[etype[i]];
      iapc = apc[etype[i]];
      for (int j = 0; j < inpe; j++) 
      for (int k = 0; k < iapc; k++) {
      vel[0] += nodev[i][j][k][0];
      vel[1] += nodev[i][j][k][1];
      vel[2] += nodev[i][j][k][2];
      }
      buf[k][2] = vel[0]/inpe/iapc;
      buf[k][3] = vel[1]/inpe/iapc;
      buf[k][4] = vel[2]/inpe/iapc;
      m++;
      }
      */
}

/*  ----------------------------------------------------
    pack node info for Atoms section in data file from write_data_elem command
    ------------------------------------------ */

void ElementVecCAC::pack_node_data_elem(double **buf, tagint tag_offset, int ghostflag)
{
  /* 
   *NEED TO FIX
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
   */
}

/*  ----------------------------------------------------------------------
    pack node velocity info for data file from write_data_elem command
    -------------------------------------------------------------------------  */

void ElementVecCAC::pack_node_vel_data_elem(double **buf, tagint tag_offset)
{
  /* 
   *NEED TO FIX
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
   */
}

/* ----------------------------------------------------
   write node info into Atoms section in data file from write_data_elem command
   ------------------------------------------ */

void ElementVecCAC::write_data_elem(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, buf[i][2], buf[i][3], buf[i][4]);
}

/*  ----------------------------------------------------------------------
    write node velocity info to data file from write_data_elem command
    -------------------------------------------------------------------------  */

void ElementVecCAC::write_vel_data_elem(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, buf[i][2], buf[i][3], buf[i][4]);
}

/*  ----------------------------------------------------
    write element info into Elements section in data file from write_data command
    ------------------------------------------ */

void ElementVecCAC::write_element_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, buf[i][2], buf[i][3], buf[i][4]);
}

/*  ----------------------------------------------------
    pack element info for Elements section in data file from write_data command
    ------------------------------------------ */

void ElementVecCAC::pack_element_data(double **buf)
{
  for (int i = 0; i < element->nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(etype[i]).d;
    buf[i][2] = x[i][0];
    buf[i][3] = x[i][1];
    buf[i][4] = x[i][2];
  }
}

/*  ----------------------------------------------------
    write node info into Nodes section in data file from write_data command
    ------------------------------------------ */

void ElementVecCAC::write_node_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %d %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, 
        (int) ubuf(buf[i][1]).i, 
        (int) ubuf(buf[i][2]).i, 
        (int) ubuf(buf[i][3]).i, 
        buf[i][4], buf[i][5], buf[i][6]);
}

/*  ----------------------------------------------------
    pack node info for Nodes section in data file from write_data command
    ------------------------------------------ */

void ElementVecCAC::pack_node_data(double **buf)
{
  int m = 0;
  int inpe, iapc;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    iapc = apc[etype[i]];

    for (int j = 0; j < iapc; j++) 
      for (int k = 0; k < inpe; k++) {
        buf[m][0] = ubuf(tag[i]).d;
        buf[m][1] = ubuf(k+1).d;
        buf[m][2] = ubuf(j+1).d;
        buf[m][3] = ubuf(ctype[i][j]).d;
        buf[m][4] = nodex[i][j][k][0];
        buf[m][5] = nodex[i][j][k][1];
        buf[m][6] = nodex[i][j][k][2];
        m++;
      }
  }
}

/*  ----------------------------------------------------------------------
    pack velocity info for data file from write_data command
    -------------------------------------------------------------------------  */

void ElementVecCAC::pack_vel_data(double **buf)
{
  int m = 0;
  int inpe, iapc;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    iapc = apc[etype[i]];

    for (int j = 0; j < iapc; j++) 
      for (int k = 0; k < inpe; k++) {
        buf[m][0] = ubuf(tag[i]).d;
        buf[m][1] = ubuf(k+1).d;
        buf[m][2] = ubuf(j+1).d;
        buf[m][3] = nodev[i][j][k][0];
        buf[m][4] = nodev[i][j][k][1];
        buf[m][5] = nodev[i][j][k][2];
        m++;
      }
  }
}

/*  ----------------------------------------------------------------------
    write velocity info to data file from write_data command
    -------------------------------------------------------------------------  */

void ElementVecCAC::write_vel_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, 
        (int) ubuf(buf[i][1]).i, 
        (int) ubuf(buf[i][2]).i, 
        buf[i][3], buf[i][4], buf[i][5]);
}

/*  ----------------------------------------------------------------------
    write element types information to Element Types section in data file from write_data command
    -------------------------------------------------------------------------  */

void ElementVecCAC::write_element_types(FILE *fp)
{
  for (int i = 1; i <= element->netypes; i++)
    fprintf(fp, "%d %s %d %d %d %d %d %d %d\n", 
        i, element_shape_names[i], apc[i], 
        ncells[i][0], ncells[i][1], ncells[i][2], 
        ngcells[i][0], ngcells[i][1], ngcells[i][2]);
}

/* ----------------------------------------------------
   write node info into Nodes section in dat file from write_tecplot command
   ------------------------------------------ */

void ElementVecCAC::write_node_tecplot(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, "%-1.16e %-1.16e %-1.16e " TAGINT_FORMAT " %d %d\n", 
        buf[i][0], buf[i][1], buf[i][2], (tagint) ubuf(buf[i][3]).i
        , (int) ubuf(buf[i][4]).i, (int) ubuf(buf[i][5]).i);
}

/* ----------------------------------------------------
   write node connectivity info into Nodes section in .dat file from write_tecplot_command
   ------------------------------------------ */

void ElementVecCAC::write_node_connect_tecplot(FILE *fp, int n, int *buf)
{
  int k = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < size_tecplot_node_connect; j++)
      fprintf(fp, TAGINT_FORMAT" ", (tagint) buf[k++]);
    fprintf(fp, "\n");
  }
}

/* ----------------------------------------------------
   pack node info for Nodes section in .dat file from write_tecplot command
   ------------------------------------------ */

void ElementVecCAC::pack_node_tecplot(double **buf, int aveflag)
{
  int inpe, iapc;
  int m = 0;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    iapc = apc[etype[i]];
    if (aveflag) {
      for (int k = 0; k < inpe; k++) {
        buf[m][0] = buf[m][1] = buf[m][2] = 0.0;
        for (int j = 0; j < iapc; j++) {
          buf[m][0] += nodex[i][j][k][0];
          buf[m][1] += nodex[i][j][k][1];
          buf[m][2] += nodex[i][j][k][2];
        }
        buf[m][0] /= iapc;
        buf[m][1] /= iapc;
        buf[m][2] /= iapc;
        buf[m][3] = ubuf(tag[i]).d;
        buf[m][4] = ubuf(etype[i]).d;
        buf[m][5] = ubuf(0).d;
        m++;
      }
    } else {
      for (int j = 0; j < iapc; j++) 
        for (int k = 0; k < inpe; k++) {
          buf[m][0] = nodex[i][j][k][0];
          buf[m][1] = nodex[i][j][k][1];
          buf[m][2] = nodex[i][j][k][2];
          buf[m][3] = ubuf(tag[i]).d;
          buf[m][4] = ubuf(etype[i]).d;
          buf[m][5] = ubuf(ctype[i][j]).d;
          m++;
        }
    }
  }
}

/* ----------------------------------------------------
   pack node info for Nodes section in plt or szplt file from write_tecplot command
   ------------------------------------------ */

void ElementVecCAC::pack_node_tecplot_binary(double *buf, int ival, int aveflag)
{
  int inpe, iapc;
  int m = 0;
  for (int i = 0; i < element->nlocal; i++) {
    inpe = npe[etype[i]];
    iapc = apc[etype[i]];
    if (aveflag) {
      for (int k = 0; k < inpe; k++) {
        buf[m] = 0.0;
        if (ival < 3 && ival >= 0) {
          for (int j = 0; j < iapc; j++) 
            buf[m] += nodex[i][j][k][ival];
          buf[m++] /= iapc;
        } else if (ival == 3) buf[m++] = tag[i];
        else if (ival == 4) buf[m++] = etype[i];
        else if (ival == 5) buf[m++] = 0;
        else error->one(FLERR, "Invalid ival in pack_tecplot_binary");
      }
    } else {
      for (int j = 0; j < iapc; j++) 
        for (int k = 0; k < inpe; k++) {
          if (ival < 3 && ival >= 0) buf[m++] = nodex[i][j][k][ival];
          else if (ival == 3) buf[m++] = tag[i];
          else if (ival == 4) buf[m++] = etype[i];
          else if (ival == 5) buf[m++] = ctype[i][j];
          else error->one(FLERR, "Invalid ival in pack_tecplot_binary");
        }
    }
  }
}

/* ----------------------------------------------------
   pack node connectivity info for Nodes section in dat file
   specifically set for each element shape
   ------------------------------------------ */

void ElementVecCAC::pack_node_connect_tecplot(int *buf, int **nodecell_ids)
{
  int itype;
  int m = 0;
  int *node_connect = new int[size_tecplot_node_connect];
  for (int i = 0; i < element->nlocal; i++) {
    itype = etype[i];
    element->node_connectivity(node_connect, element_shape_ids[itype], nodecell_ids[i]);
    for (int j = 0; j < size_tecplot_node_connect; j++) 
      buf[m++] = node_connect[j];
  }
  delete [] node_connect;
}

/* ----------------------------------------------------
   pack node connectivity info for Nodes section in dat file
   specifically set for each element shape
   ------------------------------------------ */

void ElementVecCAC::pack_node_connect_tecplot(int *buf, int ***node_ids)
{
  int itype;
  int m = 0;
  int *node_connect = new int[size_tecplot_node_connect];
  for (int i = 0; i < element->nlocal; i++) {
    itype = etype[i];
    for (int j = 0; j < apc[itype]; j++) {
      element->node_connectivity(node_connect, element_shape_ids[itype], node_ids[i][j]);
      for (int k = 0; k < size_tecplot_node_connect; k++) 
        buf[m++] = node_connect[k];
    }
  }
  delete [] node_connect;
}

/*  ---------------------------------------------
    Setup element types info (shape, interpolation, integration)
    ------------------------------- */

void ElementVecCAC::set_element_types(const char *str, int type_offset)
{
  int itype, iapc;
  int ncellx, ncelly, ncellz, ngcellx, ngcelly, ngcellz;
  char element_shape_name[256];

  // remove user input for intergration points later and have a fixed scheme

  int n = sscanf(str, "%d %s %d %d %d %d %d %d %d", &itype, element_shape_name, &iapc, &ncellx, &ncelly, &ncellz, &ngcellx, &ngcelly, &ngcellz);
  if (n != 9) error->all(FLERR, "Invalid Element Type line in data file");
  itype += type_offset;

  if (itype < 1 || itype > element->netypes) {
    printf("itype = %d netypes = %d\n", itype, element->netypes);
    error->all(FLERR, "Invalid element type for interpolate set");
  }
  if (iapc < 1) error->all(FLERR, "Invalid apc value for this element type");
  if (iapc > element->maxapc) error->all(FLERR, "apc for this type was greater than maxapc defined, increase maxapc through element_style command");

  // determine element shape id

  int element_shape_id = element->find_element_shape_id(element_shape_name);
  if (element_shape_id == -1) error->all(FLERR, "Invalid element shape name");

  // restrict number of integration to element size, if larger, then reduce to element size

  if (ngcellx > ncellx) {
    char errstr[200];
    sprintf(errstr, "The number of gaussian cells along x direction for "
        "element type #%d has been reduced to match the element size", itype);
    error->warning(FLERR, errstr);
  }
  if (ngcelly > ncelly) {
    char errstr[200];
    sprintf(errstr, "The number of gaussian cells along y direction for "
        "element type #%d has been reduced to match the element size", itype);
    error->warning(FLERR, errstr);
  }
  if (ngcellz > ncellz) {
    char errstr[200];
    sprintf(errstr, "The number of gaussian cells along z direction for "
        "element type #%d has been reduced to match the element size", itype);
    error->warning(FLERR, errstr);
  }


  ngcellx = MIN(ngcellx, ncellx);
  ngcelly = MIN(ngcelly, ncelly);
  ngcellz = MIN(ngcellz, ncellz);

  if (element_type_setflag[itype]) {

    // if identical element type is defined, skip the setup
    // else throw error;

    if (ncellx == ncells[itype][0] &&
        ncelly == ncells[itype][1] &&
        ncellz == ncells[itype][2] && 
        ngcellx == ngcells[itype][0] &&
        ngcelly == ngcells[itype][1] &&
        ngcellz == ngcells[itype][2] &&
        apc[itype] == iapc &&
        element_shape_id == element_shape_ids[itype]) 
      return;
    error->all(FLERR, "This element type has been set");
  }


  apc[itype] = iapc;
  element_shape_ids[itype] = element_shape_id;
  strcpy(element_shape_names[itype], element_shape_name);
  element_type_setflag[itype] = 1;

  ncells[itype][0] = ncellx;
  ncells[itype][1] = ncelly;
  ncells[itype][2] = ncellz;

  ngcells[itype][0] = ngcellx;
  ngcells[itype][1] = ngcelly;
  ngcells[itype][2] = ngcellz;

  // check input error for each element shapes, specify constants for this type
  // add here to include more element shapes

  if (element_shape_id == Element::QUADRILATERAL) {
    //if (ncellz != 1) error->all(FLERR, "Number of cells along third direction must be one for 2D simulation"); 
    if (ncellx < 2 || ncelly < 2) 
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell || 
        ncelly-2 > element->max_edge_ucell) 
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx-2) * (ncelly-2) > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceed max_surface_ucell");
    npe[itype] = 4;
    n = ncellx * ncelly;
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcelly < 2)
      error->all(FLERR, "Invalid integration values");
    else {
      ngcelly = MIN(ngcelly, ncelly);
      n = ngcellx * ngcelly;
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }

    ncells[itype][2] = 1;
    ngcells[itype][2] = 1;

  } else if (element_shape_id == Element::TRIANGLE) {
    error->all(FLERR,"TRIANGLE element is not fully implemented yet"); 
    if (domain->dimension == 3) error->all(FLERR, "2D element requires 2D simulation");
    //if (ncellz != 1) error->all(FLERR, "Number of cells along third direction must be one for 2D simulation"); 
    if (ncellx < 2)
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell) 
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx - 2) * (ncellx - 1)/2 > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceed max_surface_ucell");
    npe[itype] = 3;
    n = ncellx * (ncellx + 1)/2;
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcellx > ncellx)
      error->all(FLERR, "Invalid integration values");
    else {
      n = ngcellx * (ngcellx + 1)/2;
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }

  } else if (element_shape_id == Element::HEXAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR, "3D element requires 3D simulation");
    if (ncellx < 2 || ncelly < 2 || ncellz < 2) {
      if (comm->me == 0) printf("ietype = %d ncell = %d %d %d ngcell = %d %d %d str = %s\n", itype, ncellx, ncelly, ncellz, ngcellx, ngcelly, ngcellz, str);
      error->all(FLERR, "Invalid interpolate values");
    }
    if (ncellx-2 > element->max_edge_ucell ||
        ncelly-2 > element->max_edge_ucell ||
        ncellz-2 > element->max_edge_ucell) {
      if (comm->me == 0) printf("ietype = %d ncell = %d %d %d ngcell = %d %d %d str = %s\n", itype, ncellx, ncelly, ncellz, ngcellx, ngcelly, ngcellz, str);
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    }
    if ((ncellx-2) * (ncelly-2) > element->max_surface_ucell ||
        (ncellz-2) * (ncelly-2) > element->max_surface_ucell ||
        (ncellx-2) * (ncellz-2) > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceeds max_surface_ucell");
    npe[itype] = 8;
    n = ncellx * ncelly * ncellz;
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcelly < 2 || ngcellz < 2)
      error->all(FLERR, "Invalid integration values");
    else {
      n = ngcellx * ngcelly * ngcellz;
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }

  } else if (element_shape_id == Element::TETRAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR, "3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell)
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx-2) * (ncellx-1)/2 > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceeds max_surface_ucell");
    npe[itype] = 4;
    n = sum_tetrahedron(ncellx);
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcellx > 4 || ngcellx > ncellx)
      error->all(FLERR, "Invalid integration values");
    else {
      if (ncellx > 4)
        n = 32;
      else
        n = nucell[itype];
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }
  } else if (element_shape_id == Element::OCTAHEDRON) {
    if (domain->dimension == 2) error->all(FLERR, "3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell) 
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx-2) * (ncellx-1)/2 > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceeds max_surface_ucell");

    npe[itype] = 6;
    n = sum_octahedron(ncellx);
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcellx > 4 ||
        (ngcellx < 4 && ngcellx != ncellx))
      error->all(FLERR, "Invalid integration values");
    else {
      n = sum_octahedron(ngcellx);
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }
  } else if (element_shape_id == Element::PYRAMID) {
    if (domain->dimension == 2) error->all(FLERR, "3D element requires 3D simulation");
    if (ncellx < 2)
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell) 
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx-2) * (ncellx-2) > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceeds max_surface_ucell");
    npe[itype] = 5;
    n = sum_pyramid(ncellx);
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcellx < 2 || ngcellx > 4 ||
        (ngcellx < 4 && ngcellx != ncellx))
      error->all(FLERR, "Invalid integration values");
    else {

      if (ncellx > 4)
        n = 42;
      else n = nucell[itype];
      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }
  } else if (element_shape_id == Element::WEDGE) {
    if (domain->dimension == 2) error->all(FLERR, "3D element requires 3D simulation");

    npe[itype] = 6;

    if (ncellx < 2 || ncelly < 2)
      error->all(FLERR, "Invalid interpolate values");
    if (ncellx-2 > element->max_edge_ucell || ncelly-2 > element->max_edge_ucell) 
      error->all(FLERR, "Number of virtual atoms on an edge exceeds max_edge_ucell");
    if ((ncellx-2) * (ncelly-2) > element->max_surface_ucell ||
        (ncelly-1) * (ncelly-2)/2 > element->max_surface_ucell)
      error->all(FLERR, "Number of virtual atoms on a surface exceeds max_surface_ucell");

    n = ncellx * (ncelly + 1) * ncelly/2;
    if (n > element->maxucell) error->all(FLERR, "Too many virtual atoms per element, boost maxucell");
    nucell[itype] = n;

    if (ngcellx == 0 && ngcelly == 0 && ngcellz == 0) ngcell[itype] = 0;
    else if (ngcelly < 2 || ngcelly > 4 ||
        (ngcelly < 4 && ngcelly != ncelly) ||
        ngcellx < 2 || ngcellx > ncellx)
      error->all(FLERR, "Invalid integration values");
    else {
      if (ncelly > 4)
        n = 12 * ngcellx;
      else 
        n = (ngcelly + 1) * ngcelly * ngcellx/2;

      if (n > element->maxgcell) error->all(FLERR, "Too many gaussian cells per element, boost maxgcell");
      ngcell[itype] = n;
    }
  }

  nodal_weight[itype] = (static_cast<double> (nucell[itype])) 
    / (static_cast<double> (npe[itype]));
  element->max_nucell = MAX(element->max_nucell, nucell[itype]);
  element->max_ngcell = MAX(element->max_ngcell, ngcell[itype]);

  // setup interpolation, sub-element division, and integration scheme for this element type
  setup_interpolate_element(itype);
  if (ngcell[itype]) setup_integration_point(itype);
  setup_sub_element(itype);
}

/* -----------------------------------------------------------------------------------
   Calculate unit cell related arrays for an element type
   - Shape function array to interpolate unit cells
   - List of element nodes in natural coordinate system
   + Quadrilateral
   - 1: (-1, -1)
   - 2: (1, -1)
   - 3: (1, 1)
   - 4: (-1, 1)
   + Triangle
   - 1: (0, 0)
   - 2: (1, 0)
   - 3: (0, 1)
   + Hexahedron
   - 1: (-1, -1, -1)
   - 2: (1, -1, -1)
   - 3: (1, 1, -1)
   - 4: (-1, 1, -1)
   - 5: (-1, -1, 1)
   - 6: (1, -1, 1)
   - 7: (1, 1, 1)
   - 8: (-1, 1, 1)
   + Tetrahedron
   - 1: (0, 0, 0)
   - 2: (1, 0, 0)
   - 3: (0, 1, 0)
   - 4: (0, 0, 1)
   + Pyramid
   - 1: (-1, -1, 0)
   - 2: (1, -1, 0)
   - 3: (1, 1, 0)
   - 4: (-1, 1, 0)
   - 5: (0, 0, 1)
   + Wedge
   - 1: (-1, 0, 0)
   - 2: (-1, 1, 0)
   - 3: (-1, 0, 1)
   - 4: (1, 0, 0)
   - 5: (1, 1, 0)
   - 6: (1, 0, 1)

   - Lists of indices of virtual atoms on edges and surfaces
   ------------------------------------------------------------------------------------ */

void ElementVecCAC::setup_interpolate_element(int itype)
{
  int ncellx = ncells[itype][0];
  int ncelly = ncells[itype][1];
  int ncellz = ncells[itype][2];

  // virtual atom index counter;

  int n = 0;

  // outer virtual atom index counter;

  int outer = 1;

  // calculate shape function arrays based on its shape and size
  // add here to include other element shapes

  if (element_shape_ids[itype] == Element::QUADRILATERAL) {

    int nedge[4];
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = 0;
    int nsur = 0;

    double px, py;
    double dx_pc, dy_pc;
    dx_pc = 2.0/(ncellx - 1);
    dy_pc = 2.0/(ncelly - 1);
    for (int i = 0; i < ncellx; i++) {
      px = dx_pc * i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dy_pc * j - 1.0;
        shape_array[itype][n][0] = 0.25 * (1 - px) * (1 - py);
        shape_array[itype][n][1] = 0.25 * (1 + px) * (1 - py);
        shape_array[itype][n][2] = 0.25 * (1 + px) * (1 + py);
        shape_array[itype][n][3] = 0.25 * (1 - px) * (1 + py);

        if (domain->dimension == 3) {
          is_outer[itype][n] = outer++;
        } else {
          if (i == 0 || j == 0 || i == ncellx -1 || j == ncelly - 1)
            is_outer[itype][n] = outer++;
          else is_outer[itype][n] = 0;
        }

        // store indices of unit cells on edges and surfaces for this element type       

        if (i == 0 && j != 0 && j != ncelly - 1)
          edge_ucell[itype][0][nedge[0]++] = n; 
        if (i == ncellx - 1 && j != 0 && j != ncelly - 1)
          edge_ucell[itype][1][nedge[1]++] = n; 
        if (j == 0 && i != 0 && i != ncellx - 1)
          edge_ucell[itype][2][nedge[2]++] = n; 
        if (j == ncelly - 1 && i != 0 && i != ncellx - 1)
          edge_ucell[itype][3][nedge[3]++] = n; 
        if (i != 0 && i != ncellx - 1 && j != 0 && j != ncelly - 1)
          surface_ucell[itype][0][nsur++] = n;

        n++;
      }
    }
    for (int i = 0; i < 4; i++)
      nedge_ucell[itype][i] = nedge[i];
    nsurface_ucell[itype][0] = nsur;


  } else if (element_shape_ids[itype] == Element::TRIANGLE) {
    int nedge[3];
    double px, py;
    nedge[0] = nedge[1] = nedge[2] = 0;
    double dxy_pc = 1.0/(ncellx - 1);

    for (int i = 0; i < ncellx; i++) {
      px = dxy_pc * i;
      for (int j = 0; j < ncellx - i; j++) {
        py = dxy_pc * j;
        shape_array[itype][n][0] = px;
        shape_array[itype][n][1] = py;
        shape_array[itype][n][2] = 1 - px - py;

        if (i == 0 || j == 0 || i + j == ncellx- 1)
          is_outer[itype][n] = outer++;
        else is_outer[itype][n] = 0;

        // store indices of unit cells on edges and surfaces for this element type       

        if (i == 0 && j != 0 && j != ncelly - 1)
          edge_ucell[itype][0][nedge[0]++] = n; 
        if (j == 0 && i != 0 && i != ncellx - 1)
          edge_ucell[itype][1][nedge[1]++] = n; 
        if (i+j == ncellx - 1 && i != 0 && j != 0)
          edge_ucell[itype][2][nedge[2]++] = n; 

        n++;
      }
    }
    for (int i = 0; i < 3; i++)
      nedge_ucell[itype][i] = nedge[i];


  } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

    int nsur[6], nedge[12];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = nsur[5] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = nedge[9] = nedge[10] = nedge[11] = 0;
    double px, py, pz;
    double dx_pc, dy_pc, dz_pc;
    dx_pc = 2.0/(ncellx - 1);
    dy_pc = 2.0/(ncelly - 1);
    dz_pc = 2.0/(ncellz - 1);

    for (int i = 0; i < ncellx; i++) {
      px = dx_pc * i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dy_pc * j - 1.0;
        for (int k = 0; k < ncellz; k++) {
          pz = dz_pc * k - 1.0;
          shape_array[itype][n][0] = 0.125 * (1 - px) * (1 - py) * (1 - pz);
          shape_array[itype][n][1] = 0.125 * (1 + px) * (1 - py) * (1 - pz);
          shape_array[itype][n][2] = 0.125 * (1 + px) * (1 + py) * (1 - pz);
          shape_array[itype][n][3] = 0.125 * (1 - px) * (1 + py) * (1 - pz);
          shape_array[itype][n][4] = 0.125 * (1 - px) * (1 - py) * (1 + pz);
          shape_array[itype][n][5] = 0.125 * (1 + px) * (1 - py) * (1 + pz);
          shape_array[itype][n][6] = 0.125 * (1 + px) * (1 + py) * (1 + pz);
          shape_array[itype][n][7] = 0.125 * (1 - px) * (1 + py) * (1 + pz);

          if (i == 0 || i == ncellx -1 ||
              j == 0 || j == ncelly -1 ||
              k == 0 || k == ncellz -1)
            is_outer[itype][n] = outer++;
          else is_outer[itype][n] = 0;

          // store indices of unit cells on edges and surfaces for this element type       

          if (k != 0 && k != ncellz-1) {
            if (i == 0 && j == 0) {
              edge_ucell[itype][0][nedge[0]++] = n; 
            } else if (i == 0 && j == ncelly-1) {
              edge_ucell[itype][1][nedge[1]++] = n; 
            } else if (i == ncellx-1 && j == 0) {
              edge_ucell[itype][2][nedge[2]++] = n; 
            } else if (i == ncellx-1 && j == ncelly-1) {
              edge_ucell[itype][3][nedge[3]++] = n; 
            } else if (i == 0) {
              surface_ucell[itype][0][nsur[0]++] = n;
            } else if (i == ncellx-1) {
              surface_ucell[itype][1][nsur[1]++] = n;
            } else if (j == 0) {
              surface_ucell[itype][2][nsur[2]++] = n;
            } else if (j == ncelly-1) {
              surface_ucell[itype][3][nsur[3]++] = n;
            }
          } else if (k == 0) {
            if (j != 0 && j != ncelly-1) {
              if (i == 0) {
                edge_ucell[itype][4][nedge[4]++] = n; 
              } else if (i == ncellx-1) {
                edge_ucell[itype][6][nedge[6]++] = n; 
              } else {
                surface_ucell[itype][4][nsur[4]++] = n;
              }
            } else if (i != 0 && i != ncellx-1) {
              if (j == 0) {
                edge_ucell[itype][8][nedge[8]++] = n; 
              } else {
                edge_ucell[itype][10][nedge[10]++] = n; 
              }
            }
          } else {
            if (j != 0 && j != ncelly-1) {
              if (i == 0) {
                edge_ucell[itype][5][nedge[5]++] = n; 
              } else if (i == ncellx-1) {
                edge_ucell[itype][7][nedge[7]++] = n; 
              } else {
                surface_ucell[itype][5][nsur[5]++] = n;
              }
            } else if (i != 0 && i != ncellx-1) {
              if (j == 0) {
                edge_ucell[itype][9][nedge[9]++] = n; 
              } else {
                edge_ucell[itype][11][nedge[11]++] = n; 
              }
            }
          }
          n++;
        }
      }
    }
    for (int i = 0; i < 6; i++)
      nsurface_ucell[itype][i] = nsur[i];
    for (int i = 0; i < 12; i++)
      nedge_ucell[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

    int nsur[4], nedge[6];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = 0;
    nedge[0] = nedge[1] = nedge[2] = 0;
    nedge[3] = nedge[4] = nedge[5] = 0;
    double px, py, pz;
    double dxyz_pc;
    int xnum, ynum;
    dxyz_pc = 1.0/(ncellx-1);

    for (int k = 0; k < ncellx; k++) {
      pz = dxyz_pc * k;
      xnum = ncellx - k;
      for (int i = 0; i < xnum; i++) {
        px = dxyz_pc * i;
        ynum = xnum - i;
        for (int j = 0; j < ynum; j++) {
          py = dxyz_pc * j;
          shape_array[itype][n][0] = 1 - px - py - pz;
          shape_array[itype][n][1] = px;
          shape_array[itype][n][2] = py;
          shape_array[itype][n][3] = pz;

          if (i == 0 || j == 0 || k == 0 || i + j + k == ncellx- 1)
            is_outer[itype][n] = outer++;
          else is_outer[itype][n] = 0;

          // store indices of unit cells on edges and surfaces for this element type       

          if (i != 0 && j != 0 && k == 0 && i+j+k != ncellx-1)          // xy plane
            surface_ucell[itype][0][nsur[0]++] = n;
          if (i == 0 && j != 0 && k != 0 && i+j+k != ncellx-1)          // yz plane
            surface_ucell[itype][1][nsur[1]++] = n;
          if (i != 0 && j == 0 && k != 0 && i+j+k != ncellx-1)          // xz plane
            surface_ucell[itype][2][nsur[2]++] = n;
          if (i != 0 && j != 0 && k != 0 && i+j+k == ncellx-1)          // x+y+z = 1 plane
            surface_ucell[itype][3][nsur[3]++] = n;

          if (i == 0 && j == 0 && k != 0 && k != ncellx-1)              // z axis
            edge_ucell[itype][0][nedge[0]++] = n;
          if (i == 0 && k == 0 && j != 0 && j != ncellx-1)              // y axis
            edge_ucell[itype][1][nedge[1]++] = n;
          if (j == 0 && k == 0 && i != 0 && i != ncellx-1)              // x axis
            edge_ucell[itype][2][nedge[2]++] = n;
          if (i == 0 && j+k == ncellx-1 && j != 0 && k != 0)            // y+z = 1 line
            edge_ucell[itype][3][nedge[3]++] = n;
          if (j == 0 && i+k == ncellx-1 && i != 0 && k != 0)            // x+z = 1 line
            edge_ucell[itype][4][nedge[4]++] = n;
          if (k == 0 && i+j == ncellx-1 && i != 0 && j != 0)            // x+y = 1 line
            edge_ucell[itype][5][nedge[5]++] = n;

          n++;
        }
      }
    }

    for (int i = 0; i < 4; i++) 
      nsurface_ucell[itype][i] = nsur[i];
    for (int i = 0; i < 6; i++) 
      nedge_ucell[itype][i] = nedge[i];

  } else if (element_shape_ids[itype] == Element::OCTAHEDRON) {

    error->all(FLERR, "Shape function for octahedron element is not found yet");
    int nsur[8], nedge[12];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = nsur[5] = nsur[6] = nsur[7] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = nedge[9] = nedge[10] = nedge[11] = 0;
    double px, py, pz;
    double dxy_pc, dz_pc;
    double xycorner;
    int xynum;
    dxy_pc = 2.0/(ncellx - 1);
    dz_pc = 1.0/(ncellx - 1);

    for (int k = 1 - ncellx; k <= ncellx - 1; k++) {
      pz = dz_pc * k;
      xynum = ncellx - abs(k);            // size of current xy plane
      xycorner = -(xynum - 1) * dz_pc;
      for (int i = 0; i < xynum; i++) {
        px = dxy_pc * i + xycorner;
        for (int j = 0; j < xynum; j++) {
          py = dxy_pc * j + xycorner;

          // THIS IS INCORRECT, NEED TO FIND NEW SHAPE FUNCTION
          shape_array[itype][n][0] = 0.5 * pz * (pz-1);
          shape_array[itype][n][1] = 0.25 * (1-pz * pz) * (1-px) * (1-py);
          shape_array[itype][n][2] = 0.25 * (1-pz * pz) * (1+px) * (1-py);
          shape_array[itype][n][3] = 0.25 * (1-pz * pz) * (1+px) * (1+py);
          shape_array[itype][n][4] = 0.25 * (1-pz * pz) * (1-px) * (1+py);
          shape_array[itype][n][5] = 0.5 * pz * (pz+1);

          // store indices of unit cells on edges and surfaces for this element type       

          if (k > 1 - ncellx && k < 0) {          // lower half
            if (i == 0 && j == 0)                       // on edge
              edge_ucell[itype][0][nedge[0]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_ucell[itype][1][nedge[1]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_ucell[itype][2][nedge[2]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_ucell[itype][3][nedge[3]++] = n;
            else if (i == 0)                            // on surface
              surface_ucell[itype][0][nsur[0]++] = n;
            else if (j == 0)                            // on surface
              surface_ucell[itype][1][nsur[1]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_ucell[itype][2][nsur[2]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_ucell[itype][3][nsur[3]++] = n;
          } else if (k == 0) {                    // middle plane
            if (!(i == 0 && j == 0) &&       
                !(i == xynum - 1 && j == 0) &&
                !(i == 0 && j == xynum - 1) &&
                !(i == xynum - 1 && j == xynum - 1)) {  // exclude 4 corners
              if (i == 0)                               // on edge
                edge_ucell[itype][4][nedge[4]++] = n;
              else if (j == 0)                          // on edge
                edge_ucell[itype][5][nedge[5]++] = n;
              else if (i == xynum - 1)                  // on edge
                edge_ucell[itype][6][nedge[6]++] = n;
              else if (j == xynum - 1)                  // on edge
                edge_ucell[itype][7][nedge[7]++] = n;
            }
          } else if (k > 0 && k < ncellx - 1) {   // upper half
            if (i == 0 && j == 0)                       // on edge
              edge_ucell[itype][8][nedge[8]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_ucell[itype][9][nedge[9]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_ucell[itype][10][nedge[10]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_ucell[itype][11][nedge[11]++] = n;
            else if (i == 0)                            // on surface
              surface_ucell[itype][4][nsur[4]++] = n;
            else if (j == 0)                            // on surface
              surface_ucell[itype][5][nsur[5]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_ucell[itype][6][nsur[6]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_ucell[itype][7][nsur[7]++] = n;
          }
          n++;
        }
      }
    }

    for (int i = 0; i < 8; i++) 
      nsurface_ucell[itype][i] = nsur[i];
    for (int i = 0; i < 12; i++) 
      nedge_ucell[itype][i] = nedge[i];


  } else if (element_shape_ids[itype] == Element::PYRAMID) {

    int nsur[5], nedge[8];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = 0;
    nedge[0] = nedge[1] = nedge[2] = nedge[3] = 0;
    nedge[4] = nedge[5] = nedge[6] = nedge[7] = 0;
    double px, py, pz;
    double dxy_pc, dz_pc;
    double xycorner;
    int xynum;
    dxy_pc = 2.0/(ncellx-1);
    dz_pc = 1.0/(ncellx-1);

    for (int k = 0; k < ncellx-1; k++) {
      pz = dz_pc * k;
      xynum = ncellx - k;                // size of current xy plane
      xycorner = -(xynum - 1) * dz_pc;
      for (int i = 0; i < xynum; i++) {
        px = dxy_pc * i + xycorner;
        for (int j = 0; j < xynum; j++) {
          py = dxy_pc * j + xycorner;
          shape_array[itype][n][0] = 0.25 * (1.0 - px - pz) * (1.0 - py - pz)/(1 - pz);
          shape_array[itype][n][1] = 0.25 * (1.0 + px - pz) * (1.0 - py - pz)/(1 - pz);
          shape_array[itype][n][2] = 0.25 * (1.0 + px - pz) * (1.0 + py - pz)/(1 - pz);
          shape_array[itype][n][3] = 0.25 * (1.0 - px - pz) * (1.0 + py - pz)/(1 - pz);
          shape_array[itype][n][4] = pz;

          if (k == 0 || i == 0 || i == xynum - 1 || j == 0 || j == xynum - 1)
            is_outer[itype][n] = outer++;
          else is_outer[itype][n] = 0;

          // store indices of unit cells on edges and surfaces for this element type       

          if (k == 0) {                                 // base plane
            if (!(i == 0 && j == 0) &&       
                !(i == xynum - 1 && j == 0) &&
                !(i == 0 && j == xynum - 1) &&
                !(i == xynum - 1 && j == xynum - 1)) {  // exclude 4 corners
              if (i == 0)                               // on edge
                edge_ucell[itype][0][nedge[0]++] = n;
              else if (j == 0)                          // on edge
                edge_ucell[itype][1][nedge[1]++] = n;
              else if (i == xynum - 1)                  // on edge
                edge_ucell[itype][2][nedge[2]++] = n;
              else if (j == xynum - 1)                  // on edge
                edge_ucell[itype][3][nedge[3]++] = n;
              else                                      // on surface
                surface_ucell[itype][0][nsur[0]++] = n;
            }
          } else if (k < ncellx - 1) {                  
            if (i == 0 && j == 0)                       // on edge
              edge_ucell[itype][4][nedge[4]++] = n;
            else if (i == xynum - 1 && j == 0)          // on edge
              edge_ucell[itype][5][nedge[5]++] = n;
            else if (i == 0 && j == xynum - 1)          // on edge
              edge_ucell[itype][6][nedge[6]++] = n;
            else if (i == xynum - 1 && j == xynum - 1)  // on edge
              edge_ucell[itype][7][nedge[7]++] = n;
            else if (i == 0)                            // on surface
              surface_ucell[itype][1][nsur[1]++] = n;
            else if (j == 0)                            // on surface
              surface_ucell[itype][2][nsur[2]++] = n;
            else if (i == xynum - 1)                    // on surface
              surface_ucell[itype][3][nsur[3]++] = n;
            else if (j == xynum - 1)                    // on surface
              surface_ucell[itype][4][nsur[4]++] = n;
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
    is_outer[itype][n] = outer;
    n++;

    for (int i = 0; i < 5; i++) 
      nsurface_ucell[itype][i] = nsur[i];
    for (int i = 0; i < 8; i++) 
      nedge_ucell[itype][i] = nedge[i];



  } else if (element_shape_ids[itype] == Element::WEDGE) {

    int nsur[5], nedge[9];
    nsur[0] = nsur[1] = nsur[2] = nsur[3] = nsur[4] = 0;
    nedge[0] = nedge[1] = nedge[2] = 0;
    nedge[3] = nedge[4] = nedge[5] = 0;
    nedge[6] = nedge[7] = nedge[8] = 0;
    double px, py, pz;
    int znum;
    double dx_pc = 2.0/(ncellx-1);
    double dyz_pc = 1.0/(ncelly-1);

    for (int i = 0; i < ncellx; i++) {
      px = dx_pc * i - 1.0;
      for (int j = 0; j < ncelly; j++) {
        py = dyz_pc * j;
        znum = ncelly - j;
        for (int k = 0; k < znum; k++) {
          pz = dyz_pc * k;

          shape_array[itype][n][0] = 0.5 * (1.0 - px) * (1.0 - py - pz);
          shape_array[itype][n][1] = 0.5 * (1.0 - px) * py;
          shape_array[itype][n][2] = 0.5 * (1.0 - px) * pz;
          shape_array[itype][n][3] = 0.5 * (1.0 + px) * (1.0 - py - pz);
          shape_array[itype][n][4] = 0.5 * (1.0 + px) * py;
          shape_array[itype][n][5] = 0.5 * (1.0 + px) * pz;

          if (j == 0 || k == 0 || j + k == ncelly- 1 ||
              i == 0 || i == ncellx - 1)
            is_outer[itype][n] = outer++;
          else is_outer[itype][n] = 0;

          // store indices of unit cells on edges and surfaces for this element type       

          if (i == ncellx-1 && j != 0 && k != 0 && j+k != ncelly-1)             // +yz plane
            surface_ucell[itype][0][nsur[0]++] = n;
          if (i == 0 && j != 0 && k != 0 && j+k != ncelly-1)                    // -yz plane
            surface_ucell[itype][1][nsur[1]++] = n;
          if (i != 0 && i != ncellx-1 && j != 0 && k == 0 && j+k != ncelly-1)   // xy plane
            surface_ucell[itype][2][nsur[2]++] = n;
          if (i != 0 && i != ncellx-1 && j == 0 && k != 0 && j+k != ncelly-1)   // xz plane
            surface_ucell[itype][3][nsur[3]++] = n;
          if (i != 0 && i != ncellx-1 && j != 0 && k != 0 && j+k == ncelly-1)   // y + z = 1 plane
            surface_ucell[itype][4][nsur[4]++] = n;

          if (i == ncellx-1 && j == 0 && k != 0 && k != ncelly-1)               // x = 1, y = 0 line
            edge_ucell[itype][0][nedge[0]++] = n;
          if (i == ncellx-1 && k == 0 && j != 0 && j != ncelly-1)               // x = 1, z = 0 line
            edge_ucell[itype][1][nedge[1]++] = n;
          if (i == ncellx-1 && j+k == ncelly-1 && j != 0 && k != 0)             // x = 1, y+z = 1 line
            edge_ucell[itype][2][nedge[2]++] = n;
          if (i == 0 && j == 0 && k != 0 && k != ncelly-1)                      // x = -1, y = 0 line
            edge_ucell[itype][3][nedge[3]++] = n;
          if (i == 0 && k == 0 && j != 0 && j != ncelly-1)                      // x = -1, z = 0 line
            edge_ucell[itype][4][nedge[4]++] = n;
          if (i == 0 && j+k == ncelly-1 && j != 0 && k != 0)                    // x = -1, y+z = 1 line
            edge_ucell[itype][5][nedge[5]++] = n;
          if (i != 0 && i != ncellx-1 && j == 0 && k == 0)                      // y = 0, z = 0 line
            edge_ucell[itype][6][nedge[6]++] = n;
          if (i != 0 && i != ncellx-1 && j == ncelly-1)                         // y = 0 line
            edge_ucell[itype][7][nedge[7]++] = n;
          if (i != 0 && i != ncellx-1 && k == ncelly-1)                         // z = 0 line
            edge_ucell[itype][8][nedge[8]++] = n;

          n++;
        }
      }
    }

    for (int i = 0; i < 5; i++) 
      nsurface_ucell[itype][i] = nsur[i];
    for (int i = 0; i < 9; i++) 
      nedge_ucell[itype][i] = nedge[i];

  } else if (0) {
    // add new interpolation scheme for new element shape here
  }

  if (n != nucell[itype]) error->all(FLERR, "TEST"); 

  // check is_outer arrays
  for (int i = 0; i < element->nsurface[element->element_shape_ids[itype]]; i++) {
    for (int j = 0; j < nsurface_ucell[itype][i]; j++) {
      if (is_outer[itype][surface_ucell[itype][i][j]] == 0) {
        printf("i = %d j = %d nsurface = %d nsurface_ucell = %d surface_ucell = %d \n", i, j, element->nsurface[itype], nsurface_ucell[itype][i], surface_ucell[itype][i][j]);
        printf("@@@@@@@@@ Error is_outer arrays for etype %d %d %d %d %d\n", itype, apc[itype], ncellx, ncelly, ncellz);
        error->all(FLERR, "TEST"); 
      }
    }
  }

  // setup node to unit cell index mapping 

  for (int i = 0; i < n; i++)
    for (int j = 0; j < npe[itype]; j++) 
      if (fabs(shape_array[itype][i][j]-1.0) < EPSILON)
        n2u[itype][j] = i;

}

/*  -----------------------------------------------------------------------------------
    Calculate gaussian cell related arrays:
    - Weight of gaussian cells
    - g2u to convert gaussian cell index 
    to virtual atom index in an element, 
    used for interpolating gaussian cells
    - n2g to convert node index to gaussian cell index
    Must be done after setting up interpolate element
    ------------------------------------------------------------------------------------ */

void ElementVecCAC::setup_integration_point(int itype)
{
  int ngcellx, ngcelly, ngcellz;
  int ncellx, ncelly, ncellz;
  int iucell;

  int igcell = 0;                      // gaussian cell index counter

  double **xgauss = element->xgauss;
  double **wgauss = element->wgauss;
  int weight_scale_flag = element->weight_scale_flag;
  double snode = element->weight_scale[0];
  double sedge = element->weight_scale[1];
  double sface = element->weight_scale[2];
  double sinner = element->weight_scale[3];

  // number of virtual atoms along each edge

  ncellx = ncells[itype][0];
  ncelly = ncells[itype][1];
  ncellz = ncells[itype][2];

  // number of gaussian cells along each edge

  ngcellx = ngcells[itype][0];
  ngcelly = ngcells[itype][1];
  ngcellz = ngcells[itype][2];

  // determine the indices of gaussian cells and their weights for this element type

  if (element_shape_ids[itype] == Element::QUADRILATERAL) {

    // calculate weight fraction for gaussian cells

    int a, b, w;
    int px, py;
    double w_edge_x;
    double w_edge_y;
    double w_inner;
    double scale = ncellx * ncelly;
    w_edge_x = static_cast<double> (ncellx - 2)/2.0;
    w_edge_y = static_cast<double> (ncelly - 2)/2.0;
    w_inner = w_edge_x * w_edge_y;
    for (int i = 0; i < ngcellx; i++) {
      if (i == 0) px = 0;
      else if (i == ngcellx-1) px = ncellx-1;
      else px = (int) ((1.0 + xgauss[ngcellx-2][i])/2.0 * (ncellx-2)+1.0);
      for (int j = 0; j < ngcelly; j++) {
        if (j == 0) py = 0;
        else if (j == ngcelly-1) py = ncelly-1;
        else py = (int) ((1.0 + xgauss[ngcelly-2][j])/2.0 * (ncelly-2)+1.0);

        a = (i == 0 || i == (ngcellx - 1));
        b = (j == 0 || j == (ngcelly - 1));
        w = a + b;

        // gaussian cell is a node

        if (w == 2) weight[itype][igcell] = 1;

        // gaussian cell is on an edge

        else if (w == 1) {

          // gaussian cell is on edge in x direction

          if (!a) weight[itype][igcell] = w_edge_x * wgauss[ngcellx-2][i];

          // gaussian cell is on edge in y direction

          else if (!b) weight[itype][igcell] = w_edge_y * wgauss[ngcelly-2][j];

          // gaussian cell is inside element

        } else weight[itype][igcell] = w_inner * wgauss[ngcellx-2][i] * wgauss[ngcelly-2][j];


        if (element->debug_mode) debug_gcell_type[itype][igcell] = w;

        // list to convert gaussian cell index to 
        // virtual atom index

        iucell = px * ncelly + py;
        g2u[itype][igcell] = iucell;
        u2g[itype][iucell] = igcell;

        igcell++;
      }
    }

  } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

    // calculate weight fraction for gaussian cells

    int a, b, c, w;
    int px, py, pz;
    double w_edge_x, w_face_xy;
    double w_edge_y, w_face_yz;
    double w_edge_z, w_face_xz;
    double w_inner, w_node;

    //printf("scale_flag = %d\n", weight_scale_flag);

    if (weight_scale_flag == Element::UNIFORM) {

      int nedge;
      int nface;
      int ninner;
      if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 0) {
        nedge = 24;
        nface = 24;
        ninner = 8;
      } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 1) {
        nedge = 16;
        nface = 8;
        ninner = 0;
      } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 2) {
        nedge = 8;
        nface = 0;
        ninner = 0;
      } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 3) {
        nedge = 0;
        nface = 0;
        ninner = 0;
      }
      w_node = static_cast<double> (nucell[itype])/(snode * 8+sedge * nedge+sface * nface+sinner * ninner) * snode;
      //printf("w_node = %g nedge = %d nface = %d ninner = %d s = %g %g %g %g\n", w_node, nedge, nface, ninner, snode, sedge, sface, sinner);
      //printf("ngcell = %d %d %d condition = %d\n", ngcellx, ngcelly, ngcellz, ngcellx==2, ngcelly==2, ngcellz==2);
      w_edge_x = w_edge_y = w_edge_z = w_node * sedge / snode;
      w_face_xy = w_face_yz = w_face_xz = w_node * sface / snode;
      w_inner = w_node * sinner / snode;
    } else if (weight_scale_flag == Element::NONUNIFORM) {
      error->all(FLERR, "Non-uniform weigth scale style not supported yet");
      /* 
         int nedge;
         int nface;
         int ninner;
         if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 0) {
         nedge = 24;
         nface = 24;
         ninner = 8;
         } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 1) {
         nedge = 16;
         nface = 8;
         ninner = 0;
         } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 2) {
         nedge = 8;
         nface = 0;
         ninner = 0;
         } else if (((ngcellx == 2) + (ngcelly == 2) + (ngcellz == 2)) == 3) {
         nedge = 0;
         nface = 0;
         ninner = 0;
         }
         w_node = static_cast<double> (nucell[itype])/(snode * 8+sedge * nedge+sface * nface+sinner * ninner) * snode;
         printf("w_node = %g nedge = %d nface = %d ninner = %d s = %g %g %g %g\n", w_node, nedge, nface, ninner, snode, sedge, sface, sinner);
         printf("ngcell = %d %d %d condition = %d\n", ngcellx, ngcelly, ngcellz, ngcellx==2, ngcelly==2, ngcellz==2);
         w_edge_x = w_edge_y = w_edge_z = w_node * sedge / snode;
         w_face_xy = w_face_yz = w_face_xz = w_node * sface / snode;
         w_inner = w_node * sinner / snode;
         */ 

    } else {
      w_edge_x = static_cast<double> (ncellx - 2)/2.0;
      w_edge_y = static_cast<double> (ncelly - 2)/2.0;
      w_edge_z = static_cast<double> (ncellz - 2)/2.0;
      w_face_xy = w_edge_x * w_edge_y;
      w_face_yz = w_edge_y * w_edge_z;
      w_face_xz = w_edge_x * w_edge_z;
      w_inner = w_edge_x * w_edge_y * w_edge_z;
      w_node = 1.0;
    }

    for (int i = 0; i < ngcellx; i++) {
      if (i == 0) px = 0;
      else if (i == ngcellx-1) px = ncellx-1;
      else px = (int) ((1.0 + xgauss[ngcellx-2][i])/2.0 * (ncellx-2)+1.0);
      for (int j = 0; j < ngcelly; j++) {
        if (j == 0) py = 0;
        else if (j == ngcelly-1) py = ncelly-1;
        else py = (int) ((1.0 + xgauss[ngcelly-2][j])/2.0 * (ncelly-2)+1.0);
        for (int k = 0; k < ngcellz; k++) {
          if (k == 0) pz = 0;
          else if (k == ngcellz-1) pz = ncellz-1;
          else pz = (int) ((1.0 + xgauss[ngcellz-2][k])/2.0 * (ncellz-2)+1.0);

          // calculate weight of gaussian cells

          //printf("P = %d %d %d\n", px, py, pz);
          a = (i == 0 || i == (ngcellx - 1));
          b = (j == 0 || j == (ngcelly - 1));
          c = (k == 0 || k == (ngcellz - 1));
          w = a + b + c;

          // gaussian cell is a node

          if (w == 3) weight[itype][igcell] = w_node;

          // gaussian cell is on an edge

          else if (w == 2) {

            // gaussian cell is on edge in x direction

            if (!a) weight[itype][igcell] = w_edge_x * wgauss[ngcellx-2][i];

            // gaussian cell is on edge in y direction

            else if (!b) weight[itype][igcell] = w_edge_y * wgauss[ngcelly-2][j];

            // gaussian cell is on edge in z direction

            else if (!c) weight[itype][igcell] = w_edge_z * wgauss[ngcellz-2][k];

            // gaussian cell is on a face

          } else if (w == 1) {

            // gaussian cell is on yz face 

            if (a) weight[itype][igcell] = w_face_yz * wgauss[ngcelly-2][j] * wgauss[ngcellz-2][k];

            // gaussian cell is on xz face 

            else if (b) weight[itype][igcell] = w_face_xz * wgauss[ngcellx-2][i] * wgauss[ngcellz-2][k];

            // gaussian cell is on xy face 

            else if (c) weight[itype][igcell] = w_face_xy * wgauss[ngcellx-2][i] * wgauss[ngcelly-2][j];

            // gaussian cell is inside element

          } else {
            weight[itype][igcell] = w_inner * wgauss[ngcellx-2][i] * wgauss[ngcelly-2][j] * wgauss[ngcellz-2][k];
          }

          if (element->debug_mode) debug_gcell_type[itype][igcell] = w;

          // list to convert gaussian cell index to 
          // virtual atom index

          iucell = px * ncelly * ncellz + py * ncellz + pz;
          g2u[itype][igcell] = iucell;
          u2g[itype][iucell] = igcell;
          igcell++;
        }
      }
    }

  } else if (element_shape_ids[itype] == Element::OCTAHEDRON) {

    // calculate weight fraction for gaussian cells
    // for element size <= 4 unit cells, all atoms are gaussian cells

    if (ncellx <= 4) {
      for (igcell = 0; igcell < nucell[itype]; igcell++) {

        // virtual atom index
        iucell = igcell; 
        g2u[itype][igcell] = igcell;
        u2g[itype][igcell] = igcell;
        weight[itype][igcell] = 1.0;
      }
    } else {
      int xynum;
      double dxy_pc = 2.0/(ncellx-1);
      double dz_pc = 1.0/(ncellx-1);
      double w_inner = static_cast<double> (sum_octahedron(ncellx - 1))/6.0;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iucell = 0;


      // position of gaussian cells on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0 * (ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0 * (ncellx - 2) + 1.0);

      // position of gaussian cells on surface

      int isur = (int) (1.0/6.0 * (ncellx - 4.0) + 1);

      for (int k = 1 - ncellx; k <= ncellx - 1; k++) {
        xynum = ncellx - abs(k);            // size of current xy plane
        for (int i = 0; i < xynum; i++) {
          for (int j = 0; j < xynum; j++) {

            // nodal gaussian cell

            if (k == 1 - ncellx || k == ncellx - 1 || 
                (k == 0 && (i == 0 || i == xynum - 1) && (j == 0 || j == xynum -1))) {
              weight[itype][igcell] = 1.0;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              igcell++;
            }

            // edge gaussian cell 

            if (k != 0 && k != 1 - ncellx && k != ncellx - 1) {     // upper and lower half edges
              if ((abs(k) == i1 || abs(k) == i2) &&
                  (i == 0 || i == xynum - 1) &&
                  (j == 0 || j == xynum - 1)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                igcell++;
              }
            } else if (k == 0) {                                    // middle plane edges
              if ((i == 0 || i == xynum - 1) && (j == i1 || j == i2) ||
                  (j == 0 || j == xynum - 1) && (i == i1 || i == i2)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                igcell++;
              }
            }

            // surface gaussian cell

            if (abs(k) == ncellx - 1 - isur * 2) {
              if (((i == 0 || i == xynum - 1) && j == (xynum-1)/2) ||
                  ((j == 0 || j == xynum - 1) && i == (xynum-1)/2)) {
                weight[itype][igcell] = w_face;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                igcell++;
              }
            } else if (abs(k) == isur) {
              if (((i == 0 || i == xynum - 1) && (j == isur || j == xynum - 1 - isur)) ||
                  ((j == 0 || j == xynum - 1) && (i == isur || i == xynum - 1 - isur))) {
                weight[itype][igcell] = w_face;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                igcell++;
              }
            }

            // inner gaussian cells

            if ((k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) ||                 // on either inner half
                (abs(k) == ncellx - 1 - i1 * 2 && i == (xynum-1)/2 && j == (xynum-1)/2)) {    // on middle plane
              weight[itype][igcell] = w_inner;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              igcell++;
            }

            iucell++;
          }
        }
      }
    }

  } else if (element_shape_ids[itype] == Element::PYRAMID) {

    // for element size <= 4 unit cells, all atoms are gaussian cells

    if (ncellx <= 4) {
      for (igcell = 0; igcell < nucell[itype]; igcell++) {

        // virtual atom index
        iucell = igcell; 
        g2u[itype][igcell] = igcell;
        u2g[itype][igcell] = igcell;
        weight[itype][igcell] = 1.0;

      }
    } else {

      double px, py, pz;
      double xycorner;
      int xynum;
      double dxy_pc = 2.0/(ncellx - 1.0);
      double dz_pc = 1.0/(ncellx - 1.0);
      double w_inner_total = static_cast<double> (sum_pyramid(ncellx - 3));
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face_tri = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      double w_face_quad = w_edge * w_edge;
      iucell = 0;

      // position of gaussian cells on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0 * (ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0 * (ncellx - 2) + 1.0);

      // position of gaussian cells on surface

      int isur = (int) ((ncellx - 4.0)/6.0 + 1.0);

      // for inner gaussian cells, locate the closest atoms in natural coords
      // parameters from Liping Liu - Numerical Integration over Pyramids
      // remap the coordinates to the inner pyramid

      // position parameters of lower 4 inner gaussian cells

      double z1 = (35.0 - sqrt(140.0))/140.0;
      z1 = (z1 * (ncellx - 4.0) + 1.0)/(ncellx - 1.0);
      double a = sqrt(5.0/21.0);
      a = a * (ncellx - 4.0)/(ncellx - 1.0);
      double w1 = 0.21;

      // position parameter of upper inner gaussian cells

      double z0 = (70.0 + 21.0 * sqrt(35.0))/280.0;
      z0 = (z0 * (ncellx - 4.0) + 1.0)/(ncellx - 1.0);
      double w0 = 0.16;

      // list of minimum distance and indexes of closest atoms for 5 inner gaussian cells

      double mindistsq[5];
      int minucell[5];
      mindistsq[0] = mindistsq[1] = mindistsq[2] = mindistsq[3] = mindistsq[4] = BIG; 
      double rsq;

      for (int k = 0; k < ncellx; k++) {
        xynum = ncellx - k;            // size of current xy plane
        pz = k * dz_pc;
        xycorner = -(xynum - 1) * dz_pc;
        for (int i = 0; i < xynum; i++) {
          px = dxy_pc * i + xycorner;
          for (int j = 0; j < xynum; j++) {
            py = dxy_pc * j + xycorner;

            // nodal gaussian cell

            if (k == ncellx - 1 || 
                (k == 0 && (i == 0 || i == xynum - 1) && (j == 0 || j == xynum -1))) {
              weight[itype][igcell] = 1.0;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              if (element->debug_mode) debug_gcell_type[itype][igcell] = 3;
              igcell++;
            }

            // edge gaussian cell 

            if (k != 0 && k != ncellx - 1) {     // upper edges
              if ((k == i1 || k == i2) &&
                  (i == 0 || i == xynum - 1) &&
                  (j == 0 || j == xynum - 1)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            } else if (k == 0) {                 // base plane edges
              if ((i == 0 || i == xynum - 1) && (j == i1 || j == i2) ||
                  (j == 0 || j == xynum - 1) && (i == i1 || i == i2)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            }

            // surface gaussian cell

            if (k == ncellx - 1 - isur * 2) {
              if (((i == 0 || i == xynum - 1) && j == (xynum-1)/2) ||
                  ((j == 0 || j == xynum - 1) && i == (xynum-1)/2)) {
                weight[itype][igcell] = w_face_tri;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            } else if (k == isur) {
              if (((i == 0 || i == xynum - 1) && (j == isur || j == xynum - 1 - isur)) ||
                  ((j == 0 || j == xynum - 1) && (i == isur || i == xynum - 1 - isur))) {
                weight[itype][igcell] = w_face_tri;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            } else if (k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) {
              weight[itype][igcell] = w_face_quad;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
              igcell++;
            }

            // inner gaussian cells

            // upper inner gaussian cell

            rsq = px * px + py * py + (pz - z0) * (pz - z0);
            if (rsq < mindistsq[0]) {
              mindistsq[0] = rsq;
              minucell[0] = iucell; 
            }

            // lower 4 inner gaussian cells

            rsq = (px + a) * (px + a) + (py + a) * (py + a) + (pz - z1) * (pz - z1);    // (-a, -a, z1)
            if (rsq < mindistsq[1]) {
              mindistsq[1] = rsq;
              minucell[1] = iucell; 
            }
            rsq = (px - a) * (px - a) + (py + a) * (py + a) + (pz - z1) * (pz - z1);    // (+a, -a, z1)
            if (rsq < mindistsq[2]) {
              mindistsq[2] = rsq;
              minucell[2] = iucell; 
            }
            rsq = (px - a) * (px - a) + (py - a) * (py - a) + (pz - z1) * (pz - z1);    // (+a, +a, z1)
            if (rsq < mindistsq[3]) {
              mindistsq[3] = rsq;
              minucell[3] = iucell; 
            }
            rsq = (px + a) * (px + a) + (py - a) * (py - a) + (pz - z1) * (pz - z1);    // (-a, +a, z1)
            if (rsq < mindistsq[4]) {
              mindistsq[4] = rsq;
              minucell[4] = iucell; 
            }

            iucell++;
          }
        }
      }

      // save indices and weights of 5 inner gaussian cells

      iucell = minucell[0];
      weight[itype][igcell] = w_inner_total * w0;
      g2u[itype][igcell] = iucell;
      u2g[itype][iucell] = igcell;
      if (element->debug_mode) debug_gcell_type[itype][igcell] = 0;
      igcell++;

      for (int i = 1; i < 5; i++) {
        iucell = minucell[i];
        weight[itype][igcell] = w_inner_total * w1;
        g2u[itype][igcell] = iucell;
        u2g[itype][iucell] = igcell;
        if (element->debug_mode) debug_gcell_type[itype][igcell] = 0;
        igcell++;
      }
    }

  } else if (element_shape_ids[itype] == Element::TRIANGLE) {

    // for element size <= 4 unit cells, all atoms are gaussian cells

    if (ncellx <= 4) {
      for (igcell = 0; igcell < nucell[itype]; igcell++) {

        // virtual atom index

        iucell = igcell; 
        g2u[itype][igcell] = igcell;
        u2g[itype][igcell] = igcell;
        weight[itype][igcell] = 1.0;
      }
    } else {

      int ynum;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iucell = 0;

      // position of gaussian cells on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0 * (ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0 * (ncellx - 2) + 1.0);

      // position of gaussian cells on surface

      int isur = (int) (1.0/6.0 * (ncellx - 4.0) + 1.0);


      for (int i = 0; i < ncellx; i++) {
        ynum = ncellx - i;
        for (int j = 0; j < ynum; j++) {

          // nodal gaussian cells

          if (i == ncellx - 1 ||
              j == ncellx - 1 ||
              (i == 0 && j == 0)) {
            weight[itype][igcell] = 1.0;
            g2u[itype][igcell] = iucell;
            u2g[itype][iucell] = igcell;
            if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
            igcell++;
          }

          // edge gaussian cells

          if (i == 0 && (j == i1 || j == i2) ||
              (j == 0 || j == ynum - 1) && (i == i1 || i == i2)) {
            weight[itype][igcell] = w_edge;
            g2u[itype][igcell] = iucell;
            u2g[itype][iucell] = igcell;
            if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
            igcell++;
          }

          // surface gaussian cell

          if ((i == isur || i == i2) && (j == i1 || j == i2)) {
            weight[itype][igcell] = w_face;
            g2u[itype][igcell] = iucell;
            u2g[itype][iucell] = igcell;
            if (element->debug_mode) debug_gcell_type[itype][igcell] = 0;
            igcell++;
          }

          iucell++;
        }
      }
    }

  } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

    // for element size <= 4 unit cells, all atoms are gaussian cells

    if (ncellx <= 4) {
      for (igcell = 0; igcell < nucell[itype]; igcell++) {

        // virtual atom index
        iucell = igcell; 
        g2u[itype][igcell] = igcell;
        u2g[itype][igcell] = igcell;
        weight[itype][igcell] = 1.0;
      }
    } else {

      int xnum, ynum;
      double px, py, pz;
      double dxyz_pc = 1.0/(ncellx-1);
      double w_inner = static_cast<double> (sum_tetrahedron(ncellx - 4))/4.0;
      double w_edge = static_cast<double> (ncellx - 2)/2.0;
      double w_face = static_cast<double> (sum_triangle(ncellx - 3))/3.0;
      iucell = 0;

      // position of gaussian cells on edge

      int i1 = (int) ((1.0 + xgauss[2][1])/2.0 * (ncellx - 2) + 1.0);
      int i2 = (int) ((1.0 + xgauss[2][2])/2.0 * (ncellx - 2) + 1.0);

      // position of gaussian cells on surface

      int isur = (int) (1.0/6.0 * (ncellx - 4.0) + 1.0);


      // list of minimum distance and indexes of closest atoms for 4 inner gaussian cells

      double mindistsq[4];
      int minucell[4];
      mindistsq[0] = mindistsq[1] = mindistsq[2] = mindistsq[3] = BIG; 
      double rsq;

      // coordinates of inner gaussian cells
      // map coordinates to the inner tetrahedron domain

      double a = 0.1381966011250110;
      a = (a * (ncellx - 5.0) + 1.0)/(ncellx - 1.0);
      double b = 0.5854101966249680;
      b = (b * (ncellx - 5.0) + 1.0)/(ncellx - 1.0);

      for (int k = 0; k < ncellx; k++) {
        xnum = ncellx - k;            // size of current xy plane
        pz = k * dxyz_pc;
        for (int i = 0; i < xnum; i++) {
          ynum = xnum - i;
          px = dxyz_pc * i;
          for (int j = 0; j < ynum; j++) {
            py = dxyz_pc * j;

            // nodal gaussian cells

            if (k == ncellx - 1 || 
                i == ncellx - 1 ||
                j == ncellx - 1 ||
                (k == 0 && i == 0 && j == 0)) {
              weight[itype][igcell] = 1.0;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              if (element->debug_mode) debug_gcell_type[itype][igcell] = 3;
              igcell++;
            }

            // edge gaussian cells

            if (k != 0 && k != ncellx - 1) {     // upper edges
              if ((k == i1 || k == i2) &&
                  (i == 0 || i == xnum - 1) &&
                  (j == 0 || j == ynum - 1)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            } else if (k == 0) {                 // base plane edges
              if (i == 0 && (j == i1 || j == i2) ||
                  (j == 0 || j == ynum - 1) && (i == i1 || i == i2)) {
                weight[itype][igcell] = w_edge;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            }

            // surface gaussian cell

            if (k == ncellx - 1 - isur * 2) {
              if ((i == 0 && j == (xnum-1)/2) ||
                  (j == 0 && i == (xnum-1)/2) ||
                  (j == (ynum-1) && i == j)) {
                weight[itype][igcell] = w_face;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            } else if (k == isur) {
              if ((i == 0 && (j == isur || j == ynum - 1 - isur)) ||
                  ((j == 0 || j == ynum - 1) && (i == isur || i == xnum - 1 - isur))) {
                weight[itype][igcell] = w_face;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                igcell++;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
              }
            } else if (k == 0 && (i == i1 || i == i2) && (j == i1 || j == i2)) {
              weight[itype][igcell] = w_face;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
              igcell++;
            }


            // inner gaussian cells

            rsq = (px - a) * (px - a) + (py - a) * (py - a) + (pz - a) * (pz - a);    // (a, a, a)
            if (rsq < mindistsq[0]) {
              mindistsq[0] = rsq;
              minucell[0] = iucell; 
            }
            rsq = (px - b) * (px - b) + (py - a) * (py - a) + (pz - a) * (pz - a);    // (b, a, a)
            if (rsq < mindistsq[1]) {
              mindistsq[1] = rsq;
              minucell[1] = iucell; 
            }
            rsq = (px - a) * (px - a) + (py - b) * (py - b) + (pz - a) * (pz - a);    // (a, b, a)
            if (rsq < mindistsq[2]) {
              mindistsq[2] = rsq;
              minucell[2] = iucell; 
            }
            rsq = (px - a) * (px - a) + (py - a) * (py - a) + (pz - b) * (pz - b);    // (a, a, b)
            if (rsq < mindistsq[3]) {
              mindistsq[3] = rsq;
              minucell[3] = iucell; 
            }

            iucell++;
          }
        }
      }

      // save indices and weights of 4 inner gaussian cells

      for (int i = 0; i < 4; i++) {
        iucell = minucell[i];
        weight[itype][igcell] = w_inner;
        g2u[itype][igcell] = iucell;
        u2g[itype][iucell] = igcell;
        if (element->debug_mode) debug_gcell_type[itype][igcell] = 0;
        igcell++;
      }
    }

  } else if (element_shape_ids[itype] == Element::WEDGE) {

    // position of gaussian cells on edge

    int i1x = (int) ((1.0 + xgauss[2][1])/2.0 * (ncellx - 2) + 1.0);
    int i2x = (int) ((1.0 + xgauss[2][2])/2.0 * (ncellx - 2) + 1.0);
    int i1y = (int) ((1.0 + xgauss[2][1])/2.0 * (ncelly - 2) + 1.0);
    int i2y = (int) ((1.0 + xgauss[2][2])/2.0 * (ncelly - 2) + 1.0);

    // for element size in yz <= 4 unit cells, all atoms on yz plane are gaussian cells

    if (ncelly <= 4) {
      int znum;
      double dx_pc = 2.0/(ncellx-1);
      double dyz_pc = 1.0/(ncelly-1);
      double w_edge_x = (ncellx-2.0)/(ngcellx-2.0);

      iucell = 0;
      for (int i = 0; i < ncellx; i++) {
        for (int j = 0; j < ncelly; j++) {
          znum = ncelly - j;
          for (int k = 0; k < znum; k++) {
            if (i == 0 || i == ncellx-1) {
              weight[itype][igcell] = 1;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              igcell++;
            } else if (ncellx >= 3 && (i == i1x || i == i2x)) {
              weight[itype][igcell] = w_edge_x;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              igcell++;
            }
            iucell++;
          }
        }
      }
    } else {

      int znum;
      double px, py, pz;
      double dxyz_pc = 1.0/(ncellx-1);
      double w_inner = ((ncellx - 2.0) * (ncelly - 3.0) * (ncelly - 2.0))
        /(ngcellx - 2.0)/6.0;
      double w_edge_x = (ncellx - 2.0)/(ngcellx - 2.0);
      double w_edge_yz = (ncelly - 2.0)/2.0;
      double w_face_tri = static_cast<double> (sum_triangle(ncelly - 3))/3.0;
      double w_face_quad = w_edge_x * w_edge_yz;

      iucell = 0;

      // position of gaussian cells on surface

      int isur = (int) (1.0/6.0 * (ncelly - 4.0) + 1.0);

      for (int i = 0; i < ncellx; i++) {
        for (int j = 0; j < ncelly; j++) {
          znum = ncelly - j;
          for (int k = 0; k < znum; k++) {

            // nodal gaussian cells

            if ((i == ncellx - 1 || i == 0) &&
                (j == ncelly - 1 || j == 0) &&
                (k == ncelly - 1 || k == 0)) {
              weight[itype][igcell] = 1.0;
              g2u[itype][igcell] = iucell;
              u2g[itype][iucell] = igcell;
              if (element->debug_mode) debug_gcell_type[itype][igcell] = 3;
              igcell++;
            }

            // edge gaussian cells

            if (ncellx >= 3 && (i == i1x || i == i2x)) { // x edge
              if ((j == 0 && k == 0) ||
                  j == ncelly-1 || k == ncelly-1) {
                weight[itype][igcell] = w_edge_x;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            } else if (i == 0 || i == ncellx-1) {        // y, z edges
              if ((j == 0 && (k == i1y || k == i2y)) ||
                  ((j == i1y || j == i2y) && (k == 0 || k == znum-1))) {
                weight[itype][igcell] = w_edge_yz;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 2;
                igcell++;
              }
            }

            // surface gaussian cell

            if (i == ncellx - 1 || i == 0) {             // x = +-1 planes
              if ((j == isur && k == isur) ||
                  (j == isur && k == znum-1-isur) ||
                  (j == ncelly-1-2 * isur && k == znum/2)) {
                weight[itype][igcell] = w_face_tri;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            } else if (k == 0 || k == znum-1) {          // z = 0, y + z = 1 planes
              if (ncellx >= 3 && (i == i1x || i == i2x) && (j == i1y || j == i2y)) { 
                weight[itype][igcell] = w_face_quad;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            } else if (j == 0) {                         // y = 0 plane
              if (ncellx >= 3 && (i == i1x || i == i2x) && (k == i1y || k == i2y)) { 
                weight[itype][igcell] = w_face_quad;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 1;
                igcell++;
              }
            }

            // inner gaussian cells

            if (ncellx >= 3 && (i == i1x || i == i2x)) {      
              if ((j == isur && k == isur) ||
                  (j == isur && k == znum-1-isur) ||
                  (j == ncelly-1-2 * isur && k == znum/2)) {
                weight[itype][igcell] = w_inner;
                g2u[itype][igcell] = iucell;
                u2g[itype][iucell] = igcell;
                if (element->debug_mode) debug_gcell_type[itype][igcell] = 0;
                igcell++;
              }
            }

            iucell++;
          }
        }
      }

    }
  } else {

    // add new integration scheme for new element shape here

  }

  // Error check

  if (igcell != ngcell[itype]) {
    printf("itype = %d igcell = %d ngcell = %d\n", itype, igcell, ngcell[itype]);
    error->all(FLERR, "Wrong integration setup");
  }

  // setup weighted shape array and node to gaussian cell mapping

  for (igcell = 0; igcell < ngcell[itype]; igcell++) {
    iucell = g2u[itype][igcell];
    for (int inode = 0; inode < npe[itype]; inode++) {

      // weighted shape array to distribute force from gaussian cells to nodes

      if (element->mass_style == Element::LUMPED)
        weighted_shape_array[itype][igcell][inode] = 
          shape_array[itype][iucell][inode] * weight[itype][igcell] / nodal_weight[itype];
      else 
        weighted_shape_array[itype][igcell][inode] = 
          shape_array[itype][iucell][inode] * weight[itype][igcell] / nucell[itype];

      // list to convert node index to gaussian cell index

      if (fabs(shape_array[itype][iucell][inode]-1.0) < EPSILON) {
        if (n2g[itype][inode] < 0) {
          n2g[itype][inode] = igcell;
          g2n[itype][igcell] = inode;
        } else error->all(FLERR, "Node to gaussian cell mapping is not correct");

      }
    }
  }

  // check if all node to gaussian cell mapping is set

  for (int inode = 0; inode < npe[itype]; inode++) {
    if (n2g[itype][inode] < 0) 
      error->all(FLERR, "Node to gaussian cell mapping is not complete");
    if (g2n[itype][n2g[itype][inode]] < 0)
      error->all(FLERR, "Integration point to node mapping is not complete");
  }

  // check if weights are assigned correctly

  double total_weight = 0.0;
  for (int l = 0; l < ngcell[itype]; l++)  {
    total_weight += weight[itype][l];
    //printf("weight[%d] = %g\n", l, weight[itype][l]);
  }

  if (fabs(total_weight - nucell[itype]) > EPSILON) {
    printf("total_weight = %g nucell = %d ncells = %d %d %d ngcells = %d %d %d\n", total_weight, nucell[itype], ncellx, ncelly, ncellz, ngcellx, ngcelly, ngcellz);
    error->all(FLERR, "Integration weight is not assigned correctly");
  }
}

/* -----------------------------------------------------------------------------------
   Calculate sub-element related arrays for an element type
   - us2u to map atom index inside sub-element to unit cell index inside parent element
   - nucell_subelem: number of unit cells inside each sub element
   ------------------------------------------------------------------------------------ */

void ElementVecCAC::setup_sub_element(int itype)
{
  int m, n, ntotal;
  int isub = 0;
  int ncellx, ncelly, ncellz;

  ncellx = ncells[itype][0];
  ncelly = ncells[itype][1];
  ncellz = ncells[itype][2];

  nsubelem[itype] = 0;
  int iapc = apc[itype];

  // total virtual atoms counter

  ntotal = 0;

  bigint checksum = 0;

  // makes sure atoms on element boundaries are included
  // atoms might sit right on sub-element boundaries
  // so subtract a small number from the boundaries to make sure
  // all atoms are assigned to only one sub-element
  if (element_shape_ids[itype] == Element::QUADRILATERAL) {

    double px, py;
    double xlo, xhi, ylo, yhi;
    int esplit = round(sqrt(sqrt(ncellx * ncelly)));
    subsplit[itype] = esplit;
    double ds = 2.0/esplit;
    double dx_pc = 2.0/(ncellx-1);
    double dy_pc = 2.0/(ncelly-1);

    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix * ds - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds;

      for (int iy = 0; iy < esplit; iy++) {

        if (nsubelem[itype] >= element->maxsubelem) error->all(FLERR, "Too many sub-elements, boost maxsubelem through element_modify command");

        ylo = -1.0 + iy * ds - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds;

        // calculate sub-element center shape function array

        px = (ix + 0.5) * ds - 1.0;
        py = (iy + 0.5) * ds - 1.0;

        shape_array_center_subelem[itype][isub][0] = 0.25 * (1 - px) * (1 - py) / iapc;
        shape_array_center_subelem[itype][isub][1] = 0.25 * (1 + px) * (1 - py) / iapc;
        shape_array_center_subelem[itype][isub][2] = 0.25 * (1 + px) * (1 + py) / iapc;
        shape_array_center_subelem[itype][isub][3] = 0.25 * (1 - px) * (1 + py) / iapc;

        // place virtual atoms into sub element isub

        n = 0; // number of atoms inside a sub element counter
        m = 0; // index of virtual atom in element

        for (int i = 0; i < ncellx; i++) {
          px = dx_pc * i - 1.0;
          for (int j = 0; j < ncelly; j++) {
            py = dy_pc * j - 1.0;
            if (px >= xlo && px < xhi &&
                py >= ylo && py < yhi) {
              if (n >= element->maxsubucell) error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
              us2u[itype][isub][n++] = m;
              checksum += m;
            }
            m++;
          }
        }

        // don't store info if sub-element contains no atom

        if (n == 0) continue;

        nsubelem[itype]++;
        nucell_subelem[itype][isub++] = n;
        ntotal += n;
      }
    }

  } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

    double px, py, pz;
    double dx_pc, dy_pc, dz_pc;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    int esplit = round(sqrt(cbrt(ncellx * ncelly * ncellz)));
    subsplit[itype] = esplit;
    double ds = 2.0/esplit;
    dx_pc = 2.0/(ncellx-1);
    dy_pc = 2.0/(ncelly-1);
    dz_pc = 2.0/(ncellz-1);
    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 + ix * ds - EPSILON;
      if (ix == (esplit-1)) xhi = 1.0 + EPSILON;
      else xhi = xlo + ds;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = -1.0 + iy * ds - EPSILON;
        if (iy == (esplit-1)) yhi = 1.0 + EPSILON;
        else yhi = ylo + ds;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = -1.0 + iz * ds - EPSILON;
          if (iz == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds;

          if (nsubelem[itype] >= element->maxsubelem) error->all(FLERR, "Too many sub-elements, boost maxsubelem through element_modify command");

          // calculate sub-element center shape function array

          px = (ix + 0.5) * ds - 1.0;
          py = (iy + 0.5) * ds - 1.0;
          pz = (iz + 0.5) * ds - 1.0;

          shape_array_center_subelem[itype][isub][0] = 0.125 * (1 - px) * (1 - py) * (1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][1] = 0.125 * (1 + px) * (1 - py) * (1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][2] = 0.125 * (1 + px) * (1 + py) * (1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][3] = 0.125 * (1 - px) * (1 + py) * (1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][4] = 0.125 * (1 - px) * (1 - py) * (1 + pz) / iapc;
          shape_array_center_subelem[itype][isub][5] = 0.125 * (1 + px) * (1 - py) * (1 + pz) / iapc;
          shape_array_center_subelem[itype][isub][6] = 0.125 * (1 + px) * (1 + py) * (1 + pz) / iapc;
          shape_array_center_subelem[itype][isub][7] = 0.125 * (1 - px) * (1 + py) * (1 + pz) / iapc;

          // place virtual atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of virtual atom in element

          for (int i = 0; i < ncellx; i++) {
            px = dx_pc * i - 1.0;
            for (int j = 0; j < ncelly; j++) {
              py = dy_pc * j - 1.0;
              for (int k = 0; k < ncellz; k++) {
                pz = dz_pc * k - 1.0;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  if (n >= element->maxsubucell) error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
                  us2u[itype][isub][n++] = m;
                  checksum += m;
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          nsubelem[itype]++;

          nucell_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }

  } else if (element_shape_ids[itype] == Element::PYRAMID) {
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double px, py, pz, px_center, py_center, pz_center;
    double xycorner;
    int xynum;
    int esplit = round(sqrt(ncellx));
    subsplit[itype] = esplit;
    double dxy_pc = 2.0/(ncellx-1);
    double dz_pc = 1.0/(ncellx-1);
    double ds_xy = (2.0+dxy_pc)/esplit;
    double ds_z = (1.0+dz_pc)/esplit;
    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 - dxy_pc/2 + ix * ds_xy;
      xhi = -1.0 - dxy_pc/2 + (ix+1) * ds_xy;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = -1.0 - dxy_pc/2 + iy * ds_xy;
        yhi = -1.0 - dxy_pc/2 + (iy+1) * ds_xy;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = -dz_pc/2 + iz * ds_z;
          zhi = -dz_pc/2 + (iz+1) * ds_z;

          if (nsubelem[itype] >= element->maxsubelem) error->all(FLERR, "Too many sub-elements, boost maxsubelem through element_modify command");


          // place virtual atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of virtual atom in element

          px_center = py_center = pz_center = 0.0;

          for (int k = 0; k < ncellx; k++) {
            pz = dz_pc * k;
            xynum = ncellx - k;            // size of current xy plane
            xycorner = -(xynum - 1) * dz_pc;
            for (int i = 0; i < xynum; i++) {
              px = dxy_pc * i + xycorner;
              for (int j = 0; j < xynum; j++) {
                py = dxy_pc * j + xycorner;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  if (n >= element->maxsubucell) error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
                  us2u[itype][isub][n++] = m;
                  checksum += m;
                  px_center += px;
                  py_center += py;
                  pz_center += pz;
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          // calculate sub-element center shape function array
          // sub-element center is centroid of atoms

          px = px_center/n;
          py = py_center/n;
          pz = pz_center/n;

          shape_array_center_subelem[itype][isub][0] = 0.25 * (1.0 - px - pz) * (1.0 - py - pz)/(1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][1] = 0.25 * (1.0 + px - pz) * (1.0 - py - pz)/(1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][2] = 0.25 * (1.0 + px - pz) * (1.0 + py - pz)/(1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][3] = 0.25 * (1.0 - px - pz) * (1.0 + py - pz)/(1 - pz) / iapc;
          shape_array_center_subelem[itype][isub][4] = pz / iapc;

          nsubelem[itype]++;
          nucell_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }
  } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double px, py, pz, px_center, py_center, pz_center;
    int xnum, ynum;
    int esplit = round(sqrt(ncellx));
    subsplit[itype] = esplit;
    double dxyz_pc = 1.0/(ncellx-1);
    double ds_xyz = (1.0+dxyz_pc)/esplit;

    for (int ix = 0; ix < esplit; ix++) {
      xlo = ix * ds_xyz - dxyz_pc/2;
      xhi = (ix+1) * ds_xyz - dxyz_pc/2;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = iy * ds_xyz - dxyz_pc/2;
        yhi = (iy+1) * ds_xyz - dxyz_pc/2;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = iz * ds_xyz - dxyz_pc/2;
          zhi = (iz+1) * ds_xyz - dxyz_pc/2;

          if (nsubelem[itype] >= element->maxsubelem) error->all(FLERR, "Too many sub-elements, boost maxsubelem through element_modify command");


          // place virtual atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of virtual atom in element

          px_center = py_center = pz_center = 0.0;

          for (int k = 0; k < ncellx; k++) {
            pz = dxyz_pc * k;
            xnum = ncellx - k;
            for (int i = 0; i < xnum; i++) {
              px = dxyz_pc * i;
              ynum = xnum - i;
              for (int j = 0; j < ynum; j++) {
                py = dxyz_pc * j;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  if (n >= element->maxsubucell) error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
                  us2u[itype][isub][n++] = m;
                  checksum += m;
                  px_center += px;
                  py_center += py;
                  pz_center += pz;
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          // calculate sub-element center shape function array
          // sub-element center is centroid of atoms

          px = px_center/n;
          py = py_center/n;
          pz = pz_center/n;

          shape_array_center_subelem[itype][isub][0] = (1.0 - px - py - pz) / iapc;
          shape_array_center_subelem[itype][isub][1] = px / iapc;
          shape_array_center_subelem[itype][isub][2] = py / iapc;
          shape_array_center_subelem[itype][isub][3] = pz / iapc;

          nsubelem[itype]++;
          nucell_subelem[itype][isub++] = n;
          ntotal += n;
        }
      }
    }

  } else if (element_shape_ids[itype] == Element::WEDGE) {

    double px, py, pz, px_center, py_center, pz_center;
    int znum;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    int esplit = round(cbrt(sqrt(ncellx * ncelly * (ncelly+1)/2)));
    subsplit[itype] = esplit;
    double dx_pc = 2.0/(ncellx-1);
    double dyz_pc = 1.0/(ncelly-1);
    double ds_x = (2.0+dx_pc)/esplit;
    double ds_yz = (1.0+dyz_pc)/esplit;

    for (int ix = 0; ix < esplit; ix++) {
      xlo = -1.0 - dx_pc/2 + ix * ds_x;
      xhi = -1.0 - dx_pc/2 + (ix+1) * ds_x;

      for (int iy = 0; iy < esplit; iy++) {
        ylo = -dyz_pc/2 + iy * ds_yz;
        yhi = -dyz_pc/2 + (iy+1) * ds_yz;

        for (int iz = 0; iz < esplit; iz++) {
          zlo = -dyz_pc/2 + iz * ds_yz;
          zhi = -dyz_pc/2 + (iz+1) * ds_yz;


          if (nsubelem[itype] >= element->maxsubelem) 
            error->all(FLERR, "Too many sub-elements, boost maxsubelem through element_modify command");

          // place virtual atoms into sub element isub

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of virtual atom in element

          px_center = py_center = pz_center = 0.0;
          for (int i = 0; i < ncellx; i++) {
            px = dx_pc * i - 1.0;
            for (int j = 0; j < ncelly; j++) {
              py = dyz_pc * j;
              znum = ncelly - j;
              for (int k = 0; k < znum; k++) {
                pz = dyz_pc * k;
                if (px >= xlo && px < xhi &&
                    py >= ylo && py < yhi &&
                    pz >= zlo && pz < zhi) {
                  if (n >= element->maxsubucell) 
                    error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
                  us2u[itype][isub][n++] = m;
                  checksum += m;
                  px_center += px;
                  py_center += py;
                  pz_center += pz;
                }
                m++;
              }
            }
          }

          // don't store info if sub-element contains no atom

          if (n == 0) continue;

          // calculate sub-element center shape function array
          // sub-element center is centroid of atoms

          px = px_center/n;
          py = py_center/n;
          pz = pz_center/n;

          shape_array_center_subelem[itype][isub][0] = 0.5 * (1 - px) * (1 - py - pz) / iapc;
          shape_array_center_subelem[itype][isub][1] = 0.5 * (1 - px) * py / iapc;
          shape_array_center_subelem[itype][isub][2] = 0.5 * (1 - px) * pz / iapc;
          shape_array_center_subelem[itype][isub][3] = 0.5 * (1 + px) * (1 - py - pz) / iapc;
          shape_array_center_subelem[itype][isub][4] = 0.5 * (1 + px) * py / iapc;
          shape_array_center_subelem[itype][isub][5] = 0.5 * (1 + px) * pz / iapc;

          nsubelem[itype]++;
          nucell_subelem[itype][isub++] = n;
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
    nucell_subelem[itype][0] = nucell[itype];
    if (nucell[itype] > element->maxsubucell) 
      error->all(FLERR, "Too many atoms in a sub-element, boost maxsubucell");
    for (int i = 0; i < nucell[itype]; i++)
      us2u[itype][0][i] = i;

  }

  bigint correctsum = (nucell[itype]-1);
  correctsum *= nucell[itype];
  correctsum /= 2;
  if (checksum != correctsum)
    error->all(FLERR, "Atoms are not assigned to sub-elements correctly");
  if (ntotal != nucell[itype]) 
    error->all(FLERR, "Total atom count in sub-elements do not match");
}

/*  -----------------------------------
    update all local element center coordinates from node coordinates
    to be called by Fix after integrating node positions
    ------------------------------------ */

void ElementVecCAC::update_center_coord()
{

  // zero center coords of elements

  size_t nbytes = sizeof(double) * element->nlocal;
  if (nbytes) {
    memset(&x[0][0], 0, 3 * nbytes);
  }

  // update center coords of elements

  for (int i = 0; i < element->nlocal; i++) {
    int inpe = npe[etype[i]];
    int iapc = apc[etype[i]];
    for (int j = 0; j < iapc; j++) 
      for (int k = 0; k < inpe; k++) {
        x[i][0] += nodex[i][j][k][0];
        x[i][1] += nodex[i][j][k][1];
        x[i][2] += nodex[i][j][k][2];
      }
    x[i][0] = x[i][0]/inpe/iapc;
    x[i][1] = x[i][1]/inpe/iapc;
    x[i][2] = x[i][2]/inpe/iapc;
  }
}

/*  -----------------------------------
    update local element I node idim coordinate 
    from center idim coordinate
    ------------------------------------ */

void ElementVecCAC::update_node_coord(int i, int idim)
{
  double xtmp = 0.0;
  int inpe = npe[etype[i]];
  int iapc = apc[etype[i]];
  for (int j = 0; j < iapc; j++)
    for (int k = 0; k < inpe; k++)
      xtmp += nodex[i][j][k][idim];
  xtmp = x[i][idim] - xtmp/inpe/iapc;
  for (int j = 0; j < iapc; j++) 
    for (int k = 0; k < inpe; k++)
      nodex[i][j][k][idim] += xtmp;
}

/*  -----------------------------------
    update all local element node coordinates 
    from center coordinates
    ------------------------------------ */

void ElementVecCAC::update_node_coord()
{
  double xtmp, ytmp, ztmp;
  for (int i = 0; i < element->nlocal; i++) {
    int inpe = npe[etype[i]];
    int iapc = apc[etype[i]];
    xtmp = ytmp = ztmp = 0.0;
    for (int j = 0; j < iapc; j++) 
      for (int k = 0; k < inpe; k++) {
        xtmp += nodex[i][j][k][0];
        ytmp += nodex[i][j][k][1];
        ztmp += nodex[i][j][k][2];
      }
    xtmp = x[i][0] - xtmp/inpe/iapc;
    ytmp = x[i][1] - ytmp/inpe/iapc;
    ztmp = x[i][2] - ztmp/inpe/iapc;
    for (int j = 0; j < iapc; j++) 
      for (int k = 0; k < inpe; k++) {
        nodex[i][j][k][0] += xtmp;
        nodex[i][j][k][1] += ytmp;
        nodex[i][j][k][2] += ztmp;
      }
  }
}

/*  --------------------------------------------------------------
    discritize element I into discrete atoms
    to be called by disc_element command
    return number of atoms created
    --------------------------------------------------------------- */

int ElementVecCAC::element2atom(int i) 
{

  // get the coordinates of virtual atom k
  // move last element to element i

  int ietype = etype[i];
  double coord[3];
  for (int ibasis = 0; ibasis < apc[ietype]; ibasis++) 
    for (int iucell = 0; iucell < nucell[ietype]; iucell++) {
      interpolate(coord, nodex, i, ibasis, iucell, 3);
      atom->avec->create_atom(coord, ctype[i][ibasis], 0); 
    }


  copy(--element->nlocal, i, 1);
  return nucell[ietype] * apc[ietype];
}

/*  --------------------------------------------------------------
    discritize element I into smaller elements with new etype and same ctype
    to be called by disc_element command
    return number of elements added
 ****NEED TO BE FIXED FOR NEW DATA STRUCTURE ****
 ****NO LONGER USED, MAY BE DELETED LATED ****
 --------------------------------------------------------------- */

int ElementVecCAC::element2element(int i, int ietype_new, tagint id_offset) 
{
  /*  
  // get the coordinates of virtual atom k
  // move last element to element i

  int ietype_old = etype[i];
  int ictype = ctype[i];
  if (apc[ietype_old] != apc[ietype_new]) 
  error->all(FLERR, "New and old etype must have the same apc");

  // number of new elements along each direction

  int ne[3];
  tagint itag = tag[i];
  tagint mytag;
  int iucell;

  int n = 0;

  if (element_shape_ids[ietype_old] == Element::HEXAHEDRON) {

  double **nodecoord;
  memory->create(nodecoord, 8, 3, "evec:nodecoord");
  int d[8][3] = {
  {0, 0, 0}, 
  {1, 0, 0}, 
  {1, 1, 0}, 
  {0, 1, 0}, 
  {0, 0, 1}, 
  {1, 0, 1}, 
  {1, 1, 1}, 
  {0, 1, 1}};

  if (element_shape_ids[ietype_new] != Element::HEXAHEDRON) 
  error->one(FLERR, "Hex elements can only be divided into Hex elements");

  // check if element i can be divided evenly

  if (ncells[ietype_old][0] % ncells[ietype_new][0] ||
  ncells[ietype_old][1] % ncells[ietype_new][1] ||
  ncells[ietype_old][2] % ncells[ietype_new][2]) 
  error->one(FLERR, "Cannot divide elements evenly");

  for (int ix = 0; ix < ncells[ietype_old][0]; ix += ncells[ietype_new][0])
  for (int iy = 0; iy < ncells[ietype_old][1]; iy += ncells[ietype_new][1])
  for (int iz = 0; iz < ncells[ietype_old][2]; iz += ncells[ietype_new][2]) {

  // interpolate new node coords from element I

  for (int new_node = 0; new_node < 8; new_node++) { 

  iucell = (ix + d[new_node][0] * (ncells[ietype_new][0]-1)) * ncells[ietype_old][1] * ncells[ietype_old][2] 
  + (iy + d[new_node][1] * (ncells[ietype_new][1]-1)) * ncells[ietype_old][2] 
  + (iz + d[new_node][2] * (ncells[ietype_new][2]-1));

  nodecoord[new_node][0] = 0.0;
  nodecoord[new_node][1] = 0.0;
  nodecoord[new_node][2] = 0.0;
  for (int inode = 0; inode < 8; inode++) {
  nodecoord[new_node][0] += nodex[i][inode][0] * shape_array[ietype_old][iucell][inode];
  nodecoord[new_node][1] += nodex[i][inode][1] * shape_array[ietype_old][iucell][inode];
  nodecoord[new_node][2] += nodex[i][inode][2] * shape_array[ietype_old][iucell][inode];
  }
  } 
  // first element will have the same tag as old element
  // other elements will start from id_offset and increment by apc
  // if id_offset = 0, set mytag = 0;

  if (id_offset) {
  if (n) mytag = id_offset + (n-1) * apc[ietype_old];
  else mytag = itag;
  } else mytag = 0;

  create_element(nodecoord, ietype_new, ictype, mytag);
  n++;
}

memory->destroy(nodecoord);
} else if (element_shape_ids[ietype_old] == Element::QUADRILATERAL) {

  double **nodecoord;
  memory->create(nodecoord, 4, 3, "evec:nodecoord");
  int d[4][2] = {
    {0, 0}, 
    {1, 0}, 
    {1, 1}, 
    {0, 1}};

  if (element_shape_ids[ietype_new] != Element::QUADRILATERAL)
    error->one(FLERR, "Quad elements can only be divided into Quad elements");

  // check if element i can be divided evenly

  if (ncells[ietype_old][0] % ncells[ietype_new][0] ||
      ncells[ietype_old][1] % ncells[ietype_new][1])
    error->one(FLERR, "Cannot divide elements evenly");

  // interpolate new node coords from element I
  //
  for (int ix = 0; ix < ncells[ietype_old][0]; ix += ncells[ietype_new][0])
    for (int iy = 0; iy < ncells[ietype_old][1]; iy += ncells[ietype_new][1]) {

      for (int new_node = 0; new_node < 4; new_node++) { 
        nodecoord[new_node][0] = 0.0;
        nodecoord[new_node][1] = 0.0;
        nodecoord[new_node][2] = 0.0;

        // interpolated index of this new node in old element

        iucell = (ix + d[new_node][0] * (ncells[ietype_new][0]-1)) * ncells[ietype_old][1]
          + (iy + d[new_node][1] * (ncells[ietype_new][1]-1));

        for (int inode = 0; inode < 4; inode++) {
          nodecoord[new_node][0] += nodex[i][inode][0] * shape_array[ietype_old][iucell][inode];
          nodecoord[new_node][1] += nodex[i][inode][1] * shape_array[ietype_old][iucell][inode];
        }
        nodecoord[new_node][3] = 0.0;
      }

      // first element will have the same tag as old element
      // other elements will start from id_offset and increment by apc
      // if id_offset = 0, set mytag = 0;

      if (id_offset) {
        if (n) mytag = id_offset + (n-1) * apc[ietype_old];
        else mytag = itag;
      } else mytag = 0;

      create_element(nodecoord, ietype_new, ictype, mytag);
      n++;
    }
  memory->destroy(nodecoord);

} else {
  // add new element shape here
  error->one(FLERR, "Element split for this type has not been set yet");
}

// init per-element fix/compute/variable values for created elements

element->data_fix_compute_dump_variable(element->nlocal-n, element->nlocal);

// remove element I by replacing it with last element in list

copy(--element->nlocal, i, 1);

// clean up and return number of additional elements

return n;
*/
return 0;
}


/*  ----------------------------------------------------------------------
    return # of bytes of allocated memory
    -------------------------------------------------------------------------  */

bigint ElementVecCAC::memory_usage()
{
  bigint bytes = 0;

  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;
  int max_ngcell = element->max_ngcell;
  if (element->memcheck("tag")) bytes += memory->usage(tag, nmax);
  if (element->memcheck("ctype")) bytes += memory->usage(ctype, nmax, maxapc);
  if (element->memcheck("etype")) bytes += memory->usage(etype, nmax);
  if (element->memcheck("mask")) bytes += memory->usage(mask, nmax);
  if (element->memcheck("image")) bytes += memory->usage(image, nmax);
  if (element->memcheck("nodemask")) bytes += memory->usage(nodemask, nmax, maxapc, maxnpe);
  if (element->memcheck("x")) bytes += memory->usage(x, nmax, 3);
  if (element->memcheck("surface_plane")) bytes += memory->usage(surface_plane, nmax, 3, 3);
  if (element->memcheck("cell_size")) bytes += memory->usage(cell_size, nmax, 3);
  if (element->memcheck("nodex")) bytes += memory->usage(nodex, nmax, maxapc, maxnpe, 3);
  if (element->memcheck("nodev")) bytes += memory->usage(nodev, nmax, maxapc, maxnpe, 3);
  if (element->memcheck("nodef")) bytes += memory->usage(nodef, nmax, maxapc, maxnpe, 3);
  if (element->memcheck("gaussf")) bytes += memory->usage(gaussf, nmax, maxapc, max_ngcell, 3);
  if (element->memcheck("subelem_size")) bytes += memory->usage(subelem_size, nmax);
  if (element->memcheck("initial_box_size")) bytes += memory->usage(initial_box_size, nmax, 3);
  if (element->memcheck("element_bound_box")) bytes += memory->usage(element_bound_box, nmax, 6);
  return bytes;
}

/*  ----------------------------------------------------------------------
    update internal surface_plane (axes) of all elements
    surface_plane is defined as unit normal vector of the plane
    -------------------------------------------------------------------------  */

void ElementVecCAC::update_surface_plane()
{
  double xnorm, ynorm, znorm;
  double xdir[3], ydir[3], zdir[3];
  for (int i = 0; i < element->nlocal+element->nghost; i++) {
    int itype = etype[i];
    if (element_shape_ids[itype] == Element::QUADRILATERAL) {
      if (domain->dimension == 2) {
        for (int j = 0; j < 2; j++) {
          for (int k = 0; k < apc[itype]; k++) {
            xdir[j] = 
              (nodex[i][k][1][j] - nodex[i][k][0][j] + 
               nodex[i][k][2][j] - nodex[i][k][3][j]);
            ydir[j] = 
              (nodex[i][k][3][j] - nodex[i][k][0][j] +
               nodex[i][k][2][j] - nodex[i][k][1][j]);
            zdir[j] = 0.0;
          }
          xdir[j] /= 2 * apc[itype];
          ydir[j] /= 2 * apc[itype];
          zdir[j] /= 2 * apc[itype];
        }
        xdir[2] = 0.0;
        ydir[2] = 0.0;
        zdir[2] = 1.0;
      } else {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < apc[itype]; k++) {
            xdir[j] = 
              (nodex[i][k][1][j] - nodex[i][k][0][j] + 
               nodex[i][k][2][j] - nodex[i][k][3][j]);
            ydir[j] = 
              (nodex[i][k][3][j] - nodex[i][k][0][j] +
               nodex[i][k][2][j] - nodex[i][k][1][j]);
          }
          xdir[j] /= 2 * apc[itype];
          ydir[j] /= 2 * apc[itype];
          zdir[j] /= 2 * apc[itype];
        }
        cross3(xdir, ydir, zdir);
        norm3(zdir);
      }
    } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < apc[itype]; k++) {
          xdir[j] = 
            (nodex[i][k][1][j] - nodex[i][k][0][j] + 
             nodex[i][k][2][j] - nodex[i][k][3][j] +
             nodex[i][k][5][j] - nodex[i][k][4][j] +
             nodex[i][k][6][j] - nodex[i][k][7][j]);
          ydir[j] = 
            (nodex[i][k][3][j] - nodex[i][k][0][j] +
             nodex[i][k][2][j] - nodex[i][k][1][j] +
             nodex[i][k][7][j] - nodex[i][k][4][j] +
             nodex[i][k][6][j] - nodex[i][k][5][j]);
          zdir[j] = 
            (nodex[i][k][4][j] - nodex[i][k][0][j] +
             nodex[i][k][7][j] - nodex[i][k][3][j] +
             nodex[i][k][5][j] - nodex[i][k][1][j] +
             nodex[i][k][6][j] - nodex[i][k][2][j]);
        }
        xdir[j] /= 4 * apc[itype];
        ydir[j] /= 4 * apc[itype];
        zdir[j] /= 4 * apc[itype];
      }
    } else {
      // add new element shape here
      // error->one(FLERR, "Slip plane for this element shape has not been set");
      continue;
    }

    // calculate cell sizes and
    // normalize vectors

    xnorm = norm3(xdir);
    ynorm = norm3(ydir);
    znorm = norm3(zdir);
    cross3(ydir, zdir, surface_plane[i][0]);
    cross3(zdir, xdir, surface_plane[i][1]);
    cross3(xdir, ydir, surface_plane[i][2]);
    norm3(surface_plane[i][0]);
    norm3(surface_plane[i][1]);
    norm3(surface_plane[i][2]);
    if (dot3(xdir, surface_plane[i][0]) < 0) negate3(surface_plane[i][0]);
    if (dot3(ydir, surface_plane[i][1]) < 0) negate3(surface_plane[i][1]);
    if (dot3(zdir, surface_plane[i][2]) < 0) negate3(surface_plane[i][2]);
    cell_size[i][0] = xnorm * fabs(dot3(xdir, surface_plane[i][0]));
    cell_size[i][1] = ynorm * fabs(dot3(ydir, surface_plane[i][1]));

    if (element_shape_ids[itype] == Element::QUADRILATERAL) {
      cell_size[i][2] = 0;
    } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {
      cell_size[i][2] = znorm * fabs(dot3(zdir, surface_plane[i][2]));
    } else {

    }
  }
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at an virtual atom 
    n = number of values
    -------------------------------------------------------------------------  */

void ElementVecCAC::interpolate(double *value, double ****nodevalue, int i, int ibasis, int iucell, int n) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  for (int k = 0; k < n; k++) {
    value[k] = 0.0;
    for (int l = 0; l < inpe; l++) 
      value[k] += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l][k];
  }
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at an virtual atom 
    n = number of values
    -------------------------------------------------------------------------  */

void ElementVecCAC::interpolate(double *value, int ****nodevalue, int i, int ibasis, int iucell, int n) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  for (int k = 0; k < n; k++) {
    value[k] = 0.0;
    for (int l = 0; l < inpe; l++) 
      value[k] += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l][k];
  }
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at virtual atom
    k = index of nodal value
    -------------------------------------------------------------------------  */

double ElementVecCAC::interpolate(double ****nodevalue, int i, int ibasis, int iucell, int index) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  double value = 0.0;
  for (int l = 0; l < inpe; l++) 
    value += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l][index];
  return value;
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at virtual atom
    -------------------------------------------------------------------------  */

double ElementVecCAC::interpolate(double ***nodevalue, int i, int ibasis, int iucell) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  double value = 0.0;
  for (int l = 0; l < inpe; l++) 
    value += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l];
  return value;
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at virtual atom
    k = index of nodal value
    -------------------------------------------------------------------------  */

double ElementVecCAC::interpolate(int ****nodevalue, int i, int ibasis, int iucell, int index) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  double value = 0.0;
  for (int l = 0; l < inpe; l++) 
    value += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l][index];
  return value;
}

/*  ----------------------------------------------------------------------
    interpolate values from nodevalues in element I at virtual atom
    -------------------------------------------------------------------------  */

double ElementVecCAC::interpolate(int ***nodevalue, int i, int ibasis, int iucell) 
{
  int ietype = etype[i];
  int inpe = npe[ietype];
  int iapc = apc[ietype];
  if (iucell >= nucell[ietype] || iucell < 0 ||
      ibasis >= iapc || ibasis < 0) {
    printf("Etype %d me = %d i = %d iucell = %d ibasis = %d nucell = %d ucell %d %d %d gcell %d %d %d\n", 
        ietype, comm->me, i, iucell, ibasis, nucell[ietype], 
        ncells[ietype][0], ncells[ietype][1], ncells[ietype][2], 
        ngcells[ietype][0], ngcells[ietype][1], ngcells[ietype][2]);
    error->one(FLERR, "Virtual atom index out of bound");
  }
  double value = 0.0;
  for (int l = 0; l < inpe; l++) 
    value += shape_array[ietype][iucell][l] * nodevalue[i][ibasis][l];
  return value;
}

/*  ----------------------------------------------------------------------
    add atom plane next to elements
    Called from input add_atoms command
    -------------------------------------------------------------------------  */

void ElementVecCAC::add_atoms(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR, "Illegal add_atoms command");

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR, "Could not find add_atoms group ID");
  int groupbit = group->bitmask[igroup];

  // check which type to modify

  int n1list[4], n2list[4], dim0, dim1, dim2;
  if (strcmp(arg[1], "+x") == 0) {
    n1list[0] = 1; n2list[0] = 0;
    n1list[1] = 2; n2list[1] = 3;
    n1list[2] = 6; n2list[2] = 7;
    n1list[3] = 5; n2list[3] = 4;
    dim0 = 0; dim1 = 1; dim2 = 2;
  } else if (strcmp(arg[1], "-x") == 0) { 
    n1list[0] = 0; n2list[0] = 1;
    n1list[1] = 3; n2list[1] = 2;
    n1list[2] = 7; n2list[2] = 6;
    n1list[3] = 4; n2list[3] = 5;
    dim0 = 0; dim1 = 1; dim2 = 2;
  } else if (strcmp(arg[1], "+y") == 0) { 
    n1list[0] = 2; n2list[0] = 1;
    n1list[1] = 3; n2list[1] = 0;
    n1list[2] = 7; n2list[2] = 4;
    n1list[3] = 6; n2list[3] = 5;
    dim0 = 1; dim1 = 0; dim2 = 2;
  } else if (strcmp(arg[1], "-y") == 0) { 
    n1list[0] = 1; n2list[0] = 2;
    n1list[1] = 0; n2list[1] = 3;
    n1list[2] = 4; n2list[2] = 7;
    n1list[3] = 5; n2list[3] = 6;
    dim0 = 1; dim1 = 0; dim2 = 2;
  } else if (strcmp(arg[1], "+z") == 0) { 
    n1list[0] = 4; n2list[0] = 0;
    n1list[1] = 5; n2list[1] = 1;
    n1list[2] = 6; n2list[2] = 2;
    n1list[3] = 7; n2list[3] = 3;
    dim0 = 2; dim1 = 0; dim2 = 1;
  } else if (strcmp(arg[1], "-z") == 0) { 
    n1list[0] = 0; n2list[0] = 4;
    n1list[1] = 1; n2list[1] = 5;
    n1list[2] = 2; n2list[2] = 6;
    n1list[3] = 3; n2list[3] = 7;
    dim0 = 2; dim1 = 0; dim2 = 1;
  } else error->all(FLERR, "Illegal add_atoms command");

  int nlayers = universe->inumeric(FLERR, arg[2]);

  // modify elements and add atoms

  double coord[3];
  int iucell, itype, ncellx, ncelly, ncellz, iapc;
  int natoms_added = 0;
  int nlocal_previous = atom->nlocal;
  for (int i = 0; i < element->nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    itype = etype[i];
    if (element_shape_ids[itype] != Element::HEXAHEDRON)
      error->all(FLERR, "add_atoms command only works for Hex element for now");

    ncellx = ncells[itype][0];
    ncelly = ncells[itype][1];
    ncellz = ncells[itype][2];
    for (int n = 0; n < nlayers; n++) {
      for (int ibasis = 0; ibasis < iapc; ibasis++) {
        for (int j = 0; j < ncells[itype][dim1]; j++)  {
          for (int k = 0; k < ncells[itype][dim2]; k++) {
            // crazy!!! 
            iucell = (j * (dim1 == 2) + k * (dim2 == 2)  + (ncellz-1) * ((dim1!=2)&&(dim2!=2)&&(n1list[0]==4)))
              + ncellz * (j * (dim1 == 1) + k * (dim2 == 1) + (ncelly-1) * ((dim1!=1)&&(dim2!=1)&&(n1list[0]==2)))
              + ncelly * ncellz * (j * (dim1 == 0) + k * (dim2 == 0) + (ncellx-1) * ((dim1!=0)&&(dim2!=0)&&(n1list[0]==1)));
            interpolate(coord, nodex, i, ibasis, iucell, 3);
            atom->avec->create_atom(coord, ctype[i][ibasis], 0); 
            natoms_added++;
          }
        }
        for (int j = 0; j < 4; j++) 
          for (int k = 0; k < 3; k++) 
            nodex[i][ibasis][n1list[j]][k] += 
              (nodex[i][ibasis][n2list[j]][k] - nodex[i][ibasis][n1list[j]][k])
              / (ncells[itype][dim0] + nlayers - 1.0);
      }
    }
  }
  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_dump_variable(nlocal_previous, atom->nlocal);

  // add new tags to atoms

  atom->tag_extend();
  int total_atoms_added;
  MPI_Allreduce(&natoms_added, &total_atoms_added, 1, MPI_INT, MPI_SUM, world);
  atom->natoms += total_atoms_added;
  if (comm->me == 0) {
    if (screen) fprintf(screen, " %d atoms added by add_atoms command.\n", total_atoms_added);
    if (logfile) fprintf(logfile, " %d atoms added by add_atoms command.\n", total_atoms_added);
  }
}

/*  ----------------------------------------------------------------------
    called from fix_adaptive to split element I along 'idim' direction 
    flag = 0: first pass, check if element need to be splitted and request new etype. Element is not splitted yet
    flag = 1: second pass, new etype has been created and element can now be splitted
    return 0 if element does not need splitted or discretized to all atoms
    1 if 1 new element created
    2 if 2 new elements created
NOTE: OLD METHOD
-------------------------------------------------------------------------  */
/* 
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
   int index_offset = (idim == 2) + ncellz * (idim == 1) + ncelly * ncellz * (idim == 0);
   int type_info[8];
   type_info[0] = ncellx;
   type_info[1] = ncelly;
   type_info[2] = ncellz;
   type_info[3] = ngcells[ietype][0];
   type_info[4] = ngcells[ietype][1];
   type_info[5] = ngcells[ietype][2];
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

// limit the number of gaussian cells to the number of unit cells

type_info[3+idim] = MIN(type_info[idim], ngcells[ietype][idim]);

// request new etype

request_new_etype(type_info);
}

// upper element

if (split+1 < chi) {

type_info[idim] = chi - split;

// limit the number of gaussian cells to the number of unit cells

type_info[3+idim] = MIN(type_info[idim], ngcells[ietype][idim]);

// request new etype

request_new_etype(type_info);
}

// if flag != 0, split elements

} else {
  double **nodecoord;
  memory->create(nodecoord, element->maxnpe, 3, "evec:nodecoord");
  int *iucell_list = new int[element->maxnpe];
  int *n2u = element->n2u[ietype];

  // lower element
  // if has only one plane, discretize into atoms

  if (split == clo) {
    int hi[3], iucell;
    hi[0] = ncellx;
    hi[1] = ncelly;
    hi[2] = ncellz;
    hi[idim] = 1;
    for (int ii = 0; ii < hi[0]; ii ++)
      for (int jj = 0; jj < hi[1]; jj ++)
        for (int kk = 0; kk < hi[2]; kk ++) {
          iucell = kk + jj * ncellz + ii * ncellz * ncelly;
          interpolate(coord, nodex, i, iucell, 3);
          atom->avec->create_atom(coord, ictype, 0);

          // pass on stored data for new atom J from fix/compute/variable and velocity

          int j = atom->nlocal-1;
          atom->data_pass_on_fix_compute_variable(j, i, iucell);
          interpolate(atom->v[j], nodev, i, iucell, 3);
        } 
  } else {
    split_flag++;
    type_info[idim] = split - clo + 1;

    // limit the number of gaussian cells to the number of unit cells

    type_info[3+idim] = MIN(type_info[idim], ngcells[ietype][idim]);

    int myetype = find_etype(type_info);
    if (myetype == 0) error->one(FLERR, "Etype not defined yet");

    if (domain->dimension == 3) {
      for (int ii = 0; ii < 4; ii++) {

        // lower set of nodes stays are the same as the origional element
        // calculate upper set of nodes by scaling along edges

        int lnode = node_set_3D[idim][0][ii];
        int unode = node_set_3D[idim][1][ii]; 

        iucell_list[lnode] = n2u[lnode];
        iucell_list[unode] = n2u[lnode] + (type_info[idim]-1) * index_offset;

        for (int jj = 0; jj < 3; jj++) {
          nodecoord[lnode][jj] = nodex[i][lnode][jj];
          nodecoord[unode][jj] = nodecoord[lnode][jj] + 
            (nodex[i][unode][jj]-nodex[i][lnode][jj])
            *(type_info[idim]-1)/(ncells[ietype][idim]-1);
        }
      }
      create_element(nodecoord, myetype, ictype, 0);

      // pass on stored data for new element in fix/compute/variable and velocity

      int j = element->nlocal-1;
      element->data_pass_on_fix_compute_variable(j, i, iucell_list);
      for (int ii = 0; ii < 8; ii++)
        interpolate(nodev[j][ii], nodev, i, iucell_list[ii], 3);

    } else {
      for (int ii = 0; ii < 2; ii++) {

        // lower set of nodes stays are the same as the origional element
        // calculate upper set of nodes by scaling along edges

        int lnode = node_set_2D[idim][0][ii];
        int unode = node_set_2D[idim][1][ii]; 

        iucell_list[lnode] = n2u[lnode];
        iucell_list[unode] = n2u[lnode] + (type_info[idim]-1) * index_offset;

        for (int jj = 0; jj < 2; jj++) {
          nodecoord[lnode][jj] = nodex[i][lnode][jj];
          nodecoord[unode][jj] = nodecoord[lnode][jj] + 
            (nodex[i][unode][jj]-nodex[i][lnode][jj])
            *(type_info[idim]-1)/(ncells[ietype][idim]-1);
          coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
        }
      }

      // z direction coord stay the same

      for (int ii = 0; ii < 4; ii++) 
        nodecoord[ii][2] = nodex[i][ii][2];
      create_element(nodecoord, ietype, ictype, 0);
      int j = element->nlocal-1;
      element->data_pass_on_fix_compute_variable(j, i, iucell_list);
      for (int ii = 0; ii < 4; ii++)
        interpolate(nodev[j][ii], nodev, i, iucell_list[ii], 3);

    }

    // assign new plane_id for new element

    plane_id[element->nlocal-1] = static_cast<double> (split+clo);
    plane_id[element->nlocal-1] /= 2;
  }

  // upper element
  // if has only one plane, discretize into atoms

  if (split+1 == chi) {
    int lo[3], iucell;
    lo[0] = lo[1] = lo[2] = 0;
    lo[idim] = ncells[ietype][idim]-1;
    for (int ii = lo[0]; ii < ncellx; ii ++)
      for (int jj = lo[1]; jj < ncelly; jj ++)
        for (int kk = lo[2]; kk < ncellz; kk ++) {
          iucell = kk + jj * ncellz + ii * ncellz * ncelly;
          interpolate(coord, nodex, i, iucell, 3);
          atom->avec->create_atom(coord, ictype, 0);
          int j = atom->nlocal-1;
          atom->data_pass_on_fix_compute_variable(j, i, iucell);
          interpolate(atom->v[j], nodev, i, iucell, 3);
        } 

  } else {

    split_flag++;
    type_info[idim] = chi - split;

    // limit the number of gaussian cells to the number of unit cells

    type_info[3+idim] = MIN(type_info[idim], ngcells[ietype][idim]);

    int myetype = find_etype(type_info);
    if (myetype == 0) error->one(FLERR, "Etype not defined yet");

    if (domain->dimension == 3) {
      for (int ii = 0; ii < 4; ii++) {

        // upper set of nodes stays are the same as the origional element
        // calculate lower set of nodes by scaling along edges

        int lnode = node_set_3D[idim][0][ii];
        int unode = node_set_3D[idim][1][ii]; 

        iucell_list[unode] = n2u[unode];
        iucell_list[lnode] = n2u[unode] - (type_info[idim]-1) * index_offset;

        for (int jj = 0; jj < 3; jj++) {
          nodecoord[unode][jj] = nodex[i][unode][jj];
          nodecoord[lnode][jj] = nodecoord[unode][jj] + 
            (nodex[i][lnode][jj]-nodex[i][unode][jj])
            *(type_info[idim]-1)/(ncells[ietype][idim]-1);
          coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
        }
      }
      create_element(nodecoord, myetype, ictype, 0);
      int j = element->nlocal-1;
      element->data_pass_on_fix_compute_variable(j, i, iucell_list);
      for (int ii = 0; ii < 8; ii++)
        interpolate(nodev[j][ii], nodev, i, iucell_list[ii], 3);

    } else {
      for (int ii = 0; ii < 2; ii++) {

        // lower set of nodes stays are the same as the origional element
        // calculate upper set of nodes by scaling along edges

        int lnode = node_set_2D[idim][0][ii];
        int unode = node_set_2D[idim][1][ii]; 

        iucell_list[unode] = n2u[unode];
        iucell_list[lnode] = n2u[unode] - (type_info[idim]-1) * index_offset;

        for (int jj = 0; jj < 2; jj++) {
          nodecoord[unode][jj] = nodex[i][unode][jj];
          nodecoord[lnode][jj] = nodecoord[unode][jj] + 
            (nodex[i][lnode][jj]-nodex[i][unode][jj])
            *(type_info[idim]-1)/(ncells[ietype][idim]-1);
          coord[jj] += nodecoord[lnode][jj] + nodecoord[unode][jj];
        }
      }

      // z direction coord stay the same

      for (int ii = 0; ii < 4; ii++) 
        nodecoord[ii][2] = nodex[i][ii][2];

      create_element(nodecoord, ietype, ictype, 0);
      int j = element->nlocal-1;
      element->data_pass_on_fix_compute_variable(j, i, iucell_list);
      for (int ii = 0; ii < 4; ii++)
        interpolate(nodev[j][ii], nodev, i, iucell_list[ii], 3);

    }

    // assign new plane_id for new element

    plane_id[element->nlocal-1] = static_cast<double> (split+chi+1);
    plane_id[element->nlocal-1] /= 2;

  }

  // delete element I

  copy(element->nlocal-1, i, 1);
  element->nlocal--;
  delete [] iucell_list;
  memory->destroy(nodecoord);
}
return split_flag;
}
*/
/*  ----------------------------------------------------------------------
    check if a plane with normal vector <n1, n2, n3> that
    go through point P intersects with element I
    -------------------------------------------------------------------------  */

int ElementVecCAC::check_split_element(int i, double n[3], double P[3], 
    double tol, double split_width, int *cut, int &cutdim)
{
  /*  NEED TO FIX
      int ietype = etype[i];

  // check if cut plane is parallel to any surface plane of I
  // if not, return 0

  cutdim = -1;
  for (int dim = 0; dim < 3; dim++) {
  double tmp = 1.0 - fabs(dot3(surface_plane[i][dim], n));
  if (tmp < tol) {
  cutdim = dim;
  break;
  }
  }
  //if (element->tag[i] == 481 || element->tag[i] == 482)
  //  printf("Element %d step = %d P = %g %g %g normal = %g %g %g\n   surface = %g %g %g dot = %g tol = %g cutdim = %d\n", 
  //      tag[i], update->ntimestep, 
  //      P[0], P[1], P[2], 
  //      n[0], n[1], n[2], 
  //      surface_plane[i][2][0], 
  //      surface_plane[i][2][1], 
  //      surface_plane[i][2][2], 
  //      dot3(surface_plane[i][2], n), tol, 
  //      cutdim
  //      );

  if (cutdim < 0) return 0;
  if (ncells[ietype][cutdim] < 3) return 0;
  //if (element->tag[i] == 481 || element->tag[i] == 482)
  //  printf("Element = %d step = %d surface_plane = %g %g %g\n", tag[i], update->ntimestep, 
  //      surface_plane[i][cutdim][0], 
  //      surface_plane[i][cutdim][1], 
  //      surface_plane[i][cutdim][2]);

  // project point P onto slip plane normal line of element I
  // if projection falls inside element --> need splitting


  double dist = dot3(surface_plane[i][cutdim], P) 
  - dot3(surface_plane[i][cutdim], x[i]);
  double a = cell_size[i][cutdim]/(ncells[ietype][cutdim]-1);
  double hi = cell_size[i][cutdim]/2;
  double lo = -hi;

  //if (element->tag[i] == 481 || element->tag[i] == 482)
  //  printf("Element %d step = %d dist = %g lo/hi = %g/%g x = %g %g %g\n", 
  //      tag[i], update->ntimestep, dist, lo, hi, x[i][0], x[i][1], x[i][2]);
  if (dist < lo || dist > hi) return 0;

  double cutd = (dist-lo)/a;
  cut[0] = static_cast<int> (cutd);
  //if (element->tag[i] == 481 || element->tag[i] == 482)
  //  printf("Element %d step = %d cutd = %g cut = %d split_width = %g\n", 
  //      tag[i], update->ntimestep, cutd, cut, split_width);

  if (cutd-cut[0] < (1-split_width)/2 || cutd-cut[0] > (1+split_width)/2) return 0;
  cut[0]++;

  //printf("cut = %d cutdim = %d surface_plane = %g %g %g \nx= %g %g %g \nP-x = %g %g %g \ndist = %g lo = %g a = %g \n"
  //    , cut, cutdim, 
  //    surface_plane[i][cutdim][0], 
  //    surface_plane[i][cutdim][1], 
  //    surface_plane[i][cutdim][2], 
  //    x[i][0], x[i][1], x[i][2], 
  //    P[0] - x[i][0], 
  //    P[1] - x[i][1], 
  //    P[2] - x[i][2], 
  //    dist, lo, a);
  //printf("nodex = \n");
  //for (int iii = 0; iii < 8; iii++)
  //  printf("%d %g %g %g\n", iii+1, 
  //      element->nodex[i][iii][0], 
  //      element->nodex[i][iii][1], 
  //      element->nodex[i][iii][2]);
  //printf("Cut element %d normal = %g %g %g \nP = %g %g %g \ncellsize = %g timestep = %d\n\n", tag[i], n[0], n[1], n[2], P[0], P[1], P[2], cell_size[i][cutdim], update->ntimestep);

  if (cut[0] > 0 && cut[0] < ncells[ietype][cutdim]) 
    return 1;

  */
    return 0;

}

/*  ----------------------------------------------------------------------
    check if element J's surface K cut through element I
    -------------------------------------------------------------------------  */

int ElementVecCAC::check_split_element(int i, int j, int surface, double tol, int *cut, int &cutdim, int maxcut)
{
  /*  NEED TO FIX
      int jetype = etype[j];
      int j_shape_id = element_shape_ids[jetype];
  // check if cut surface of J parallel to any surface plane of I
  // if not, return 0

  int jcutdim, dir;

  if (j_shape_id == Element::HEXAHEDRON) {
  jcutdim = surface/2;
  dir = surface%2;
  } else if (j_shape_id == Element::QUADRILATERAL) {
  if (domain->dimension == 3) {
  jcutdim = 2;
  dir = 0;
  } else {
  jcutdim = surface/2;
  dir = surface%2;
  }
  }

  cutdim = -1;
  for (int dim = 0; dim < 3; dim++) {
  double tmp = 1.0 - fabs(dot3(surface_plane[i][dim], surface_plane[j][jcutdim]));
  if (tmp < tol) {
  cutdim = dim;
  break;
  }
  }

  if (cutdim < 0) return 0;

  // look for node of element J closest to element I

  int jnode_closest, jnode, jj;
  double delx, dely, delz, rsq;
  double rsqmin = BIG;
  double ix = x[i][0];
  double iy = x[i][1];
  double iz = x[i][2];
  int ietype = etype[i];

  if (j_shape_id == Element::HEXAHEDRON) {
  for (jj = 0; jj < 4; jj++) {
  jnode = node_set_3D[jcutdim][dir][jj];
  delx = ix - nodex[j][jnode][0];
  dely = iy - nodex[j][jnode][1];
  delz = iz - nodex[j][jnode][2];
  rsq = delx * delx + dely * dely + delz * delz;
  if (rsq < rsqmin) { 
  rsqmin = rsq;
  jnode_closest = jnode;
  }
  }

  // project closest node coord of element J onto slip plane normal line of element I
  // if projection falls inside element --> need splitting

  //if (ncells[ietype][cutdim] < 3) return 0;

  double a = cell_size[i][cutdim]/(ncells[ietype][cutdim]-1);
  double dist = (dot3(surface_plane[i][cutdim], nodex[j][jnode_closest]) 
  - dot3(surface_plane[i][cutdim], x[i]))/a;
  double orient = dot3(surface_plane[i][cutdim], nodex[j][jnode_closest]) 
  - dot3(surface_plane[i][cutdim], x[j]);
  double hi = ncells[ietype][cutdim]/2;
  double lo = -hi;

  if (orient > 0) {
  hi -= 1.0;
  } else {
  lo += 1.0;
}

if (dist < lo || dist > hi) return 0;

cut[0] = static_cast<int> (dist-lo);
cut[0]++;


//if (cut == 1) cut++;
//else if (cut == ncells[ietype][cutdim]-1) cut--;

//printf(" ******me = %d\ndist = %g a = %g orient = %g hi = %g lo = %g cut = %d\n ****closest node[%d] = %g %g %g ix = %g %g %g\n", comm->me, 
//dist, a, orient, hi, lo, cut, jnode_closest
//, nodex[j][jnode_closest][0]
//, nodex[j][jnode_closest][1]
//, nodex[j][jnode_closest][2]
//, element->x[i][0]
//, element->x[i][1]
//, element->x[i][2]
//);


//printf("me = %d Cut element %d from element %d surface = %d timestep = %d cut = %d cutdim = %d ncells = %d\n", comm->me, tag[i], tag[j], surface, update->ntimestep, cut, cutdim, ncells[ietype][cutdim]);
//error->one(FLERR, "TEST"); 

if (cut[0] > 0 && cut[0] < ncells[ietype][cutdim]) 
  return 1;

  } else if (j_shape_id == Element::QUADRILATERAL) {
    if (domain->dimension == 3) {
      for (jnode = 0; jnode < 4; jnode++) {
        delx = ix - nodex[j][jnode][0];
        dely = iy - nodex[j][jnode][1];
        delz = iz - nodex[j][jnode][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < rsqmin) { 
          rsqmin = rsq;
          jnode_closest = jnode;
        }
      }
      double a = cell_size[i][cutdim]/(ncells[ietype][cutdim]-1);
      double dist = (dot3(surface_plane[i][cutdim], nodex[j][jnode_closest]) 
          - dot3(surface_plane[i][cutdim], x[i]))/a;
      double hi = ncells[ietype][cutdim]/2;
      double lo = -hi;

      if (dist < lo || dist > hi) return 0;

      cut[0] = static_cast<int> (dist-lo);

      if (cut[0] == 0) {
        cut[0] = 1;
        return 1;
      } else if (cut[0] == ncells[ietype][cutdim]-1) {
        return 1;
      } else if (cut[0] > 0 && cut[0] < ncells[ietype][cutdim]-1) {
        cut[1] = cut[0] + 1;
        return 2;
      }
    } else {
      // Add 2d here
    }

  } else {
    // add new element types here
  }
*/
return 0;
}

/*  ----------------------------------------------------------------------
    called from fix_adaptive to split elements in list
    return number of elements being splitted
    -------------------------------------------------------------------------  */

int ElementVecCAC::split_element(int **splitlist, int nsplit)
{
  return 0;
  /*  NEED TO FIX
      int ietype, ictype;
      int type_info[8];
      int *nelem, **cellsize;
      int ii = 0;
      int i, icut, icutdim;
      int *ilist, **nelem_array, ***cellsize_array;

  // counter for number of elements being splitted

  int n = 0;

  if (nsplit > 0) {
  memory->create(ilist, nsplit, "evec:ilist");
  memory->create(nelem_array, nsplit, 3, "evec:nelem_array");
  memory->create(cellsize_array, nsplit, 3, MAXSUB, "evec:cellsize_array");

  // first pass: check for new etypes
  // loop through all splits in list
  // each element might have several splits and
  // they are all next to each other in splitlist

  while (ii < nsplit) {
  nelem = nelem_array[n];
  cellsize = cellsize_array[n];
  nelem[0] = nelem[1] = nelem[2] = 1;

  i = splitlist[ii][0];

  ietype = etype[i];
  ietype = 1;
  type_info[6] = element->element_shape_ids[ietype];
  type_info[7] = element->apc[ietype];

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

  // apply each additional split to element cellsize list 
  // if next split still cuts element I

  if (ii < nsplit-1) {
  while (i == splitlist[ii+1][0]) {

  icut = splitlist[ii+1][1];
  icutdim = splitlist[ii+1][2];

  if (nelem[icutdim] == MAXSUB) 
  error->one(FLERR, "Too many split in one direction");
  int sum = 0;

  for (int l = 0; l < nelem[icutdim]; l++) {
  if (icut > sum + cellsize[icutdim][l])
  sum += cellsize[icutdim][l];
  else if (icut < sum + cellsize[icutdim][l]) {
  for (int m = nelem[icutdim]-1; m >= l; m--) 
  cellsize[icutdim][m+1] = cellsize[icutdim][m];
  cellsize[icutdim][l] = icut-sum;
  cellsize[icutdim][l+1] -= icut-sum;
  if (cellsize[icutdim][l] < 1 || cellsize[icutdim][l] > ncells[ietype][icutdim]
  || cellsize[icutdim][l+1] < 1 || cellsize[icutdim][l+1] > ncells[ietype][icutdim])
  {
  printf("icut = %d sum = %d cellsize[%d] = %d cellsize[%d] = %d icutdim = %d\n", icut, sum, l, cellsize[icutdim][l], l+1, cellsize[icutdim][l+1], icutdim);
  error->one(FLERR, "TEST");
}
nelem[icutdim]++;
break;
} else break;
}
ii++;
if (ii == nsplit-1) break;
}
}

// check rearrange splits and split any elements of size 2

int index;
for (int idim = 0; idim < domain->dimension; idim++) {
  index = 0; 
  while (index < nelem[idim]) {
    //printf("cellsize[%d][%d] = %d\n", idim, index, cellsize[idim][index]);
    if (cellsize[idim][index] == 2) {
      if (nelem[icutdim] == MAXSUB) 
        error->one(FLERR, "Too many split in one direction");
      for (int m = nelem[idim]; m > index+1; m--) {
        cellsize[idim][m] = cellsize[idim][m-1];
      }
      cellsize[idim][index] = 1;
      cellsize[idim][index+1] = 1;
      nelem[idim]++;
    } 
    index++;
  }
}     

// make sure the splits tally up to the original element sizes

int checksum = 0;
for (int ix = 0; ix < nelem[0]; ix++) 
checksum += cellsize[0][ix];
if (checksum > ncells[ietype][0]) {
  printf("Ncells = %d cellsize x =", ncells[ietype][0]);
  for (int ix = 0; ix < nelem[0]; ix++) 
    printf(" %d", cellsize[0][ix]);
  printf("\n");
  error->one(FLERR, "Incorrect split");
}

checksum = 0;
for (int iy = 0; iy < nelem[1]; iy++) 
checksum += cellsize[1][iy];
if (checksum > ncells[ietype][1]) {
  printf("Ncells = %d cellsize y =", ncells[ietype][1]);
  for (int iy = 0; iy < nelem[1]; iy++) 
    printf(" %d", cellsize[1][iy]);
  printf("\n");
  error->one(FLERR, "Incorrect split");
}

checksum = 0;
for (int iz = 0; iz < nelem[2]; iz++) 
checksum += cellsize[2][iz];
if (checksum > ncells[ietype][2]) {
  printf("Ncells = %d cellsize z =", ncells[ietype][2]);
  for (int iz = 0; iz < nelem[2]; iz++) 
    printf(" %d", cellsize[2][iz]);
  printf("\n");
  error->one(FLERR, "Incorrect split");
}

for (int ix = 0; ix < nelem[0]; ix++) {
  for (int iy = 0; iy < nelem[1]; iy++) {
    for (int iz = 0; iz < nelem[2]; iz++) {
      type_info[0] = cellsize[0][ix];
      type_info[3] = MIN(type_info[0], ngcells[ietype][0]);
      type_info[1] = cellsize[1][iy];
      type_info[4] = MIN(type_info[1], ngcells[ietype][1]);
      type_info[2] = cellsize[2][iz];
      type_info[5] = MIN(type_info[2], ngcells[ietype][2]);
      if (type_info[0] <= 0 || type_info[0] > ncells[ietype][0]) {
        printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
            , type_info[0]
            , type_info[1]
            , type_info[2]);
        error->one(FLERR, "TEST");
      }
      if (type_info[1] <= 0 || type_info[1] > ncells[ietype][1]) {
        printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
            , type_info[0]
            , type_info[1]
            , type_info[2]);
        error->one(FLERR, "TEST");
      }
      if (type_info[2] <= 0 || type_info[2] > ncells[ietype][2]) {
        printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
            , type_info[0]
            , type_info[1]
            , type_info[2]);
        error->one(FLERR, "TEST");
      }

      if (domain->dimension == 3) {

        // skip if 2 or more directions have size 1
        // will split to atoms

        if (((type_info[0] == 1) + (type_info[1] == 1) + (type_info[2] == 1)) > 1) 
          continue;

        // if x direction size is 1, yz -> xy
        if (type_info[0] == 1) {
          type_info[0] = type_info[1];
          type_info[1] = type_info[2];
          type_info[2] = 1;
          type_info[3] = type_info[4];
          type_info[4] = type_info[5];
          type_info[5] = 1;
          type_info[6] = Element::QUADRILATERAL;
        }

        // if y direction size is 1, xz -> xy
        if (type_info[1] == 1) {
          type_info[1] = type_info[2];
          type_info[2] = 1;
          type_info[4] = type_info[5];
          type_info[5] = 1;
          type_info[6] = Element::QUADRILATERAL;
        }

        if (type_info[2] == 1) 
          type_info[6] = Element::QUADRILATERAL;

      }
      request_new_etype(type_info);
      type_info[6] = Element::HEXAHEDRON;
    }
  }
}
n++;
ii++;
}
}

add_requested_etype(0);

// second pass: split elements

if (nsplit > 0) {
  int index_offset[3], index_node0, iucell;
  int myetype, ncellx, ncelly, ncellz;

  // iucell_list = list of ucell index in old element for each node of new element

  int *iucell_list = new int[element->maxnpe]; 

  for (ii = 0; ii < n; ii++) {

    int nalocal_previous = atom->nlocal;
    int nelocal_previous = element->nlocal;
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
      index_offset[0] = ncelly * ncellz * (cellsize[0][ix]-1);
      for (int iy = 0; iy < nelem[1]; iy++) {
        index_offset[1] = ncellz * (cellsize[1][iy]-1);
        if (domain->dimension == 3) {
          for (int iz = 0; iz < nelem[2]; iz++) {
            type_info[0] = cellsize[0][ix];
            type_info[3] = MIN(type_info[0], ngcells[ietype][0]);
            type_info[1] = cellsize[1][iy];
            type_info[4] = MIN(type_info[1], ngcells[ietype][1]);
            type_info[2] = cellsize[2][iz];
            type_info[5] = MIN(type_info[2], ngcells[ietype][2]);

            index_offset[2] = cellsize[2][iz]-1;

            //              printf("condition = %d cell = %d %d %d typeinfo = %d %d %d\n", ((type_info[0] == 1) + (type_info[1] == 1) + (type_info[2] == 1)), cellsize[0][ix], cellsize[1][iy], cellsize[2][iz], type_info[0], type_info[1], type_info[2]);
            if (((type_info[0] == 1) + (type_info[1] == 1) + (type_info[2] == 1)) > 1) {
              for (int iix = 0; iix < type_info[0]; iix++)
                for (int iiy = 0; iiy < type_info[1]; iiy++)
                  for (int iiz = 0; iiz < type_info[2]; iiz++) {
                    iucell = index_node0 
                      + iix * ncelly * ncellz 
                      + iiy * ncellz 
                      + iiz;
                    if (iucell > nucell[ietype]) {
                      printf("Split to atom iucell = %d nucell = %d ietype = %d ncells = %d %d %d\n"
                          , iucell, nucell[ietype], ietype, ncells[ietype][0], ncells[ietype][1], ncells[ietype][2]);
                      error->one(FLERR, "Wrong ucell index");
                    }
                    create_pass_on_atom(i, iucell, ictype, 0);
                  }
            } else {
              if (type_info[0] == 1) {
                type_info[0] = type_info[1];
                type_info[1] = type_info[2];
                type_info[2] = 1;
                type_info[3] = type_info[4];
                type_info[4] = type_info[5];
                type_info[5] = 1;
                type_info[6] = Element::QUADRILATERAL;
                iucell_list[0] = index_node0; 
                iucell_list[1] = index_node0 + index_offset[1];
                iucell_list[2] = index_node0 + index_offset[1] 
                  + index_offset[2];
                iucell_list[3] = index_node0 + index_offset[2];
              } else if (type_info[1] == 1) {
                type_info[1] = type_info[2];
                type_info[2] = 1;
                type_info[4] = type_info[5];
                type_info[5] = 1;
                type_info[6] = Element::QUADRILATERAL;
                iucell_list[0] = index_node0; 
                iucell_list[1] = index_node0 + index_offset[0];
                iucell_list[2] = index_node0 + index_offset[0] 
                  + index_offset[2];
                iucell_list[3] = index_node0 + index_offset[2];
              } else if (type_info[2] == 1) {
                type_info[6] = Element::QUADRILATERAL;
                iucell_list[0] = index_node0; 
                iucell_list[1] = index_node0 + index_offset[0];
                iucell_list[2] = index_node0 + index_offset[0] 
                  + index_offset[1];
                iucell_list[3] = index_node0 + index_offset[1];
              } else {

                type_info[6] = Element::HEXAHEDRON;
                // calculate ucell index in element I for each node of new element

                iucell_list[0] = index_node0; 
                iucell_list[1] = index_node0 + index_offset[0];
                iucell_list[2] = index_node0 + index_offset[0] 
                  + index_offset[1];
                iucell_list[3] = index_node0 + index_offset[1];
                iucell_list[4] = iucell_list[0] + index_offset[2];
                iucell_list[5] = iucell_list[1] + index_offset[2];
                iucell_list[6] = iucell_list[2] + index_offset[2];
                iucell_list[7] = iucell_list[3] + index_offset[2];

              }
              myetype = find_etype(type_info);
              if (myetype == 0) {
                printf("type_info = %d %d %d %d %d %d %d %d\n"
                    , type_info[0]
                    , type_info[1]
                    , type_info[2]
                    , type_info[3]
                    , type_info[4]
                    , type_info[5]
                    , type_info[6]
                    , type_info[7]);
                error->one(FLERR, "Etype not defined yet");
              } 

              for (int aa = 0; aa < element->npe[myetype]; aa++)
                if (iucell_list[aa] > nucell[ietype]) {
                  printf("Split to element iucell_list[%d] = %d nucell = %d ietype = %d ncells = %d %d %d\n"
                      , aa, iucell_list[aa], nucell[ietype], ietype, ncells[ietype][0], ncells[ietype][1], ncells[ietype][2]);
                  error->one(FLERR, "Wrong ucell_list");
                }

              //printf("myetype = %d type_info[6] = %s ncells = %d %d %d npe = %d\n", myetype, element->element_shape_list[element->element_shape_ids[myetype]], ncells[myetype][0], ncells[myetype][1], ncells[myetype][2], element->npe[myetype]);
              create_pass_on_element(i, iucell_list, myetype, ictype, 0);
            } 
            index_node0 += index_offset[2]+1;
          }
          index_node0 += index_offset[1];
        } else {
          type_info[0] = cellsize[0][ix];
          type_info[3] = MIN(type_info[0], ngcells[ietype][0]);
          type_info[1] = cellsize[1][iy];
          type_info[4] = MIN(type_info[1], ngcells[ietype][1]);

          type_info[2] = type_info[5] = 1;
          if (type_info[0] == 1) {
            for (int iiy = 0; iiy < type_info[1]; iiy++)
              create_pass_on_atom(i, index_node0+iiy, ictype, 0);
          } else if (type_info[1] == 1) {
            error->one(FLERR, "TEST"); 
            //            for (int iix = 0; iix < info[0]; iix++) 
            //              create_pass_on_atom(i, index_node0+iix * ncellyz, ictype, 0);
          } else {
            myetype = find_etype(type_info);
            if (myetype == 0) error->one(FLERR, "Etype not defined yet");
            iucell_list[0] = index_node0; 
            iucell_list[1] = index_node0 + index_offset[0];
            iucell_list[2] = index_node0 + index_offset[0] 
              + index_offset[1];
            iucell_list[3] = index_node0 + index_offset[1];
            create_pass_on_element(i, iucell_list, myetype, ictype, 0);
          }
          index_node0 += index_offset[1]+1;
        }
      }
      index_node0 += index_offset[0];
    }

    // check mass conservation 
    // number of atoms in splitted elements and atoms must be equal 
    // to the number of atoms in the original element
    int sum = atom->nlocal - nalocal_previous;
    for (int j = nelocal_previous; j < element->nlocal; j++)
      sum += nucell[etype[j]];
    if (sum != nucell[ietype]) {
      printf("split element tag i = %d i = %d ietype = %d ncells = %d %d %d sum = %d new atom = %d\n", element->tag[i], i, ietype, ncellx, ncelly, ncellz, sum, atom->nlocal-nalocal_previous);
      for (int j = 0; j < nsplit; j++)
        if (splitlist[j][0] == i)
          printf("split %d %d\n", splitlist[j][1], splitlist[j][2]);
      printf("Ncells = %d cellsize x =", ncells[ietype][0]);
      for (int iz = 0; iz < nelem[0]; iz++) 
        printf(" %d", cellsize[0][iz]);
      printf("\n");
      printf("Ncells = %d cellsize y =", ncells[ietype][1]);
      for (int iz = 0; iz < nelem[1]; iz++) 
        printf(" %d", cellsize[1][iz]);
      printf("\n");
      printf("Ncells = %d cellsize z =", ncells[ietype][2]);
      for (int iz = 0; iz < nelem[2]; iz++) 
        printf(" %d", cellsize[2][iz]);
      printf("\n");

      error->one(FLERR, "Mass not conserved during splitting");
    } 

    // delete element I

    copy(element->nlocal-1, i, 1);
    element->nlocal--;

  }
  memory->destroy(ilist); 
  memory->destroy(nelem_array); 
  memory->destroy(cellsize_array);
  delete [] iucell_list;
}
//error->all(FLERR, "TEST"); 
//error->one(FLERR, "TEST");
return n;
*/
}

/*  ----------------------------------------------------------------------
    create new element from an old element
    -------------------------------------------------------------------------  */

void ElementVecCAC::create_pass_on_element(int i, int *iucell_list, int ietype, int *ictype, tagint itag)
{
  double ***nodecoord;
  int inpe = element->npe[ietype];
  int iapc = element->apc[ietype];
  memory->create(nodecoord, iapc, inpe, 3, "evec:nodecoord");
  for (int ibasis = 0; ibasis < iapc; ibasis++)
    for (int inode = 0; inode < inpe; inode++)
      interpolate(nodecoord[ibasis][inode], nodex, i, ibasis, iucell_list[inode], 3);

  create_element(nodecoord, ietype, ictype, itag);

  // pass on stored data for new element in fix/compute/variable (if needed) and velocity
  // new elements will be in the same group as old element

  int j = element->nlocal-1;
  element->data_pass_on_fix_compute_dump_variable(j, i, iucell_list);
  mask[j] = mask[i];
  image[j] = image[i];
  for (int ibasis = 0; ibasis < iapc; ibasis++)
    for (int inode = 0; inode < inpe; inode++)
      interpolate(nodev[j][ibasis][inode], nodev, i, ibasis, iucell_list[inode], 3);
  memory->destroy(nodecoord);
}

/*  ----------------------------------------------------------------------
    create new atom from an old element
    -------------------------------------------------------------------------  */

void ElementVecCAC::create_pass_on_atom(int i, int ibasis, int iucell, int itype, tagint itag)
{
  double coord[3];
  interpolate(coord, nodex, i, ibasis, iucell, 3);
  atom->avec->create_atom(coord, itype, itag);

  // pass on stored data for new element in fix/compute/variable and velocity
  // new atoms will be in the same group as old element

  int j = atom->nlocal-1;
  int iindex = iucell * apc[etype[i]] + ibasis;
  atom->data_pass_on_fix_compute_dump_variable(j, i, iindex);
  atom->mask[j] = mask[i];
  atom->image[j] = image[i];
  interpolate(atom->v[j], nodev, i, ibasis, iucell, 3);
}

/*  ----------------------------------------------------------------------
    split hex element I to 2 hex elements
    scan_flag = 1 --> check for new etype
    scan_flag = 0 --> split element
    -------------------------------------------------------------------------  */

void ElementVecCAC::hex2hex(int i, int split_cell, int split_dim, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  int type_info[8];
  int left_etype, right_etype;
  type_info[0] = ncells[ietype][0];
  type_info[1] = ncells[ietype][1];
  type_info[2] = ncells[ietype][2];
  type_info[3] = ngcells[ietype][0];
  type_info[4] = ngcells[ietype][1];
  type_info[5] = ngcells[ietype][2];
  type_info[6] = Element::HEXAHEDRON;
  type_info[7] = iapc;

  if (split_cell <= 1 || split_cell >= ncells[ietype][split_dim]-1) 
    error->one(FLERR, "Illegal split_cell");
  type_info[split_dim] = split_cell;
  type_info[split_dim+3] = MIN(type_info[split_dim], ngcells[ietype][split_dim]);

  if (scan_flag) request_new_etype(type_info);
  else left_etype = find_etype(type_info);

  type_info[split_dim] = ncells[ietype][split_dim]-split_cell;
  type_info[split_dim+3] = MIN(type_info[split_dim], ngcells[ietype][split_dim]);

  if (scan_flag) request_new_etype(type_info);
  else right_etype = find_etype(type_info);

  if (scan_flag) return;

  double **leftnodecoord, **rightnodecoord;
  double **inodex;
  int ncell = ncells[ietype][split_dim]-1;
  tagint itag = tag[i];

  int *ictype = ctype[i];

  double ***leftnodes, ***rightnodes;
  memory->create(leftnodes, iapc, 8, 3, "evec:leftnode");
  memory->create(rightnodes, iapc, 8, 3, "evec:rightnode");

  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    leftnodecoord = leftnodes[ibasis];
    rightnodecoord = rightnodes[ibasis];
    if (split_dim == 0) {
      copy3(inodex[0], leftnodecoord[0]);
      copy3(inodex[3], leftnodecoord[3]);
      copy3(inodex[4], leftnodecoord[4]);
      copy3(inodex[7], leftnodecoord[7]);

      for (int dim = 0; dim < 3; dim++) {
        leftnodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell * (split_cell-1);
        leftnodecoord[2][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell * (split_cell-1);
        leftnodecoord[5][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell * (split_cell-1);
        leftnodecoord[6][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell * (split_cell-1);
      }

      copy3(inodex[1], rightnodecoord[1]);
      copy3(inodex[2], rightnodecoord[2]);
      copy3(inodex[5], rightnodecoord[5]);
      copy3(inodex[6], rightnodecoord[6]);
      for (int dim = 0; dim < 3; dim++) {
        rightnodecoord[0][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell * split_cell;
        rightnodecoord[3][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell * split_cell;
        rightnodecoord[4][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell * split_cell;
        rightnodecoord[7][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell * split_cell;
      }
    } else if (split_dim == 1) {
      copy3(inodex[0], leftnodecoord[0]);
      copy3(inodex[1], leftnodecoord[1]);
      copy3(inodex[4], leftnodecoord[4]);
      copy3(inodex[5], leftnodecoord[5]);
      for (int dim = 0; dim < 3; dim++) {
        leftnodecoord[3][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell * (split_cell-1);
        leftnodecoord[2][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell * (split_cell-1);
        leftnodecoord[7][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell * (split_cell-1);
        leftnodecoord[6][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell * (split_cell-1);
      }

      copy3(inodex[2], rightnodecoord[2]);
      copy3(inodex[3], rightnodecoord[3]);
      copy3(inodex[6], rightnodecoord[6]);
      copy3(inodex[7], rightnodecoord[7]);
      for (int dim = 0; dim < 3; dim++) {
        rightnodecoord[0][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell * split_cell;
        rightnodecoord[1][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell * split_cell;
        rightnodecoord[4][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell * split_cell;
        rightnodecoord[5][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell * split_cell;
      }

    } else if (split_dim == 2) {
      copy3(inodex[0], leftnodecoord[0]);
      copy3(inodex[1], leftnodecoord[1]);
      copy3(inodex[2], leftnodecoord[2]);
      copy3(inodex[3], leftnodecoord[3]);
      for (int dim = 0; dim < 3; dim++) {
        leftnodecoord[4][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell * (split_cell-1);
        leftnodecoord[5][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell * (split_cell-1);
        leftnodecoord[6][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell * (split_cell-1);
        leftnodecoord[7][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell * (split_cell-1);
      }

      copy3(inodex[4], rightnodecoord[4]);
      copy3(inodex[5], rightnodecoord[5]);
      copy3(inodex[6], rightnodecoord[6]);
      copy3(inodex[7], rightnodecoord[7]);
      for (int dim = 0; dim < 3; dim++) {
        rightnodecoord[0][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell * split_cell;
        rightnodecoord[1][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell * split_cell;
        rightnodecoord[2][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell * split_cell;
        rightnodecoord[3][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell * split_cell;
      }

    }
  }
  create_element(rightnodes, right_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(leftnodes, left_etype, ictype, 0);
  mask[nlocal+1] = mask[i];

  memory->destroy(rightnodes);
  memory->destroy(leftnodes);
}

/*  ----------------------------------------------------------------------
    split hex element I to 1 hex and 1 quad elements
    scan_flag = 1 --> check for new etype
    scan_flag = 0 --> split element
    -------------------------------------------------------------------------  */

void ElementVecCAC::hex2hexquad(int i, int split_cell, int split_dim, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  int type_info[8];
  int left_etype, right_etype;
  type_info[0] = ncells[ietype][0];
  type_info[1] = ncells[ietype][1];
  type_info[2] = ncells[ietype][2];
  type_info[3] = ngcells[ietype][0];
  type_info[4] = ngcells[ietype][1];
  type_info[5] = ngcells[ietype][2];
  type_info[6] = Element::HEXAHEDRON;
  type_info[7] = iapc;

  if (split_cell == 0) split_cell = 1;
  else if (split_cell == 1) split_cell = ncells[ietype][split_dim]-1;

  if (split_cell == 1) {
    type_info[6] = Element::QUADRILATERAL;
    if (split_dim == 0) {
      type_info[0] = ncells[ietype][1];
      type_info[1] = ncells[ietype][2];
    } else if (split_dim == 1) {
      type_info[0] = ncells[ietype][0];
      type_info[1] = ncells[ietype][2];
    } else if (split_dim == 2) {
      type_info[0] = ncells[ietype][0];
      type_info[1] = ncells[ietype][1];
    }
    type_info[2] = type_info[5] = 1;
  } else {
    type_info[split_dim] = split_cell;
    type_info[split_dim+3] = MIN(type_info[split_dim], ngcells[ietype][split_dim]);
  }

  if (scan_flag) request_new_etype(type_info);
  else left_etype = find_etype(type_info);

  type_info[0] = ncells[ietype][0];
  type_info[1] = ncells[ietype][1];
  type_info[2] = ncells[ietype][2];
  type_info[3] = ngcells[ietype][0];
  type_info[4] = ngcells[ietype][1];
  type_info[5] = ngcells[ietype][2];
  type_info[6] = Element::HEXAHEDRON;
  type_info[7] = apc[ietype];

  if (ncells[ietype][split_dim] - split_cell == 1) {
    type_info[6] = Element::QUADRILATERAL;
    if (split_dim == 0) {
      type_info[0] = ncells[ietype][1];
      type_info[1] = ncells[ietype][2];
    } else if (split_dim == 1) {
      type_info[0] = ncells[ietype][0];
      type_info[1] = ncells[ietype][2];
    } else if (split_dim == 2) {
      type_info[0] = ncells[ietype][0];
      type_info[1] = ncells[ietype][1];
    }
    type_info[2] = type_info[5] = 1;
  } else {
    type_info[split_dim] = ncells[ietype][split_dim]-split_cell;
    type_info[split_dim+3] = MIN(type_info[split_dim], ngcells[ietype][split_dim]);
  }
  if (scan_flag) request_new_etype(type_info);
  else right_etype = find_etype(type_info);

  if (scan_flag) return;

  double **leftnodecoord, **rightnodecoord;
  double **inodex;
  int ncell = ncells[ietype][split_dim]-1;
  tagint itag = tag[i];

  int *ictype = ctype[i];
  double ***leftnodes, ***rightnodes;
  memory->create(leftnodes, iapc, 8, 3, "evec:leftnode");
  memory->create(rightnodes, iapc, 8, 3, "evec:rightnode");

  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    leftnodecoord = leftnodes[ibasis];
    rightnodecoord = rightnodes[ibasis];

    if (split_dim == 0) {
      if (split_cell == 1) {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[1], leftnodecoord[3]);
        copy3(inodex[2], leftnodecoord[7]);
        copy3(inodex[3], leftnodecoord[4]);
      } else {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[3], leftnodecoord[3]);
        copy3(inodex[4], leftnodecoord[4]);
        copy3(inodex[7], leftnodecoord[7]);

        for (int dim = 0; dim < 3; dim++) {
          leftnodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell * (split_cell-1);
          leftnodecoord[2][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell * (split_cell-1);
          leftnodecoord[5][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell * (split_cell-1);
          leftnodecoord[6][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell * (split_cell-1);
        }
      }

      if (split_cell == ncell -1) {
        copy3(inodex[0], leftnodecoord[1]);
        copy3(inodex[1], leftnodecoord[2]);
        copy3(inodex[2], leftnodecoord[6]);
        copy3(inodex[3], leftnodecoord[5]);
      } else { 
        copy3(inodex[1], rightnodecoord[1]);
        copy3(inodex[2], rightnodecoord[2]);
        copy3(inodex[5], rightnodecoord[5]);
        copy3(inodex[6], rightnodecoord[6]);
        for (int dim = 0; dim < 3; dim++) {
          rightnodecoord[0][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell * split_cell;
          rightnodecoord[3][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell * split_cell;
          rightnodecoord[4][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell * split_cell;
          rightnodecoord[7][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell * split_cell;
        }
      }
    } else if (split_dim == 1) {
      if (split_cell == 1) {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[1], leftnodecoord[1]);
        copy3(inodex[2], leftnodecoord[5]);
        copy3(inodex[3], leftnodecoord[4]);
      } else {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[1], leftnodecoord[1]);
        copy3(inodex[4], leftnodecoord[4]);
        copy3(inodex[5], leftnodecoord[5]);
        for (int dim = 0; dim < 3; dim++) {
          leftnodecoord[3][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell * (split_cell-1);
          leftnodecoord[2][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell * (split_cell-1);
          leftnodecoord[7][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell * (split_cell-1);
          leftnodecoord[6][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell * (split_cell-1);
        }
      }
      if (split_cell == ncell -1) {
        copy3(inodex[0], leftnodecoord[3]);
        copy3(inodex[1], leftnodecoord[2]);
        copy3(inodex[2], leftnodecoord[6]);
        copy3(inodex[3], leftnodecoord[7]);
      } else { 
        copy3(inodex[2], rightnodecoord[2]);
        copy3(inodex[3], rightnodecoord[3]);
        copy3(inodex[6], rightnodecoord[6]);
        copy3(inodex[7], rightnodecoord[7]);
        for (int dim = 0; dim < 3; dim++) {
          rightnodecoord[0][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell * split_cell;
          rightnodecoord[1][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell * split_cell;
          rightnodecoord[4][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell * split_cell;
          rightnodecoord[5][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell * split_cell;
        }
      }
    } else if (split_dim == 2) {
      if (split_cell == 1) {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[1], leftnodecoord[1]);
        copy3(inodex[2], leftnodecoord[2]);
        copy3(inodex[3], leftnodecoord[3]);
      } else {
        copy3(inodex[0], leftnodecoord[0]);
        copy3(inodex[1], leftnodecoord[1]);
        copy3(inodex[2], leftnodecoord[2]);
        copy3(inodex[3], leftnodecoord[3]);
        for (int dim = 0; dim < 3; dim++) {
          leftnodecoord[4][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell * (split_cell-1);
          leftnodecoord[5][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell * (split_cell-1);
          leftnodecoord[6][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell * (split_cell-1);
          leftnodecoord[7][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell * (split_cell-1);
        }
      }
      if (split_cell == ncell -1) {
        copy3(inodex[0], leftnodecoord[4]);
        copy3(inodex[1], leftnodecoord[5]);
        copy3(inodex[2], leftnodecoord[6]);
        copy3(inodex[3], leftnodecoord[7]);
      } else { 
        copy3(inodex[4], rightnodecoord[4]);
        copy3(inodex[5], rightnodecoord[5]);
        copy3(inodex[6], rightnodecoord[6]);
        copy3(inodex[7], rightnodecoord[7]);
        for (int dim = 0; dim < 3; dim++) {
          rightnodecoord[0][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell * split_cell;
          rightnodecoord[1][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell * split_cell;
          rightnodecoord[2][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell * split_cell;
          rightnodecoord[3][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell * split_cell;
        }
      }
    }
  }
  create_element(rightnodes, right_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(leftnodes, left_etype, ictype, 0);
  mask[nlocal+1] = mask[i];

  memory->destroy(rightnodes);
  memory->destroy(leftnodes);
}

/*  ----------------------------------------------------------------------
    split hex element I to 2 wedge elements
    scan_flag = 1 --> check for new etype
    = 0 --> split element
    split_dir = 0 --> split across diagonal with slope = -1
    = 1 --> split across the other diagonal
    larger_side = 0 --> larger side is lower half
    = 1 --> larger side is upper half
    -------------------------------------------------------------------------  */

void ElementVecCAC::hex2wedge(int i, int split_axis, int split_dir, int larger_side, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  if (element->element_shape_ids[ietype] != Element::HEXAHEDRON)
    return;
  int larger_etype, smaller_etype;
  int type_info[8];
  if ((split_axis == 0 && ncells[ietype][1] != ncells[ietype][2]) || 
      (split_axis == 1 && ncells[ietype][0] != ncells[ietype][2]) ||
      (split_axis == 2 && ncells[ietype][1] != ncells[ietype][0]))
    error->one(FLERR, "Hex element is not square to split");

  // larger wedge

  type_info[0] = ncells[ietype][split_axis];
  type_info[1] = ncells[ietype][(split_axis+1)%3];
  type_info[2] = 0;
  type_info[3] = ngcells[ietype][split_axis];
  if (type_info[1] >= 4) type_info[4] = 4;
  else type_info[4] = type_info[1];
  type_info[5] = 0;
  type_info[6] = Element::WEDGE;
  type_info[7] = iapc;

  if (scan_flag) request_new_etype(type_info);
  else larger_etype = find_etype(type_info);

  // smaller wedge

  type_info[1] -= 1;
  if (type_info[1] >= 4) type_info[4] = 4;
  else type_info[4] = type_info[1];

  if (scan_flag) request_new_etype(type_info);
  else smaller_etype = find_etype(type_info);

  if (scan_flag) return;

  double **largernodecoord;
  double **smallernodecoord;
  double **inodex;
  tagint itag = tag[i];

  int ncell = type_info[1];
  int *ictype = ctype[i];
  double ***largernodes, ***smallernodes;
  memory->create(largernodes, iapc, 6, 3, "evec:largernodes");
  memory->create(smallernodes, iapc, 6, 3, "evec:smallernodes");

  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    largernodecoord = largernodes[ibasis];
    smallernodecoord = smallernodes[ibasis];

    if (split_axis == 0) {
      if (split_dir == 0) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[0], largernodecoord[0]);
          copy3(inodex[3], largernodecoord[1]);
          copy3(inodex[4], largernodecoord[2]);
          copy3(inodex[1], largernodecoord[3]);
          copy3(inodex[2], largernodecoord[4]);
          copy3(inodex[5], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[7], smallernodecoord[0]);
          copy3(inodex[6], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell;
            smallernodecoord[2][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell;
            smallernodecoord[4][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell;
            smallernodecoord[5][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[7], largernodecoord[0]);
          copy3(inodex[4], largernodecoord[1]);
          copy3(inodex[3], largernodecoord[2]);
          copy3(inodex[6], largernodecoord[3]);
          copy3(inodex[5], largernodecoord[4]);
          copy3(inodex[2], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[0], smallernodecoord[0]);
          copy3(inodex[1], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
            smallernodecoord[2][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;
            smallernodecoord[4][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
            smallernodecoord[5][dim] = inodex[5][dim] + (inodex[1][dim] - inodex[5][dim])/ncell;
          }
        }
      } else if (split_dir == 1) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[3], largernodecoord[0]);
          copy3(inodex[0], largernodecoord[1]);
          copy3(inodex[7], largernodecoord[2]);
          copy3(inodex[2], largernodecoord[3]);
          copy3(inodex[1], largernodecoord[4]);
          copy3(inodex[6], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[4], smallernodecoord[0]);
          copy3(inodex[5], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[7][dim] + (inodex[4][dim] - inodex[7][dim])/ncell;
            smallernodecoord[2][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell;
            smallernodecoord[4][dim] = inodex[6][dim] + (inodex[5][dim] - inodex[6][dim])/ncell;
            smallernodecoord[5][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[4], largernodecoord[0]);
          copy3(inodex[7], largernodecoord[1]);
          copy3(inodex[0], largernodecoord[2]);
          copy3(inodex[5], largernodecoord[3]);
          copy3(inodex[6], largernodecoord[4]);
          copy3(inodex[1], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[3], smallernodecoord[0]);
          copy3(inodex[2], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
            smallernodecoord[2][dim] = inodex[7][dim] + (inodex[3][dim] - inodex[7][dim])/ncell;
            smallernodecoord[4][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
            smallernodecoord[5][dim] = inodex[6][dim] + (inodex[2][dim] - inodex[6][dim])/ncell;
          }
        }

      }
    } else if (split_axis == 1) {
      if (split_dir == 0) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[0], largernodecoord[0]);
          copy3(inodex[4], largernodecoord[1]);
          copy3(inodex[1], largernodecoord[2]);
          copy3(inodex[3], largernodecoord[3]);
          copy3(inodex[7], largernodecoord[4]);
          copy3(inodex[2], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[5], smallernodecoord[0]);
          copy3(inodex[6], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell;
            smallernodecoord[2][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;
            smallernodecoord[4][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell;
            smallernodecoord[5][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[5], largernodecoord[0]);
          copy3(inodex[1], largernodecoord[1]);
          copy3(inodex[4], largernodecoord[2]);
          copy3(inodex[6], largernodecoord[3]);
          copy3(inodex[2], largernodecoord[4]);
          copy3(inodex[7], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[0], smallernodecoord[0]);
          copy3(inodex[3], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;
            smallernodecoord[2][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
            smallernodecoord[4][dim] = inodex[7][dim] + (inodex[3][dim] - inodex[7][dim])/ncell;
            smallernodecoord[5][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
          }
        }
      } else if (split_dir == 1) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[1], largernodecoord[0]);
          copy3(inodex[0], largernodecoord[1]);
          copy3(inodex[5], largernodecoord[2]);
          copy3(inodex[2], largernodecoord[3]);
          copy3(inodex[3], largernodecoord[4]);
          copy3(inodex[6], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[4], smallernodecoord[0]);
          copy3(inodex[7], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;
            smallernodecoord[2][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell;
            smallernodecoord[4][dim] = inodex[6][dim] + (inodex[7][dim] - inodex[6][dim])/ncell;
            smallernodecoord[5][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[4], largernodecoord[0]);
          copy3(inodex[5], largernodecoord[1]);
          copy3(inodex[0], largernodecoord[2]);
          copy3(inodex[7], largernodecoord[3]);
          copy3(inodex[6], largernodecoord[4]);
          copy3(inodex[3], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[1], smallernodecoord[0]);
          copy3(inodex[2], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
            smallernodecoord[2][dim] = inodex[5][dim] + (inodex[1][dim] - inodex[5][dim])/ncell;
            smallernodecoord[4][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
            smallernodecoord[5][dim] = inodex[6][dim] + (inodex[2][dim] - inodex[6][dim])/ncell;
          }
        }

      }

    } else if (split_axis == 2) {
      if (split_dir == 0) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[0], largernodecoord[0]);
          copy3(inodex[1], largernodecoord[1]);
          copy3(inodex[3], largernodecoord[2]);
          copy3(inodex[4], largernodecoord[3]);
          copy3(inodex[5], largernodecoord[4]);
          copy3(inodex[7], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[2], smallernodecoord[0]);
          copy3(inodex[6], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
            smallernodecoord[2][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
            smallernodecoord[4][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell;
            smallernodecoord[5][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[2], largernodecoord[0]);
          copy3(inodex[3], largernodecoord[1]);
          copy3(inodex[1], largernodecoord[2]);
          copy3(inodex[6], largernodecoord[3]);
          copy3(inodex[7], largernodecoord[4]);
          copy3(inodex[5], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[0], smallernodecoord[0]);
          copy3(inodex[4], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
            smallernodecoord[2][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
            smallernodecoord[4][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;
            smallernodecoord[5][dim] = inodex[7][dim] + (inodex[4][dim] - inodex[7][dim])/ncell;
          }
        }
      } else if (split_dir == 1) {
        if (larger_side == 0) {

          // larger wedge           

          copy3(inodex[1], largernodecoord[0]);
          copy3(inodex[2], largernodecoord[1]);
          copy3(inodex[0], largernodecoord[2]);
          copy3(inodex[5], largernodecoord[3]);
          copy3(inodex[6], largernodecoord[4]);
          copy3(inodex[4], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[3], smallernodecoord[0]);
          copy3(inodex[7], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
            smallernodecoord[2][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
            smallernodecoord[4][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell;
            smallernodecoord[5][dim] = inodex[6][dim] + (inodex[7][dim] - inodex[6][dim])/ncell;
          }

        } else if (larger_side == 1) {

          // larger wedge           

          copy3(inodex[3], largernodecoord[0]);
          copy3(inodex[0], largernodecoord[1]);
          copy3(inodex[2], largernodecoord[2]);
          copy3(inodex[7], largernodecoord[3]);
          copy3(inodex[4], largernodecoord[4]);
          copy3(inodex[6], largernodecoord[5]);

          // smaller wedge

          copy3(inodex[1], smallernodecoord[0]);
          copy3(inodex[5], smallernodecoord[3]);

          for (int dim = 0; dim < 3; dim++) {
            smallernodecoord[1][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
            smallernodecoord[2][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
            smallernodecoord[4][dim] = inodex[6][dim] + (inodex[5][dim] - inodex[6][dim])/ncell;
            smallernodecoord[5][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;
          }
        }
      }
    } else error->one(FLERR, "Incorrect split_axis");
  }
  create_element(largernodes, larger_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(smallernodes, smaller_etype, ictype, 0);
  mask[nlocal+1] = mask[i];

  memory->destroy(largernodes);
  memory->destroy(smallernodes);

}

/*  ----------------------------------------------------------------------
    split wedge element I to 1 pyramid and 1 tetrahedron elements
    scan_flag = 1 --> check for new etype
    = 0 --> split element
    split_dir = 0 --> split across diagonal with slope = -1
    = 1 --> split across the other diagonal
    -------------------------------------------------------------------------  */

void ElementVecCAC::wedge2pyrtet(int i, int base_axis, int split_dir, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  if (element->element_shape_ids[ietype] != Element::WEDGE)
    return;
  int pyr_etype, tet_etype;
  int type_info[8];
  if (ncells[ietype][1] != ncells[ietype][0])
    error->one(FLERR, "Wedge element is not square to split");

  // pyramid

  type_info[0] = ncells[ietype][0];
  type_info[1] = 0;
  type_info[2] = 0;
  type_info[3] = MIN(type_info[0], 4);
  type_info[4] = 0;
  type_info[5] = 0;
  type_info[6] = Element::PYRAMID;
  type_info[7] = iapc;

  if (scan_flag) request_new_etype(type_info);
  else pyr_etype = find_etype(type_info);

  // tetrahedron

  type_info[0]--;
  type_info[3] = MIN(type_info[0], 4);
  type_info[6] = Element::TETRAHEDRON;

  if (scan_flag) request_new_etype(type_info);
  else tet_etype = find_etype(type_info);

  if (scan_flag) return;

  double **pyrnodecoord;
  double **tetnodecoord;
  double **inodex;
  tagint itag = tag[i];

  int ncell = ncells[ietype][0]-1;
  int *ictype = ctype[i];
  double ***pyrnodes, ***tetnodes;
  memory->create(pyrnodes, iapc, 5, 3, "evec:pyrnodes");
  memory->create(tetnodes, iapc, 4, 3, "evec:tetnodes");

  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    pyrnodecoord = pyrnodes[ibasis];
    tetnodecoord = tetnodes[ibasis];
    if (base_axis == 0) {

      // pyramid base

      copy3(inodex[1], pyrnodecoord[0]);
      copy3(inodex[2], pyrnodecoord[1]);
      copy3(inodex[5], pyrnodecoord[2]);
      copy3(inodex[4], pyrnodecoord[3]);

      if (split_dir == 0) {

        // pyramid peak

        copy3(inodex[0], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[3], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
          tetnodecoord[2][dim] = inodex[4][dim] + (inodex[3][dim] - inodex[4][dim])/ncell;
          tetnodecoord[3][dim] = inodex[5][dim] + (inodex[3][dim] - inodex[5][dim])/ncell;
        }
      } else if (split_dir == 1) {

        // pyramid peak

        copy3(inodex[3], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[0], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
          tetnodecoord[2][dim] = inodex[2][dim] + (inodex[0][dim] - inodex[2][dim])/ncell;
          tetnodecoord[3][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
        }
      }
    } else if (base_axis == 1) {

      // pyramid base

      copy3(inodex[0], pyrnodecoord[0]);
      copy3(inodex[2], pyrnodecoord[1]);
      copy3(inodex[5], pyrnodecoord[2]);
      copy3(inodex[3], pyrnodecoord[3]);

      if (split_dir == 0) {

        // pyramid peak

        copy3(inodex[1], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[4], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[1][dim] + (inodex[4][dim] - inodex[1][dim])/ncell;
          tetnodecoord[2][dim] = inodex[3][dim] + (inodex[4][dim] - inodex[3][dim])/ncell;
          tetnodecoord[3][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;
        }
      } else if (split_dir == 1) {

        // pyramid peak

        copy3(inodex[4], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[1], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
          tetnodecoord[2][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
          tetnodecoord[3][dim] = inodex[4][dim] + (inodex[1][dim] - inodex[4][dim])/ncell;
        }
      }
    } else if (base_axis == 2) {

      // pyramid base

      copy3(inodex[0], pyrnodecoord[0]);
      copy3(inodex[1], pyrnodecoord[1]);
      copy3(inodex[4], pyrnodecoord[2]);
      copy3(inodex[3], pyrnodecoord[3]);

      if (split_dir == 0) {

        // pyramid peak

        copy3(inodex[2], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[5], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[2][dim] + (inodex[5][dim] - inodex[2][dim])/ncell;
          tetnodecoord[2][dim] = inodex[3][dim] + (inodex[5][dim] - inodex[3][dim])/ncell;
          tetnodecoord[3][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;
        }
      } else if (split_dir == 1) {

        // pyramid peak

        copy3(inodex[5], pyrnodecoord[4]);

        // tetrahedron

        copy3(inodex[2], tetnodecoord[0]);

        for (int dim = 0; dim < 3; dim++) {
          tetnodecoord[1][dim] = inodex[0][dim] + (inodex[2][dim] - inodex[0][dim])/ncell;
          tetnodecoord[2][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
          tetnodecoord[3][dim] = inodex[5][dim] + (inodex[2][dim] - inodex[5][dim])/ncell;
        }
      }
    } else error->one(FLERR, "Incorrect split_axis");
  }
  create_element(pyrnodes, pyr_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(tetnodes, tet_etype, ictype, 0);
  mask[nlocal+1] = mask[i];

  memory->destroy(pyrnodes);
  memory->destroy(tetnodes);

}

/*  ----------------------------------------------------------------------
    transform pyramid element I to tetrahedron elements
    scan_flag = 1 --> check for new etype
    = 0 --> split element
    split_dir = 0 --> split across diagonal with slope = -1
    = 1 --> split across the other diagonal
    larger_side = 0 --> larger side is lower half
    = 1 --> larger side is upper half
    -------------------------------------------------------------------------  */

void ElementVecCAC::pyr2tet(int i, int split_dir, int larger_side, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  if (element->element_shape_ids[ietype] != Element::PYRAMID)
    return;
  int larger_etype, smaller_etype;
  int type_info[8];

  // larger tetrahedron 

  type_info[0] = ncells[ietype][0];
  type_info[1] = 0;
  type_info[2] = 0;
  type_info[3] = MIN(type_info[0], 4);
  type_info[4] = 0;
  type_info[5] = 0;
  type_info[6] = Element::TETRAHEDRON;
  type_info[7] = iapc;

  if (scan_flag) request_new_etype(type_info);
  else larger_etype = find_etype(type_info);

  // smaller tetrahedron 

  type_info[0] = ncells[ietype][0]-1;
  type_info[1] = 0;
  type_info[2] = 0;
  type_info[3] = MIN(type_info[0], 4);
  type_info[4] = 0;
  type_info[5] = 0;
  type_info[6] = Element::TETRAHEDRON;
  type_info[7] = apc[ietype];

  if (scan_flag) request_new_etype(type_info);
  else smaller_etype = find_etype(type_info);

  if (scan_flag) return;

  double **smallernodecoord;
  double **largernodecoord;
  double **inodex;
  tagint itag = tag[i];

  int ncell = ncells[ietype][0]-1;
  int *ictype = ctype[i];
  double ***smallernodes;
  double ***largernodes;
  memory->create(smallernodes, iapc, 4, 3, "evec:smallertnodes");
  memory->create(largernodes, iapc, 4, 3, "evec:largertnodes");
  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    largernodecoord = largernodes[ibasis];
    smallernodecoord = smallernodes[ibasis];

    copy3(inodex[4], largernodecoord[3]);

    double delx, dely, delz;
    delx = inodex[1][0]-inodex[3][0];
    dely = inodex[1][1]-inodex[3][1];
    delz = inodex[1][2]-inodex[3][2];
    double rsq0 = delx * delx + dely * dely + delz * delz;
    delx = inodex[0][0]-inodex[2][0];
    dely = inodex[0][1]-inodex[2][1];
    delz = inodex[0][2]-inodex[2][2];
    double rsq1 = delx * delx + dely * dely + delz * delz;

    if (split_dir == 0 || (split_dir == 2 && rsq0 < rsq1)) {
      if (larger_side == 0) {

        // larger tetrahedral

        copy3(inodex[0], largernodecoord[0]);
        copy3(inodex[1], largernodecoord[1]);
        copy3(inodex[3], largernodecoord[2]);

        // smaller tetrahedral

        copy3(inodex[2], smallernodecoord[0]);
        for (int dim = 0; dim < 3; dim++) {
          smallernodecoord[1][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
          smallernodecoord[2][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
          smallernodecoord[3][dim] = inodex[4][dim] + (inodex[2][dim] - inodex[4][dim])/ncell;
        }
      } else if (larger_side == 1) {

        // larger tetrahedral

        copy3(inodex[1], largernodecoord[0]);
        copy3(inodex[2], largernodecoord[1]);
        copy3(inodex[3], largernodecoord[2]);

        // smaller tetrahedral

        copy3(inodex[0], smallernodecoord[0]);
        for (int dim = 0; dim < 3; dim++) {
          smallernodecoord[1][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
          smallernodecoord[2][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
          smallernodecoord[3][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;
        }
      }
    } else {
      if (larger_side == 0) {

        // larger tetrahedral

        copy3(inodex[0], largernodecoord[0]);
        copy3(inodex[1], largernodecoord[1]);
        copy3(inodex[2], largernodecoord[2]);

        // smaller tetrahedral

        copy3(inodex[3], smallernodecoord[0]);
        for (int dim = 0; dim < 3; dim++) {
          smallernodecoord[1][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
          smallernodecoord[2][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
          smallernodecoord[3][dim] = inodex[4][dim] + (inodex[3][dim] - inodex[4][dim])/ncell;
        }
      } else if (larger_side == 1) {

        // larger tetrahedral

        copy3(inodex[0], largernodecoord[0]);
        copy3(inodex[2], largernodecoord[1]);
        copy3(inodex[3], largernodecoord[2]);

        // smaller tetrahedral

        copy3(inodex[1], smallernodecoord[0]);
        for (int dim = 0; dim < 3; dim++) {
          smallernodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
          smallernodecoord[2][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
          smallernodecoord[3][dim] = inodex[4][dim] + (inodex[1][dim] - inodex[4][dim])/ncell;
        }
      }
    }
  }
  create_element(smallernodes, smaller_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(largernodes, larger_etype, ictype, 0);
  mask[nlocal] = mask[i];

  memory->destroy(smallernodes);
  memory->destroy(largernodes);
}



/*  ----------------------------------------------------------------------
    transform wedge element I to tetrahedron elements (removing left over atoms)
    scan_flag = 1 --> check for new etype
    = 0 --> split element
    split_dir = 0 --> split across diagonal with slope = -1
    = 1 --> split across the other diagonal
    -------------------------------------------------------------------------  */

void ElementVecCAC::wedge2tet(int i, int base_axis, int split_dir, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  if (element->element_shape_ids[ietype] != Element::WEDGE)
    return;
  int tet_etype;
  int type_info[8];
  if (ncells[ietype][1] != ncells[ietype][0]) {
    char errmsg[100];
    sprintf(errmsg, "Wedge element type %d: %dx%d is not square to split", ietype, ncells[ietype][0], ncells[ietype][1]);
    error->one(FLERR, errmsg);
  }
  // tetrahedron

  type_info[0] = ncells[ietype][0];
  type_info[1] = 0;
  type_info[2] = 0;
  type_info[3] = MIN(type_info[0], 4);
  type_info[4] = 0;
  type_info[5] = 0;
  type_info[6] = Element::TETRAHEDRON;
  type_info[7] = iapc;

  if (scan_flag) request_new_etype(type_info);
  else tet_etype = find_etype(type_info);

  if (scan_flag) return;

  double **tetnodecoord, ***tetnodes;
  double **inodex;
  tagint itag = tag[i];

  int ncell = ncells[ietype][0]-1;
  int *ictype = ctype[i];
  memory->create(tetnodes, 4, 3, "evec:tetnodes");
  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    tetnodecoord = tetnodes[ibasis];

    if (base_axis == 0) {
      if (split_dir == 0) {
        copy3(inodex[3], tetnodecoord[0]);
        copy3(inodex[0], tetnodecoord[1]);
        copy3(inodex[4], tetnodecoord[2]);
        copy3(inodex[5], tetnodecoord[3]);
      } else if (split_dir == 1) {
        copy3(inodex[0], tetnodecoord[0]);
        copy3(inodex[1], tetnodecoord[1]);
        copy3(inodex[2], tetnodecoord[2]);
        copy3(inodex[3], tetnodecoord[3]);
      }
    } else if (base_axis == 1) {
      if (split_dir == 0) {
        copy3(inodex[4], tetnodecoord[0]);
        copy3(inodex[1], tetnodecoord[1]);
        copy3(inodex[3], tetnodecoord[2]);
        copy3(inodex[5], tetnodecoord[3]);
      } else if (split_dir == 1) {
        copy3(inodex[1], tetnodecoord[0]);
        copy3(inodex[0], tetnodecoord[1]);
        copy3(inodex[2], tetnodecoord[2]);
        copy3(inodex[4], tetnodecoord[3]);
      }
    } else if (base_axis == 2) {
      if (split_dir == 0) {
        copy3(inodex[5], tetnodecoord[0]);
        copy3(inodex[2], tetnodecoord[1]);
        copy3(inodex[3], tetnodecoord[2]);
        copy3(inodex[4], tetnodecoord[3]);
      } else if (split_dir == 1) {
        copy3(inodex[2], tetnodecoord[0]);
        copy3(inodex[0], tetnodecoord[1]);
        copy3(inodex[1], tetnodecoord[2]);
        copy3(inodex[5], tetnodecoord[3]);
      }
    } else error->one(FLERR, "Incorrect split_axis");
  }
  create_element(tetnodes, tet_etype, ictype, itag);
  mask[nlocal] = mask[i];
  memory->destroy(tetnodes);
}


/*  ----------------------------------------------------------------------
    split hex element I to 2 tetrahedron elements and 2 pyramid elements
    scan_flag = 1 --> check for new etype
    = 0 --> split element
    split_dir = 0 --> split across diagonal with slope = -1
    = 1 --> split across the other diagonal
    larger_side = 0 --> larger side is lower half
    = 1 --> larger side is upper half
    -------------------------------------------------------------------------  */

void ElementVecCAC::hex2pyrtet(int i, int split_axis, int split_dir, int larger_side, int scan_flag) 
{
  int nlocal = element->nlocal;
  int ietype = etype[i];
  int iapc = apc[ietype];
  if (element->element_shape_ids[ietype] != Element::HEXAHEDRON)
    return;
  int large_pyr_etype, small_pyr_etype, tet_etype;
  int type_info[12];
  if (ncells[ietype][1] != ncells[ietype][0] ||
      ncells[ietype][2] != ncells[ietype][0] ||
      ncells[ietype][1] != ncells[ietype][2])
    error->one(FLERR, "Hex element must have same size in all direction to split");

  // large pyramid

  type_info[0] = ncells[ietype][0];
  type_info[1] = 0;
  type_info[2] = 0;
  type_info[3] = MIN(type_info[0], 4);
  type_info[4] = 0;
  type_info[5] = 0;
  type_info[6] = Element::PYRAMID;
  type_info[7] = iapc;

  if (scan_flag) request_new_etype(type_info);
  else large_pyr_etype = find_etype(type_info);

  // large pyramid

  type_info[0]--;
  type_info[3] = MIN(type_info[0], 4);

  if (scan_flag) request_new_etype(type_info);
  else small_pyr_etype = find_etype(type_info);

  // tetrahedron

  type_info[6] = Element::TETRAHEDRON;
  if (scan_flag) request_new_etype(type_info);
  else tet_etype = find_etype(type_info);

  if (scan_flag) return;

  double **largepyrnodecoord;
  double **smallpyrnodecoord;
  double **tet1nodecoord;
  double **tet2nodecoord;

  double **inodex;
  tagint itag = tag[i];

  int ncell = ncells[ietype][0]-1;
  int *ictype = ctype[i];

  double ***largepyrnodes;
  double ***smallpyrnodes;
  double ***tet1nodes;
  double ***tet2nodes;
  memory->create(largepyrnodes, iapc, 5, 3, "evec:largepyrnodes");
  memory->create(smallpyrnodes, iapc, 5, 3, "evec:smallpyrnodes");
  memory->create(tet1nodes, iapc, 4, 3, "evec:tet1nodes");
  memory->create(tet2nodes, iapc, 4, 3, "evec:tet2nodes");
  for (int ibasis; ibasis < iapc; ibasis++) {
    inodex = nodex[i][ibasis];
    largepyrnodecoord = largepyrnodes[ibasis];
    smallpyrnodecoord = smallpyrnodes[ibasis];
    tet1nodecoord = tet1nodes[ibasis];
    tet2nodecoord = tet2nodes[ibasis];

    if (split_dir == 0) {

      // tetrahedrons

      copy3(inodex[0], tet1nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet1nodecoord[1][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
        tet1nodecoord[2][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
        tet1nodecoord[3][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;
      }

      copy3(inodex[6], tet2nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet2nodecoord[1][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell;
        tet2nodecoord[2][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell;
        tet2nodecoord[3][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell;
      }

      if (split_axis == 0) {
        // larger pyramid base

        copy3(inodex[1], largepyrnodecoord[0]);
        copy3(inodex[2], largepyrnodecoord[1]);
        copy3(inodex[7], largepyrnodecoord[2]);
        copy3(inodex[4], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[5], largepyrnodecoord[4]);
          copy3(inodex[3], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[3][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[3][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[3][dim] - inodex[4][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[3], largepyrnodecoord[4]);
          copy3(inodex[5], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[5][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[5][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;

          }
        }
      } else if (split_axis == 1) {

        // larger pyramids base

        copy3(inodex[2], largepyrnodecoord[0]);
        copy3(inodex[5], largepyrnodecoord[1]);
        copy3(inodex[4], largepyrnodecoord[2]);
        copy3(inodex[3], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[7], largepyrnodecoord[4]);
          copy3(inodex[1], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[5][dim] + (inodex[1][dim] - inodex[5][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[4][dim] + (inodex[1][dim] - inodex[4][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[3][dim] + (inodex[1][dim] - inodex[3][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[1], largepyrnodecoord[4]);
          copy3(inodex[7], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[2][dim] + (inodex[7][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[5][dim] + (inodex[7][dim] - inodex[5][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell;

          }
        }
      } else if (split_axis == 2) {

        // larger pyramids base

        copy3(inodex[1], largepyrnodecoord[0]);
        copy3(inodex[3], largepyrnodecoord[1]);
        copy3(inodex[7], largepyrnodecoord[2]);
        copy3(inodex[5], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[4], largepyrnodecoord[4]);
          copy3(inodex[2], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[2][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[2][dim] - inodex[5][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[2], largepyrnodecoord[4]);
          copy3(inodex[4], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[4][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[4][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[4][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;
          }
        }
      }
    } else if (split_dir == 1) {

      // tetrahedrons

      copy3(inodex[1], tet1nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet1nodecoord[1][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
        tet1nodecoord[2][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
        tet1nodecoord[3][dim] = inodex[5][dim] + (inodex[1][dim] - inodex[5][dim])/ncell;
      }

      copy3(inodex[7], tet2nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet2nodecoord[1][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell;
        tet2nodecoord[2][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell;
        tet2nodecoord[3][dim] = inodex[6][dim] + (inodex[7][dim] - inodex[6][dim])/ncell;
      }

      if (split_axis == 0) {

        // larger pyramid base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[3], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[5], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[4], largepyrnodecoord[4]);
          copy3(inodex[2], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[2][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[2][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[2][dim] - inodex[5][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[2], largepyrnodecoord[4]);
          copy3(inodex[4], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[4][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[4][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;

          }
        }
      } else if (split_axis == 1) {

        // larger pyramids base

        copy3(inodex[2], largepyrnodecoord[0]);
        copy3(inodex[5], largepyrnodecoord[1]);
        copy3(inodex[4], largepyrnodecoord[2]);
        copy3(inodex[3], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[0], largepyrnodecoord[4]);
          copy3(inodex[6], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[4][dim] + (inodex[6][dim] - inodex[4][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[3][dim] + (inodex[6][dim] - inodex[3][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[6], largepyrnodecoord[4]);
          copy3(inodex[0], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[2][dim] + (inodex[0][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[5][dim] + (inodex[0][dim] - inodex[5][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;

          }
        }
      } else if (split_axis == 2) {

        // larger pyramids base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[2], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[4], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[3], largepyrnodecoord[4]);
          copy3(inodex[5], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[5][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[5][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[5][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[5], largepyrnodecoord[4]);
          copy3(inodex[3], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[3][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[3][dim] - inodex[4][dim])/ncell;
          }
        }
      }
    } else if (split_dir == 2) {

      // tetrahedrons

      copy3(inodex[2], tet1nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet1nodecoord[1][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;
        tet1nodecoord[2][dim] = inodex[3][dim] + (inodex[2][dim] - inodex[3][dim])/ncell;
        tet1nodecoord[3][dim] = inodex[6][dim] + (inodex[2][dim] - inodex[6][dim])/ncell;
      }

      copy3(inodex[4], tet2nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet2nodecoord[1][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell;
        tet2nodecoord[2][dim] = inodex[5][dim] + (inodex[4][dim] - inodex[5][dim])/ncell;
        tet2nodecoord[3][dim] = inodex[7][dim] + (inodex[4][dim] - inodex[7][dim])/ncell;
      }

      if (split_axis == 0) {
        // larger pyramid base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[3], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[5], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[1], largepyrnodecoord[4]);
          copy3(inodex[7], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[7][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[7][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[7][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[7][dim] - inodex[5][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[7], largepyrnodecoord[4]);
          copy3(inodex[1], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[1][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[1][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[1][dim] - inodex[5][dim])/ncell;

          }
        }
      } else if (split_axis == 1) {

        // larger pyramids base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[7], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[1], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[3], largepyrnodecoord[4]);
          copy3(inodex[5], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[5][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[7][dim] + (inodex[5][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[5][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[5], largepyrnodecoord[4]);
          copy3(inodex[3], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[7][dim] + (inodex[3][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[3][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[1][dim] + (inodex[3][dim] - inodex[1][dim])/ncell;

          }
        }
      } else if (split_axis == 2) {

        // larger pyramids base

        copy3(inodex[1], largepyrnodecoord[0]);
        copy3(inodex[3], largepyrnodecoord[1]);
        copy3(inodex[7], largepyrnodecoord[2]);
        copy3(inodex[5], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[0], largepyrnodecoord[4]);
          copy3(inodex[6], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[6][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[6][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[6][dim] - inodex[5][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[6], largepyrnodecoord[4]);
          copy3(inodex[0], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[3][dim] + (inodex[0][dim] - inodex[3][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[0][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[5][dim] + (inodex[0][dim] - inodex[5][dim])/ncell;
          }
        }
      }
    } else if (split_dir == 3) {

      // tetrahedrons

      copy3(inodex[3], tet1nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet1nodecoord[1][dim] = inodex[0][dim] + (inodex[3][dim] - inodex[0][dim])/ncell;
        tet1nodecoord[2][dim] = inodex[2][dim] + (inodex[3][dim] - inodex[2][dim])/ncell;
        tet1nodecoord[3][dim] = inodex[7][dim] + (inodex[3][dim] - inodex[7][dim])/ncell;
      }

      copy3(inodex[5], tet2nodecoord[0]);
      for (int dim = 0; dim < 3; dim++) {
        tet2nodecoord[1][dim] = inodex[1][dim] + (inodex[5][dim] - inodex[1][dim])/ncell;
        tet2nodecoord[2][dim] = inodex[4][dim] + (inodex[5][dim] - inodex[4][dim])/ncell;
        tet2nodecoord[3][dim] = inodex[6][dim] + (inodex[5][dim] - inodex[6][dim])/ncell;
      }

      if (split_axis == 0) {
        // larger pyramid base

        copy3(inodex[1], largepyrnodecoord[0]);
        copy3(inodex[2], largepyrnodecoord[1]);
        copy3(inodex[7], largepyrnodecoord[2]);
        copy3(inodex[4], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[0], largepyrnodecoord[4]);
          copy3(inodex[6], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[6][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[6][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[6][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[6][dim] - inodex[4][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[6], largepyrnodecoord[4]);
          copy3(inodex[0], smallpyrnodecoord[4]);

          // smaller pyramid base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[1][dim] + (inodex[0][dim] - inodex[1][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[0][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[7][dim] + (inodex[0][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[0][dim] - inodex[4][dim])/ncell;

          }
        }
      } else if (split_axis == 1) {

        // larger pyramids base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[7], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[1], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[4], largepyrnodecoord[4]);
          copy3(inodex[2], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[2][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[7][dim] + (inodex[2][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[2][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[1][dim] + (inodex[2][dim] - inodex[1][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[2], largepyrnodecoord[4]);
          copy3(inodex[4], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[4][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[7][dim] + (inodex[4][dim] - inodex[7][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[4][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[1][dim] + (inodex[4][dim] - inodex[1][dim])/ncell;

          }
        }
      } else if (split_axis == 2) {

        // larger pyramids base

        copy3(inodex[0], largepyrnodecoord[0]);
        copy3(inodex[2], largepyrnodecoord[1]);
        copy3(inodex[6], largepyrnodecoord[2]);
        copy3(inodex[4], largepyrnodecoord[3]);

        if (larger_side == 0) {

          // pyramids peaks

          copy3(inodex[1], largepyrnodecoord[4]);
          copy3(inodex[7], smallpyrnodecoord[4]);

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[7][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[7][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[7][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[7][dim] - inodex[4][dim])/ncell;

          }
        } else if (larger_side == 1) {

          // pyramids peaks

          copy3(inodex[7], largepyrnodecoord[4]);
          copy3(inodex[1], smallpyrnodecoord[4]);

          // smaller pyramids base

          for (int dim = 0; dim < 3; dim++) {
            smallpyrnodecoord[0][dim] = inodex[0][dim] + (inodex[1][dim] - inodex[0][dim])/ncell;
            smallpyrnodecoord[1][dim] = inodex[2][dim] + (inodex[1][dim] - inodex[2][dim])/ncell;
            smallpyrnodecoord[2][dim] = inodex[6][dim] + (inodex[1][dim] - inodex[6][dim])/ncell;
            smallpyrnodecoord[3][dim] = inodex[4][dim] + (inodex[1][dim] - inodex[4][dim])/ncell;
          }
        }
      }
    }
  }
  create_element(smallpyrnodes, small_pyr_etype, ictype, itag);
  mask[nlocal] = mask[i];
  create_element(largepyrnodes, large_pyr_etype, ictype, 0);
  mask[nlocal+1] = mask[i];
  create_element(tet1nodes, tet_etype, ictype, 0);
  mask[nlocal+2] = mask[i];
  create_element(tet2nodes, tet_etype, ictype, 0);
  mask[nlocal+3] = mask[i];

  memory->destroy(smallpyrnodes);
  memory->destroy(largepyrnodes);
  memory->destroy(tet1nodes);
  memory->destroy(tet2nodes);
}
