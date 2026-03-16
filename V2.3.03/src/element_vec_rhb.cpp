#include <stdlib.h>
#include <string.h>
#include "element_vec_rhb.h"
#include "element.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "memory.h"
#include "group.h"
#include "error.h"
#include "math_extra.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathExtra;
#define EPSILON 1e-6
#define BIG     1e30
/*----------------------------------------------------------------*/

ElementVecRHB::ElementVecRHB(CAC *cac) : ElementVec(cac)
{
  if (domain->dimension == 3) npe = 8;
  else npe = 4;
  element->npe = npe;

  linkcut = 0.0;
  apc = element->apc;
  naae = niae = NULL;
  intg = NULL;

  size_data_element = 6;
  size_data_node = 6;
  size_dat_node = 3;
  size_dat_node_connect = npe;
  size_data_vel = 5;
  size_velocity = 3*npe;
  size_forward = 3*(npe+1);
  size_reverse = 3*npe;
  size_border = 7 + 5*npe;
  xcol_data = 4;

}

/*----------------------------------------------------------------*/

ElementVecRHB::~ElementVecRHB()
{
  memory->destroy(naae);
  memory->destroy(niae);
  memory->destroy(intg);
}

/*----------------------------------------------------------------*/

void ElementVecRHB::check_type()
{
  // check if every element has defined interpolate and integration
 
  for (int i = 0; i < element->nlocal; i++) {
    if (!interpolate_setflag[etype[i]]) 
      error->all(FLERR,"Element interpolated atoms have not been set for all elements");
    if (!integration_setflag[etype[i]]) 
      error->all(FLERR,"Element integration points have not been set for all elements");
  }
}

/* ----------------------------------------------------------------------
   copy element I info to element J
   ------------------------------------------------------------------------- */

void ElementVecRHB::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  ctype[j] = ctype[i];
  etype[j] = etype[i];
  mask[j] = mask[i];
  image[j] = image[i];
  initial_size[j] = initial_size[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  for (int k = 0; k < npe; k++) {
    nodemask[j][k] = nodemask[i][k];
    nodetag[j][k] = nodetag[i][k];
    nodex[j][k][0] = nodex[i][k][0];
    nodex[j][k][1] = nodex[i][k][1];
    nodex[j][k][2] = nodex[i][k][2];
    nodev[j][k][0] = nodev[i][k][0];
    nodev[j][k][1] = nodev[i][k][1];
    nodev[j][k][2] = nodev[i][k][2];
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
   ------------------------------------------------------------------------- */

void ElementVecRHB::grow(int n)
{
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
  nodemask = memory->grow(element->nodemask,nmax,npe,"element:nodemask");
  image = memory->grow(element->image,nmax,"element:image");
  x = memory->grow(element->x,nmax,3,"element:x");
  slip_plane = memory->grow(element->slip_plane,nmax,3,3,"element:slip_plane");
  nodex = memory->grow(element->nodex,nmax,npe,3,"element:nodex");
  nodev = memory->grow(element->nodev,nmax,npe,3,"element:nodev");
  nodef = memory->grow(element->nodef,nmax,npe,3,"element:nodef");
  nodetag = memory->grow(element->nodetag,nmax,npe,"element:nodetag");
  elem_size = memory->grow(element->elem_size,nmax,"element:elem_size");
  initial_size = memory->grow(element->initial_size,nmax,"element:initial_size");
}

/* ---------------------------------------------------------
   grow etype arrays to n types
 * ------------------------------------------------------*/

void ElementVecRHB::grow_etype_arrays(int n)
{
  int old_netypes;

  if (n <= element->netypes) return;
  else {
    old_netypes = element->netypes;
    element->netypes = n;
  }

  // grow style specific arrays

  // grow integration point arrays

  memory->grow(niae,n+1,3,"evec:niae");
  memory->grow(intg,n+1,element->maxintg,3,"evec:intg");
  memory->grow(integration_setflag,n+1,"evec:integration_setflag");

  // interpolated atom arrays

  memory->grow(naae,n+1,3,"element:naae");
  memory->grow(interpolate_setflag,n+1,"element:interpolate_setflag");

  // grow element arrays

  // grow arrays common to all element style

  // grow integration point arrays

  nintg = memory->grow(element->nintg,n+1,"element:nintg");
  i2ia = memory->grow(element->i2ia,n+1,element->maxintg,"element:i2ia");
  weight = memory->grow(element->weight,n+1,element->maxintg,"element:weight");
  i2n = memory->grow(element->i2n,n+1,element->maxintg,"element:i2n");

  // grow node arrays

  n2i = memory->grow(element->n2i,n+1,npe,"element:n2i");

  // interpolated atom arrays

  nintpl = memory->grow(element->nintpl,n+1,"element:nintpl");
  shape_array = memory->grow(element->shape_array,n+1,element->maxintpl,npe,"element:shape_array");
  ia2i = memory->grow(element->ia2i,n+1,element->maxintpl,"element:ia2i");

  // sub-element array

  natom_subelem = memory->grow(element->natom_subelem,n+1,nsubelem,"element:natom_subelem");
  ias2ia = memory->grow(element->ias2ia,n+1,nsubelem,element->maxsubintpl,"element:ias2ia");

  for (int itype = old_netypes + 1; itype <= n; itype++) {
    interpolate_setflag[itype] = 0;
    integration_setflag[itype] = 0;
    nintg[itype] = 0; 
    nintpl[itype] = 0; 
    for (int i = 0; i < element->maxintpl; i++) 
      ia2i[itype][i] = -1;
    for (int i = 0; i < element->maxintg; i++) {
      i2ia[itype][i] = -1;
      i2n[itype][i] = -1;
    }
    for (int i = 0; i < npe; i++)
      n2i[itype][i] = -1;
  }
}

/* ----------------------------------------------------------------------
   unpack one line from Elements section of data file
   initialize other element quantities
   ------------------------------------------------------------------------- */

void ElementVecRHB::data_element(double *coord, imageint imagetmp, char **values)
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
  if (ctype[nlocal] <= 0 || ctype[nlocal] > element->nctypes)
    error->one(FLERR,"Invalid atom type in Elements section of data file"); 
  image[nlocal] = imagetmp;

  // 1 equal to the all group, 2 equal to the element group

  mask[nlocal] = 1|2;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  initial_size[nlocal] = 0.0;
  element->nintpls += nintpl[etype[nlocal]];
  element->nlocal++;

}

/* ----------------------------------------------------------------------
   unpack one line from Node Velocities section of data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::data_vel(int m, char **values)
{

  int i = atoi(values[0])-1; 

  if (i >= npe) error->all(FLERR,"Invalid node index value in Node section of data file");
  nodev[m][i][0] = atof(values[1]);
  nodev[m][i][1] = atof(values[2]);
  nodev[m][i][2] = atof(values[3]);
}

/* ----------------------------------------------------------------------
   unpack one line from Node of data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::data_node(int m, char **values)
{
  int i = atoi(values[0])-1; 
  if (i >= npe) error->all(FLERR,"Invalid node index value in Node section of data file");
  nodex[m][i][0] = atof(values[1]);
  nodex[m][i][1] = atof(values[2]);
  nodex[m][i][2] = atof(values[3]);
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

int ElementVecRHB::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = initial_size[i];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(ctype[i]).d;
  buf[m++] = ubuf(etype[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  for (int j = 0; j < npe; j++){
    buf[m++] = ubuf(nodemask[i][j]).d;
    buf[m++] = ubuf(nodetag[i][j]).d;
    buf[m++] = nodex[i][j][0];
    buf[m++] = nodex[i][j][1];
    buf[m++] = nodex[i][j][2];
    buf[m++] = nodev[i][j][0];
    buf[m++] = nodev[i][j][1];
    buf[m++] = nodev[i][j][2];
  }
  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------
   unpack element from exchanged proc
   -------------------------------------------------------------- */

int ElementVecRHB::unpack_exchange(double *buf)
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
  initial_size[nlocal] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  ctype[nlocal] = (int) ubuf(buf[m++]).i;
  etype[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  for (int j = 0; j < npe; j++) {
    nodemask[nlocal][j] = (int) ubuf(buf[m++]).i;
    nodetag[nlocal][j] = (tagint) ubuf(buf[m++]).i;
    nodex[nlocal][j][0] = buf[m++];
    nodex[nlocal][j][1] = buf[m++];
    nodex[nlocal][j][2] = buf[m++];
    nodev[nlocal][j][0] = buf[m++];
    nodev[nlocal][j][1] = buf[m++];
    nodev[nlocal][j][2] = buf[m++];
  }

  // update nlocal;

  element->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   pack data for elements in list for border exchange
   ------------------------------------------------------------------------- */

int ElementVecRHB::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k,m;
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
      for (k = 0; k < npe; k++) {
        buf[m++] = ubuf(nodemask[j][k]).d;
        buf[m++] = ubuf(nodetag[j][k]).d;
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
      for (k = 0; k < npe; k++) {
        buf[m++] = ubuf(nodemask[j][k]).d;
        buf[m++] = ubuf(nodetag[j][k]).d;
        buf[m++] = nodex[j][k][0]+dx;
        buf[m++] = nodex[j][k][1]+dy;
        buf[m++] = nodex[j][k][2]+dz;
      }
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for elements received from border exchange
   ------------------------------------------------------------------------- */


void ElementVecRHB::unpack_border(int n, int first, double *buf)
{
  int m,last;
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
    for (int j = 0; j < npe; j++) {
      nodemask[i][j] = (int) ubuf(buf[m++]).i;
      nodetag[i][j] = (tagint) ubuf(buf[m++]).i;
      nodex[i][j][0] = buf[m++];
      nodex[i][j][1] = buf[m++];
      nodex[i][j][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int ElementVecRHB::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k;
  double dx,dy,dz;

  int m = 0;

  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      for (k = 0; k < npe; k++) {
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
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      for (k = 0; k < npe; k++) {
        buf[m++] = nodex[j][k][0] + dx;
        buf[m++] = nodex[j][k][1] + dy;
        buf[m++] = nodex[j][k][2] + dz;
      }
    }
  }
  return m;
}


/* ---------------------------------------------------------------------- */

void ElementVecRHB::unpack_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    for (j = 0; j < npe; j++) {
      nodex[i][j][0] = buf[m++];
      nodex[i][j][1] = buf[m++];
      nodex[i][j][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int ElementVecRHB::pack_reverse(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
    for (j = 0; j < npe; j++) {
      buf[m++] = nodef[i][j][0];
      buf[m++] = nodef[i][j][1];
      buf[m++] = nodef[i][j][2];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void ElementVecRHB::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,k,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < npe; k++) {
      nodef[j][k][0] += buf[m++];
      nodef[j][k][1] += buf[m++];
      nodef[j][k][2] += buf[m++];
    }
  }
}


/* ----------------------------------------------------------------------
   add a new element
   ------------------------------------------------------------------------- */

void ElementVecRHB::create_element(double *coord, double **nodecoord, int ietype, int ictype, int itag)
{

  int nlocal = element->nlocal;
  if (nlocal == nmax) grow(0);
  while ((nlocal+1) * element->max_nintg >= nmaxintg) grow_nmaxintg();
  while ((nlocal+1) * element->max_nintpl >= nmaxintpl) grow_nmaxintpl();

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  for (int i = 0; i < npe; i++) {
    nodex[nlocal][i][0] = nodecoord[i][0];
    nodex[nlocal][i][1] = nodecoord[i][1];
    nodex[nlocal][i][2] = nodecoord[i][2];
    nodev[nlocal][i][0] = 0.0;
    nodev[nlocal][i][1] = 0.0;
    nodev[nlocal][i][2] = 0.0;
  }

  etype[nlocal] = ietype;
  ctype[nlocal] = ictype;
  tag[nlocal] = itag;
  mask[nlocal] = 1|2; 
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  element->nlocal++;
}

/* ----------------------------------------------------------------------
   add a new element, node coords are calculated from lattice
   ------------------------------------------------------------------------- */

void ElementVecRHB::create_element(double *coord, int ietype, int ictype, int itag)
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

  domain->lattice->box2lattice(xtmp,ytmp,ztmp);

  double *a1 = domain->lattice->a1;
  double *a2 = domain->lattice->a2;
  double *a3 = domain->lattice->a3;

  int c1,c2,c3;

  c1 = naae[ietype][0]-1; 
  c2 = naae[ietype][1]-1; 
  c3 = naae[ietype][2]-1;

  double xnode,ynode,znode;

  if (domain->dimension == 3) {

    double d[8][3] = {{-1,-1,-1},
      { 1,-1,-1},
      { 1, 1,-1},
      {-1, 1,-1},
      {-1,-1, 1},
      { 1,-1, 1},
      { 1, 1, 1},
      {-1, 1, 1}};

    for (int i = 0; i < npe; i++) {
      xnode = xtmp + c1*d[i][0]/2.0;
      ynode = ytmp + c2*d[i][1]/2.0;
      znode = ztmp + c3*d[i][2]/2.0;

      domain->lattice->lattice2box(xnode,ynode,znode);

      nodex[nlocal][i][0] = xnode;
      nodex[nlocal][i][1] = ynode;
      nodex[nlocal][i][2] = znode;

      nodev[nlocal][i][0] = 0.0;
      nodev[nlocal][i][1] = 0.0;
      nodev[nlocal][i][2] = 0.0;
      nodetag[nlocal][i] = 0;
      nodemask[nlocal][i] = 1|2;
    }
  } else {
    double d[4][2] = {{-1,-1},
      { 1,-1},
      { 1, 1},
      {-1, 1}};

    for (int i = 0; i < npe; i++) {
      xnode = xtmp + c1*d[i][0]/2.0;
      ynode = ytmp + c2*d[i][1]/2.0;
      znode = ztmp;

      domain->lattice->lattice2box(xnode,ynode,znode);

      nodex[nlocal][i][0] = xnode;
      nodex[nlocal][i][1] = ynode;
      nodex[nlocal][i][2] = znode;

      nodev[nlocal][i][0] = 0.0;
      nodev[nlocal][i][1] = 0.0;
      nodev[nlocal][i][2] = 0.0;
      nodetag[nlocal][i] = 0;
      nodemask[nlocal][i] = 1|2;
    }
  }
  etype[nlocal] = ietype;
  ctype[nlocal] = ictype;
  initial_size[nlocal] = 0.0;
  tag[nlocal] = itag;
  mask[nlocal] = 1|2; 
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  element->nlocal++;
}


/*----------------------------------------------------
  write element info into Elements section in data file
  ------------------------------------------*/

void ElementVecRHB::write_element_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %d %-1.16e %-1.16e %-1.16e \n",
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, (int) ubuf(buf[i][2]).i,buf[i][3],buf[i][4],buf[i][5]);
}

/*----------------------------------------------------
  pack element info for Elements section in data file
  ------------------------------------------*/

void ElementVecRHB::pack_element_data(double **buf)
{
  for (int i = 0; i < element->nlocal; i++){
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(etype[i]).d;
    buf[i][2] = ubuf(ctype[i]).d;
    buf[i][3] = x[i][0];
    buf[i][4] = x[i][1];
    buf[i][5] = x[i][2];
  }
}

/*----------------------------------------------------
  write node info into Nodes section in data file
  ------------------------------------------*/

void ElementVecRHB::write_node_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e "TAGINT_FORMAT" \n",
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4],(tagint) ubuf(buf[i][5]).i);
}

/*----------------------------------------------------
  write node info into Nodes section in dat file
  ------------------------------------------*/

void ElementVecRHB::write_node_dat(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",
        buf[i][0],buf[i][1],buf[i][2]);
}

/*----------------------------------------------------
  write node connectivity info into Nodes section in dat file
  ------------------------------------------*/

void ElementVecRHB::write_node_connect_dat(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < npe; j++)
      fprintf(fp,TAGINT_FORMAT" ",(tagint) ubuf(buf[i][j]).i);
    fprintf(fp,"\n");
  }
}


/*----------------------------------------------------
  pack node info for Nodes section in data file
  ------------------------------------------*/

void ElementVecRHB::pack_node_data(double **buf)
{
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    for (int j = 0; j < npe; j++) {
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

/*----------------------------------------------------
  pack node info for Nodes section in dat file
  ------------------------------------------*/

void ElementVecRHB::pack_node_dat(double **buf)
{
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    for (int j = 0; j < npe; j++) {
      buf[k][0] = nodex[i][j][0];
      buf[k][1] = nodex[i][j][1];
      buf[k][2] = nodex[i][j][2];
      k++;
    }
  }
}

/*----------------------------------------------------
  pack node connectivity info for Nodes section in dat file
  ------------------------------------------*/

void ElementVecRHB::pack_node_connect_dat(double **buf)
{
  for (int i = 0; i < element->nlocal; i++) {
    for (int j = 0; j < npe; j++) {
      buf[i][j] = ubuf(nodetag[i][j]).d;
    }
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::pack_vel(double **buf)
{
  int k = 0;
  for (int i = 0; i < element->nlocal; i++) {
    for (int j = 0; j < npe; j++) {
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
   write velocity info to data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n",
        (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,buf[i][2],buf[i][3],buf[i][4]);
}

/* ---------------------------------------------
   Setup interpolated atoms
   -------------------------------*/

void ElementVecRHB::set_interpolate(const char *str, int type_offset)
{
  int itype;
  int nx,ny,nz;
  int n = sscanf(str,"%d %d %d %d", &itype,&nx,&ny,&nz);
  if (n != 4) error->all(FLERR, "Invalid interpolate line in data file");
  if (interpolate_setflag[itype]) error->all(FLERR,"Interpolated atoms for this element type has been set");
  interpolate_setflag[itype] = 1;

  itype += type_offset;

  if (itype < 1 || itype > element->netypes)
    error->all(FLERR,"Invalid element type for interpolate set");
  naae[itype][0] = nx;
  naae[itype][1] = ny;
  naae[itype][2] = nz;

  if (nx < 2 || ny < 2 || nz < 2) error->all(FLERR,"Invalid interpolate values");

  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Number of atoms along third direction must be one for 2D"); 

  n = nx*ny*nz;
  if (n > element->maxintpl) error->all(FLERR,"Too many interpolated atoms per element, boost maxintpl");
  element->max_nintpl = MAX(element->max_nintpl,n);
  nintpl[itype] = n;
  setup_interpolate_element(itype);
  if (subelemflag == 0) {
    setup_sub_element();
    subelemflag = 1;
  }
  setup_sub_element(itype);
}

/* ------------------------------------------
   Setup integration points from data file
 * ----------------------------------*/

void ElementVecRHB::set_integration(const char *str, int type_offset)
{
  int itype;
  int nx,ny,nz;
  int n = sscanf(str,"%d %d %d %d", &itype,&nx,&ny,&nz);
  if (n != 4) error->all(FLERR, "Invalid integration line in data file");
  if (integration_setflag[itype]) 
    error->all(FLERR,"Integration points for this element type has been set");
  integration_setflag[itype] = 1;

  itype += type_offset;
  if (!interpolate_setflag[itype]) {
    char *erstr = new char[100];
    sprintf(erstr,"Number of interpolated atoms for element type %d has not been defined yet",itype);
    error->all(FLERR,erstr);
    delete [] erstr;
  }
  if (itype < 1 || itype > element->netypes)
    error->all(FLERR,"Invalid element type for integration");
  niae[itype][0] = nx;
  niae[itype][1] = ny;
  niae[itype][2] = nz;

  if (nx < 2 || ny < 2 || nz < 2) error->all(FLERR, "Invalid integration values");

  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Number of integration point along third direction must be 1 for 2D"); 
  n = nx * ny * nz;
  if (n > element->maxintg) error->all(FLERR,"Too many integration points per element, boost maxintg");
  element->max_nintg = MAX(element->max_nintg,n);
  nintg[itype] = n;

  setup_integration_point(itype);
}

/* ------------------------------------------
   Read number of integration points from Integration List 
   section in data file
 * ----------------------------------*/

void ElementVecRHB::set_nintg_user(char *str, int type_offset)
{
  int itype;
  int n,num;
  n = sscanf(str,"%d %d", &itype,&num);
  if (n != 2) error->all(FLERR, "Invalid Integration points line in data file");
  itype += type_offset;
  if (itype < 1 || itype > element->netypes)
    error->all(FLERR,"Invalid type for Integration points set");
  if (!interpolate_setflag[itype]) {
    char *erstr = new char[100];
    sprintf(erstr,"Number of interpolated atoms for element type %d has not been defined yet",itype);
    error->all(FLERR,erstr);
    delete [] erstr;
  }
  if (integration_setflag[itype]) 
    error->all(FLERR,"Integration points for this element type has already been set");

  nintg[itype] = num;
  if (nintg[itype] <= 0 || nintg[itype] > element->maxintg)
    error->all(FLERR,"Invalid number of integration points value in Integration List section");

}

/* ------------------------------------------
   Read integration points from Integration List 
   section in data file
 * ----------------------------------*/

void ElementVecRHB::set_intg_user(char *str, int type_offset)
{
  int itype,iintg,i;
  int nx,ny,nz;
  int niax,niay,niaz;
  char **values = new char*[6];

  values[0] = strtok(str," \t\n\r\f");

  for (i = 1; i < 6; i++) 
    values[i] = strtok(NULL," \t\n\r\f");

  if (values[0] == NULL) 
    error->all(FLERR,"Invalid integration point line in data file");
  itype = atoi(values[0]);
  if (!interpolate_setflag[itype]) {
    char *errstr = new char[100];
    sprintf(errstr,"Number of interpolated atoms for element type %d has not been defined yet",itype);
    error->all(FLERR,errstr);
    delete [] errstr;
  }

  iintg = atoi(values[1])-1;
  nx = atoi(values[2]);
  ny = atoi(values[3]);
  nz = atoi(values[4]);

  niax = naae[itype][0];
  niay = naae[itype][1];
  niaz = naae[itype][2];
  if (nx >= niax || ny >= niay || nz >= niaz ||
      nx < 0 || ny < 0 || nz < 0) 
    error->all(FLERR,"Invalid Integration points line in data file");

  intg[itype][iintg][0] = nx;
  intg[itype][iintg][1] = ny;
  intg[itype][iintg][2] = nz;
  weight[itype][iintg] = atof(values[5]); 

  // mapping between integration point index and interpolated atom index

  int iintpl = nx*niay*niaz + ny*niaz + nz;
  if (i2ia[itype][iintg] < 0) {
    i2ia[itype][iintg] = iintpl;
    ia2i[itype][iintpl] = iintg;
  } else error->all(FLERR,"Integration point to interpolated atom mapping is not correct");

  // mapping between node index and integration point index

  for (i = 0; i < npe; i++)
    if (fabs(shape_array[itype][iintpl][i]-1.0) < EPSILON) {
      if (n2i[itype][i] < 0) {
        n2i[itype][i] = iintg;
        i2n[itype][iintg] = i;
      } else error->all(FLERR,"Node to integration point mapping is not correct");
    }

  delete [] values;

} 

/* -----------------------------------------------------------------------------------
   Calculate integration point related arrays:
   - Weight of integration points
   - i2ia to convert integration point index 
   to interpolated atom index in an element,
   used for interpolating integration points
   - n2i to convert node index to integration point index
   ------------------------------------------------------------------------------------*/

void ElementVecRHB::setup_integration_point(int itype)
{
  double dx,dy,dz;
  double px,py,pz;
  int x,y,z,w;
  double w_edge_x,w_face_xy;
  double w_edge_y,w_face_yz;
  double w_edge_z,w_face_xz;
  double w_inner;
  double total_weight;
  int niax,niay,niaz;
  int nix,niy,niz;
  int nx,ny,nz;
  int i,j,k;
  int inode,iintg,iintpl;

  double **xgauss = element->xgauss;
  double **wgauss = element->wgauss;

  // number of interpolated atoms along each edge

  niax = naae[itype][0];
  niay = naae[itype][1];
  niaz = naae[itype][2];

  // number of integration points along each edge

  nix = niae[itype][0];
  niy = niae[itype][1];
  niz = niae[itype][2];

  // calculate weight for integration points

  w_edge_x = static_cast<double> (niax - 2)/2.0;
  w_edge_y = static_cast<double> (niay - 2)/2.0;
  if (domain->dimension == 3) {
    w_edge_z = static_cast<double> (niaz - 2)/2.0;
    w_face_xy = w_edge_x * w_edge_y;
    w_face_yz = w_edge_y * w_edge_z;
    w_face_xz = w_edge_x * w_edge_z;
    w_inner = w_edge_x * w_edge_y * w_edge_z;
  } else w_inner = w_edge_x * w_edge_y;

  // loop over integration points
  // nx,ny,nz are indices of interpolated atoms
  // that are integration points (rounded to nearest whole integer)

  dx = (niax-1) / static_cast<double>(nix-1);
  dy = (niay-1) / static_cast<double>(niy-1);
  if (domain->dimension == 3)
    dz = (niaz-1.0) / static_cast<double>(niz-1.0);
  iintg = 0; 
  for (i = 0; i < nix; i++){
    if (i == 0) nx = 0;
    else if (i == nix-1) nx = niax-1;
    else nx = (int) ((1.0 + xgauss[nix-2][i])/2.0*(niax-2)+1.0);
    for (j = 0; j < niy; j++){
      if (j == 0) ny = 0;
      else if (j == niy-1) ny = niay-1;
      else ny = (int) ((1.0 + xgauss[niy-2][j])/2.0*(niay-2)+1.0);
      if (domain->dimension == 3) {
        for (k = 0; k < niz; k++){
          if (k == 0) nz = 0;
          else if (k == niz-1) nz = niaz-1;
          else nz = (int) ((1.0 + xgauss[niz-2][k])/2.0*(niaz-2)+1.0);

          // calculate weight of integration points

          x = (i == 0 || i == (nix - 1));
          y = (j == 0 || j == (niy - 1));
          z = (k == 0 || k == (niz - 1));
          w = x + y + z;

          // integration point is a node

          if (w == 3) weight[itype][iintg] = 1.0;

          // integration point is on an edge

          else if (w == 2) {
            // integration point is on edge in x direction
            if (!x) weight[itype][iintg] = w_edge_x*wgauss[nix-2][i];
            // integration point is on edge in y direction
            else if (!y) weight[itype][iintg] = w_edge_y*wgauss[niy-2][j];
            // integration point is on edge in z direction
            else if (!z) weight[itype][iintg] = w_edge_z*wgauss[niy-2][k];

            // integration point is on a face

          } else if (w == 1) {
            // integration point is on yz face 
            if (x) weight[itype][iintg] = w_face_yz*wgauss[niy-2][j]*wgauss[niz-2][k];
            // integration point is on xz face 
            else if (y) weight[itype][iintg] = w_face_xz*wgauss[nix-2][i]*wgauss[niz-2][k];
            // integration point is on xy face 
            else if (z) weight[itype][iintg] = w_face_xy*wgauss[nix-2][i]*wgauss[niy-2][j];

            // integration point is inside element

          } else weight[itype][iintg] = w_inner*wgauss[nix-2][i]*wgauss[niy-2][j]*wgauss[niz-2][k];

          intg[itype][iintg][0] = nx;
          intg[itype][iintg][1] = ny;
          intg[itype][iintg][2] = nz;

          // list to convert integration point index to 
          // interpolated atom index
          iintpl = nx*niay*niaz + ny*niaz + nz;
          i2ia[itype][iintg] = iintpl;
          ia2i[itype][iintpl] = iintg;

          // list to convert node index to integration point index

          for (inode = 0; inode < npe; inode++)
            if (fabs(shape_array[itype][i2ia[itype][iintg]][inode]-1.0) < EPSILON) {
              if (n2i[itype][inode] < 0) {
                n2i[itype][inode] = iintg;
                i2n[itype][iintg] = inode;
              } else error->all(FLERR,"Node to integration point mapping is not correct");
            }
          iintg++;
        }
      } else {

        x = (i == 0 || i == (nix - 1));
        y = (j == 0 || j == (niy - 1));
        w = x + y;

        // integration point is a node

        if (w == 2) weight[itype][iintg] = 1.0;

        // integration point is on an edge

        else if (w == 1) {
          // integration point is on edge in x direction
          if (!x) weight[itype][iintg] = w_edge_x;
          // integration point is on edge in y direction
          else if (!y) weight[itype][iintg] = w_edge_y;

          // integration point is inside element

        } else weight[itype][iintg] = w_inner;

        intg[itype][iintg][0] = nx;
        intg[itype][iintg][1] = ny;
        intg[itype][iintg][2] = 0;

        // list to convert integration point index to 
        // interpolated atom index

        iintpl = nx*niay + ny;
        i2ia[itype][iintg] = iintpl;
        ia2i[itype][iintpl] = iintg;

        // list to convert node index to integration point index

        for (inode = 0; inode < npe; inode++)
          if (fabs(shape_array[itype][i2ia[itype][iintg]][inode]-1.0) < EPSILON) {
            if (n2i[itype][inode] < 0) {
              n2i[itype][inode] = iintg;
              i2n[itype][iintg] = inode;
            } else error->all(FLERR,"Node to integration point mapping is not correct");
          }
        iintg++;
      }
    }
  }

  // check if all node to integration point mapping is set

  for (inode = 0; inode < npe; inode++) {
    if (n2i[itype][inode] < 0) 
      error->all(FLERR,"Node to integration point mapping is not complete");
    if (i2n[itype][n2i[itype][inode]] < 0)
      error->all(FLERR,"Integration point to node mapping is not complete");
  }
  // check if weight is assigned correctly

  total_weight = 0.0;
  for (int l = 0; l < nintg[itype]; l++) total_weight += weight[itype][l];
  if (fabs(total_weight - nintpl[itype]) > EPSILON) error->all(FLERR,"Integraion weight is not assigned correctly");
}

/*-----------------------------------------------------------------------------------
  Calculate interpolated atom related arrays for etype itype
  - Shape function array to interpolate atoms
  - Use linear interpolation function
  x direction: node 1 -> 2
  y direction: node 1 -> 4 
  z direction: node 1 -> 5 
  ------------------------------------------------------------------------------------*/

void ElementVecRHB::setup_interpolate_element(int itype)
{
  int i,j,k,n;
  double px,py,pz;
  double dx_pc,dy_pc,dz_pc;
  int nx,ny,nz;

  nx = naae[itype][0];
  ny = naae[itype][1];
  nz = naae[itype][2];
  dx_pc = 2.0/(nx-1);
  dy_pc = 2.0/(ny-1);
  if (domain->dimension == 3) dz_pc = 2.0/(nz-1);

  // calculate shape function arrays

  n = 0;

  for (i = 0; i < nx; i++) {
    px = dx_pc*i - 1.0;
    for (j = 0; j < ny; j++) {
      py = dy_pc*j - 1.0;
      if (domain->dimension == 3) {
        for (k = 0; k < nz; k++) {
          pz = dz_pc*k - 1.0;
          shape_array[itype][n][0] = 0.125*(1-px)*(1-py)*(1-pz);
          shape_array[itype][n][1] = 0.125*(1+px)*(1-py)*(1-pz);
          shape_array[itype][n][2] = 0.125*(1+px)*(1+py)*(1-pz);
          shape_array[itype][n][3] = 0.125*(1-px)*(1+py)*(1-pz);
          shape_array[itype][n][4] = 0.125*(1-px)*(1-py)*(1+pz);
          shape_array[itype][n][5] = 0.125*(1+px)*(1-py)*(1+pz);
          shape_array[itype][n][6] = 0.125*(1+px)*(1+py)*(1+pz);
          shape_array[itype][n][7] = 0.125*(1-px)*(1+py)*(1+pz);
          n++;
        }
      } else {
        shape_array[itype][n][0] = 0.25*(1-px)*(1-py);
        shape_array[itype][n][1] = 0.25*(1+px)*(1-py);
        shape_array[itype][n][2] = 0.25*(1+px)*(1+py);
        shape_array[itype][n][3] = 0.25*(1-px)*(1+py);
        n++;
      }
    }
  }
}
/*-----------------------------------------------------------------------------------
  Calculate sub-element related arrays for etype itype
  - ias2ia to map atom index inside sub-element to atom index inside parent element
  - natom_subelem: number of atoms inside each sub element
  ------------------------------------------------------------------------------------*/

void ElementVecRHB::setup_sub_element(int itype)
{

  double xlo,xhi,ylo,yhi,zlo,zhi;
  int i,j,k,m,n,ntotal;
  int x,y,z;
  int isub = 0;
  double px,py,pz;
  double dx_pc,dy_pc,dz_pc;
  int nx,ny,nz;
  int esplit = element->esplit;
  double ds = 2.0/esplit;

  nx = naae[itype][0];
  ny = naae[itype][1];
  nz = naae[itype][2];

  dx_pc = 2.0/(nx-1);
  dy_pc = 2.0/(ny-1);

  if (domain->dimension == 3)
    dz_pc = 2.0/(nz-1);

  ntotal = 0;

  // makes sure atoms on element boundaries are included
  // atoms might sit right on sub-element boundaries
  // so subtract a small number from the boundaries to make sure
  // all atoms are assigned to only one sub-element

  for (x = 0; x < esplit; x++) {

    xlo = -1.0 + x*ds - EPSILON;
    if (x == (esplit-1)) xhi = 1.0 + EPSILON;
    else xhi = xlo + ds - EPSILON;

    for (y = 0; y < esplit; y++) {

      ylo = -1.0 + y*ds - EPSILON;
      if (y == (esplit-1)) yhi = 1.0 + EPSILON;
      else yhi = ylo + ds - EPSILON;

      if (domain->dimension == 3) {
        for (z = 0; z < esplit; z++) {

          zlo = -1.0 + z*ds - EPSILON;
          if (z == (esplit-1)) zhi = 1.0 + EPSILON;
          else zhi = zlo + ds - EPSILON;

          n = 0; // number of atoms inside a sub element counter
          m = 0; // index of interpolated atom in element

          // place interpolated atoms into sub element isub

          for (i = 0; i < nx; i++) {
            px = dx_pc*i - 1.0;
            for (j = 0; j < ny; j++) {
              py = dy_pc*j - 1.0;
              for (k = 0; k < nz; k++) {
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
          natom_subelem[itype][isub++] = n;
          ntotal += n;
        }
      } else {

        n = 0; // number of atoms inside a sub element counter
        m = 0; // index of interpolated atom in element

        // place interpolated atoms into sub element isub

        for (i = 0; i < nx; i++) {
          px = dx_pc*i - 1.0;
          for (j = 0; j < ny; j++) {
            py = dy_pc*j - 1.0;
            if (px >= xlo && px < xhi &&
                py >= ylo && py < yhi) {
              ias2ia[itype][isub][n++] = m;
              if (n > element->maxsubintpl) error->all(FLERR,"Too many atoms in a sub-element, boost maxsubintpl");
            }
            m++;
          }
        }
        natom_subelem[itype][isub++] = n;
        ntotal += n;
      }
    }
  }

  if (ntotal != nintpl[itype])
    error->all(FLERR,"Total atom count in sub-elements do not match");
}
/* -------------------------------------------------
   Setup sub-element related arrays in common for all element types
   Called whenever esplit is changed
   ------------------------------------------------*/

void ElementVecRHB::setup_sub_element()
{
  int esplit = element->esplit;
  if (domain->dimension == 3) nsubelem = esplit*esplit*esplit;
  else nsubelem = esplit*esplit;
  element->nsubelem = nsubelem;
  memory->destroy(element->shape_array_center_subelem);
  shape_array_center_subelem = memory->create(element->shape_array_center_subelem,
      nsubelem,npe,"element:shape_array_center_subelem");

  double ds = 2.0/esplit;
  // calculate sub-element center shape function array

  int isub = 0;
  int i,j,k; 
  double px,py,pz;

  for (i = 0; i < esplit; i++){
    px = (i + 0.5) * ds - 1.0;
    for (j = 0; j < esplit; j++){
      py = (j + 0.5) * ds - 1.0;
      if (domain->dimension == 3) {
        for (k = 0; k < esplit; k++) {
          pz = (k + 0.5) * ds - 1.0;
          shape_array_center_subelem[isub][0] = 0.125*(1-px)*(1-py)*(1-pz);
          shape_array_center_subelem[isub][1] = 0.125*(1+px)*(1-py)*(1-pz);
          shape_array_center_subelem[isub][2] = 0.125*(1+px)*(1+py)*(1-pz);
          shape_array_center_subelem[isub][3] = 0.125*(1-px)*(1+py)*(1-pz);
          shape_array_center_subelem[isub][4] = 0.125*(1-px)*(1-py)*(1+pz);
          shape_array_center_subelem[isub][5] = 0.125*(1+px)*(1-py)*(1+pz);
          shape_array_center_subelem[isub][6] = 0.125*(1+px)*(1+py)*(1+pz);
          shape_array_center_subelem[isub][7] = 0.125*(1-px)*(1+py)*(1+pz);
          isub++;
        }
      } else {
        shape_array_center_subelem[isub][0] = 0.25*(1-px)*(1-py);
        shape_array_center_subelem[isub][1] = 0.25*(1+px)*(1-py);
        shape_array_center_subelem[isub][2] = 0.25*(1+px)*(1+py);
        shape_array_center_subelem[isub][3] = 0.25*(1-px)*(1+py);
        isub++;
      }
    }
  }
}


/* -----------------------------------
   update all local element center coordinates from node coordinates
   to be called by Fix after integrating node positions
   ------------------------------------*/

void ElementVecRHB::update_center_coord()
{

  // zero center coords of elements

  size_t nbytes = sizeof(double) * element->nlocal;
  if (nbytes) {
    memset(&x[0][0],0,3*nbytes);
  }

  // update center coords of elements

  for (int i = 0; i < element->nlocal; i++) {
    for (int j = 0; j < npe; j++) {
      x[i][0] += nodex[i][j][0];
      x[i][1] += nodex[i][j][1];
      x[i][2] += nodex[i][j][2];
    }
    x[i][0] = x[i][0]/npe;
    x[i][1] = x[i][1]/npe;
    x[i][2] = x[i][2]/npe;
  }
}

/* -----------------------------------
   update local element I node idim coordinate 
   from center idim coordinate
   ------------------------------------*/

void ElementVecRHB::update_node_coord(int i, int idim)
{
  double oldx;
  oldx = 0.0;
  for (int j = 0; j < npe; j++)
    oldx += nodex[i][j][idim];
  oldx /= npe;
  for (int j = 0; j < npe; j++) 
    nodex[i][j][idim] += x[i][idim] - oldx;
}

/* -----------------------------------
   update all local element node coordinates 
   from center coordinates
   ------------------------------------*/

void ElementVecRHB::update_node_coord()
{
  double oldx[3];
  for (int i = 0; i < element->nlocal; i++) {
    oldx[0] = 0.0;
    oldx[1] = 0.0;
    oldx[2] = 0.0;
    for (int j = 0; j < npe; j++) {
      oldx[0] += nodex[i][j][0];
      oldx[1] += nodex[i][j][1];
      oldx[2] += nodex[i][j][2];
    }
    oldx[0] = oldx[0]/npe;
    oldx[1] = oldx[1]/npe;
    oldx[2] = oldx[2]/npe;
    for (int j = 0; j < npe; j++) {
      nodex[i][j][0] += x[i][0] - oldx[0];
      nodex[i][j][1] += x[i][1] - oldx[1];
      nodex[i][j][2] += x[i][2] - oldx[2];
    }
  }
}

/* --------------------------------------------------------------
   discritize element I into discrete atoms
   to be called by disc_element command
   return number of atoms created
   ---------------------------------------------------------------*/

int ElementVecRHB::element2atom(int i){

  // get the coordinates of interpolated atom k
  // move last element to element i

  int ietype = etype[i];
  int ictype = ctype[i];
  double coord[3];
  for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++){
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    for (int inode = 0; inode < npe; inode++) 
      for (int idim = 0; idim < 3; idim++)
        coord[idim] += shape_array[ietype][iintpl][inode]*nodex[i][inode][idim];
    atom->avec->create_atom(coord,ictype,0); 
  }

  copy(--element->nlocal,i,1);
  return nintpl[ietype];
}


/*-------------------------------------------------------------*/

void ElementVecRHB::check_element_size()
{
  double max = 0.0;
  int i,j,k;
  int nlocal = element->nlocal;
  double maxcoord[3],mincoord[3];
  int ndim = domain->dimension;

  for (i = 0; i < nlocal; i++) elem_size[i] = 0.0;

  // calculate the size of the element

  for (i = 0; i < nlocal; i++) {

    maxcoord[0] = maxcoord[1] = maxcoord[2] = -BIG;
    mincoord[0] = mincoord[1] = mincoord[2] = BIG;
    for (k = 0; k < ndim; k++) {
      for (j = 0; j < npe; j++) {
        maxcoord[k] = MAX(maxcoord[k],nodex[i][j][k]);
        mincoord[k] = MIN(mincoord[k],nodex[i][j][k]);
      }
      elem_size[i] = MAX(elem_size[i],maxcoord[k]-mincoord[k]);
    }
    if (elem_size[i] <= 0.0) error->one(FLERR,"Negative element size");
    if (initial_size[i] == 0.0) initial_size[i] = elem_size[i];
    max = MAX(max,elem_size[i]);
  }

  // update element->max_size

  MPI_Allreduce(&max,&element->max_size,1,MPI_DOUBLE,MPI_MAX,world);

  // check element deformation

  for (i = 0; i < nlocal; i++) {
    if (elem_size[i]/initial_size[i] > element->maxelemchange || initial_size[i]/elem_size[i] > element->maxelemchange) {
      char *erstr = new char[100];
      sprintf(erstr,"Element ID %d deformed by more than %3.1f times",tag[i],element->maxelemchange);
      error->one(FLERR,erstr);
      delete [] erstr;
    }
  } 
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint ElementVecRHB::memory_usage()
{
  bigint bytes = 0;

  if (element->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (element->memcheck("nodetag")) bytes += memory->usage(nodetag,nmax,npe);
  if (element->memcheck("etype")) bytes += memory->usage(etype,nmax);
  if (element->memcheck("ctype")) bytes += memory->usage(ctype,nmax);
  if (element->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (element->memcheck("nodemask")) bytes += memory->usage(nodemask,nmax,npe);
  if (element->memcheck("image")) bytes += memory->usage(image,nmax);
  if (element->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (element->memcheck("slip_plane")) bytes += memory->usage(slip_plane,3,3);
  if (element->memcheck("nodex")) bytes += memory->usage(nodex,nmax,npe,3);
  if (element->memcheck("nodev")) bytes += memory->usage(nodev,nmax,npe,3);
  if (element->memcheck("nodef")) bytes += memory->usage(nodef,nmax,npe,3);

  bytes += memory->usage(intg,element->netypes+1,element->maxintg,3);

  return bytes;
}

/* ----------------------------------------------------------------------
   write interpolate information to Interpolate section in data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::write_interpolate(FILE *fp)
{
  for (int i = 1; i <= element->netypes; i++) 
    fprintf(fp,"%d %d %d %d\n",i,naae[i][0],naae[i][1],naae[i][2]);
}

/* ----------------------------------------------------------------------
   write integration information to Integration section in data file
   ------------------------------------------------------------------------- */

void ElementVecRHB::write_integration(FILE *fp)
{
  if (element->intg_input_style == 2) {
    for (int i = 1; i <= element->netypes; i++)
      fprintf(fp,"%d %d\n",i,nintg[i]);
    for (int i = 1; i <= element->netypes; i++)
      for (int j = 0; j < nintg[i]; j++)
        fprintf(fp,"%d %d %d %d %d %f\n",i,j+1,
            intg[i][j][0],intg[i][j][1],intg[i][j][2],weight[i][j]);
  } else {
    for (int i = 1; i <= element->netypes; i++) 
      fprintf(fp,"%d %d %d %d\n",i,niae[i][0],niae[i][1],niae[i][2]);
  }
}

/* ----------------------------------------------------------------------
   update internal slip_plane (axes) of all elements
   slip_plane is defined as unit normal vector of the plane
   ------------------------------------------------------------------------- */

void ElementVecRHB::update_slip_plane()
{
  double norm;
  for (int i = 0; i < element->nlocal; i++) {
    if (domain->dimension == 3) {
      for (int j = 0; j < 3; j++) {

        // yz plane (x axis)

        slip_plane[i][0][j] = 
          nodex[i][2][j] - nodex[i][1][j] + 
          nodex[i][3][j] - nodex[i][4][j] +
          nodex[i][6][j] - nodex[i][5][j] +
          nodex[i][7][j] - nodex[i][8][j];

        // xz plane (y axis)

        slip_plane[i][1][j] = 
          nodex[i][4][j] - nodex[i][1][j] +
          nodex[i][3][j] - nodex[i][2][j] +
          nodex[i][8][j] - nodex[i][5][j] +
          nodex[i][7][j] - nodex[i][6][j];

        // xy plane (z axis)

        slip_plane[i][2][j] = 
          nodex[i][5][j] - nodex[i][1][j] +
          nodex[i][8][j] - nodex[i][4][j] +
          nodex[i][6][j] - nodex[i][2][j] +
          nodex[i][7][j] - nodex[i][3][j];
      }
    } else {
      for (int j = 0; j < 2; j++) {

        // yz plane (x axis)

        slip_plane[i][0][j] = 
          nodex[i][2][j] - nodex[i][1][j] + 
          nodex[i][3][j] - nodex[i][4][j];

        // xz plane (y axis)

        slip_plane[i][1][j] = 
          nodex[i][4][j] - nodex[i][1][j] +
          nodex[i][3][j] - nodex[i][2][j];

        // xy plane (z axis)

        slip_plane[i][2][j] = 0.0;
      }
      slip_plane[i][0][2] = 0.0;
      slip_plane[i][1][2] = 0.0;
      slip_plane[i][2][2] = 1.0;
    }

    // normalize vectors

    norm3(slip_plane[i][0]);
    norm3(slip_plane[i][1]);
    norm3(slip_plane[i][2]);
  }
}

/* ----------------------------------------------------------------------
   interpolate values from nodevalues in element I at interpolated atom J
   n = number of values
   ------------------------------------------------------------------------- */

void ElementVecRHB::interpolate(double *value, double ***nodevalue, int i, int j, int n) 
{
  int ietype = etype[i];
  for (int k = 0; k < n; k++) {
    value[k] = 0.0;
    for (int l = 0; l < npe; l++) 
      value[k] += shape_array[ietype][j][l]*nodevalue[i][l][k];
  }
}

/* ----------------------------------------------------------------------
   add atom plane next to elements
   Called from input add_atoms command
   ------------------------------------------------------------------------- */

void ElementVecRHB::add_atoms(int narg, char **arg)
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
  int iintpl,nx,ny,nz;
  int natoms_added = 0;
  for (int i = 0; i < element->nlocal; i++) 
    if (mask[i] & groupbit) {
      nx = naae[etype[i]][0];
      ny = naae[etype[i]][1];
      nz = naae[etype[i]][2];
      for (int n = 0; n < nlayers; n++) {
        for (int j = 0; j < naae[etype[i]][dim1]; j++) {
          for (int k = 0; k < naae[etype[i]][dim2]; k++) {
            iintpl = (j*(dim1 == 2) + k*(dim2 == 2) + (nz-1)*((dim1!=2)&&(dim2!=2)&&(n1list[0]==4)))
              + nz * (j*(dim1 == 1) + k*(dim2 == 1) + (ny-1)*((dim1!=1)&&(dim2!=1)&&(n1list[0]==2)))
              + ny * nz * (j*(dim1 == 0) + k*(dim2 == 0) + (nx-1)*((dim1!=0)&&(dim2!=0)&&(n1list[0]==1)));
            coord[0] = coord[1] = coord[2] = 0.0;
            for (int l = 0; l < npe; l++) {
              coord[0] += shape_array[etype[i]][iintpl][l]*nodex[i][l][0];
              coord[1] += shape_array[etype[i]][iintpl][l]*nodex[i][l][1];
              coord[2] += shape_array[etype[i]][iintpl][l]*nodex[i][l][2];
            }
            atom->avec->create_atom(coord,ctype[i],0); 
            natoms_added++;
          }
        }

        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 3; k++) 
            nodex[i][n1list[j]][k] += (nodex[i][n2list[j]][k] - nodex[i][n1list[j]][k])/naae[etype[i]][dim0];
        }
      }
    }

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

bigint ElementVecRHB::create_element_template(
    int klo, int khi, int jlo, int jhi, int ilo, int ihi, 
    int ietype, int *basisctype, double **xlist, int *clist)
{
  bigint n = 0;

  bigint maxlist = (khi-klo+1)*(jhi-jlo+1)*(ihi-ilo+1);
  int nbasis = domain->lattice->nbasis;
  double **basis = domain->lattice->basis;
  int i,j,k,m;
  for (k = klo; k <= khi; k += naae[ietype][2]) 
    for (j = jlo; j <= jhi; j += naae[ietype][1]) 
      for (i = ilo; i <= ihi; i += naae[ietype][0]) 
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


