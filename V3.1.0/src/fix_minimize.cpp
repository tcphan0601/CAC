#include "fix_minimize.h"
#include "atom.h"
#include "element.h"
#include "domain.h"
#include "memory.h"

#include "comm.h"
using namespace CAC_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMinimize::FixMinimize(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  nvector(0), peratom(nullptr), atomvectors(nullptr), nodevectors(nullptr)
{
  // register callback to this fix from Atom class
  // don't perform initial allocation here, must wait until add_vector()

  atom->add_callback(0);
  element->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixMinimize::~FixMinimize()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, 0);
  element->delete_callback(id, 0);

  // delete locally stored data

  memory->destroy(peratom);
  if (atomvectors) {
    for (int m = 0; m < nvector; m++) memory->destroy(atomvectors[m]);
    memory->sfree(atomvectors);
  }
  if (nodevectors) {
    for (int m = 0; m < nvector; m++) memory->destroy(nodevectors[m]);
    memory->sfree(nodevectors);
  }

}

/* ---------------------------------------------------------------------- */

int FixMinimize::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixMinimize::add_vector(int n)
{
  int ntotal;
  memory->grow(peratom, nvector + 1, "minimize:peratom");
  peratom[nvector] = n;

  atomvectors = (double **)
    memory->srealloc(atomvectors, (nvector + 1) * sizeof(double *), "minimize:atomvectors");
  memory->create(atomvectors[nvector], atom->nmax * n, "minimize:atomvector");
  ntotal = n * atom->nlocal;
  for (int i = 0; i < ntotal; i++) atomvectors[nvector][i] = 0.0;

  nodevectors = (double **)
    memory->srealloc(nodevectors, (nvector + 1) * sizeof(double *), "minimize:nodevectors");
  memory->create(nodevectors[nvector], element->nmax * element->maxapc * element->maxnpe * n, "minimize:nodevector");
  ntotal = n * element->nlocal * element->maxapc * element->maxnpe;
  for (int i = 0; i < ntotal; i++) nodevectors[nvector][i] = 0.0;

  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double *FixMinimize::request_vector(int m, int which)
{
  if (which == 0) return atomvectors[m];
  else return nodevectors[m];
}

/* ----------------------------------------------------------------------
   store box size at beginning of line search
------------------------------------------------------------------------- */

void FixMinimize::store_box()
{
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
  boxhi[0] = domain->boxhi[0];
  boxhi[1] = domain->boxhi[1];
  boxhi[2] = domain->boxhi[2];
}

/* ----------------------------------------------------------------------
   reset x0 for atoms that moved across PBC via reneighboring in line search
   x0 = 1st vector
   must do minimum_image using original box stored at beginning of line search
   swap & set_global_box() change to original box, then restore current box
------------------------------------------------------------------------- */

void FixMinimize::reset_coords()
{
  box_swap();
  domain->set_global_box();

  double **x = atom->x;
  double ****nodex = element->nodex;
  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  double dx, dy, dz, dx0, dy0, dz0, *x0;
  int ietype, n;
  int *apc = element->apc;
  int *npe = element->npe;
  int *etype = element->etype;

  n = 0;
  x0 = atomvectors[0];
  for (int i = 0; i < nalocal; i++) {
    dx = dx0 = x[i][0] - x0[n];
    dy = dy0 = x[i][1] - x0[n + 1];
    dz = dz0 = x[i][2] - x0[n + 2];
    domain->minimum_image(dx, dy, dz);
    if (dx != dx0) x0[n] = x[i][0] - dx;
    if (dy != dy0) x0[n + 1] = x[i][1] - dy;
    if (dz != dz0) x0[n + 2] = x[i][2] - dz;
    n += 3;
  }
 
  n = 0;
  x0 = nodevectors[0];
  int nper_elem = 3 * element->maxapc * element->maxnpe;
  for (int i = 0; i < nelocal; i++) {
    // use first node to compute dx, dy, dz
    // since it will be the same for all other nodes
    
    ietype = etype[i];
    dx = dx0 = nodex[i][0][0][0] - x0[n];
    dy = dy0 = nodex[i][0][0][1] - x0[n + 1];
    dz = dz0 = nodex[i][0][0][2] - x0[n + 2];
    domain->minimum_image(dx, dy, dz);
    int m = 0;
    for (int ibasis = 0; ibasis < apc[ietype]; ibasis++)
      for (int inode = 0; inode < npe[ietype]; inode++) {
        if (dx != dx0) x0[n + m] = nodex[i][ibasis][inode][0] - dx;
        if (dy != dy0) x0[n + m + 1] = nodex[i][ibasis][inode][1] - dy;
        if (dz != dz0) x0[n + m + 2] = nodex[i][ibasis][inode][2] - dz;
        m += 3;
      }
    n += nper_elem;
  }

  box_swap();
  domain->set_global_box();
}

/* ----------------------------------------------------------------------
   swap current box size with stored box size
   ------------------------------------------------------------------------- */

void FixMinimize::box_swap()
{
  double tmp;

  tmp = boxlo[0];
  boxlo[0] = domain->boxlo[0];
  domain->boxlo[0] = tmp;
  tmp = boxlo[1];
  boxlo[1] = domain->boxlo[1];
  domain->boxlo[1] = tmp;
  tmp = boxlo[2];
  boxlo[2] = domain->boxlo[2];
  domain->boxlo[2] = tmp;

  tmp = boxhi[0];
  boxhi[0] = domain->boxhi[0];
  domain->boxhi[0] = tmp;
  tmp = boxhi[1];
  boxhi[1] = domain->boxhi[1];
  domain->boxhi[1] = tmp;
  tmp = boxhi[2];
  boxhi[2] = domain->boxhi[2];
  domain->boxhi[2] = tmp;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
   ------------------------------------------------------------------------- */

double FixMinimize::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++) {
    bytes += atom->nmax * peratom[m] * sizeof(double);
    bytes += element->nmax * element->maxapc * element->maxnpe * 
      peratom[m] * sizeof(double);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
   ------------------------------------------------------------------------- */

void FixMinimize::grow_atom_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
    memory->grow(atomvectors[m], peratom[m] * nmax, "minimize:atomvector");
}

/* ----------------------------------------------------------------------
   allocate local element-based arrays
   ------------------------------------------------------------------------- */

void FixMinimize::grow_elem_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
    memory->grow(nodevectors[m], peratom[m] * nmax, "minimize:elemvector");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
   ------------------------------------------------------------------------- */

void FixMinimize::copy_atom_arrays(int i, int j, int /*delflag*/)
{
  int m, iper, nper, ni, nj;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper * i;
    nj = nper * j;
    for (iper = 0; iper < nper; iper++) 
      atomvectors[m][nj++] = atomvectors[m][ni++];
  }
}

/* ----------------------------------------------------------------------
   copy values within local element-based arrays
   ------------------------------------------------------------------------- */

void FixMinimize::copy_elem_arrays(int i, int j, int /*delflag*/)
{
  int m, iper, nper, ni, nj;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m] * element->maxapc * element->maxnpe;
    ni = nper * i;
    nj = nper * j;
    for (iper = 0; iper < nper; iper++) 
      nodevectors[m][nj++] = nodevectors[m][ni++];
  }
}
/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
   ------------------------------------------------------------------------- */

int FixMinimize::pack_atom_exchange(int i, double *buf)
{
  int m, iper, nper, ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper * i;
    for (iper = 0; iper < nper; iper++) buf[n++] = atomvectors[m][ni++];
  }
  return n;
}


/* ----------------------------------------------------------------------
   pack values in local element-based arrays for exchange with another proc
   ------------------------------------------------------------------------- */

int FixMinimize::pack_elem_exchange(int i, double *buf)
{
  int m, iper, nper, ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m] * element->maxapc * element->maxnpe;
    ni = nper * i;
    for (iper = 0; iper < nper; iper++) buf[n++] = nodevectors[m][ni++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
   ------------------------------------------------------------------------- */

int FixMinimize::unpack_atom_exchange(int nlocal, double *buf)
{
  int m, iper, nper, ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper * nlocal;
    for (iper = 0; iper < nper; iper++) atomvectors[m][ni++] = buf[n++];
  }
  return n;
}
/* ----------------------------------------------------------------------
   unpack values in local elem-based arrays from exchange with another proc
   ------------------------------------------------------------------------- */

int FixMinimize::unpack_elem_exchange(int nlocal, double *buf)
{
  int m, iper, nper, ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m] * element->maxapc * element->maxnpe;
    ni = nper * nlocal;
    for (iper = 0; iper < nper; iper++) nodevectors[m][ni++] = buf[n++];
  }
  return n;
}
