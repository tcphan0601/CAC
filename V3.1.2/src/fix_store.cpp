#include <cstdlib>
#include <cstring>
#include "fix_store.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

enum{GLOBAL, PERATOM};
enum{ELEMENT, NODE, VATOM};

/*  ----------------------------------------------------------------------  */

FixStore::FixStore(CAC *cac, int narg, char **arg) : Fix(cac, narg, arg), 
vstore(nullptr), astore(nullptr), a3store(nullptr), a4store(nullptr), rbuf(nullptr)
{
  if (narg != 7) error->all(FLERR, "Illegal fix store command");

  // 4th arg determines GLOBAL vs PERATOM values
  // 5th arg determines ELEMENT, NODE or VATOM values for elements
  // syntax: id group style global nrow ncol
  //   GLOBAL: Nrow by Ncol array of global values
  //   PERATOM: Ncol = 1 is vector, use vstore for atoms, astore for elements with Element style, for NODE & VATOM, use a3store
  //            Ncol > 1 is array, use astore for atoms, a3store for elements with ELEMENT style, for NODE & VATOM, use a4store
  // syntax: id group style peratom 0/1 nvalues
  //   0/1 flag = not-store or store peratom values in restart file
  //   nvalues = # of peratom values
 
  disable = 0;
  nvalues = vecflag = 0;

  if (strcmp(arg[3], "global") == 0) flavor = GLOBAL;
  else if (strcmp(arg[3], "peratom") == 0) flavor = PERATOM;
  else error->all(FLERR, "Illegal fix store command");

  if (strcmp(arg[4], "element") == 0) elem_flavor = ELEMENT;
  else if (strcmp(arg[4], "node") == 0) elem_flavor = NODE;
  else if (strcmp(arg[4], "vatom") == 0) elem_flavor = VATOM;
  else error->all(FLERR, "Illegal fix store command");


  // GLOBAL values are always written to restart file
  // PERATOM restart_peratom is set by caller

  if (flavor == GLOBAL) {
    restart_global = 1;
    nrow = universe->inumeric(FLERR, arg[5]);
    ncol = universe->inumeric(FLERR, arg[6]);
    if (nrow <= 0 || ncol <= 0)
      error->all(FLERR, "Illegal fix store command");
    vecflag = 0;
    if (ncol == 1) vecflag = 1;
  } else {
    restart_peratom = universe->inumeric(FLERR, arg[5]);
    nvalues = universe->inumeric(FLERR, arg[6]);
    if (restart_peratom < 0 or restart_peratom > 1 || nvalues <= 0)
      error->all(FLERR, "Illegal fix store command");
    vecflag = 0;
    if (nvalues == 1) vecflag = 1;
  }

  vstore = nullptr;
  astore = nullptr;
  a3store = nullptr;
  a4store = nullptr;

  nbasis = element->maxapc;
  if (elem_flavor == ELEMENT) ncell = 1;
  else if (elem_flavor == NODE) ncell = element->maxnpe;
  else ncell = element->maxucell;

  // allocate vector or array and restart buffer rbuf
  // for PERATOM, register with Atom & Element class

  if (flavor == GLOBAL) {
    if (vecflag) memory->create(vstore, nrow, "fix/store:vstore");
    else memory->create(astore, nrow, ncol, "fix/store:astore");
    memory->create(rbuf, nrow * ncol+2, "fix/store:rbuf");
  } else {
    if (atom->nmax) grow_atom_arrays(atom->nmax);
    else grow_atom_arrays(1);
    if (element->nmax) grow_elem_arrays(element->nmax);
    else grow_elem_arrays(1);

    atom->add_callback(0);
    element->add_callback(0);
    if (restart_peratom) {
      atom->add_callback(1);
      element->add_callback(1);
    }
    rbuf = nullptr;
  }

  // zero the storage
  // PERATOM may be comm->exchanged before filled by caller

  if (flavor == GLOBAL) {
    if (vecflag)
      for (int i = 0; i < nrow; i++) vstore[i] = 0.0;
    else
      for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
          astore[i][j] = 0.0;
  } else {
    int nalocal = atom->nlocal;
    int nelocal = element->nlocal;
    if (vecflag) {
      if (elem_flavor == ELEMENT) {
        for (int i = 0; i < nelocal; i++) 
          for (int j = 0; j < ncell; j++) 
            astore[i][j] = 0.0;
      } else {
        for (int i = 0; i < nelocal; i++) 
          for (int j = 0; j < nbasis; j++) 
            for (int k = 0; k < ncell; k++) 
              a3store[i][j][k] = 0.0;
      }
    } else {
      if (elem_flavor == ELEMENT) {
        for (int i = 0; i < nelocal; i++)
          for (int j = 0; j < ncell; j++)
            for (int k = 0; k < nvalues; k++)
              a3store[i][j][k] = 0.0;
      } else {
        for (int i = 0; i < nelocal; i++)
          for (int j = 0; j < nbasis; j++)
            for (int k = 0; k < ncell; k++)
              for (int l = 0; l < nvalues; l++) {
                a4store[i][j][k][l] = 0.0;
              }
      }
    }
  }
}

/*  ----------------------------------------------------------------------  */

FixStore::~FixStore()
{
  // unregister callbacks to this fix from Atom class

  if (flavor == PERATOM) {
    atom->delete_callback(id, 0);
    element->delete_callback(id, 0);
    if (restart_peratom) {
      atom->delete_callback(id, 1);
      element->delete_callback(id, 1);
    }
  }

  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(a3store);
  memory->destroy(a4store);
  memory->destroy(rbuf);
}

/*  ----------------------------------------------------------------------  */

int FixStore::setmask()
{
  int mask = 0;
  return mask;
}

/*  ----------------------------------------------------------------------
    reset size of global vector/array
    invoked by caller if size is unknown at time this fix is instantiated
    caller will do subsequent initialization
    -------------------------------------------------------------------------  */

void FixStore::reset_global(int nrow_caller, int ncol_caller)
{
  memory->destroy(vstore);
  memory->destroy(astore);
  memory->destroy(rbuf);
  vstore = nullptr;
  astore = nullptr;

  vecflag = 0;
  if (ncol_caller == 1) vecflag = 1;
  nrow = nrow_caller;
  ncol = ncol_caller;
  if (vecflag) memory->create(vstore, nrow, "fix/store:vstore");
  else memory->create(astore, nrow, ncol, "fix/store:astore");
  memory->create(rbuf, nrow * ncol+2, "fix/store:rbuf");
}

/*  ----------------------------------------------------------------------
    write global array to restart file
    -------------------------------------------------------------------------  */

//void FixStore::write_restart(FILE *fp)
//{
//  // fill rbuf with size and vec/array values
//
//  rbuf[0] = nrow;
//  rbuf[1] = ncol;
//  if (vecflag) memcpy(&rbuf[2], vstore, nrow * sizeof(double));
//  else memcpy(&rbuf[2], &astore[0][0], nrow * ncol * sizeof(double));
//
//  int n = nrow * ncol + 2;
//  if (comm->me == 0) {
//    int size = n * sizeof(double);
//    fwrite(&size, sizeof(int), 1, fp);
//    fwrite(rbuf, sizeof(double), n, fp);
//  }
//}

/*  ----------------------------------------------------------------------
    use global array from restart file to restart the Fix
    -------------------------------------------------------------------------  */

//void FixStore::restart(char *buf)
//{
//  // first 2 values in buf are vec/array sizes
//
//  double *dbuf = (double *) buf;
//  int nrow_restart = dbuf[0];
//  int ncol_restart = dbuf[1];
//
//  // if size of vec/array has changed, 
//  //   means the restart file is setting size of vec or array and doing init
//  //   because caller did not know size at time this fix was instantiated
//  // reallocate vstore or astore accordingly
//
//  if (nrow != nrow_restart || ncol != ncol_restart) {
//    memory->destroy(vstore);
//    memory->destroy(astore);
//    memory->destroy(rbuf);
//    vstore = nullptr;
//    astore = nullptr;
//
//    vecflag = 0;
//    if (ncol_restart == 1) vecflag = 1;
//    nrow = nrow_restart;
//    ncol = ncol_restart;
//    if (vecflag) memory->create(vstore, nrow, "fix/store:vstore");
//    else memory->create(astore, nrow, ncol, "fix/store:astore");
//    memory->create(rbuf, nrow * ncol+2, "fix/store:rbuf");
//  }
//
//  int n = nrow * ncol;
//  if (vecflag) memcpy(vstore, &dbuf[2], n * sizeof(double));
//  else memcpy(&astore[0][0], &dbuf[2], n * sizeof(double));
//}

/*  ----------------------------------------------------------------------
    allocate atom-based array
    -------------------------------------------------------------------------  */

void FixStore::grow_atom_arrays(int nmax)
{
  if (vecflag) memory->grow(vstore, nmax, "store:vstore");
  else memory->grow(astore, nmax, nvalues, "store:astore");
}

/*  ----------------------------------------------------------------------
    allocate elem-based array
    -------------------------------------------------------------------------  */

void FixStore::grow_elem_arrays(int nmax)
{
  if (vecflag) {
    if (elem_flavor == ELEMENT)
      memory->grow(astore, nmax, ncell, "store:astore");
    else 
      memory->grow(a3store, nmax, nbasis, ncell, "store:a3store");
  } else {
    if (elem_flavor == ELEMENT)
      memory->grow(a3store, nmax, ncell, nvalues, "store:a3store");
    else 
      memory->grow(a4store, nmax, nbasis, ncell, nvalues, "store:a4store");

  }
}


/*  ----------------------------------------------------------------------
    copy values within local atom-based array
    -------------------------------------------------------------------------  */

void FixStore::copy_atom_arrays(int i, int j, int /* delflag */)
{
  if (disable) return;

  if (vecflag) vstore[j] = vstore[i];
  else
    for (int m = 0; m < nvalues; m++)
      astore[j][m] = astore[i][m];
}

/*  ----------------------------------------------------------------------
    copy values within local elem-based array
    -------------------------------------------------------------------------  */

void FixStore::copy_elem_arrays(int i, int j, int /* delflag */)
{
  if (disable) return;

  if (vecflag) {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        astore[j][k] = astore[i][k];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          a3store[j][k][l] = a3store[i][k][l];
    }
  } else {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        for (int m = 0; m < nvalues; m++)
          a3store[j][k][m] = a3store[i][k][m];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          for (int m = 0; m < nvalues; m++)
            a4store[j][k][l][m] = a4store[i][k][l][m];
    }
  }
}

/*  ----------------------------------------------------------------------
    pack values in local atom-based array for exchange with another proc
    -------------------------------------------------------------------------  */

int FixStore::pack_atom_exchange(int i, double *buf)
{
  if (disable) return 0;

  if (vecflag) buf[0] = vstore[i];
  else
    for (int m = 0; m < nvalues; m++)
      buf[m] = astore[i][m];
  return nvalues;
}

/*  ----------------------------------------------------------------------
    pack values in local elem-based array for exchange with another proc
    -------------------------------------------------------------------------  */

int FixStore::pack_elem_exchange(int i, double *buf)
{
  if (disable) return 0;

  int j = 0;

  if (vecflag) {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        buf[j++] = astore[i][k];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          buf[j++] = a3store[i][k][l];
    }
  } else {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        for (int m = 0; m < nvalues; m++)
          buf[j++] = a3store[i][k][m];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          for (int m = 0; m < nvalues; m++)
            buf[j++] = a4store[i][k][l][m];
    }
  }

  return j;
}

/*  ----------------------------------------------------------------------
    unpack values in local atom-based array from exchange with another proc
    -------------------------------------------------------------------------  */

int FixStore::unpack_atom_exchange(int nlocal, double *buf)
{
  if (disable) return 0;

  if (vecflag) vstore[nlocal] = buf[0];
  else
    for (int m = 0; m < nvalues; m++)
      astore[nlocal][m] = buf[m];
  return nvalues;
}

/*  ----------------------------------------------------------------------
    unpack values in local elem-based array from exchange with another proc
    -------------------------------------------------------------------------  */

int FixStore::unpack_elem_exchange(int nlocal, double *buf)
{
  if (disable) return 0;

  int m = 0;

  if (vecflag) {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        astore[nlocal][k] = buf[m++];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          a3store[nlocal][k][l] = buf[m++];
    }
  } else {
    if (elem_flavor == ELEMENT) {
      for (int k = 0; k < ncell; k++)
        for (int n = 0; n < nvalues; n++)
          a3store[nlocal][k][n] = buf[m++];
    } else {
      for (int k = 0; k < nbasis; k++)
        for (int l = 0; l < ncell; l++)
          for (int n = 0; n < nvalues; n++)
            a4store[nlocal][k][l][n] = buf[m++];
    }
  }

  return m;
}

/*  ----------------------------------------------------------------------
    pack values in local atom-based arrays for restart file
    -------------------------------------------------------------------------  */

//int FixStore::pack_restart(int i, double *buf)
//{
//  if (disable) {
//    buf[0] = 0;
//    return 1;
//  }
//
//  buf[0] = nvalues+1;
//  if (vecflag) buf[1] = vstore[i];
//  else
//    for (int m = 0; m < nvalues; m++)
//      buf[m+1] = astore[i][m];
//  return nvalues+1;
//}

/*  ----------------------------------------------------------------------
    unpack values from atom->extra array to restart the fix
    -------------------------------------------------------------------------  */

//void FixStore::unpack_restart(int nlocal, int nth)
//{
//  if (disable) return;
//
//  double **extra = atom->extra;
//
//  // skip to Nth set of extra values
//
//  int m = 0;
//  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
//  m++;
//
//  if (vecflag) vstore[nlocal] = extra[nlocal][m];
//  else
//    for (int i = 0; i < nvalues; i++)
//      astore[nlocal][i] = extra[nlocal][m++];
//}

/*  ----------------------------------------------------------------------
    maxsize of any atom's restart data
    -------------------------------------------------------------------------  */

//int FixStore::maxsize_restart()
//{
//  if (disable) return 1;
//  return nvalues+1;
//}

/*  ----------------------------------------------------------------------
    size of atom nlocal's restart data
    -------------------------------------------------------------------------  */

//int FixStore::size_restart(int /* nlocal */)
//{
//  if (disable) return 1;
//  return nvalues+1;
//}

/*  ----------------------------------------------------------------------
    memory usage of local atom-based/elem-based arrays
    -------------------------------------------------------------------------  */

double FixStore::memory_usage()
{
  double bytes = 0.0;

  if (flavor == GLOBAL) bytes += nrow * ncol * sizeof(double);
  else {
    if (elem_flavor == ELEMENT) 
      bytes += (atom->nmax + element->nmax * ncell) * nvalues * sizeof(double);
    else 
      bytes += (atom->nmax + element->nmax * ncell * nbasis) * nvalues * sizeof(double);
  }
  return bytes;
}
