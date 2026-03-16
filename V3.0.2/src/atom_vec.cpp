#include <string.h>
#include <stdlib.h>
#include "atom_vec.h"
#include "atom.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;

#define DELTA 16384

/* ---------------------------------------------------------------------- */

AtomVec::AtomVec(CAC *cac) : Pointers(cac)
{
  nmax = 0;
  mass_type = 0;
  forceclearflag = 0;
  size_data_bonus = 0;

  nargcopy = 0;
  argcopy = NULL;
}

/* ----------------------------------------------------------------------- */

AtomVec::~AtomVec()
{
  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;
}

/*  ----------------------------------------------------------------------
     make copy of args for use by restart & replicate
    ------------------------------------------------------------------------- */

void AtomVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  argcopy = new char*[nargcopy];
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(arg[i]) + 1;
    argcopy[i] = new char[n];
    strcpy(argcopy[i], arg[i]);
  }
}

/*  ----------------------------------------------------------------------
    no additional args by default
  -------------------------------------------------------------------------  */

void AtomVec::process_args(int narg, char **arg)
{
  if (narg) error->all(FLERR, "Invalid atom_style command");
}

/*  ----------------------------------------------------------------------
  grow nmax so it is a multiple of DELTA
 -------------------------------------------------------------------------  */

void AtomVec::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
}

/*  ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
-------------------------------------------------------------------------  */

void AtomVec::data_vel(int m, char **values)
{
  double **v = atom->v;
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
}

/*  ----------------------------------------------------------------------
   pack velocity info for data file
-------------------------------------------------------------------------  */

void AtomVec::pack_vel(double **buf)
{
  double **v = atom->v;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
  }
}

/*  ----------------------------------------------------------------------
   write velocity info to data file
-------------------------------------------------------------------------  */

void AtomVec::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e\n", 
            (tagint) ubuf(buf[i][0]).i, buf[i][1], buf[i][2], buf[i][3]);
}

/*  ----------------------------------------------------------------------
   copy of velocity remap settings from Domain
-------------------------------------------------------------------------  */

void AtomVec::init()
{
  deform_vremap = domain->deform_vremap;
  deform_groupbit = domain->deform_groupbit;
  h_rate = domain->h_rate;
}

/*  ----------------------------------------------------------------------
   roundup N so it is a multiple of DELTA
   error if N exceeds 32-bit int, since will be used as arg to grow()
-------------------------------------------------------------------------  */

bigint AtomVec::roundup(bigint n)
{
  if (n % DELTA) n = n/DELTA * DELTA + DELTA;
  if (n > MAXSMALLINT)
    error->one(FLERR, "Too many atoms created on one or more procs");
  return n;
}


