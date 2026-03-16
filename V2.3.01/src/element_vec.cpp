#include <string.h>
#include <stdlib.h>
#include "element_vec.h"
#include "element.h"
#include "domain.h"
#include "error.h"

using namespace CAC_NS;

#define BIGBIGDELTA 16384       // for nmaxintpl increment
#define BIGDELTA 4096           // for nmaxintg increment
#define DELTA 1024              // for nmax increment

/*----------------------------------------------------------------------*/

ElementVec::ElementVec(CAC *cac) : Pointers(cac)
{
  nmax = nmaxintg = 0;
  nmaxintpl = 0;
  size_data_bonus = 0;
  nargcopy = 0;
  argcopy = NULL;
  subelemflag = 0;
}

/*-----------------------------------------------------------------------*/

ElementVec::~ElementVec()
{
  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;
}

/* ----------------------------------------------------------------------
     make copy of args for use by restart & replicate
    -------------------------------------------------------------------------*/

void ElementVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  argcopy = new char*[nargcopy];
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(arg[i]) + 1;
    argcopy[i] = new char[n];
    strcpy(argcopy[i],arg[i]);
  }
}

/* ----------------------------------------------------------------------
    no additional args by default
  ------------------------------------------------------------------------- */

void ElementVec::process_args(int narg, char **arg)
{
  if (narg) error->all(FLERR,"Invalid element_style command");
}

/* ----------------------------------------------------------------------
  grow nmax so it is a multiple of DELTA
 ------------------------------------------------------------------------- */

void ElementVec::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
  element->nmax = nmax;
}

/*----------------------------------------------------------------------------*/
void ElementVec::grow_nmaxintg()
{
  nmaxintg = nmaxintg/BIGDELTA *BIGDELTA;
  nmaxintg += BIGDELTA;
  element->nmaxintg = nmaxintg;
}

/*----------------------------------------------------------------------------*/
void ElementVec::grow_nmaxintpl()
{
  nmaxintpl = nmaxintpl/BIGBIGDELTA *BIGBIGDELTA;
  nmaxintpl += BIGBIGDELTA;
  element->nmaxintpl = nmaxintpl;
}
/* ---------------------------------------------------------------------- */

void ElementVec::init()
{
  check_element_size();
  if (subelemflag == 0) {
    setup_sub_element();
    subelemflag = 1;
  }
  check_type();
}



