#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "compute.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "force.h"

#include "update.h"

using namespace CAC_NS;

#define DELTA 4
#define BIG MAXTAGINT

// allocate space for static class instance variable and initialize it

int Compute::instance_total = 0;

/*  ----------------------------------------------------------------------  */

Compute::Compute(CAC *cac, int narg, char **arg) : Pointers(cac), 
  id(nullptr), style(nullptr), 
  vector(nullptr), array(nullptr), 
  vector_atom(nullptr), array_atom(nullptr), 
  vector_node(nullptr), array_node(nullptr), 
  vector_vatom(nullptr), array_vatom(nullptr), 
  vector_local(nullptr), array_local(nullptr), 
  extlist(nullptr), tlist(nullptr), vbiasall(nullptr)

{
  instance_me = instance_total++;

  if (narg < 3) error->all(FLERR, "Illegal compute command");

  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id, arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR, 
                 "Compute ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR, "Could not find compute group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style, arg[2]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  local_flag = 0;
  group_flag = Group::ELEMENT;
  peratom_flag = local_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;
  size_vector_local_variable = size_array_local_rows_variable = 0;

  preneighflag = 0;
  tempflag = pressflag = peflag = 0;
  pressatomflag = mechatomflag = peatomflag = 0;
  ghostskinflag = 0;
  create_attribute = 0;
  tempbias = 0;

  timeflag = 0;
  comm_atom_forward = comm_atom_reverse = 0;
  comm_elem_forward = comm_elem_reverse = 0;
  dynamic = 0;
  dynamic_group_allow = 1;
  cudable = 0;
  ghostskin = 0;

  invoked_peratom = invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_scalar_local = invoked_vector_local = invoked_array_local = -1;
  invoked_flag = 0;

  // set modify defaults

  extra_dof = domain->dimension;
  dynamic_user = 0;
  fix_dof = 0;

  // setup list of timesteps

  ntime = maxtime = 0;


  // data masks

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

  // force init to zero in case these are used as logicals

}

/*  ----------------------------------------------------------------------  */

Compute::~Compute()
{
  delete [] id;
  delete [] style;
  memory->destroy(tlist);
}

/*  ----------------------------------------------------------------------
  reset extra_dof to its default value
 -------------------------------------------------------------------------  */

void Compute::reset_extra_dof()
{
  extra_dof = domain->dimension;
}

/*  ----------------------------------------------------------------------
   add ntimestep to list of timesteps the compute will be called on
   do not add if already in list
   search from top downward, since list of times is in decreasing order
-------------------------------------------------------------------------  */

void Compute::addstep(bigint ntimestep)
{
  // i = location in list to insert ntimestep

  int i;
  for (i = ntime-1; i >= 0; i--) {
    if (ntimestep == tlist[i]) return;
    if (ntimestep < tlist[i]) break;
  }
  i++;

  // extend list as needed

  if (ntime == maxtime) {
    maxtime += DELTA;
    memory->grow(tlist, maxtime, "compute:tlist");
  }

  // move remainder of list upward and insert ntimestep

  for (int j = ntime-1; j >= i; j--) tlist[j+1] = tlist[j];
  tlist[i] = ntimestep;
  ntime++;

}

/*  ----------------------------------------------------------------------
   return 1/0 if ntimestep is or is not in list of calling timesteps
   if value(s) on top of list are less than ntimestep, delete them
   search from top downward, since list of times is in decreasing order
-------------------------------------------------------------------------  */

int Compute::matchstep(bigint ntimestep)
{
  for (int i = ntime-1; i >= 0; i--) {
    if (ntimestep < tlist[i]) return 0;
    if (ntimestep == tlist[i]) return 1;
    if (ntimestep > tlist[i]) ntime--;
  }
  return 0;
}
/*  ----------------------------------------------------------------------
   clean out list of timesteps to call the compute on
-------------------------------------------------------------------------  */

void Compute::clearstep()
{
  ntime = 0;
}

/*  ----------------------------------------------------------------------  */

void Compute::reset_extra_compute_fix(const char *)
{
  error->all(FLERR, 
             "Compute does not allow an extra compute or fix to be reset");
}


