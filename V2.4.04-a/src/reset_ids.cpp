#include "reset_ids.h"
#include <mpi.h>
#include "atom.h"
#include "element.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ResetIDs::ResetIDs(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void ResetIDs::command(int narg, char ** /* arg */)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_ids command before simulation box is defined");
  if (narg != 0) error->all(FLERR,"Illegal reset_ids command");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_ids unless atoms have IDs");
  if (element->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_ids unless elements have IDs");

  // NOTE: check if any fixes exist which store atom IDs?
  // if so, this operation will mess up the fix

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Resetting atom/element IDs ...\n");
    if (logfile) fprintf(logfile,"Resetting atom/element IDs ...\n");
  }

  // reset IDs to 0 for all atoms/elements 
  // assign IDs through tag_extends()
  //
  for (int i = 0; i < atom->nlocal; i++)
    atom->tag[i] = 0;
  for (int i = 0; i < element->nlocal; i++)
    element->tag[i] = 0;

  atom->tag_extend();
  element->tag_extend();

}
