#include "reset_ids.h"
#include <stdlib.h>
#include <string.h>
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

void ResetIDs::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_ids command before simulation box is defined");
  if (narg != 1) error->all(FLERR,"Illegal reset_ids command");
  int atom_flag = 0;
  int element_flag = 0;
  if (strcmp(arg[0],"atom") == 0) atom_flag = 1;
  else if (strcmp(arg[0],"element") == 0) element_flag = 1;
  else if (strcmp(arg[0],"all") == 0) 
    atom_flag = element_flag = 1;
  if (atom_flag && atom->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_ids unless atoms have IDs");
  if (element_flag && element->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_ids unless elements have IDs");

  // NOTE: check if any fixes exist which store atom IDs?
  // if so, this operation will mess up the fix

  if (comm->me == 0) {
    if (atom_flag && element_flag) {
      if (screen) fprintf(screen,"Resetting atom and element IDs ...\n");
      if (logfile) fprintf(logfile,"Resetting atom and element IDs ...\n");
    } else if (atom_flag) {
      if (screen) fprintf(screen,"Resetting atom IDs ...\n");
      if (logfile) fprintf(logfile,"Resetting atom IDs ...\n");
    } else if (element_flag) {
      if (screen) fprintf(screen,"Resetting element IDs ...\n");
      if (logfile) fprintf(logfile,"Resetting element IDs ...\n");

    }
  }

  // reset IDs to 0 for all atoms/elements 
  // assign IDs through tag_extends()

  if (atom_flag) {
    for (int i = 0; i < atom->nlocal; i++)
      atom->tag[i] = 0;
    atom->tag_extend();
  }

  if (element_flag) {
    for (int i = 0; i < element->nlocal; i++)
      element->tag[i] = 0;
    element->tag_extend();
  }

}
