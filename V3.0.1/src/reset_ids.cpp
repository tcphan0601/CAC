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

/*  ----------------------------------------------------------------------  */

ResetIDs::ResetIDs(CAC *cac) : Pointers(cac) {}

/*  ----------------------------------------------------------------------  */

void ResetIDs::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Reset_ids command before simulation box is defined");
  if (narg < 1) error->all(FLERR, "Illegal reset_ids command");
  int atom_flag = 0;
  int element_flag = 0;
  if (strcmp(arg[0], "atom") == 0) atom_flag = 1;
  else if (strcmp(arg[0], "element") == 0) element_flag = 1;
  else if (strcmp(arg[0], "all") == 0) 
    atom_flag = element_flag = 1;

  order_flag = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "order") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal reset_ids command");
      if (strcmp(arg[iarg+1], "yes") == 0) order_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) order_flag = 0;
      else error->all(FLERR, "Illegal reset_ids command");
      iarg += 2;
    } else error->all(FLERR, "Illegal reset_ids command");
  }

  if (order_flag && comm->nprocs != 1) error->all(FLERR, "Must use 1 proc only for retaining IDs order");
  if (atom_flag && atom->tag_enable == 0)
    error->all(FLERR, "Cannot use reset_ids unless atoms have IDs");
  if (element_flag && element->tag_enable == 0)
    error->all(FLERR, "Cannot use reset_ids unless elements have IDs");

  // NOTE: check if any fixes exist which store atom IDs?
  // if so, this operation will mess up the fix

  if (comm->me == 0) {
    if (atom_flag && element_flag) {
      if (screen) fprintf(screen, "Resetting atom and element IDs ...\n");
      if (logfile) fprintf(logfile, "Resetting atom and element IDs ...\n");
    } else if (atom_flag) {
      if (screen) fprintf(screen, "Resetting atom IDs ...\n");
      if (logfile) fprintf(logfile, "Resetting atom IDs ...\n");
    } else if (element_flag) {
      if (screen) fprintf(screen, "Resetting element IDs ...\n");
      if (logfile) fprintf(logfile, "Resetting element IDs ...\n");

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
    if (order_flag) {

      tagint *tag = element->tag;

      int *list_index = new int[element->nlocal];
      for (int i = 0; i < element->nlocal; i++)
        list_index[i] = i;
      
      // after this for loop, the list_index array contains the index of elements with tag corresponding to its index + 1 in list_index array
      
      for (int i = 0; i < element->nlocal-1; i++) {
        int min_tag = tag[list_index[i]];
        int min_index = i;

        // find index of list_index array with the element index of smallest tag
        
        for (int j = i+1; j < element->nlocal; j++) {
          if (tag[list_index[j]] < min_tag) {
            min_tag = tag[list_index[j]];
            min_index = j;
          }
        }

        // swap indices in list_index
        
        if (min_index != i) {
          int tmp = list_index[i];
          list_index[i] = list_index[min_index];
          list_index[min_index] = tmp;
        }
      } 
      for (int i = 0; i < element->nlocal; i++)
        tag[list_index[i]] = i + 1;
    } else {
      for (int i = 0; i < element->nlocal; i++)
        element->tag[i] = 0;
      element->tag_extend();
    }
  }

}

