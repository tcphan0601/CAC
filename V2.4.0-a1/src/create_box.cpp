#include <cstdlib>
#include <cstring>
#include "create_box.h"
#include "atom.h"
#include "atom_vec.h"
#include "universe.h"
#include "domain.h"
#include "region.h"
#include "region_prism.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal create_box command");

  if (domain->box_exist)
    error->all(FLERR,"Cannot create_box after simulation box is defined");
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

  domain->box_exist = 1;

  // region check

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Create_box region ID does not exist");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Create_box region does not support a bounding box");

  domain->regions[iregion]->init();

  // if region not prism:
  //   setup orthogonal domain
  //   set simulation domain from region extent
  // if region is prism:
  //   seutp triclinic domain
  //   set simulation domain params from prism params

  if (strcmp(domain->regions[iregion]->style,"prism") != 0) {
    domain->triclinic = 0;
    domain->boxlo[0] = domain->regions[iregion]->extent_xlo;
    domain->boxhi[0] = domain->regions[iregion]->extent_xhi;
    domain->boxlo[1] = domain->regions[iregion]->extent_ylo;
    domain->boxhi[1] = domain->regions[iregion]->extent_yhi;
    domain->boxlo[2] = domain->regions[iregion]->extent_zlo;
    domain->boxhi[2] = domain->regions[iregion]->extent_zhi;

  } else {
    domain->triclinic = 1;
    RegPrism *region = (RegPrism *) domain->regions[iregion];
    domain->boxlo[0] = region->xlo;
    domain->boxhi[0] = region->xhi;
    domain->boxlo[1] = region->ylo;
    domain->boxhi[1] = region->yhi;
    domain->boxlo[2] = region->zlo;
    domain->boxhi[2] = region->zhi;
    domain->xy = region->xy;
    domain->xz = region->xz;
    domain->yz = region->yz;
  }

  // set atom type quantities

  atom->ntypes = universe->inumeric(FLERR,arg[0]);

  // problem setup using info from header
  // deallocate/grow insures any extra settings are used for topology arrays
  // necessary in case no create_atoms is performed

  update->ntimestep = 0;

  atom->allocate_type_arrays();
  atom->avec->grow(1);

  domain->print_box("Created ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();
}
