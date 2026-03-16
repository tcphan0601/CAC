#include <mpi.h>
#include <string.h>
#include "disc_elements.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "irregular_comm.h"
#include "universe.h"

#include <map>

using namespace CAC_NS;

enum{ATOM,ELEMENT};
enum{HEX2HEX,HEX2WEDGE,HEX2PYRTET,WEDGE2WEDGE,WEDGE2PYRTET,WEDGE2TET,PYR2TET,HEX2HEXQUAD};

/* ---------------------------------------------------------------------- */

DiscElements::DiscElements(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
}

/* ----------------------------------------------------------------------
   Called as disc_elements command in input script.
   New atoms created will have ID starting immediately 
   following the largest atom ID existing before the
   disc_elements command was evoked.
   Call IrregularComm to migrate atoms since atoms might move far away
   ------------------------------------------------------------------------- */

void DiscElements::command(int narg, char **arg)
{
  disc_style = ATOM;
  if (domain->box_exist == 0)
    error->all(FLERR,"Disc_elements command before simulation box is defined");

  if (element->nelements == 0) {
    if (me == 0) error->warning(FLERR,"No element in simulation box");
    return;
  }

  if (narg < 1) error->all(FLERR,"Illegal disc_elements command");

  // flag elements for discretizing

  if (strcmp(arg[0],"region") == 0) list_disc_region(narg,arg);
  else if (strcmp(arg[0],"group") == 0) list_disc_group(narg,arg);
  else error->all(FLERR,"Illegal disc_elements command");

  element->init(); 

  ElementVec *evec = element->evec;

  // discretize elements in list

  if (disc_style == ATOM) {
    atom->init();
    int nnewatoms = 0;
    int nlocal_previous = atom->nlocal;
    int ndel = 0;
    int i = 0;
    while (i < element->nlocal)
      if (disc_flag_list[i]) {
        ndel++;
        nnewatoms += evec->element2atom(i);
        disc_flag_list[i] = disc_flag_list[element->nlocal];
      } else i++;

    // init per-atom fix/compute/variable values for created atoms

    atom->data_fix_compute_dump_variable(nlocal_previous,atom->nlocal);

    // assign tag to new atoms

    if (atom->tag_enable) atom->tag_extend();
    atom->tag_check();

    // update new global natoms and nelements

    int totalnewatoms;
    int totaldel;
    MPI_Allreduce(&nnewatoms,&totalnewatoms,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&ndel,&totaldel,1,MPI_INT,MPI_SUM,world);
    atom->natoms += totalnewatoms;
    element->nelements -= totaldel;
    element->nintpls -= totalnewatoms;

    FILE *out;
    if (me == 0) 
      for (int m = 0; m < 2; m++) {
        if (m == 0) out = logfile; 
        else out = screen;
        if (out) {
          fprintf(out,"Discretized %d elements into %d atoms\n",totaldel,totalnewatoms);
          fprintf(out,"  New total atoms = " BIGINT_FORMAT "\n",atom->natoms);
          fprintf(out,"  New total elements = " BIGINT_FORMAT "\n",element->nelements);
        }
      }

  } else {
    int nnewelems = 0;
    int nlocal_previous = element->nlocal;
    int ndel = 0;

    // old elements have the same local index after discretize so just use for loop

    // first loop scan for new element etypes
    // second loop split elements
    // third loop delete old elements

    for (int i = 0; i < nlocal_previous; i++)
      if (disc_flag_list[i]) {
        if (split_style == HEX2HEX)
          evec->hex2hex(i,split_cell,split_dim,1);
        else if (split_style == HEX2HEXQUAD)
          evec->hex2hexquad(i,split_cell,split_dim,1);
        else if (split_style == HEX2WEDGE)
          evec->hex2wedge(i,split_axis,split_dir,larger_side,1);
        else if (split_style == HEX2PYRTET)
          evec->hex2pyrtet(i,split_axis,split_dir,larger_side,1);
        else if (split_style == WEDGE2WEDGE)
          error->all(FLERR,"Split scheme not implemented yet");
        else if (split_style == WEDGE2PYRTET)
          evec->wedge2pyrtet(i,split_axis,split_dir,1);
        else if (split_style == WEDGE2TET)
          evec->wedge2tet(i,split_axis,split_dir,1);
        else if (split_style == PYR2TET)
          evec->pyr2tet(i,split_dir,larger_side,1);
      }

    element->evec->add_requested_etype(0);
    for (int i = 0; i < nlocal_previous; i++)
      if (disc_flag_list[i]) {
        ndel++;
        if (split_style == HEX2HEX)
          evec->hex2hex(i,split_cell,split_dim,0);
        else if (split_style == HEX2HEXQUAD)
          evec->hex2hexquad(i,split_cell,split_dim,0);
        else if (split_style == HEX2WEDGE)
          evec->hex2wedge(i,split_axis,split_dir,larger_side,0);
        else if (split_style == HEX2PYRTET)
          evec->hex2pyrtet(i,split_axis,split_dir,larger_side,0);
        else if (split_style == WEDGE2WEDGE)
          error->all(FLERR,"Split scheme not implemented yet");
        else if (split_style == WEDGE2PYRTET)
          evec->wedge2pyrtet(i,split_axis,split_dir,0);
        else if (split_style == WEDGE2TET)
          evec->wedge2tet(i,split_axis,split_dir,0);
        else if (split_style == PYR2TET)
          evec->pyr2tet(i,split_dir,larger_side,0);
      }

    // init per-atom fix/compute/variable values for created atoms

    element->data_fix_compute_dump_variable(nlocal_previous,element->nlocal);

    for (int i = 0; i < nlocal_previous; i++)
      if (disc_flag_list[i]) {
        evec->copy(element->nlocal-1,i,1);
        element->nlocal--;
      }

    // assign tags for created elements

    if (element->tag_enable) element->tag_extend();
    element->tag_check();

    // update new global nelements, global nnodes, global nintpls if needed

    int totaldel;
    int increase_factor;
    MPI_Allreduce(&ndel,&totaldel,1,MPI_INT,MPI_SUM,world);
    if (split_style == HEX2HEX) {
      element->nelements += totaldel;
      element->nnodes += totaldel*8;
      increase_factor = 2;
    } else if (split_style == HEX2HEXQUAD) {
      element->nelements += totaldel;
      element->nnodes += totaldel*4;
      increase_factor = 2;
    } else if (split_style == HEX2WEDGE) {
      element->nelements += totaldel;
      element->nnodes += totaldel*4;
      increase_factor = 2;
    } else if (split_style == HEX2PYRTET) {
      element->nelements += totaldel*3;
      element->nnodes += totaldel*10;
      increase_factor = 4;
    } else if (split_style == WEDGE2WEDGE) {
      element->nelements += totaldel;
      element->nnodes += totaldel*6;
      increase_factor = 2;
    } else if (split_style == WEDGE2PYRTET) {
      element->nelements += totaldel;
      element->nnodes += totaldel*3;
      increase_factor = 2;
    } else if (split_style == WEDGE2TET) {
      element->nnodes -= totaldel*2;
      element->nintpls = element->count_intpl();
      increase_factor = 1;
    } else if (split_style == PYR2TET) {
      element->nelements += totaldel;
      element->nnodes += totaldel*3;
      increase_factor = 2;
    }

    FILE *out;
    if (me == 0) 
      for (int m = 0; m < 2; m++) {
        if (m == 0) out = logfile; 
        else out = screen;
        if (out) {
          fprintf(out,"Discretized %d elements into %d smaller elements\n",totaldel,totaldel*increase_factor);
          fprintf(out,"  New total elements = " BIGINT_FORMAT "\n",element->nelements);
        }
      }
  }

  // if compress flag set
  // reset element tags to be continuous
  // set all element IDs to 0, call tag_extend()

  if (compress_flag) {
    tagint *tag = element->tag;
    int nlocal = element->nlocal;
    for (int i = 0; i < nlocal; i++) tag[i] = 0;
    element->tag_extend();
  }

  // Invoke new irregular comm, migrate atoms/elements and destroy it
  // must be done after updating new global natoms and nelements

  double **ax = atom->x;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

  double **ex = element->x;
  double ***nodex = element->nodex;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i],npe[etype[i]]);

  if (migrate_flag) {
    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal,atom->x);
      domain->x2lamda(element->nlocal,element->x);
      domain->nodex2lamda(element->nlocal,element->nodex);
    }
    domain->reset_box();

    IrregularComm *irrcomm = new IrregularComm(cac);
    irrcomm->migrate(1);
    delete irrcomm;

    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal,atom->x);
      domain->lamda2x(element->nlocal,element->x);
      domain->lamda2nodex(element->nlocal,element->nodex);
    }
  }
  // reset element->map and atom->map if it exists
  // set nghost to 0 so old ghosts of deleted elements and new atoms won't be mapped

  if (element->map_style) {
    element->nghost = 0;
    element->map_init();
    element->map_set();
  }

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // clean up

  memory->destroy(disc_flag_list);
}

/* ----------------------------------------------------------------------
   flag all elements in group to be discretize, group will still exist
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_group(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find disc_elements group ID");

  options(narg-2,&arg[2]);

  // allocate discretization list

  int nlocal = element->nlocal;
  memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");


  int *mask = element->mask;
  int groupbit = group->bitmask[igroup];
  int *element_shape_ids = element->element_shape_ids;
  int *etype = element->etype;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (disc_style == ATOM) 
        disc_flag_list[i] = 1;
      else if (element_shape_ids[etype[i]] == origin_shape_id) 
        disc_flag_list[i] = 1;
      else disc_flag_list[i] = 0;
    } else disc_flag_list[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   flag all elements in region to discretize
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

  int iregion = domain->find_region(arg[1]);
  if (iregion == -1) error->all(FLERR,"Could not find disc_elements region ID");
  domain->regions[iregion]->prematch();

  options(narg-2,&arg[2]);

  // allocate discretization list

  int nlocal = element->nlocal;
  memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");

  double **x = element->x;
  int *apc = element->apc;
  int *element_shape_ids = element->element_shape_ids;
  int *etype = element->etype;
  //  if (apc == 1) {
  for (int i = 0; i < nlocal; i++) {
    if (domain->regions[iregion]->match(x[i])) {
      if (disc_style == ATOM) 
        disc_flag_list[i] = 1;
      else if (element_shape_ids[etype[i]] == Element::HEXAHEDRON) 
        disc_flag_list[i] = 1;
      else disc_flag_list[i] = 0;
    } else disc_flag_list[i] = 0;
  }
  //  } else {
  //    int mapflag = 0;
  //    if (element->map_style == 0) {
  //      mapflag = 1;
  //      element->map_init();
  //      element->map_set();
  //    }
  //
  //    tagint itag;
  //    tagint *tag = element->tag; 
  //    for (int i = 0; i < nlocal; i++) {
  //      disc_flag_list[i] = 0;
  //      itag = (tag[i]-1)/apc*apc+1;
  //      for (int j = 0; j < apc; j++) {
  //        int ii = element->map(itag+j);
  //        if (ii < 0) error->one(FLERR,"Ghost element not correct for link");
  //        if (domain->regions[iregion]->match(x[ii])) 
  //          disc_flag_list[i] = 1;
  //      }
  //    }
  //
  //    if (mapflag) {
  //      element->map_delete();
  //      element->map_style = 0;
  //    }
  //
  //  }
}

/* ----------------------------------------------------------------------
   process command options
   ------------------------------------------------------------------------- */

void DiscElements::options(int narg, char **arg)
{
  compress_flag = 1;
  migrate_flag = 1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
      if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
      else error->all(FLERR,"Illegal disc_elements command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"migrate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
      if (strcmp(arg[iarg+1],"yes") == 0) migrate_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) migrate_flag = 0;
      else error->all(FLERR,"Illegal disc_elements command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"style") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
      if (strcmp(arg[iarg+1],"atom") == 0) {
        disc_style = ATOM;
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"element") == 0) {
        disc_style = ELEMENT;
        if (iarg+3 > narg) error->all(FLERR,"Illegal disc_elements command");
        origin_shape_id = element->find_element_shape_id(arg[iarg+2]);
        if (origin_shape_id < 0)
          error->all(FLERR,"Illegal origin element shape in disc_elements command");
        if (origin_shape_id == Element::HEXAHEDRON) {
          if (iarg+4 > narg) error->all(FLERR,"Illegal disc_elements command");
          if (strcmp(arg[iarg+3],"Hex") == 0) {
            split_style = HEX2HEX;
            if (iarg+6 > narg) error->all(FLERR,"Illegal disc_elements command");
            if (strcmp(arg[iarg+4],"x") == 0) split_dim = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_dim = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_dim = 2;
            else error->all(FLERR,"Illegal disc_elements command");
            split_cell = universe->inumeric(FLERR,arg[iarg+5]);
            iarg += 6;
          } else if (strcmp(arg[iarg+3],"HexQuad") == 0) {
            split_style = HEX2HEXQUAD;
            if (iarg+6 > narg) error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+4],"x") == 0) split_dim = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_dim = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_dim = 2;
            else error->all(FLERR,"Illegal disc_elements command");
            
            if (strcmp(arg[iarg+5],"lower") == 0) split_cell = 0;
            else if (strcmp(arg[iarg+5],"upper") == 0) split_cell = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            iarg += 6;
          } else if (strcmp(arg[iarg+3],"Wedge") == 0) {
            split_style = HEX2WEDGE;
            if (iarg+7 > narg) error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+4],"x") == 0) split_axis = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_axis = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_axis = 2;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+5],"-") == 0) split_dir = 0;
            else if (strcmp(arg[iarg+5],"+") == 0) split_dir = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+6],"lower") == 0) larger_side = 0;
            else if (strcmp(arg[iarg+6],"upper") == 0) larger_side = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            iarg += 7;
          } else if (strcmp(arg[iarg+3],"PyrTet") == 0) {
            split_style = HEX2PYRTET;
            if (iarg+7 > narg) error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+4],"x") == 0) split_axis = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_axis = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_axis = 2;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+5],"1") == 0) split_dir = 0;
            else if (strcmp(arg[iarg+5],"2") == 0) split_dir = 1;
            else if (strcmp(arg[iarg+5],"3") == 0) split_dir = 2;
            else if (strcmp(arg[iarg+5],"4") == 0) split_dir = 3;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+6],"lower") == 0) larger_side = 0;
            else if (strcmp(arg[iarg+6],"upper") == 0) larger_side = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            iarg += 7;

          } else {
            error->all(FLERR,"Unsupported target element shape");
          }
        } else if (origin_shape_id == Element::WEDGE) {
          if (iarg+4 > narg) error->all(FLERR,"Illegal disc_elements command");
          if (strcmp(arg[iarg+3],"Wedge") == 0) {
            split_style = WEDGE2WEDGE;
            if (iarg+5 > narg) error->all(FLERR,"Illegal disc_elements command");
            split_cell = universe->inumeric(FLERR,arg[iarg+4]);
            iarg += 5;
          } else if (strcmp(arg[iarg+3],"PyrTet") == 0) {
            split_style = WEDGE2PYRTET;
            if (iarg+6 > narg) error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+4],"x") == 0) split_axis = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_axis = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_axis = 2;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+5],"-") == 0) split_dir = 0;
            else if (strcmp(arg[iarg+5],"+") == 0) split_dir = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            iarg += 6;
          } else if (strcmp(arg[iarg+3],"Tet") == 0) {
            split_style = WEDGE2TET;
            if (iarg+6 > narg) error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+4],"x") == 0) split_axis = 0;
            else if (strcmp(arg[iarg+4],"y") == 0) split_axis = 1;
            else if (strcmp(arg[iarg+4],"z") == 0) split_axis = 2;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+5],"-") == 0) split_dir = 0;
            else if (strcmp(arg[iarg+5],"+") == 0) split_dir = 1;
            else error->all(FLERR,"Illegal disc_elements command");

            iarg += 6;

          } else {
            error->all(FLERR,"Unsupported target element shape");
          }

        } else if (origin_shape_id == Element::PYRAMID) {
          if (iarg+4 > narg) error->all(FLERR,"Illegal disc_elements command");
          if (strcmp(arg[iarg+3],"Tet") == 0) {
            split_style = PYR2TET;
            if (iarg+6 > narg) error->all(FLERR,"Illegal disc_elements command");
            if (strcmp(arg[iarg+4],"-") == 0) split_dir = 0;
            else if (strcmp(arg[iarg+4],"+") == 0) split_dir = 1;
            else if (strcmp(arg[iarg+4],"auto") == 0) split_dir = 2;
            else error->all(FLERR,"Illegal disc_elements command");

            if (strcmp(arg[iarg+5],"lower") == 0) larger_side = 0;
            else if (strcmp(arg[iarg+5],"upper") == 0) larger_side = 1;
            else error->all(FLERR,"Illegal disc_elements command");
            iarg += 6;
          } else {
            error->all(FLERR,"Unsupported target element shape");
          }
        } else {
          error->all(FLERR,"Unsupported origin element shape");
        }
      } else error->all(FLERR,"Illegal disc_elements command");
    } else error->all(FLERR,"Illegal disc_elements command");
  } 
  // if more than one atom per node cluster, no compress element tags 

  if (element->element_cluster_flag) compress_flag = 0;
}

