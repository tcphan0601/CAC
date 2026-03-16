#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fill_surface.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "comm.h"
#include "memory.h"
#include "universe.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "group.h"
#include "variable.h"
#include "math_extra.h"
#include "math_const.h"
#include "error.h"

using namespace CAC_NS;
using namespace MathConst;
using namespace MathExtra;

#define EPSILON 0.1
#define DELTA 512

/* ---------------------------------------------------------------------- */

FillSurface::FillSurface(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void FillSurface::command(int narg, char **arg)
{
  
  if (domain->triclinic) error->all(FLERR,"Fill_surface command does not work with triclinic box yet");
  
  if (domain->box_exist == 0)
    error->all(FLERR,"Fill_surface command before simulation box is defined");
  if (domain->lattice->nbasis == 0)
    error->all(FLERR,"Cannot fill surface with undefined lattice");

  //if (modify->nfix_restart_peratom)
  //  error->all(FLERR,"Cannot create_atoms after "
  //             "reading restart file with per-atom info");

  // parse arguments

  if (narg < 2) error->all(FLERR,"Illegal fill_surface command");
  
  igroup = group->find(arg[0]);
  groupbit = group->bitmask[igroup];
  if (strcmp(arg[1],"-x") == 0) {
    direction = -1;
    dim = 0;
  } else if (strcmp(arg[1],"-y") == 0) {
    direction = -1;
    dim = 1;
  } else if (strcmp(arg[1],"-z") == 0) {
    direction = -1;
    dim = 2;
  } else if (strcmp(arg[1],"+x") == 0) {
    direction = 1;
    dim = 0;
  } else if (strcmp(arg[1],"+y") == 0) {
    direction = 1;
    dim = 1;
  } else if (strcmp(arg[1],"+z") == 0) {
    direction = 1;
    dim = 2;
  } else error->all(FLERR,"Illegal fill_surface command");

  int iarg = 2;
  width = EPSILON;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"width") == 0) {
      if (iarg+1 < narg) {
        width = universe->numeric(FLERR,arg[iarg+1]);
        if (width <= 0.0) error->all(FLERR,"Illegal fill_surface command");
      } else error->all(FLERR,"Illegal fill_surface command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fill_surface command");
  }

  double minmax[6];
  group->bounds(igroup,minmax);

  if (direction == -1) {
    lo = minmax[dim*2] - width;
    hi = minmax[dim*2] + width;
  } else {
    lo = minmax[dim*2+1] - width;
    hi = minmax[dim*2+1] + width;
  }

  char **newarg = new char*[10];

  char *idgroup = (char *) "FILL_SURFACE_GROUP";
  char *idgrouptemp = (char *) "FILL_SURFACE_GROUP_TEMP";
  char *idregion = (char *) "FILL_SURFACE_REGION";

  // create region to select element on surface

  newarg[0] = idregion;
  newarg[1] = (char *) "block";
  if (dim == 0) {
    newarg[2] = new char[256];
    newarg[3] = new char[256];
    sprintf(newarg[2],"%g",lo);
    sprintf(newarg[3],"%g",hi);
  } else {
    newarg[2] = (char *) "INF";
    newarg[3] = (char *) "INF";
  }
  if (dim == 1) {
    newarg[4] = new char[256];
    newarg[5] = new char[256];
    sprintf(newarg[4],"%g",lo);
    sprintf(newarg[5],"%g",hi);
  } else {
    newarg[4] = (char *) "INF";
    newarg[5] = (char *) "INF";
  }
  if (dim == 2) {
    newarg[6] = new char[256];
    newarg[7] = new char[256];
    sprintf(newarg[6],"%g",lo);
    sprintf(newarg[7],"%g",hi);
  } else {
    newarg[6] = (char *) "INF";
    newarg[7] = (char *) "INF";
  }

  newarg[8] = (char *) "units";
  newarg[9] = (char *) "box";

  domain->add_region(10,newarg);

  if (dim == 0) {
    delete [] newarg[2];
    delete [] newarg[3];
  } else if (dim == 1) {
    delete [] newarg[4];
    delete [] newarg[5];
  } else if (dim == 2) {
    delete [] newarg[6];
    delete [] newarg[7];
  }

  iregion = domain->find_region(idregion);
  domain->regions[iregion]->init();
  domain->regions[iregion]->prematch();

  // clear ghost count 
  // do it now b/c creating elements will overwrite ghost elements

  element->nghost = 0;

  bigint nelements_previous = element->nelements;
  int nlocal_previous = element->nlocal;

  double x[3];
  int i,ii,j;
  int nlocal = element->nlocal;
  int *element_shape_ids = element->element_shape_ids;
  double ***nodex = element->nodex;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int ncellx,ncelly,ncellz,nintgx;
  int ietype,ietype_new,ictype,shapeid;
  int type_info[8];
  double **nodecoord;
  memory->create(nodecoord,8,3,"fill_surface:nodecoord");

  max_surface_element = DELTA;
  memory->create(surface_element_list,DELTA,"fill_surface:surface_element_list");

  // assuming the elements in group are created from lattice using create_elements command
  // elements can be displaced, as long as the orientation stay the same as when they are created
  // otherwise, this will not work

  // first loop check for surface elements and new etypes

  num_surface_element = 0;
  int surface_flag;
  for (i = 0; i < nlocal; i++) {
    ictype = ctype[i];
    ietype = etype[i];
    if (element_shape_ids[ietype] != Element::HEXAHEDRON)
      continue;

    // if at least one node is in the region, add element to list
    // grow list if necessary

    surface_flag = 0;
    for (j = 0; j < 8; j++) 
      if (domain->regions[iregion]->match(nodex[i][j])) {
        surface_flag = 1;
        if (num_surface_element == max_surface_element) {
          max_surface_element += DELTA;
          memory->grow(surface_element_list,max_surface_element,"fill_surface:surface_element_list");
        }
        surface_element_list[num_surface_element++] = i;
        break;
      }

    if (surface_flag == 0) 
      continue;

    // add new etype if needed

    ncellx = element->evec->ncells[ietype][0];
    ncelly = element->evec->ncells[ietype][1];
    ncellz = element->evec->ncells[ietype][2];
    nintgx = element->evec->nintgs[ietype][0];

    type_info[0] = ncellx;
    type_info[1] = ncellz;
    type_info[2] = 0;
    type_info[3] = nintgx;
    if (ncellz >= 4) type_info[4] = 4;
    else type_info[4] = 3;
    type_info[5] = 0;
    type_info[6] = Element::WEDGE;
    type_info[7] = element->apc[ietype];

    element->evec->request_new_etype(type_info);

  }

  printf("num_surface_element = %d\n",num_surface_element);
  element->evec->add_requested_etype(0);

  for (int ii = 0; ii < num_surface_element; ii++) {
    i = surface_element_list[ii];
    ietype = element->etype[i];
    ictype = element->ctype[i];

    // find etype of new elements

    ncellx = element->evec->ncells[ietype][0];
    ncelly = element->evec->ncells[ietype][1];
    ncellz = element->evec->ncells[ietype][2];
    nintgx = element->evec->nintgs[ietype][0];

    type_info[0] = ncellx;
    type_info[1] = ncellz;
    type_info[2] = 0;
    type_info[3] = nintgx;
    if (ncellz >= 4) type_info[4] = 4;
    else type_info[4] = 3;
    type_info[5] = 0;
    type_info[6] = Element::WEDGE;
    type_info[7] = element->apc[ietype];

    ietype_new = element->evec->find_etype(type_info);

    copy3(element->x[i],x);

    domain->lattice->box2lattice(x[0],x[1],x[2]);

    if (domain->lattice->style == Lattice::PRIBCC) {
      if (dim == 0) {

      }

      // fill in y direction

      else if (dim == 1) {

        // first element

        nodecoord[0][0] = x[0] - ((double) ncellx - 1.0)/2.0;
        nodecoord[0][1] = x[1] + direction*((double) ncelly - 1.0)/2.0 + direction;
        nodecoord[0][2] = x[2] + direction*((double) ncellz - 1.0)/2.0;

        nodecoord[1][0] = nodecoord[0][0];
        nodecoord[1][1] = nodecoord[0][1];
        nodecoord[1][2] = nodecoord[0][2] - direction*(ncellz - 1);

        nodecoord[2][0] = nodecoord[0][0];
        nodecoord[2][1] = nodecoord[0][1] + direction*(ncellz - 1);
        nodecoord[2][2] = nodecoord[0][2] - direction*(2*ncellz - 2);

        nodecoord[3][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[3][1] = nodecoord[0][1];
        nodecoord[3][2] = nodecoord[0][2];

        nodecoord[4][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[4][1] = nodecoord[0][1];
        nodecoord[4][2] = nodecoord[0][2] - direction*(ncellz - 1);

        nodecoord[5][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[5][1] = nodecoord[0][1] + direction*(ncellz - 1);
        nodecoord[5][2] = nodecoord[0][2] - direction*(2*ncellz - 2);

        // convert all nodecoords from lattice coord to box coord
        // create new element and add it to user-specified group

        for (int inode = 0; inode < 6; inode++) 
          domain->lattice->lattice2box(nodecoord[inode][0],
              nodecoord[inode][1],
              nodecoord[inode][2]);
        element->evec->create_element(nodecoord,ietype_new,ictype,0);

        // second element

        nodecoord[0][0] = x[0] - ((double) ncellx - 1.0)/2.0;
        nodecoord[0][1] = x[1] + direction*((double) ncelly - 1.0)/2.0 + direction;
        nodecoord[0][2] = x[2] - direction*((double) ncellz + 1.0)/2.0;

        nodecoord[1][0] = nodecoord[0][0];
        nodecoord[1][1] = nodecoord[0][1];
        nodecoord[1][2] = nodecoord[0][2] - direction*(ncellz - 1);

        nodecoord[2][0] = nodecoord[0][0];
        nodecoord[2][1] = nodecoord[0][1] + direction*(ncellz - 1);
        nodecoord[2][2] = nodecoord[0][2] - direction*(ncellz - 1);

        nodecoord[3][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[3][1] = nodecoord[0][1];
        nodecoord[3][2] = nodecoord[0][2];

        nodecoord[4][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[4][1] = nodecoord[0][1];
        nodecoord[4][2] = nodecoord[0][2] - direction*(ncellz - 1);

        nodecoord[5][0] = nodecoord[0][0] + ncellx - 1.0;
        nodecoord[5][1] = nodecoord[0][1] + direction*(ncellz - 1);
        nodecoord[5][2] = nodecoord[0][2] - direction*(ncellz - 1);

        // convert all nodecoords from lattice coord to box coord
        // create new element and add it to user-specified group

        for (int inode = 0; inode < 6; inode++) 
          domain->lattice->lattice2box(nodecoord[inode][0],
              nodecoord[inode][1],
              nodecoord[inode][2]);
        element->evec->create_element(nodecoord,ietype_new,ictype,0);

      } 

    } else if (domain->lattice->style == Lattice::PRIFCC) {

    }
  }

  // init per-element fix/compute/dump/variable values for created atoms

  element->data_fix_compute_dump_variable(nlocal_previous,element->nlocal);

  // set new total # of elements and error check

  bigint nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&element->nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (element->nelements < 0 || element->nelements >= MAXBIGINT)
    error->all(FLERR,"Too many total elements");

  bigint nelements_added = element->nelements - nelements_previous;

  bigint nintpls_added_local,nintpls_added;
  nintpls_added_local = 0;

  etype = element->etype;
  int *nintpl = element->nintpl;
  int *mask = element->mask;
  for (int i = nlocal_previous; i < nblocal; i++) {
    mask[i] |= groupbit;
    nintpls_added_local += nintpl[etype[i]];
  }
  MPI_Allreduce(&nintpls_added_local,&nintpls_added,1,MPI_CAC_BIGINT,MPI_SUM,world);

  element->nintpls += nintpls_added;

  // add IDs for newly created elements
  // check that element IDs are valid

  if (element->tag_enable) element->tag_extend();
  element->tag_check();

  // assign element ids to clusters 
  // NOTE: Need to fix!!! This is not correct

  if (element->element_cluster_flag) {
    error->all(FLERR,"Fill_surface command does not work with element cluster yet");
    //tagint **element_clusters = element->element_clusters;
    //tagint *tag = element->tag;
    //for (int i = nlocal_previous; i < element->nlocal; i += nbasis) {
    //  for (int j = 0; j < nbasis; j++) {
    //    for (int k = 0; k < nbasis; k++)
    //      element_clusters[i+j][k] = tag[i+k];
    //    for (int k = nbasis; k < element->max_apc; k++)
    //      element_clusters[i+j][k] = 0;
    //  }
    //}
  }

  // if global map exists, reset it
  // invoke map_init() b/c element count has grown
  // set nghost to 0 so old ghosts won't be mapped

  if (element->map_style) {
    element->map_init();
    element->map_set();
  }

  FILE *out;
  if (comm->me == 0) {
    for (int i = 0; i < 2; i++) {
      if (i == 0) out = screen;
      else out = logfile;
      if (out) {
        fprintf(out,"Created " BIGINT_FORMAT " elements\n",nelements_added);
        fprintf(out,"  Total " BIGINT_FORMAT " elements\n",element->nelements);
      }
    }
  }

  // clean up, delete region

  newarg[0] = idregion;
  newarg[1] = (char *) "delete";
  domain->delete_region(2,newarg);
  delete [] newarg;

  memory->destroy(nodecoord);
  memory->destroy(surface_element_list);

}

