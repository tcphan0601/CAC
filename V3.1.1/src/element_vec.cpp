#include <string.h>
#include <stdlib.h>
#include "element_vec.h"
#include "element.h"
#include "memory.h"
#include "domain.h"
#include "error.h"
#include "comm.h"
#include "universe.h"

using namespace CAC_NS;

#define BIGBIGDELTA 16384       // for nmaxucell increment
#define BIGDELTA 4096           // for nmaxgcell increment
#define DELTA 1024              // for nmax increment

/* ---------------------------------------------------------------------- */

ElementVec::ElementVec(CAC *cac) : Pointers(cac)
{
  nmax = nmaxnode = nmaxgcell = nmaxucell = 0;
  size_data_bonus = 0;
  nargcopy = 0;
  argcopy = nullptr;
  element_type_setflag = nullptr;
  nrequested_etype = 0;
  maxrequested_etype = 4;
  memory->create(requested_etype, 4 * 8, "evec:requested_etype");
}

/* ----------------------------------------------------------------------- */

ElementVec::~ElementVec()
{
  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;
  memory->destroy(element_type_setflag);
  memory->destroy(ncells);
  memory->destroy(ngcells);
  memory->destroy(requested_etype);
}

/*  ----------------------------------------------------------------------
     make copy of args for use by restart & replicate
    ------------------------------------------------------------------------- */

void ElementVec::store_args(int narg, char **arg)
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
    set maxnpe and maxapc
  -------------------------------------------------------------------------  */

void ElementVec::process_args(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal element_style command");
  element->maxnpe = universe->inumeric(FLERR, arg[0]);
  element->maxapc = universe->inumeric(FLERR, arg[1]);
  set_data_size();

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "wscale") == 0) {
      if (iarg+5 > narg) error->all(FLERR, "Invalid element_style command");
      element->weight_scale[0] = universe->numeric(FLERR, arg[iarg+1]);
      element->weight_scale[1] = universe->numeric(FLERR, arg[iarg+2]);
      element->weight_scale[2] = universe->numeric(FLERR, arg[iarg+3]);
      element->weight_scale[3] = universe->numeric(FLERR, arg[iarg+4]);
      if (element->weight_scale[0] <= 0 || 
          element->weight_scale[1] <= 0 ||
          element->weight_scale[2] <= 0 ||
          element->weight_scale[3] <= 0)
        error->all(FLERR, "Weight scale must be positive");
      element->weight_scale_flag = MAX(element->weight_scale_flag, Element::UNIFORM);
      iarg += 5;
    } else if (strcmp(arg[iarg], "wstyle") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Invalid element_style command");
      if (strcmp(arg[iarg+1], "uniform") == 0) 
        element->weight_scale_flag = Element::UNIFORM;
      else if (strcmp(arg[iarg+1], "nonuniform") == 0) 
        element->weight_scale_flag = Element::NONUNIFORM;
      else error->all(FLERR, "Invalid element_style command");
      iarg += 2;
    } else error->all(FLERR, "Invalid element_style command");
  }
  if (element->weight_scale_flag && element->weight_scale[0] < 0)
    error->all(FLERR, "Weight scale has not been set");
}

/*  ----------------------------------------------------------------------
  grow nmax so it is a multiple of DELTA
 -------------------------------------------------------------------------  */

void ElementVec::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
  element->nmax = nmax;
}

/*  ----------------------------------------------------------------------  */

void ElementVec::init()
{

  // check if every element has defined interpolate and integration

  for (int i = 0; i < element->nlocal; i++) {
    if (!element_type_setflag[element->etype[i]]) 
      error->all(FLERR, "Element type have not been set for all elements");
  }

}

/*  ----------------------------------------------------------------------
   find existing etype
   return the index of etype if match
   return 0 if no etype found
   -------------------------------------------------------------------------  */

int ElementVec::find_etype(int *type_info)
{
  for (int i = 1; i <= element->netypes; i++)
    if (type_info[0] == ncells[i][0] && 
        type_info[1] == ncells[i][1] && 
        type_info[2] == ncells[i][2] && 
        type_info[3] == ngcells[i][0] && 
        type_info[4] == ngcells[i][1] && 
        type_info[5] == ngcells[i][2] &&
        type_info[6] == element_shape_ids[i] &&
        type_info[7] == apc[i])
      return i;
  return 0;
}

/*  ----------------------------------------------------------------------
   request new etype
   return if etype already exist or requested, 
   otherwise add to list of requested etype
   -------------------------------------------------------------------------  */

void ElementVec::request_new_etype(int *type_info)
{
  if (find_etype(type_info)) return;

  int n = nrequested_etype;

  // check if requested etype exists

  for (int i = 0; i < n; i++) 
    if (requested_etype[i * 8] == type_info[0] &&
        requested_etype[i * 8+1] == type_info[1] &&
        requested_etype[i * 8+2] == type_info[2] &&
        requested_etype[i * 8+3] == type_info[3] &&
        requested_etype[i * 8+4] == type_info[4] &&
        requested_etype[i * 8+5] == type_info[5] &&
        requested_etype[i * 8+6] == type_info[6] &&
        requested_etype[i * 8+7] == type_info[7])
      return;

  // grow requested_etype array if needed

  if (n == maxrequested_etype) {
    maxrequested_etype *= 2;
    memory->grow(requested_etype, maxrequested_etype * 8, "evec:requested_etype");
  }

  // add requested etype to the list

  for (int i = 0; i < 8; i++)
    requested_etype[n * 8+i] = type_info[i];

  nrequested_etype++;

}

/*  ----------------------------------------------------------------------
   union all requested etype from all procs
   add the requested etype
Note: this function must be called by all procs
-------------------------------------------------------------------------  */

void ElementVec::add_requested_etype(int flag)
{
  int nprocs = comm->nprocs;
  int *recvcount = new int[nprocs];
  int *recvbuf = new int[nprocs];
  int *displs = new int[nprocs];
  for (int i = 0; i < nprocs; i++) {
    recvcount[i] = 1;
    displs[i] = i;
  }

  MPI_Allgatherv(&nrequested_etype, 1, MPI_INT, recvbuf, recvcount, displs, MPI_INT, world);

  int sum = 0;
  for (int i = 0; i < nprocs; i++) {
    recvcount[i] = recvbuf[i] * 8;
    displs[i] = sum;
    sum += recvcount[i];
  }

  int nrequests_all = sum/8;
  int *requests_all = new int[sum];

  MPI_Allgatherv(requested_etype, nrequested_etype * 8, MPI_INT, 
      requests_all, recvcount, displs, MPI_INT, world);
  char *new_type_string = new char[256];
  for (int i = 0; i < nrequests_all; i++) {
    int *type_info = requests_all + i * 8;
    if (!find_etype(type_info)) {
      int offset = 0;
      offset += sprintf(&new_type_string[offset], "%d ", element->netypes+1);
      offset += sprintf(&new_type_string[offset], "%s ", element->element_shape_list[type_info[6]]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[7]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[0]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[1]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[2]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[3]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[4]);
      offset += sprintf(&new_type_string[offset], "%d ", type_info[5]);

      grow_etype_arrays(element->netypes + 1);
      set_element_types(new_type_string, 0);
      if (comm->me == 0 && flag == 0) {
        if (screen) fprintf(screen, "Add new etype # %d %s %d %d %d %d %d %d %d\n", 
            element->netypes, element->element_shape_list[type_info[6]], 
            type_info[7], type_info[0], type_info[1], type_info[2], type_info[3], type_info[4], type_info[5]);
        if (logfile) fprintf(logfile, "Add new etype # %d %s %d %d %d %d %d %d %d\n", 
            element->netypes, element->element_shape_list[type_info[6]], 
            type_info[7], type_info[0], type_info[1], type_info[2], type_info[3], type_info[4], type_info[5]);
      }

    }
  }

  // reset number of requested etype

  nrequested_etype = 0;

  // clean up

  delete [] new_type_string; 
  delete [] requests_all;
  delete [] displs;
  delete [] recvbuf;
  delete [] recvcount;
}


void ElementVec::visual_nodex(double *coord, int i, int iapc, int inode)
{
  /*
  int ietype = etype[i];
  int shape_id = element_shape_ids[ietype];
  double shape_column[8];
  double px, py, pz;
  double unitx, unity, unitz;
  double nodal_visual_position = element->nodal_visual_position;
  coord[0] = coord[1] = coord[2] = 0.0;
  if (shape_id == Element::QUADRILATERAL) {
    // natural coordinate of each node
    double d[4][2] = {
      {-1, -1}, 
      { 1, -1}, 
      { 1, 1}, 
      {-1, 1}, 
    };
    px = d[inode][0];
    py = d[inode][1];
    unitx = 2.0 / (ncells[ietype][0] - 1);
    unity = 2.0 / (ncells[ietype][1] - 1);
    px = d[inode][0] * (1 + unitx * (nodal_visual_position - 0.5));
    py = d[inode][1] * (1 + unity * (nodal_visual_position - 0.5));
    element->shape_function(shape_column, px, py, 0, shape_id);
  } else if (shape_id == Element::TRIANGLE) {

  } else if (shape_id == Element::HEXAHEDRON) {
    // natural coordinate of each node
    double d[8][3] = {
      {-1, -1, -1}, 
      { 1, -1, -1}, 
      { 1, 1, -1}, 
      {-1, 1, -1}, 
      {-1, -1, 1}, 
      { 1, -1, 1}, 
      { 1, 1, 1}, 
      {-1, 1, 1}
    };
    unitx = 2.0 / (ncells[ietype][0] - 1);
    unity = 2.0 / (ncells[ietype][1] - 1);
    unitz = 2.0 / (ncells[ietype][2] - 1);
    px = d[inode][0] * (1 + unitx * (nodal_visual_position - 0.5));
    py = d[inode][1] * (1 + unity * (nodal_visual_position - 0.5));
    pz = d[inode][2] * (1 + unitz * (nodal_visual_position - 0.5));
    element->shape_function(shape_column, px, py, pz, shape_id);
  } else if (shape_id == Element::TETRAHEDRON) {
    // natural coordinate of each node
    double d[4][3] = {
      {0, 0, 0}, 
      {1, 0, 0}, 
      {0, 1, 0}, 
      {0, 0, 1}, 
    };
    unitx = 1.0 / (ncells[ietype][0] - 1);
    px = d[inode][0] + unitx * (nodal_visual_position - 0.5) * ((d[inode][0] > 0) ? 1 : -1);
    py = d[inode][1] + unitx * (nodal_visual_position - 0.5) * ((d[inode][1] > 0) ? 1 : -1);
    pz = d[inode][2] + unitx * (nodal_visual_position - 0.5) * ((d[inode][2] > 0) ? 1 : -1);
    element->shape_function(shape_column, px, py, pz, shape_id);

  } else if (shape_id == Element::OCTAHEDRON) {
    error->all(FLERR, "Shape function for octahedron element is not working yet");
  } else if (shape_id == Element::PYRAMID) {
    // natural coordinate of each node
    double d[4][3] = {
      {-1, -1, 0}, 
      { 1, -1, 0}, 
      { 1, 1, 0}, 
      {-1, 1, 0}, 
      {0, 0, 1}, 
    };
    unitx = 1.0 / (ncells[ietype][0] - 1);
    px = d[inode][0] * (1 + unitx * (nodal_visual_position - 0.5));
    py = d[inode][1] * (1 + unitx * (nodal_visual_position - 0.5));
    pz = d[inode][2] + unitx * (nodal_visual_position - 0.5) * ((d[inode][2] > 0) ? 1 : -1);
    element->shape_function(shape_column, px, py, pz, shape_id);
  } else if (shape_id == Element::WEDGE) {
    shape_column[0] = 0.5 * (1.0 - px) * (1.0 - py - pz);
    shape_column[1] = 0.5 * (1.0 - px) * py;
    shape_column[2] = 0.5 * (1.0 - px) * pz;
    shape_column[3] = 0.5 * (1.0 + px) * (1.0 - py - pz);
    shape_column[4] = 0.5 * (1.0 + px) * py;
    shape_column[5] = 0.5 * (1.0 + px) * pz;
  }

  */
}
