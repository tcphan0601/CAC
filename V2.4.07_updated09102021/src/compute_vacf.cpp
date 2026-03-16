#include <cstring>
#include "compute_vacf.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "fix_store.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeVACF::ComputeVACF(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  id_fix(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute vacf command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  create_attribute = 1;

  // options

  interpolate_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"interpolate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute ke command");
      if (strcmp(arg[iarg+1],"yes") == 0) interpolate_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) interpolate_flag = 0;
      else error->all(FLERR,"Illegal compute ke command");
      iarg += 2;
    }
  }
  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "node";
  newarg[5] = (char *) "1";
  newarg[6] = (char *) "3";
  modify->add_fix(7,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // store current velocities in fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **voriginal = fix->astore;
    double **v = atom->v;
    int *amask = atom->mask;
    int nalocal = atom->nlocal;

    double ***nodevoriginal = fix->a3store;
    double ***nodev = element->nodev;
    int *emask = element->mask;
    int nelocal = element->nlocal;
    int *npe = element->npe;
    int *etype = element->etype;

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        voriginal[i][0] = v[i][0];
        voriginal[i][1] = v[i][1];
        voriginal[i][2] = v[i][2];
      } else voriginal[i][0] = voriginal[i][1] = voriginal[i][2] = 0.0;

    for (int i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) {
        for (int j = 0; j < npe[etype[i]]; j++) {
          nodevoriginal[i][j][0] = nodev[i][j][0];
          nodevoriginal[i][j][1] = nodev[i][j][1];
          nodevoriginal[i][j][2] = nodev[i][j][2];
        }
      }
  }

  // displacement vector

  vector = new double[4];
}

/* ---------------------------------------------------------------------- */

ComputeVACF::~ComputeVACF()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeVACF::init()
{
  // set fix which stores original atom velocities

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute vacf fix ID");
  fix = (FixStore *) modify->fix[ifix];

  // nvacf = # of atoms in group

  nvacf = group->count_atom(igroup) + group->count_node(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeVACF::compute_vector()
{
  invoked_vector = update->ntimestep;

  double **voriginal = fix->astore;
  double ***nodevoriginal = fix->a3store;

  double **v = atom->v;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;

  double ***nodev = element->nodev;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;
  int *nintpl = element->nintpl;
  double ***shape_array = element->shape_array;

  double vxsq,vysq,vzsq;
  double vx,vy,vz,vxoriginal,vyoriginal,vzoriginal;
  double vacf[4];
  vacf[0] = vacf[1] = vacf[2] = vacf[3] = 0.0;

  for (int i = 0; i < nalocal; i++)
    if (amask[i] & groupbit) {
      vxsq = v[i][0] * voriginal[i][0];
      vysq = v[i][1] * voriginal[i][1];
      vzsq = v[i][2] * voriginal[i][2];
      vacf[0] += vxsq;
      vacf[1] += vysq;
      vacf[2] += vzsq;
      vacf[3] += vxsq + vysq + vzsq;
    }

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit) {
      int ietype = etype[i];
      vxsq = vysq = vzsq = 0.0;
      if (interpolate_flag) {
        for (int j = 0; j < nintpl[etype[i]]; j++) {
          vx = vxoriginal = 0;
          vy = vyoriginal = 0;
          vz = vzoriginal = 0;
          for (int k = 0; k < npe[ietype]; k++) {
            vx += shape_array[ietype][j][k]*nodev[i][k][0];
            vy += shape_array[ietype][j][k]*nodev[i][k][1];
            vz += shape_array[ietype][j][k]*nodev[i][k][2];
            vxoriginal += shape_array[ietype][j][k]*nodevoriginal[i][k][0];
            vyoriginal += shape_array[ietype][j][k]*nodevoriginal[i][k][1];
            vzoriginal += shape_array[ietype][j][k]*nodevoriginal[i][k][2];
          }
          vxsq += vx*vxoriginal;
          vysq += vy*vyoriginal;
          vzsq += vz*vzoriginal;
        }
      } else {
        for (int j = 0; j < npe[ietype]; j++) {
          vxsq += nodev[i][j][0] * nodevoriginal[i][j][0];
          vysq += nodev[i][j][1] * nodevoriginal[i][j][1];
          vzsq += nodev[i][j][2] * nodevoriginal[i][j][2];
        }
      }
      vacf[0] += vxsq;
      vacf[1] += vysq;
      vacf[2] += vzsq;
      vacf[3] += vxsq + vysq + vzsq;
    }

  MPI_Allreduce(vacf,vector,4,MPI_DOUBLE,MPI_SUM,world);
  if (nvacf) {
    vector[0] /= nvacf;
    vector[1] /= nvacf;
    vector[2] /= nvacf;
    vector[3] /= nvacf;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
   ------------------------------------------------------------------------- */

void ComputeVACF::set_atom_arrays(int i)
{
  double **voriginal = fix->astore;
  double **v = atom->v;

  voriginal[i][0] = v[i][0];
  voriginal[i][1] = v[i][1];
  voriginal[i][2] = v[i][2];

}

/* ----------------------------------------------------------------------
   initialize one element's storage values, called when element is created
   ------------------------------------------------------------------------- */

void ComputeVACF::set_elem_arrays(int i)
{
  double ***nodevoriginal = fix->a3store;
  double ***nodev = element->nodev;
  int *npe = element->npe;
  int *etype = element->etype;

  for (int j = 0; j < npe[etype[i]]; j++) {
    nodevoriginal[i][j][0] = nodev[i][j][0];
    nodevoriginal[i][j][1] = nodev[i][j][1];
    nodevoriginal[i][j][2] = nodev[i][j][2];
  }
}
