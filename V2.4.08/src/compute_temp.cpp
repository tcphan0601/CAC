#include <mpi.h>
#include <cstring>
#include "compute_temp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  
  // optional args

  accurate_flag = 0;
  mass_style = element->mass_style;
  group_flag = Group::ELEMENT;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp command");
      if (strcmp(arg[iarg+1],"atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1],"node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR,"Illegal compute temp command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"accurate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp command");
      if (strcmp(arg[iarg+1],"yes") == 0) accurate_flag = 1;
      if (strcmp(arg[iarg+1],"no") == 0) accurate_flag = 0;
      else error->all(FLERR,"Illegal compute temp command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp command");
      if (strcmp(arg[iarg+1],"consistent") == 0) mass_style = Element::CONSISTENT;
      if (strcmp(arg[iarg+1],"lumped") == 0) mass_style = Element::LUMPED;
      else error->all(FLERR,"Illegal compute temp command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp command");
  }

  if (group_flag == Group::NODE && mass_style == Element::CONSISTENT)
    error->all(FLERR,"Cannot use group node option with consistent mass style for compute temp command");
 
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTemp::~ComputeTemp()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::init()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::dof_compute()
{
  //adjust_dof_fix();
  
  natoms_temp = group->count_atom(igroup);

  int *etype = element->etype;
  int *npe = element->npe;
  int *nintg = element->nintg;
  int *mask = element->mask;
  int **nodemask = element->nodemask;

  for (int i = 0; i < element->nlocal; i++)
    if (nintg[etype[i]])
      if (group_flag == Group::NODE) {
        for (int j = 0; j < npe[etype[i]]; j++) 
          if (nodemask[i][j] & groupbit) natoms_temp++;
      } else if (group_flag == Group::ELEMENT) {
        if (mask[i] & groupbit) natoms_temp += npe[etype[i]];
      }

  dof = domain->dimension * natoms_temp;
  //dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTemp::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int i,j,k,ietype;

  double t = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      t += mass[type[i]] *
        (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);

  double ***nodev = element->nodev;
  int **nodemask = element->nodemask;
  int *etype  = element->etype;
  int *ctype = element->ctype;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *nintpl = element->nintpl;
  double *nodal_weight = element->nodal_weight;
  double ***mass_matrices = element->mass_matrices;
  int *element_shape_ids = element->element_shape_ids;
  double **mass_matrix;
  double sum;

  if (accurate_flag) {
    error->all(FLERR,"Accurate temp has not been implemented yet");
  } else {
    for (i = 0; i < nelocal; i++) {
      if (group_flag == Group::NODE) {
        sum = 0.0;
        for (j = 0; j < npe[etype[i]]; j++)
          if (nodemask[i][j] & groupbit)
            sum += nodev[i][j][0] * nodev[i][j][0] + 
                   nodev[i][j][1] * nodev[i][j][1] + 
                   nodev[i][j][2] * nodev[i][j][2];
        if (sum != 0.0)
          t += sum * mass[ctype[i]] * nodal_weight[etype[i]];

      } else if (group_flag == Group::ELEMENT) {
        if (emask[i] & groupbit) {
          if (mass_style == Element::CONSISTENT) {
            mass_matrix = mass_matrices[element_shape_ids[etype[i]]];
            sum = 0.0;
            for (j = 0; j < npe[etype[i]]; j++)
              for (k = 0; k < npe[etype[i]]; k++)
                sum += mass_matrix[j][k] * (
                       nodev[i][j][0] * nodev[i][k][0] + 
                       nodev[i][j][1] * nodev[i][k][1] + 
                       nodev[i][j][2] * nodev[i][k][2]);
            t += sum * mass[ctype[i]] * nintpl[etype[i]];
          } else {
            sum = 0.0;
            for (j = 0; j < npe[etype[i]]; j++)
              sum += nodev[i][j][0] * nodev[i][j][0] + 
                     nodev[i][j][1] * nodev[i][j][1] + 
                     nodev[i][j][2] * nodev[i][j][2];
            t += sum * mass[ctype[i]] * nodal_weight[etype[i]];
          }
        } 
      }
    }
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::compute_vector()
{
  //if (comm->me == 0) error->warning(FLERR,"Compute temp only output KE tensors for atomistic region for now");
  int i,j,k;

  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }

  double ***nodev = element->nodev;
  int **nodemask = element->nodemask;
  int *etype  = element->etype;
  int *ctype = element->ctype;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *nintpl = element->nintpl;
  double *nodal_weight = element->nodal_weight;
  double ***mass_matrices = element->mass_matrices;
  int *element_shape_ids = element->element_shape_ids;
  double **mass_matrix;

  if (accurate_flag) {
    error->all(FLERR,"Accurate temp has not been implemented yet");
  } else {
    for (i = 0; i < nelocal; i++) {
      if (group_flag == Group::NODE) {
        massone = mass[ctype[i]] * nodal_weight[etype[i]];
        for (j = 0; j < npe[etype[i]]; j++)
          if (nodemask[i][j] & groupbit) {
            t[0] += massone * nodev[i][j][0] * nodev[i][j][0];
            t[1] += massone * nodev[i][j][1] * nodev[i][j][1];
            t[2] += massone * nodev[i][j][2] * nodev[i][j][2];
            t[3] += massone * nodev[i][j][0] * nodev[i][j][1];
            t[4] += massone * nodev[i][j][0] * nodev[i][j][2];
            t[5] += massone * nodev[i][j][1] * nodev[i][j][2];
          }
      } else if (group_flag == Group::ELEMENT) {
        if (emask[i] & groupbit) {
          if (mass_style == Element::CONSISTENT) {
            massone = mass[ctype[i]] * nintpl[etype[i]];
            mass_matrix = mass_matrices[element_shape_ids[etype[i]]];
            for (j = 0; j < npe[etype[i]]; j++)
              for (k = 0; k < npe[etype[i]]; k++) {
                t[0] += massone * mass_matrix[j][k] * nodev[i][j][0] * nodev[i][k][0];
                t[1] += massone * mass_matrix[j][k] * nodev[i][j][1] * nodev[i][k][1];
                t[2] += massone * mass_matrix[j][k] * nodev[i][j][2] * nodev[i][k][2];
                t[3] += massone * mass_matrix[j][k] * nodev[i][j][0] * nodev[i][k][1];
                t[4] += massone * mass_matrix[j][k] * nodev[i][j][0] * nodev[i][k][2];
                t[5] += massone * mass_matrix[j][k] * nodev[i][j][1] * nodev[i][k][2];
              }
          } else {
            massone = mass[ctype[i]] * nodal_weight[etype[i]];
            for (j = 0; j < npe[etype[i]]; j++) {
              t[0] += massone * nodev[i][j][0] * nodev[i][j][0];
              t[1] += massone * nodev[i][j][1] * nodev[i][j][1];
              t[2] += massone * nodev[i][j][2] * nodev[i][j][2];
              t[3] += massone * nodev[i][j][0] * nodev[i][j][1];
              t[4] += massone * nodev[i][j][0] * nodev[i][j][2];
              t[5] += massone * nodev[i][j][1] * nodev[i][j][2];
            }
          }
        } 
      }
    }
  }


  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;

}
