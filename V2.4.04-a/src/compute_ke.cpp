#include <mpi.h>
#include <string.h>
#include "compute_ke.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "element.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeKE::ComputeKE(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ke command");
  if (igroup) error->all(FLERR,"Compute ke must use group all");

  scalar_flag = 1;
  extscalar = 1;

  pairflag = 1;
  fixflag = 1;
  accurate_flag = 0;
  mass_style = element->mass_style;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute ke command");
      if (strcmp(arg[iarg+1],"atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1],"node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR,"Illegal compute ke command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"accurate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute ke command");
      if (strcmp(arg[iarg+1],"yes") == 0) accurate_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) accurate_flag = 0;
      else error->all(FLERR,"Illegal compute ke command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mass") == 0) {
    } else if (strcmp(arg[iarg],"accurate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute ke command");
      if (strcmp(arg[iarg+1],"consistent") == 0) mass_style = Element::CONSISTENT;
      if (strcmp(arg[iarg+1],"lumped") == 0) mass_style = Element::LUMPED;
      else error->all(FLERR,"Illegal compute ke command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute ke command");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeKE::init()
{
  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeKE::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int i,j,k;

  double ke = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      ke += mass[type[i]] *
        (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);

  double ***nodev = element->nodev;
  int *emask = element->mask;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int *etype  = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;
  double *nodal_weight = element->nodal_weight;
  int *nintpl = element->nintpl;
  double ***mass_matrices = element->mass_matrices;
  int *element_shape_ids = element->element_shape_ids;
  double **mass_matrix;
  double sum;

  if (accurate_flag) {
    error->all(FLERR,"Accurate compute ke has not been implemented yet");
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
          ke += sum * mass[ctype[i]] * nodal_weight[etype[i]];

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
            ke += sum * mass[ctype[i]] * nintpl[etype[i]];
          } else {
            sum = 0.0;
            for (j = 0; j < npe[etype[i]]; j++)
              sum += nodev[i][j][0] * nodev[i][j][0] + 
                nodev[i][j][1] * nodev[i][j][1] + 
                nodev[i][j][2] * nodev[i][j][2];
            ke += sum * mass[ctype[i]] * nodal_weight[etype[i]];
          }
        } 
      }
    }
  }

  MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}
