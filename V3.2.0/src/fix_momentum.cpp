#include <cstdlib>
#include <cstring>
#include "fix_momentum.h"
#include "atom.h"
#include "element.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "force.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

/*  ----------------------------------------------------------------------  */

FixMomentum::FixMomentum(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal fix momentum command");
  nevery = universe->inumeric(FLERR, arg[3]);
  if (nevery <= 0) error->all(FLERR, "Illegal fix momentum command");

  dynamic = linear = angular = rescale = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal fix momentum command");
      linear = 1;
      xflag = universe->inumeric(FLERR, arg[iarg+1]);
      yflag = universe->inumeric(FLERR, arg[iarg+2]);
      zflag = universe->inumeric(FLERR, arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg], "angular") == 0) {
      angular = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "rescale") == 0) {
      rescale = 1;
      iarg += 1;
    } else error->all(FLERR, "Illegal fix momentum command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR, "Illegal fix momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR, "Illegal fix momentum command");

  dynamic_group_allow = 1;
}

/*  ----------------------------------------------------------------------  */

int FixMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixMomentum::init()
{
  if (group->dynamic[igroup]) {
    dynamic = 1;
  } else {
   if (group->count_atom(igroup) + group->count_elem(igroup) == 0)
     error->all(FLERR, "Fix momentum group has no atoms/elements");
  }

  masstotal = group->mass(igroup);
}

/*  ----------------------------------------------------------------------  */

void FixMomentum::end_of_step()
{
  double **v = atom->v;
  int *amask = atom->mask;
  const int nalocal = atom->nlocal;

  double ****nodev = element->nodev;
  int *emask = element->mask;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  const int nelocal = element->nlocal;
  double ekin_old, ekin_new;
  ekin_old = ekin_new = 0.0;

  if (dynamic)
    masstotal = group->mass(igroup);

  // do nothing is group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;

  // compute kinetic energy before momentum removal, if needed

  if (rescale) {

    double *mass = atom->mass;
    int *atype = atom->type;
    int **ctype = element->ctype;
    double ke = 0.0;
    double sum, **mass_matrix;
    double *nodal_weight = element->nodal_weight;
    double ***mass_matrices = element->mass_matrices;
    int *nucell = element->nucell;
    int *element_shape_ids = element->element_shape_ids;

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit)
        ke += mass[atype[i]] *
          (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

    for (int i = 0; i < nelocal; i++) {
      if (emask[i] & groupbit) {
        if (element->nodal_force_style == Element::CONSISTENT) {
          mass_matrix = mass_matrices[element_shape_ids[etype[i]]];
          sum = 0.0;
          for (int j = 0; j < apc[etype[i]]; j++) {
            for (int k = 0; k < npe[etype[i]]; k++)
              for (int l = 0; l < npe[etype[i]]; l++)
                sum += mass_matrix[k][l] * (
                    nodev[i][j][k][0] * nodev[i][j][l][0] + 
                    nodev[i][j][k][1] * nodev[i][j][l][1] + 
                    nodev[i][j][k][2] * nodev[i][j][l][2]);
            ke += sum * mass[ctype[i][j]] * nucell[etype[i]];
          }
        } else if (element->nodal_force_style == Element::LUMPED) {
          sum = 0.0;
          for (int j = 0; j < apc[etype[i]]; j++) {
            for (int k = 0; k < npe[etype[i]]; k++)
              sum += nodev[i][j][k][0] * nodev[i][j][k][0] + 
                nodev[i][j][k][1] * nodev[i][j][k][1] + 
                nodev[i][j][k][2] * nodev[i][j][k][2];
            ke += sum * mass[ctype[i][j]] * nodal_weight[etype[i]];
          }
        } else {

        }
      } 
    }

    MPI_Allreduce(&ke, &ekin_old, 1, MPI_DOUBLE, MPI_SUM, world);
  }

  if (linear) {
    double vcm[3];
    group->vcm(igroup, masstotal, vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        if (xflag) v[i][0] -= vcm[0];
        if (yflag) v[i][1] -= vcm[1];
        if (zflag) v[i][2] -= vcm[2];
      }
    for (int i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) 
        for (int j = 0; j < apc[etype[i]]; j++) 
          for (int k = 0; k < npe[etype[i]]; k++) {
            if (xflag) nodev[i][j][k][0] -= vcm[0];
            if (yflag) nodev[i][j][k][1] -= vcm[1];
            if (zflag) nodev[i][j][k][2] -= vcm[2];
          }

  }

  if (angular) {
    double xcm[3], angmom[3], inertia[3][3], omega[3];
    group->xcm(igroup, masstotal, xcm);
    group->angmom(igroup, xcm, angmom);
    group->inertia(igroup, xcm, inertia);
    group->omega(angmom, inertia, omega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    double **x = atom->x;
    imageint *image = atom->image;

    double dx, dy, dz;
    double unwrap[3];

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        domain->unmap(x[i], image[i], unwrap);
        dx = unwrap[0] - xcm[0];
        dy = unwrap[1] - xcm[1];
        dz = unwrap[2] - xcm[2];
        v[i][0] -= omega[1] * dz - omega[2] * dy;
        v[i][1] -= omega[2] * dx - omega[0] * dz;
        v[i][2] -= omega[0] * dy - omega[1] * dx;
      }

    x = element->x;
    double ****nodex = element->nodex;
    image = element->image;

    for (int i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) {
        domain->unmap(x[i], image[i], unwrap);
        for (int j = 0; j < apc[etype[i]]; j++) 
          for (int k = 0; k < npe[etype[i]]; k++) {
            dx = nodex[i][j][k][0] - x[i][0] + unwrap[0] - xcm[0];
            dy = nodex[i][j][k][1] - x[i][1] + unwrap[1] - xcm[1];
            dz = nodex[i][j][k][2] - x[i][2] + unwrap[2] - xcm[2];
            nodev[i][j][k][0] -= omega[1] * dz - omega[2] * dy;
            nodev[i][j][k][1] -= omega[2] * dx - omega[0] * dz;
            nodev[i][j][k][2] -= omega[0] * dy - omega[1] * dx;
          }
      } 
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    double ke = 0.0, factor = 1.0;
    double *mass = atom->mass;
    int *atype = atom->type;
    int **ctype = element->ctype;
    double sum, **mass_matrix;
    double *nodal_weight = element->nodal_weight;
    double ***mass_matrices = element->mass_matrices;
    int *element_shape_ids = element->element_shape_ids;
    int *nucell = element->nucell;

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit)
        ke += mass[atype[i]] *
          (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

    for (int i = 0; i < nelocal; i++) {
      if (emask[i] & groupbit) {
        if (element->nodal_force_style == Element::CONSISTENT) {
          mass_matrix = mass_matrices[element_shape_ids[etype[i]]];
          sum = 0.0;
          for (int j = 0; j < apc[etype[i]]; j++) {
            for (int k = 0; k < npe[etype[i]]; k++) 
              for (int l = 0; l < npe[etype[i]]; l++)
                sum += mass_matrix[k][l] * (
                    nodev[i][j][k][0] * nodev[i][j][l][0] + 
                    nodev[i][j][k][1] * nodev[i][j][l][1] + 
                    nodev[i][j][k][2] * nodev[i][j][l][2]);
            ke += sum * mass[ctype[i][j]] * nucell[etype[i]];
          }
        } else if (element->nodal_force_style == Element::LUMPED) {
          sum = 0.0;
          for (int j = 0; j < apc[etype[i]]; j++) {
            for (int k = 0; k < npe[etype[i]]; k++) 
              sum += nodev[i][j][k][0] * nodev[i][j][k][0] + 
                nodev[i][j][k][1] * nodev[i][j][k][1] + 
                nodev[i][j][k][2] * nodev[i][j][k][2];
            ke += sum * mass[ctype[i][j]] * nodal_weight[etype[i]];
          }
        } else {

        }
      } 
    }

    MPI_Allreduce(&ke, &ekin_new, 1, MPI_DOUBLE, MPI_SUM, world);

    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);

    for (int i = 0; i < nalocal; i++) 
      if (amask[i] & groupbit) {
        v[i][0] *= factor;
        v[i][1] *= factor;
        v[i][2] *= factor;
      }

    for (int i = 0; i < nelocal; i++) 
      if (emask[i] & groupbit) 
        for (int j = 0; j < apc[etype[i]]; j++) 
          for (int k = 0; k < npe[etype[i]]; k++) {
            nodev[i][j][k][0] *= factor;
            nodev[i][j][k][1] *= factor;
            nodev[i][j][k][2] *= factor;
          }
  }
}
