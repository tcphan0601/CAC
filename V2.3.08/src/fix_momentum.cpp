#include <cstdlib>
#include <cstring>
#include "fix_momentum.h"
#include "atom.h"
#include "element.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

FixMomentum::FixMomentum(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix momentum command");
  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix momentum command");

  dynamic = linear = angular = rescale = 0;
  vector_flag = 1;
  size_vector = 6;
  
  lcm[0] = lcm[1] = lcm[2] = 
  lcm[3] = lcm[4] = lcm[5] = 0.0;
  
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix momentum command");
      linear = 1;
      xflag = universe->inumeric(FLERR,arg[iarg+1]);
      yflag = universe->inumeric(FLERR,arg[iarg+2]);
      zflag = universe->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angular") == 0) {
      angular = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"rescale") == 0) {
      rescale = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix momentum command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR,"Illegal fix momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR,"Illegal fix momentum command");

  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

int FixMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMomentum::init()
{
  if (group->dynamic[igroup]) {
    dynamic = 1;
  } else {
   if ((group->count_atom(igroup) + group->count_elem(igroup)) == 0)
     error->all(FLERR,"Fix momentum group has no atoms or elements");
  }

  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixMomentum::end_of_step()
{
  double **v = atom->v;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  
  double ***nodev = element->nodev;
  int **nodemask = element->nodemask;
  const int nelocal = element->nlocal;
  const int npe = element->npe;
  int *nintpl = element->nintpl;
  int i,j;
	
  double ekin_old,ekin_new;
  ekin_old = ekin_new = 0.0;

  if (dynamic)
    masstotal = group->mass(igroup);

  // do nothing is group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;

  // compute kinetic energy before momentum removal, if needed

  if (rescale) {

    double *mass = atom->mass;
    int *type = atom->type;
    int *etype  = element->etype;
    int *ectype = element->ctype;
  
    double ke = 0.0;
		
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke +=  mass[type[i]] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
		  
    double weight;
    
    for (i = 0; i < nelocal; i++) {
      weight =  double(nintpl[etype[i]]) / double(npe);
	  for (j = 0; j < npe; j++) {
		if (nodemask[i][j] & groupbit) {
	      ke += mass[ectype[i]] * weight *
		    (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]);
	    }
	  } 
    }
		  
    MPI_Allreduce(&ke,&ekin_old,1,MPI_DOUBLE,MPI_SUM,world);
  }

  if (linear) {
    double vcm[3];
    group->vcm(igroup,masstotal,vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set
    
	lcm[0] = vcm[0];
	lcm[1] = vcm[1];
	lcm[2] = vcm[2];
	
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xflag) v[i][0] -= vcm[0];
        if (yflag) v[i][1] -= vcm[1];
        if (zflag) v[i][2] -= vcm[2];
      }
	  
	for (i = 0; i < nelocal; i++) {
	  for (j = 0; j < npe; j++) {
        if (nodemask[i][j] & groupbit) {
		  if (xflag) nodev[i][j][0] -= vcm[0];
		  if (yflag) nodev[i][j][1] -= vcm[1];
		  if (zflag) nodev[i][j][2] -= vcm[2];
	    }
	  }
    }
  }

  if (angular) {
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    group->xcm(igroup,masstotal,xcm);
    group->angmom(igroup,xcm,angmom);
    group->inertia(igroup,xcm,inertia);
    group->omega(angmom,inertia,omega);
    
	lcm[3] = angmom[0];
	lcm[4] = angmom[1];
	lcm[5] = angmom[2];
	
    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    double **x = atom->x;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    double dx,dy,dz;
    double unwrap[3];

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - xcm[0];
        dy = unwrap[1] - xcm[1];
        dz = unwrap[2] - xcm[2];
        v[i][0] -= omega[1]*dz - omega[2]*dy;
        v[i][1] -= omega[2]*dx - omega[0]*dz;
        v[i][2] -= omega[0]*dy - omega[1]*dx;
      }
	  
    double ***nodex = element->nodex;
    imageint *image_ele = element->image;
  
    for (i = 0; i < nelocal; i++) {
	  for (j = 0; j < npe; j++) {
        if (nodemask[i][j] & groupbit) {
		  domain->unmap(nodex[i][j],image_ele[i],unwrap);
		  dx = unwrap[0] - xcm[0];
          dy = unwrap[1] - xcm[1];
          dz = unwrap[2] - xcm[2];
		  nodev[i][j][0] -= omega[1]*dz - omega[2]*dy;;
		  nodev[i][j][1] -= omega[2]*dx - omega[0]*dz;
		  nodev[i][j][2] -= omega[0]*dy - omega[1]*dx;
	    }
      }
    }
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    double ke=0.0, factor=1.0;
    double *mass = atom->mass;
	int *etype  = element->etype;
    int *ectype = element->ctype;
    int *type = atom->type;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke +=  mass[type[i]] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
	
	double weight;
    
    for (i = 0; i < nelocal; i++) {
      weight =  double(nintpl[etype[i]]) / double(npe);
	  for (j = 0; j < npe; j++) {
		if (nodemask[i][j] & groupbit) {
	      ke += mass[ectype[i]] * weight *
		    (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]);
	    }
	  } 
    }
	
    MPI_Allreduce(&ke,&ekin_new,1,MPI_DOUBLE,MPI_SUM,world);

    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] *= factor;
        v[i][1] *= factor;
        v[i][2] *= factor;
      }
    }
	
	for (i = 0; i < nelocal; i++) {
	  for (j = 0; j < npe; j++) {
		if (nodemask[i][j] & groupbit) {
		  nodev[i][j][0] *= factor;
		  nodev[i][j][1] *= factor;
		  nodev[i][j][2] *= factor;
		}
	  }
	}
  }
}

double FixMomentum::compute_vector(int n)
{

  return lcm[n];
}
