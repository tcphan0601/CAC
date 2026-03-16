#include <cmath>
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "fix_thermal_conductivity.h"
#include "atom.h"
#include "element.h"
#include "force.h"
#include "universe.h"
#include "domain.h"
#include "modify.h"
#include "error.h"

using namespace CAC_NS;
using namespace FixConst;

#define BIG 1.0e10

/* ---------------------------------------------------------------------- */

FixThermalConductivity::FixThermalConductivity(CAC *cac,
                                               int narg, char **arg) :
  Fix(cac, narg, arg),
  index_lo(NULL), index_hi(NULL), ke_lo(NULL), ke_hi(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix thermal/conductivity command");

  if (atom->natoms > 0) error->all(FLERR,"Fix thermal/conductivity command is only available for pure CG right now");
  
  MPI_Comm_rank(world,&me);

  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix thermal/conductivity command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 0;

  if (strcmp(arg[4],"x") == 0) edim = 0;
  else if (strcmp(arg[4],"y") == 0) edim = 1;
  else if (strcmp(arg[4],"z") == 0) edim = 2;
  else error->all(FLERR,"Illegal fix thermal/conductivity command");

  nbin = universe->inumeric(FLERR,arg[5]);
  if (nbin % 2 || nbin <= 2)
    error->all(FLERR,"Illegal fix thermal/conductivity command, nbin must be even");

  // optional keywords

  nswap = 1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"swap") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix thermal/conductivity command");
      nswap = universe->inumeric(FLERR,arg[iarg+1]);
      if (nswap <= 0)
        error->all(FLERR,
                   "Fix thermal/conductivity swap value must be positive");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix thermal/conductivity command");
  }

  // initialize array sizes to nswap+1 so have space to shift values down

  index_lo = new int[nswap+1][2];
  index_hi = new int[nswap+1][2];
  ke_lo = new double[nswap+1];
  ke_hi = new double[nswap+1];

  e_exchange = 0.0;
}

/* ---------------------------------------------------------------------- */

FixThermalConductivity::~FixThermalConductivity()
{
  delete [] index_lo;
  delete [] index_hi;
  delete [] ke_lo;
  delete [] ke_hi;
}

/* ---------------------------------------------------------------------- */

int FixThermalConductivity::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixThermalConductivity::init()
{
  // warn if any fix ave/spatial comes after this fix
  // can cause glitch in averaging since ave will happen after swap
  /*
  int foundme = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i] == this) foundme = 1;
    if (foundme && strcmp(modify->fix[i]->style,"ave/spatial") == 0 && me == 0)
      error->warning(FLERR,
                     "Fix thermal/conductivity comes before fix ave/spatial");
  }
  */
  // set bounds of 2 slabs in edim
  // only necessary for static box, else re-computed in end_of_step()
  // lo bin is always bottom bin
  // hi bin is just above half height
  if (domain->box_change == 0) {
    prd = domain->prd[edim];
    boxlo = domain->boxlo[edim];
    boxhi = domain->boxhi[edim];
    double binsize = (boxhi-boxlo) / nbin;
    slablo_lo = boxlo;
    slablo_hi = boxlo + binsize;
    slabhi_lo = boxlo + (nbin/2)*binsize;
    slabhi_hi = boxlo + (nbin/2+1)*binsize;
  }
  periodicity = domain->periodicity[edim];
}

/* ---------------------------------------------------------------------- */

void FixThermalConductivity::end_of_step()
{
  int i,j,m,insert,k,l;
  double coord,ke;
  struct {
    double value;
    int proc;
  } mine[2],all[2];

  // if box changes, recompute bounds of 2 slabs in edim

  if (domain->box_change) {
    prd = domain->prd[edim];
    boxlo = domain->boxlo[edim];
    boxhi = domain->boxhi[edim];
    double binsize = (boxhi-boxlo) / nbin;
    slablo_lo = boxlo;
    slablo_hi = boxlo + binsize;
    slabhi_lo = boxlo + (nbin/2)*binsize;
    slabhi_hi = boxlo + (nbin/2+1)*binsize;
  }
  
  // make 2 lists of up to nswap atoms
  // hottest atoms in lo slab, coldest atoms in hi slab (really mid slab)
  // lo slab list is sorted by hottest, hi slab is sorted by coldest
  // map atoms back into periodic box if necessary
  // insert = location in list to insert new atom

  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double *mass = atom->mass;
  int *type = atom->type;
  int *etype = element->etype;
  int *ectype = element->ctype;
  int **nodemask = element->nodemask;
  int *nintpl = element->nintpl;
  int nelocal = element->nlocal;
  int npe = element->npe;

  nhi = nlo = 0;

  for (i = 0; i < nelocal; i++)
	for (int j = 0; j < npe; j++) {
	  if (nodemask[i][j] & groupbit) {
        coord = nodex[i][j][edim];
        if (coord < boxlo && periodicity) coord += prd;
        else if (coord >= boxhi && periodicity) coord -= prd;

        if (coord >= slablo_lo && coord < slablo_hi) {
          ke = nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2];
          ke *= 0.5*mass[ectype[i]]/double(npe)*double(nintpl[etype[i]]);
          if (nlo < nswap || ke > ke_lo[nswap-1]) {
            for (insert = nlo-1; insert >= 0; insert--)
              if (ke < ke_lo[insert]) break;
            insert++;
            for (m = nlo-1; m >= insert; m--) {
              ke_lo[m+1] = ke_lo[m];
              index_lo[m+1][0] = index_lo[m][0];
			  index_lo[m+1][1] = index_lo[m][1];
            }
            ke_lo[insert] = ke;
            index_lo[insert][0] = i;
			index_lo[insert][1] = j;
            if (nlo < nswap) nlo++;
          }
        }
		
        if (coord >= slabhi_lo && coord < slabhi_hi) {
          ke = nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2];
          ke *= 0.5*mass[ectype[i]]/double(npe)*double(nintpl[etype[i]]);
          if (nhi < nswap || ke < ke_hi[nswap-1]) {
            for (insert = nhi-1; insert >= 0; insert--)
              if (ke > ke_hi[insert]) break;
            insert++;
            for (m = nhi-1; m >= insert; m--) {
              ke_hi[m+1] = ke_hi[m];
              index_hi[m+1][0] = index_hi[m][0];
			  index_hi[m+1][1] = index_hi[m][1];
            }
            ke_hi[insert] = ke;
            index_hi[insert][0] = i;
			index_hi[insert][1] = j;
            if (nhi < nswap) nhi++;
          }
        }
      }
    }
  
  // loop over nswap pairs
  // pair 2 global atoms at beginning of sorted lo/hi slab lists via Allreduce
  // BIG values are for procs with no atom to contribute
  // use negative of hottest KE since is doing a MINLOC
  // MINLOC also communicates which procs own them
  // exchange kinetic energy between the 2 particles
  // if I own both particles just swap, else point2point comm of velocities

  double sbuf[4],rbuf[4],vcm[3];
  double eswap = 0.0;

  mine[0].proc = mine[1].proc = me;
  int ilo = 0;
  int ihi = 0;
  
  for (m = 0; m < nswap; m++) {
    if (ilo < nlo) mine[0].value = -ke_lo[ilo];
    else mine[0].value = BIG;
    if (ihi < nhi) mine[1].value = ke_hi[ihi];
    else mine[1].value = BIG;

    MPI_Allreduce(mine,all,2,MPI_DOUBLE_INT,MPI_MINLOC,world);
    if (all[0].value == BIG || all[1].value == BIG) continue;

    if (me == all[0].proc && me == all[1].proc) {
      i = index_lo[ilo++][0];
      k = index_hi[ihi++][0];
	  j = index_lo[ilo++][1];
	  l = index_hi[ihi++][1];
      sbuf[0] = nodev[k][l][0];
      sbuf[1] = nodev[k][l][1];
      sbuf[2] = nodev[k][l][2];
      sbuf[3] = mass[ectype[k]]/double(npe)*double(nintpl[etype[k]]);
      rbuf[0] = nodev[i][j][0];
      rbuf[1] = nodev[i][j][1];
      rbuf[2] = nodev[i][j][2];
      rbuf[3] = mass[ectype[i]]/double(npe)*double(nintpl[etype[i]]);
      vcm[0] = (sbuf[3]*sbuf[0] + rbuf[3]*rbuf[0]) / (sbuf[3] + rbuf[3]);
      vcm[1] = (sbuf[3]*sbuf[1] + rbuf[3]*rbuf[1]) / (sbuf[3] + rbuf[3]);
      vcm[2] = (sbuf[3]*sbuf[2] + rbuf[3]*rbuf[2]) / (sbuf[3] + rbuf[3]);
      nodev[k][l][0] = 2.0 * vcm[0] - sbuf[0];
      nodev[k][l][1] = 2.0 * vcm[1] - sbuf[1];
      nodev[k][l][2] = 2.0 * vcm[2] - sbuf[2];
      eswap += sbuf[3] * (vcm[0] * (vcm[0] - sbuf[0]) +
                          vcm[1] * (vcm[1] - sbuf[1]) +
                          vcm[2] * (vcm[2] - sbuf[2]));
      nodev[i][j][0] = 2.0 * vcm[0] - rbuf[0];
      nodev[i][j][1] = 2.0 * vcm[1] - rbuf[1];
      nodev[i][j][2] = 2.0 * vcm[2] - rbuf[2];
      eswap -= rbuf[3] * (vcm[0] * (vcm[0] - rbuf[0]) +
                          vcm[1] * (vcm[1] - rbuf[1]) +
                          vcm[2] * (vcm[2] - rbuf[2]));

    } else if (me == all[0].proc) {
      k = index_lo[ilo++][0];
	  l = index_lo[ilo++][1];
      sbuf[0] = nodev[k][l][0];
      sbuf[1] = nodev[k][l][1];
      sbuf[2] = nodev[k][l][2];
      sbuf[3] = mass[ectype[k]]/double(npe)*double(nintpl[etype[k]]);
      MPI_Sendrecv(sbuf,4,MPI_DOUBLE,all[1].proc,0,
                   rbuf,4,MPI_DOUBLE,all[1].proc,0,world,MPI_STATUS_IGNORE);
      vcm[0] = (sbuf[3]*sbuf[0] + rbuf[3]*rbuf[0]) / (sbuf[3] + rbuf[3]);
      vcm[1] = (sbuf[3]*sbuf[1] + rbuf[3]*rbuf[1]) / (sbuf[3] + rbuf[3]);
      vcm[2] = (sbuf[3]*sbuf[2] + rbuf[3]*rbuf[2]) / (sbuf[3] + rbuf[3]);
      nodev[k][l][0] = 2.0 * vcm[0] - sbuf[0];
      nodev[k][l][1] = 2.0 * vcm[1] - sbuf[1];
      nodev[k][l][2] = 2.0 * vcm[2] - sbuf[2];
      eswap -= sbuf[3] * (vcm[0] * (vcm[0] - sbuf[0]) +
                          vcm[1] * (vcm[1] - sbuf[1]) +
                          vcm[2] * (vcm[2] - sbuf[2]));
						  
    } else if (me == all[1].proc) {
      k = index_hi[ihi++][0];
	  l = index_hi[ihi++][1];
      sbuf[0] = nodev[k][l][0];
      sbuf[1] = nodev[k][l][1];
      sbuf[2] = nodev[k][l][2];
      sbuf[3] = mass[ectype[k]]/double(npe)*double(nintpl[etype[k]]);
      MPI_Sendrecv(sbuf,4,MPI_DOUBLE,all[0].proc,0,
                   rbuf,4,MPI_DOUBLE,all[0].proc,0,world,MPI_STATUS_IGNORE);
      vcm[0] = (sbuf[3]*sbuf[0] + rbuf[3]*rbuf[0]) / (sbuf[3] + rbuf[3]);
      vcm[1] = (sbuf[3]*sbuf[1] + rbuf[3]*rbuf[1]) / (sbuf[3] + rbuf[3]);
      vcm[2] = (sbuf[3]*sbuf[2] + rbuf[3]*rbuf[2]) / (sbuf[3] + rbuf[3]);
      nodev[k][l][0] = 2.0 * vcm[0] - sbuf[0];
      nodev[k][l][1] = 2.0 * vcm[1] - sbuf[1];
      nodev[k][l][2] = 2.0 * vcm[2] - sbuf[2];
      eswap += sbuf[3] * (vcm[0] * (vcm[0] - sbuf[0]) +
                          vcm[1] * (vcm[1] - sbuf[1]) +
                          vcm[2] * (vcm[2] - sbuf[2]));
    }
  }

  // tally energy exchange from all swaps

  double eswap_all;
  MPI_Allreduce(&eswap,&eswap_all,1,MPI_DOUBLE,MPI_SUM,world);
  e_exchange += force->mvv2e * eswap_all;
}

/* ---------------------------------------------------------------------- */

double FixThermalConductivity::compute_scalar()
{
  return e_exchange;
}
