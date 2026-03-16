#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_temppress.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "force.h"
#include "math_extra.h"
#include "random_park.h"
#include "math_const.h"
#include "comm.h"

using namespace CAC_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE,CONSTANT,EQUAL,ATOM,RANDOM};
enum{NODE,ELEMENT};

/* ---------------------------------------------------------------------- */

FixTempPress::FixTempPress(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  xstr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg < 4) error->all(FLERR,"Illegal fix temppress command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 6;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  xstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    tvalue = universe->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }

  // optional args
  nevery = 1;
  iregion = -1;
  arandom = 0.0;
  adynamic = 0.0;
  seed = 1;
  const_flag = 1;
  rand_flag = 0;
  dynamic_flag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix temppress command");
      arandom = universe->numeric(FLERR,arg[iarg+1]);
	  seed = universe->inumeric(FLERR,arg[iarg+2]);
      if (arandom <= 0.0) error->all(FLERR,"Illegal fix temppress command");
	  if (seed < 1) error->all(FLERR,"Illegal fix temppress command");
	  rand_flag = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix temppress command");
	  adynamic = universe->numeric(FLERR,arg[iarg+1]);
	  seed = universe->inumeric(FLERR,arg[iarg+2]);
	  if (adynamic <= 0.0) error->all(FLERR,"Illegal fix temppress command");
	  if (seed < 1) error->all(FLERR,"Illegal fix temppress command");
	  dynamic_flag = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix temppress command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix temppress command");
      iarg += 2;
	} else if (strcmp(arg[iarg],"constant") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) const_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) const_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix addforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix temppress command");
  }

  force_flag = 0;
  net_force[0] = net_force[1] = net_force[2] =
  net_force[3] = net_force[4] = net_force[5] = net_force[6] = 0.0;
  
  maxatom = 1;
  memory->create(sforce,maxatom,3,"temppress:sforce");
}

/* ---------------------------------------------------------------------- */

FixTempPress::~FixTempPress()
{
  delete [] xstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixTempPress::setmask()
{
  //datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  //mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempPress::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix temppress does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix temppress is invalid style");
  }

  //if (estr) {
  //  evar = input->variable->find(estr);
  //  if (evar < 0)
  //    error->all(FLERR,"Variable name for fix temppress does not exist");
  //  if (input->variable->atomstyle(evar)) estyle = ATOM;
  //  else error->all(FLERR,"Variable for fix temppress is invalid style");
  //} else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix temppress does not exist");
  }

  if (xstyle == ATOM)
    //varflag = ATOM;
    error->all(FLERR,"Atom style variable not ready yet");
  else if (xstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  //if (varflag == CONSTANT && estyle != NONE)
  //  error->all(FLERR,"Cannot use variable energy with "
  //      "constant force in fix temppress");
  //if ((varflag == EQUAL || varflag == ATOM) &&
  //    update->whichflag == 2 && estyle == NONE)
  //  error->all(FLERR,"Must use variable energy with fix temppress");
  
  if (dynamic_flag == 1){
	  element->evec->update_ele_rand_amp(adynamic, seed);
  }

}

/* ---------------------------------------------------------------------- */

void FixTempPress::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTempPress::post_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"temppress:sforce");
  }

  // foriginal = force on atoms before extra force added
  net_force[0] = net_force[1] = net_force[2] = net_force[3] =
  net_force[4] = net_force[5] = net_force[6] = 0.0;
  force_flag = 0;

  double ***nodef = element->nodef;
  double ***nodex = element->nodex;
  double ***surface_normvec = element->surface_normvec;
  double **ele_surf_area = element->ele_surf_area;
  double **ex = element->x;
  double *ele_rand_amp = element->ele_rand_amp;
  double freq, factor;
  double dist[3],moment[3];
  
  int *emask = element->mask;
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int *nintpl = element->nintpl;
  int nelocal = element->nlocal;
  int npe = element->npe;
  int i,j;
  int surface_node[6][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {3, 4, 7, 8},
                 {1, 2, 5, 6}, {2, 3, 6, 7}, {1, 4, 5, 8}};
  double vec[3];
  double inv_volume;
  
  // For test
  double freq_temp, amp_temp;
  // For test
  
  RanPark *random = new RanPark(cac,seed + comm->me);
  
  // random phase for different freqencies
  int num_freq;
  num_freq = 36;
  double rand_phase;
  //double rand_phase[nelocal][num_freq];
  //double rand_phase[nelocal];
  //if (dynamic_flag == 1) {
  //    for (i = 0; i < nelocal; i++) {
		  //rand_phase[i] = random->uniform() * MY_PI * 2.0;
	      //for (j = 0; j < num_freq; j ++) {
		      //rand_phase[i][j] = random->uniform() * MY_PI * 2.0;
	      //}
      //}
  //}
  // random phase for different freqencies
  
  boltz = force-> boltz;
  pvalue = 0.0;
  pnorm[0] = pnorm[1] = pnorm[2] = 0.0;
  freq = 5.0; // THz for highest freqency of iron along [001] direction. Only for test.
  
  // constant force
  // potential energy = - x dot f in unwrapped coords

  element->evec->update_ele_surf_area();
  element->evec->update_surface_normal_vector();
  
  if (domain->dimension == 3)
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  else
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
  pvalue_c = tvalue/2.0*3.0*boltz*inv_volume*2*(element->max_nintpl-8)*element->nelements;  

  if (varflag == CONSTANT) {
	//Temperature force will only apply on element.
    for (i = 0; i < nelocal; i++) 
        if (emask[i] & groupbit) {
          if (region && !region->match(ex[i])) continue;
		  
		  factor = double(npe)/4.0/double(nintpl[etype[i]]) / double(element->apc);
		  pvalue = 0.0;
		  if (const_flag == 1) {
			pvalue = pvalue_c * factor;
		  }
		  if (rand_flag == 1){
			pvalue += pvalue_c * arandom * (random->uniform()-0.5)*2.0 * factor;
		  } 
		  if (dynamic_flag == 1){
			//amp_temp = pvalue_c * ele_rand_amp[i] * factor;
			amp_temp = pvalue_c * adynamic * factor;
			//rand_phase = random->uniform() * MY_PI * 2.0;
			//pvalue += amp_temp * sin(update->ntimestep * update->dt * 5.0 * 2 * MY_PI + rand_phase);
		    //rand_phase = random->uniform() * MY_PI * 2.0;
			//pvalue += 1.0 * amp_temp * sin(update->ntimestep * update->dt * 10.0 * 2 * MY_PI  + rand_phase);
		    //pvalue += 1.0 * amp_temp * sin(update->ntimestep * update->dt * 10.0 * 2 * MY_PI  + rand_phase[i][1]);
			//pvalue += 1.0 * amp_temp * sin(update->ntimestep * update->dt * 10.0 * 2 * MY_PI  + rand_phase[i][2]);
			for (int ii = 0; ii < num_freq; ii++ ) {
			  freq_temp = double(ii + 5);
			  rand_phase = random->uniform() * MY_PI * 2.0;
			  pvalue += 2.0 * amp_temp * random->uniform() * sin(update->ntimestep * update->dt * freq_temp * 2 * MY_PI + rand_phase);
			}
		  }
		  moment[0] = moment[1] = moment[2] = 0.0;
		  net_force[6] = pvalue_c;

		  for (j=0; j < 6; j++) {
			vec[0] = pvalue * ele_surf_area[i][j] * surface_normvec[i][j][0];
			vec[1] = pvalue * ele_surf_area[i][j] * surface_normvec[i][j][1];
			vec[2] = pvalue * ele_surf_area[i][j] * surface_normvec[i][j][2];
			for (int k = 0; k < 4; k++) {
			  if (xstyle) nodef[i][surface_node[j][k]-1][0] += vec[0];
			  if (xstyle) nodef[i][surface_node[j][k]-1][1] += vec[1];
			  if (xstyle) nodef[i][surface_node[j][k]-1][2] += vec[2];
			  dist[0] = nodex[i][surface_node[j][k]-1][0]-ex[i][0];
			  dist[1] = nodex[i][surface_node[j][k]-1][1]-ex[i][1];
			  dist[2] = nodex[i][surface_node[j][k]-1][2]-ex[i][2];
			  moment[0] += dist[0] * vec[0];
			  moment[1] += dist[1] * vec[1];
			  moment[2] += dist[2] * vec[2];
			  net_force[3] += moment[0];
			  net_force[4] += moment[1];
			  net_force[5] += moment[2];
			}
			  net_force[0] += vec[0] * 4.0;
			  net_force[1] += vec[1] * 4.0;
			  net_force[2] += vec[2] * 4.0;
		  }
        }

    // variable force, wrap with clear/add
    // potential energy = evar if defined, else 0.0
    // wrap with clear/add

  } else {
	error->all(FLERR,"Variable has not been added yet.");
	/*
    modify->clearstep_compute();

    if (xstyle == EQUAL) tvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);

    modify->addstep_compute(update->ntimestep + 1);

	//Temperature force will only apply on element.

	
    for (i = 0; i < nelocal; i++) 
      if (group_flag) {
        if (emask[i] & groupbit) {
          if (region && !region->match(ex[i])) continue;
          for (j = 0; j < npe; j++) {
            foriginal[1] += nodef[i][j][0];
            foriginal[2] += nodef[i][j][1];
            foriginal[3] += nodef[i][j][2];
            if (xstyle) nodef[i][j][0] += tvalue;
          }
        }
      } else error->all(FLERR,"Temperature pressure could only be added on element.");
	  */
  }
  // clean up
  delete random;
}
/* ----------------------------------------------------------------------
   potential energy of added force
   ------------------------------------------------------------------------- */

double FixTempPress::compute_scalar()
{
  // only sum across procs one time

  //if (force_flag == 0) {
  //  MPI_Allreduce(net_force,net_force_all,4,MPI_DOUBLE,MPI_SUM,world);
  //  force_flag = 1;
  //}
  return net_force[6];
  //return 0.0;
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
   ------------------------------------------------------------------------- */

double FixTempPress::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(net_force,net_force_all,6,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return net_force_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double FixTempPress::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
