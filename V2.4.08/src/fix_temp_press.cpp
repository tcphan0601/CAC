#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fix_temp_press.h"
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
using namespace MathExtra;

enum{NONE,CONSTANT,EQUAL,ATOM,RANDOM};

#define MAXLINE 1024


/* ---------------------------------------------------------------------- */

FixTempPress::FixTempPress(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  xstr(NULL), idregion(NULL)//, sforce(NULL)

{
  if (narg < 4) error->all(FLERR,"Illegal fix temp/press command");
  if (domain->dimension == 2) error->all(FLERR,"Fix temp/press command does not support 2D yet");

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
  shear_direc = 1;
  const_flag = 1;
  rand_flag = 0;
  srand_flag = 0;
  dynamic_flag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix temp/press command");
      arandom = universe->numeric(FLERR,arg[iarg+1]);
      seed = universe->inumeric(FLERR,arg[iarg+2]);
      if (arandom <= 0.0) error->all(FLERR,"Illegal fix temp/press command");
      if (seed < 1) error->all(FLERR,"Illegal fix temp/press command");
      rand_flag = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"random_shear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix temp/press command");
      arandom = universe->numeric(FLERR,arg[iarg+1]);
      shear_direc = universe->inumeric(FLERR,arg[iarg+2]);
      seed = universe->inumeric(FLERR,arg[iarg+3]);
      if (arandom <= 0.0) error->all(FLERR,"Illegal fix temp/press command");
      if (shear_direc != 1 && shear_direc != 2 && shear_direc != 3) 
        error->all(FLERR,"Illegal fix temp/press command");
      if (seed < 1) error->all(FLERR,"Illegal fix temp/press command");
      srand_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix temp/press command");
      read_file(arg[iarg+1]);
      adynamic = universe->numeric(FLERR,arg[iarg+2]);
      seed = universe->inumeric(FLERR,arg[iarg+3]);
      if (adynamic <= 0.0) error->all(FLERR,"Illegal fix temp/press command");
      if (seed < 1) error->all(FLERR,"Illegal fix temp/press command");
      dynamic_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix temp/press command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix temp/press command");
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
    } else error->all(FLERR,"Illegal fix temp/press command");
  }

  force_flag = 0;
  net_force[0] = net_force[1] = net_force[2] =
    net_force[3] = net_force[4] = net_force[5] = net_force[6] = 0.0;

  //maxatom = 1;
  //memory->create(sforce,maxatom,3,"temp/press:sforce");
}

/* ---------------------------------------------------------------------- */

FixTempPress::~FixTempPress()
{
  delete [] xstr;
  delete [] idregion;
  //memory->destroy(sforce);
  memory->destroy(ffreq);
  memory->destroy(vfreq);
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
      error->all(FLERR,"Variable name for fix temp/press does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix temp/press is invalid style");
  }

  //if (estr) {
  //  evar = input->variable->find(estr);
  //  if (evar < 0)
  //    error->all(FLERR,"Variable name for fix temp/press does not exist");
  //  if (input->variable->atomstyle(evar)) estyle = ATOM;
  //  else error->all(FLERR,"Variable for fix temp/press is invalid style");
  //} else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix temp/press does not exist");
  }

  if (xstyle == ATOM)
    //varflag = ATOM;
    error->all(FLERR,"Atom style variable not ready yet");
  else if (xstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  //if (varflag == CONSTANT && estyle != NONE)
  //  error->all(FLERR,"Cannot use variable energy with "
  //      "constant force in fix temp/press");
  //if ((varflag == EQUAL || varflag == ATOM) &&
  //    update->whichflag == 2 && estyle == NONE)
  //  error->all(FLERR,"Must use variable energy with fix temp/press");
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

  //if ((varflag == ATOM) && atom->nmax > maxatom) {
  //  maxatom = atom->nmax;
  //  memory->destroy(sforce);
  //  memory->create(sforce,maxatom,4,"temp/press:sforce");
  //}

  
  net_force[0] = net_force[1] = net_force[2] = net_force[3] =
    net_force[4] = net_force[5] = net_force[6] = 0.0;
  force_flag = 0;

  double ***nodef = element->nodef;
  double ***nodex = element->nodex;
  double **ex = element->x;
  double *nodal_weight = element->nodal_weight;
  double factor;
  double dist[3],moment[3];
  double surface_area,surface_normal_vector[3];

  int *emask = element->mask;
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int *nintpl = element->nintpl;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  int *nsurface = element->nsurface;
  int ***surface_node_list = element->surface_node_list;
  int **surface_num_node = element->surface_num_node;
  int *element_shape_ids = element->element_shape_ids;
  int i,j,k,ietype;
  int *surface_node,num_node;
  //int surface_node[6][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {3, 4, 7, 8},
  //  {1, 2, 5, 6}, {2, 3, 6, 7}, {1, 4, 5, 8}};
  double vec[3],v[3],v1[3],v2[3];
  int p1,p2,p3;
  double inv_volume;

  // For test
  double freq_temp, amp_temp;
  // For test

  int seed_time;
  seed_time = seed + comm->me + update->ntimestep;
  if (seed_time > 900000000) seed_time = seed_time - floor(seed_time / 900000000);
  RanPark *random = new RanPark(cac,seed_time);

  RanPark *random_p = new RanPark(cac,seed + comm->me);

  // random phase for different freqencies
  
  double rand_phase;

  boltz = force->boltz;
  pvalue = 0.0;
  pnorm[0] = pnorm[1] = pnorm[2] = 0.0;

  // constant force

  if (domain->dimension == 3)
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  else
    inv_volume = 1.0 / (domain->xprd * domain->yprd);

  pvalue_c = tvalue/2.0*3.0*boltz*inv_volume*2*(element->nintpls - element->nnodes);  

  if (varflag == CONSTANT) {

    // Temperature force will only apply on element.

    for (i = 0; i < nelocal; i++) 
      if (emask[i] & groupbit) {
        if (region && !region->match(ex[i])) continue;

        ietype = etype[i];

        factor = 0.25/nodal_weight[ietype]; //weight for nodal force calculation
        pvalue = 0.0;
        if (const_flag) 
          pvalue = pvalue_c * factor;

        if (rand_flag || srand_flag)
          pvalue += pvalue_c * arandom * (random->uniform()-0.5)*2.0 * factor;

        if (dynamic_flag == 1) {
          amp_temp = factor * adynamic * pvalue_c;
          for (k = 0; k < nfreq; k++) {
            rand_phase = random_p->uniform() * MY_2PI;
            pvalue += amp_temp * ffreq[k] * sin(update->ntimestep * update->dt * vfreq[k] * MY_2PI + rand_phase);
          }
        }

        moment[0] = moment[1] = moment[2] = 0.0;
        net_force[6] = pvalue_c;

        for (j = 0; j < nsurface[element_shape_ids[ietype]]; j++) {

          surface_node = surface_node_list[element_shape_ids[ietype]][j];
          num_node = surface_num_node[element_shape_ids[ietype]][j];

          // calculate the surface area for element I

          sub3(nodex[i][surface_node[0]],nodex[i][surface_node[1]],v1);
          sub3(nodex[i][surface_node[2]],nodex[i][surface_node[1]],v2);
          cross3(v1,v2,v);
          surface_area = len3(v)/2.0;

          if (num_node == 4) {
            sub3(nodex[i][surface_node[0]],nodex[i][surface_node[3]],v1);
            sub3(nodex[i][surface_node[2]],nodex[i][surface_node[3]],v2);
            cross3(v1,v2,v);
            surface_area += len3(v)/2.0;
          }
          
          // calculate the surface normal vector for element I
          // normal verctor pointing outward
          
          surface_normal_vector[0] = 0;
          surface_normal_vector[1] = 0;
          surface_normal_vector[2] = 0;
          for (k = 0; k < num_node; k++) {
            p1 = k;
            p2 = (k + 1) % num_node;
            p3 = (k + 2) % num_node;
            sub3(nodex[i][surface_node[p3]],nodex[i][surface_node[p2]],v1);
            sub3(nodex[i][surface_node[p1]],nodex[i][surface_node[p2]],v2);
            cross3(v1,v2,v);
            norm3(v);
            add3(v,surface_normal_vector,surface_normal_vector);
          }
          surface_normal_vector[0] /= num_node;
          surface_normal_vector[1] /= num_node;
          surface_normal_vector[2] /= num_node;

          vec[0] = pvalue * surface_area * surface_normal_vector[0];
          vec[1] = pvalue * surface_area * surface_normal_vector[1];
          vec[2] = pvalue * surface_area * surface_normal_vector[2];
          for (k = 0; k < num_node; k++) {
            if (srand_flag) {
              if (xstyle) nodef[i][surface_node[k]][shear_direc - 1] += vec[shear_direc - 1];
            }  else {
              if (xstyle) nodef[i][surface_node[k]][0] += vec[0];
              if (xstyle) nodef[i][surface_node[k]][1] += vec[1];
              if (xstyle) nodef[i][surface_node[k]][2] += vec[2];
            }
            dist[0] = nodex[i][surface_node[k]][0] - ex[i][0];
            dist[1] = nodex[i][surface_node[k]][1] - ex[i][1];
            dist[2] = nodex[i][surface_node[k]][2] - ex[i][2];

            //??
            moment[0] += dist[0] * vec[0];
            moment[1] += dist[1] * vec[1];
            moment[2] += dist[2] * vec[2];
            net_force[3] += moment[0];
            net_force[4] += moment[1];
            net_force[5] += moment[2];
          }
          net_force[0] += vec[0] * num_node;
          net_force[1] += vec[1] * num_node;
          net_force[2] += vec[2] * num_node;
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
  delete random_p;
}

/*---------------------------------------------------------------------------------
  read PDOS values from file
  -----------------------------------------------------------------------------------*/

void FixTempPress::read_file(char *filename)
{
  int me = comm->me;
  FILE *fpdos;
  int i;
  char *dos;
  char line[MAXLINE];

  if (me == 0) {
    fpdos = force->open_potential(filename);
    if (fpdos == NULL) {
      char str[128];
      sprintf(str,"Cannot open PDOS force file %s",filename);
      error->one(FLERR,str);
    }
  }

  if (me == 0) {
    fgets(line,MAXLINE,fpdos);
    fgets(line,MAXLINE,fpdos);
    sscanf(line,"%d",&nfreq);
  }

  MPI_Bcast(&nfreq,1,MPI_INT,0,world);

  memory->create(ffreq,nfreq,"temp/press:ffreq");
  memory->create(vfreq,nfreq,"temp/press:vfreq");

  double tmp1, tmp2;

  if (me == 0) {
    for (i = 0; i < nfreq; i++) {
      fgets(line,MAXLINE,fpdos);
      sscanf(line,"%lf %lf",&tmp1, &tmp2);
      vfreq[i] = tmp1;
      ffreq[i] = tmp2;
    }
  }

  MPI_Bcast(vfreq,nfreq,MPI_DOUBLE,0,world);
  MPI_Bcast(ffreq,nfreq,MPI_DOUBLE,0,world);

  if (me == 0) fclose(fpdos); 
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
  //if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
