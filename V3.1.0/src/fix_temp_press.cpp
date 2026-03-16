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
#include "random_mars.h"
#include "math_const.h"
#include "comm.h"

using namespace CAC_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

enum{NONE, CONSTANT, EQUAL, ATOM, RANDOM};
enum{CONST, RAND, PHONON_ENRICH};
enum{USER, BOX, DYNAMIC};
#define MAXLINE 1024


/*  ----------------------------------------------------------------------  */

FixTempPress::FixTempPress(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  tstr(nullptr)

{
  if (narg < 8) error->all(FLERR, "Illegal fix temp/press command");
  if (domain->dimension == 2) error->all(FLERR, "Fix temp/press command does not support 2D yet");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 6;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  if (strstr(arg[3], "v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr = new char[n];
    strcpy(tstr, &arg[3][2]);
  } else {
    t_start = universe->numeric(FLERR, arg[3]);
    t_target = t_start;
    tstyle = CONSTANT;
  }
  t_stop = universe->numeric(FLERR, arg[4]);


  direction_flag[0] = direction_flag[1] = direction_flag[2] = 0;
  seed = universe->inumeric(FLERR, arg[5]);
  if (seed <= 0) error->all(FLERR, "Illegal fix temp/press command, seed must be >= 1");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(cac, seed + comm->me);

  if (strcmp(arg[6], "constant") == 0) {
    style = CONST;
    if (strcmp(arg[7], "x") == 0) 
      direction_flag[0] = 1;
    else if (strcmp(arg[7], "y") == 0) 
      direction_flag[1] = 1;
    else if (strcmp(arg[7], "z") == 0) 
      direction_flag[2] = 1;
    else if (strcmp(arg[7], "xy") == 0 || strcmp(arg[7], "yx") == 0)
      direction_flag[0] = direction_flag[1] = 1;
    else if (strcmp(arg[7], "xz") == 0 || strcmp(arg[7], "zx") == 0)
      direction_flag[0] = direction_flag[2] = 1;
    else if (strcmp(arg[7], "zy") == 0 || strcmp(arg[7], "yz") == 0)
      direction_flag[1] = direction_flag[2] = 1;
    else if (strcmp(arg[7], "all") == 0) 
      direction_flag[0] = direction_flag[1] = direction_flag[2] = 1;
    else error->all(FLERR, "Illegal fix temp/press command");
  } else if (strcmp(arg[6], "random") == 0) {
    style = RAND;
    if (strcmp(arg[7], "x") == 0) 
      direction_flag[0] = 1;
    else if (strcmp(arg[7], "y") == 0) 
      direction_flag[1] = 1;
    else if (strcmp(arg[7], "z") == 0) 
      direction_flag[2] = 1;
    else if (strcmp(arg[7], "xy") == 0 || strcmp(arg[7], "yx") == 0)
      direction_flag[0] = direction_flag[1] = 1;
    else if (strcmp(arg[7], "xz") == 0 || strcmp(arg[7], "zx") == 0)
      direction_flag[0] = direction_flag[2] = 1;
    else if (strcmp(arg[7], "zy") == 0 || strcmp(arg[7], "yz") == 0)
      direction_flag[1] = direction_flag[2] = 1;
    else if (strcmp(arg[7], "all") == 0) 
      direction_flag[0] = direction_flag[1] = direction_flag[2] = 1;
    else error->all(FLERR, "Illegal fix temp/press command");
  } else if (strcmp(arg[6], "phonon/enrich") == 0) {
    style = PHONON_ENRICH;
    direction_flag[0] = direction_flag[1] = direction_flag[2] = 1;
    read_file(arg[7]);
  }

  press_fluct_amplitude = universe->numeric(FLERR, arg[8]);;
  if (press_fluct_amplitude < 0.0) error->all(FLERR, "Illegal fix temp/press command, pressure fluctuation amplitude must non-negative");

  // optional args

  nevery = 1;
  zero_mean_flag = 0;
  ucell_volume_user = 0.0;
  ucell_volume_style = BOX;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "every") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix temp/press command");
      nevery = universe->inumeric(FLERR, arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR, "Illegal fix temp/press command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ucell/volume") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix temp/press command");
      if (strcmp(arg[iarg+1], "box") == 0) ucell_volume_style = BOX;
      else if (strcmp(arg[iarg+1], "dynamic") == 0) ucell_volume_style = DYNAMIC;
      else {
        ucell_volume_style = USER;
        ucell_volume_user = universe->numeric(FLERR, arg[iarg+1]);
        if (ucell_volume_user <= 0) error->all(FLERR, "Illegal fix temp/press command, unitcell volume must be positive");
      }
      iarg += 2;

    } else if (strcmp(arg[iarg], "zero/mean") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal temp/press command");
      if (strcmp(arg[iarg+1], "no") == 0) zero_mean_flag = 0;
      else if (strcmp(arg[iarg+1], "yes") == 0) zero_mean_flag = 1;
      else error->all(FLERR, "Illegal fix temp/press command");
      iarg += 2;
    } else error->all(FLERR, "Illegal fix temp/press command");
  }

  force_flag = 0;
  net_force[0] = net_force[1] = net_force[2] =
    net_force[3] = net_force[4] = net_force[5] = net_force[6] = 0.0;
}

/*  ----------------------------------------------------------------------  */

FixTempPress::~FixTempPress()
{
  delete [] tstr;
  delete [] random;
  memory->destroy(ffreq);
  memory->destroy(vfreq);
}

/*  ----------------------------------------------------------------------  */

int FixTempPress::setmask()
{
  //datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  //mask |= THERMO_ENERGY;
  return mask;
}

/*  ----------------------------------------------------------------------  */

void FixTempPress::init()
{
  // check variables

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR, "Variable name for fix temp/press does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else if (input->variable->atomstyle(tvar)) tstyle = ATOM;
    else error->all(FLERR, "Variable for fix temp/press is invalid style");
  }

}

/*  ----------------------------------------------------------------------  */

void FixTempPress::setup(int vflag)
{
  if (strstr(update->integrate_style, "verlet"))
    post_force(vflag);
}

/*  ----------------------------------------------------------------------  */

void FixTempPress::post_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  // reallocate sforce array if necessary

  net_force[0] = net_force[1] = net_force[2] = net_force[3] =
    net_force[4] = net_force[5] = net_force[6] = 0.0;
  force_flag = 0;

  double ****nodef = element->nodef;
  double ****nodex = element->nodex;
  double **ex = element->x;
  double *nodal_weight = element->nodal_weight;
  double *ucell_vol = element->ucell_vol;
  double weight_factor, t;
  double **inodex;
  double del[3];
  double surface_normal_vector[3], surface_force[3];

  int *emask = element->mask;
  int *etype = element->etype;
  int *nucell = element->nucell;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  int *nsurface = element->nsurface;
  int ***surface_node_list = element->surface_node_list;
  int **surface_num_node = element->surface_num_node;
  int *element_shape_ids = element->element_shape_ids;
  int i, ietype, ibasis, isurface, inode;
  int *surface_node, num_node;
  double v[3], v1[3], v2[3];
  int p1, p2, p3;
  double inv_volume;



  compute_target();


  boltz = force->boltz;
  t = update->ntimestep * update->dt;

  // mean force
  // if user define unit cell volumes, the (N_i - n * N_e) / Omega part is multiplied within the loop for each element
  // otherwise, assuming the simulation cell is fully filled and use it to estimate the (N_i - n * N_e) / Omega part

  mean_press = t_target * 3.0 * boltz;
  if (ucell_volume_style == BOX) {
    element->count_vatoms();
    element->count_nodes();
    double CG_volume = domain->xprd * domain->yprd * domain->zprd * element->nvatoms / (element->nvatoms + atom->natoms);
    mean_press *= (element->nvatoms - element->nnodes) / CG_volume;
  }


  // Temperature force will only apply on element.

  for (i = 0; i < nelocal; i++) {
    if (emask[i] & groupbit) {
      ietype = etype[i];

      weight_factor = 0.25 / nodal_weight[ietype]; //weight for nodal force calculation
      if (zero_mean_flag) 
        press = 0.0;
      else 
        press = mean_press * weight_factor;

      if (style == RAND)
        press += mean_press * press_fluct_amplitude * (random->uniform() - 0.5) * 2.0 * weight_factor;

      else if (style == PHONON_ENRICH) {
        double amp_temp = weight_factor * press_fluct_amplitude * mean_press;

        // random phases for different freqencies
        
        for (int ifreq = 0; ifreq < nfreq; ifreq++) 
          press += amp_temp * ffreq[ifreq] * sin((t * vfreq[ifreq] + random->uniform()) * MY_2PI);
      }

      if (ucell_volume_style != BOX) {
        double ucell_volume;
        if (ucell_volume_style == USER)
          ucell_volume = ucell_volume_user;
        else if (ucell_volume_style == DYNAMIC)
          ucell_volume = ucell_vol[i];
        press *= (nucell[ietype] - npe[ietype]) / nucell[ietype] / ucell_volume;
        net_force[6] += mean_press * (nucell[ietype] - npe[ietype]) / nucell[ietype] / ucell_volume;
      } else net_force[6] += mean_press;



      for (isurface = 0; isurface < nsurface[element_shape_ids[ietype]]; isurface++) {
        surface_node = surface_node_list[element_shape_ids[ietype]][isurface];
        num_node = surface_num_node[element_shape_ids[ietype]][isurface];

        for (ibasis = 0; ibasis < apc[ietype]; ibasis++) {

          inodex = nodex[i][ibasis];

          // calculate the surface normal vector for element I 
          // with magnitude scaled to average area
          // normal verctor pointing outward

          surface_normal_vector[0] = 0.0;
          surface_normal_vector[1] = 0.0;
          surface_normal_vector[2] = 0.0;

          if (num_node == 4) {
            for (inode = 0; inode < 4; inode++) {
              p1 = inode;
              p2 = (inode + 1) % num_node;
              p3 = (inode + 2) % num_node;
              sub3(inodex[surface_node[p3]], inodex[surface_node[p2]], v1);
              sub3(inodex[surface_node[p1]], inodex[surface_node[p2]], v2);
              surface_normal_vector[0] += v1[1] * v2[2] - v1[2] * v2[1];
              surface_normal_vector[1] += v1[2] * v2[0] - v1[0] * v2[2];
              surface_normal_vector[2] += v1[0] * v2[1] - v1[1] * v2[0];
            }
            surface_normal_vector[0] /= 4.0;
            surface_normal_vector[1] /= 4.0;
            surface_normal_vector[2] /= 4.0;
          } else if (num_node == 3) {
            sub3(inodex[surface_node[2]], inodex[surface_node[1]], v1);
            sub3(inodex[surface_node[0]], inodex[surface_node[1]], v2);
            surface_normal_vector[0] += (v1[1] * v2[2] - v1[2] * v2[1]) / 2.0;
            surface_normal_vector[1] += (v1[2] * v2[0] - v1[0] * v2[2]) / 2.0;
            surface_normal_vector[2] += (v1[0] * v2[1] - v1[1] * v2[0]) / 2.0;
          } else error->one(FLERR,"TEST"); // sanity check

          surface_force[0] = press * surface_normal_vector[0] * direction_flag[0];
          surface_force[1] = press * surface_normal_vector[1] * direction_flag[1];
          surface_force[2] = press * surface_normal_vector[2] * direction_flag[2];

          for (inode = 0; inode < num_node; inode++) {
            nodef[i][ibasis][surface_node[inode]][0] += surface_force[0];
            nodef[i][ibasis][surface_node[inode]][1] += surface_force[1];
            nodef[i][ibasis][surface_node[inode]][2] += surface_force[2];

            del[0] = inodex[surface_node[inode]][0] - ex[i][0];
            del[1] = inodex[surface_node[inode]][1] - ex[i][1];
            del[2] = inodex[surface_node[inode]][2] - ex[i][2];
            double moment[3];
            cross3(surface_force, del, moment);
            net_force[3] += moment[0];
            net_force[4] += moment[1];
            net_force[5] += moment[2];
          }
          net_force[0] += surface_force[0] * num_node;
          net_force[1] += surface_force[1] * num_node;
          net_force[2] += surface_force[2] * num_node;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------------------
   read PDOS values from file
   ----------------------------------------------------------------------------------- */

void FixTempPress::read_file(char *filename)
{
  int me = comm->me;
  FILE *fpdos;
  int i;
  char *dos;
  char line[MAXLINE];

  if (me == 0) {
    fpdos = force->open_potential(filename);
    if (fpdos == nullptr) {
      char str[128];
      sprintf(str, "Cannot open PDOS force file %s", filename);
      error->one(FLERR, str);
    }
  }

  if (me == 0) {
    fgets(line, MAXLINE, fpdos);
    fgets(line, MAXLINE, fpdos);
    sscanf(line, "%d", &nfreq);
  }

  MPI_Bcast(&nfreq, 1, MPI_INT, 0, world);

  memory->create(ffreq, nfreq, "temp/press:ffreq");
  memory->create(vfreq, nfreq, "temp/press:vfreq");

  double tmp1, tmp2;

  if (me == 0) {
    for (i = 0; i < nfreq; i++) {
      fgets(line, MAXLINE, fpdos);
      sscanf(line, "%lf %lf", &tmp1, &tmp2);
      vfreq[i] = tmp1;
      ffreq[i] = tmp2;
    }
  }

  MPI_Bcast(vfreq, nfreq, MPI_DOUBLE, 0, world);
  MPI_Bcast(ffreq, nfreq, MPI_DOUBLE, 0, world);

  if (me == 0) fclose(fpdos); 
}

/*  ----------------------------------------------------------------------
    set current t_target and t_sqrt
    -------------------------------------------------------------------------  */

void FixTempPress::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // if variable temp, evaluate variable, wrap with clear/add
  // reallocate tforce array if necessary

  if (style == CONSTANT) {
    t_target = t_start + delta * (t_stop-t_start);
  } else {
    modify->clearstep_compute();
    if (tstyle == EQUAL) {
      t_target = input->variable->compute_equal(tvar);
      if (t_target < 0.0)
        error->one(FLERR, "Fix langevin variable returned negative temperature");
    } else {
    }
  }
}

/*  ----------------------------------------------------------------------
    potential energy of added force
    -------------------------------------------------------------------------  */

double FixTempPress::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(net_force, net_force_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return net_force[6];
}

/*  ----------------------------------------------------------------------
    return components of total force on fix group before force was changed
    -------------------------------------------------------------------------  */

double FixTempPress::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(net_force, net_force_all, 6, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return net_force_all[n];
}

/*  ----------------------------------------------------------------------
    memory usage of local atom-based array
    -------------------------------------------------------------------------  */

double FixTempPress::memory_usage()
{
  double bytes = nfreq * 2 * sizeof(double);
  return bytes;
}
