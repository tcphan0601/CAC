#include <stdlib.h>
#include <string.h>
#include "compute_stress_atom.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

ComputeStressAtom::ComputeStressAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  id_temp(NULL), atom_stress(NULL), node_stress(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute stress/atom command");

  peratom_flag = 1;
  size_peratom_cols = 6;
  pressatomflag = 1;
  timeflag = 1;
  npe = element->npe;
  comm_atom_reverse = 6;
  comm_elem_reverse = 6*npe;


  // store temperature ID used by stress computation
  // insure it is valid for temperature computation

  if (strcmp(arg[3],"NULL") == 0) id_temp = NULL;
  else {
    //int n = strlen(arg[3]) + 1;
    //id_temp = new char[n];
    //strcpy(id_temp,arg[3]);
    error->all(FLERR,"Temperature calculation is not available yet");
    //int icompute = modify->find_compute(id_temp);
    //if (icompute < 0)
    //  error->all(FLERR,"Could not find compute stress/atom temperature ID");
    //if (modify->compute[icompute]->tempflag == 0)
    //  error->all(FLERR,
		// "Compute stress/atom temperature ID does not "
    //             "compute temperature");
  }

  // process optional args

  if (narg == 4) {
    keflag = 0;
    pairflag = 1;
    fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      //if (strcmp(arg[iarg],"ke") == 0) keflag = 1;
      if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        fixflag = 1;
      } else error->all(FLERR,"Illegal compute stress/atom command");
      iarg++;
    }
  }

  atom_nmax = 0;
  elem_nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStressAtom::~ComputeStressAtom()
{
  delete [] id_temp;
  memory->destroy(atom_stress);
  memory->destroy(node_stress);
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::init()
{
  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (id_temp) {
    error->all(FLERR,"Temperature calculation is not available yet");
    //int icompute = modify->find_compute(id_temp);
    //if (icompute < 0)
    //  error->all(FLERR,"Could not find compute stress/atom temperature ID");
    //temperature = modify->compute[icompute];
    //if (temperature->tempbias) biasflag = BIAS;
    //else biasflag = NOBIAS;
  } else biasflag = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::compute_peratom()
{
  int i,j;
  double onemass;

  invoked_peratom = update->ntimestep;
  if (update->vflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom virial was not tallied on needed timestep");

  // grow local stress arrays if necessary
  // needs to be atom->nmax and element->nmax in length

  if (atom->nmax > atom_nmax) {
    memory->destroy(atom_stress);
    atom_nmax = atom->nmax;
    array_atom = memory->create(atom_stress,atom_nmax,6,"stress/atom:atom_stress");
  }

  if (element->nmax > elem_nmax) {
    memory->destroy(node_stress);
    elem_nmax = element->nmax;
    array_node = memory->create(node_stress,elem_nmax,npe,6,"stress/atom:atom_stress");
  }

  // npair includes ghosts if either newton flag is set

  int atom_nlocal = atom->nlocal;
  int atom_ntotal = atom_nlocal;
  int elem_nlocal = element->nlocal;
  int elem_ntotal = elem_nlocal;

  if (force->newton) {
    atom_ntotal += atom->nghost;
    elem_ntotal += element->nghost;
  }

  // clear local stress arrays

  for (i = 0; i < atom_ntotal; i++) {
      atom_stress[i][0] = 0.0;
      atom_stress[i][1] = 0.0;
      atom_stress[i][2] = 0.0;
      atom_stress[i][3] = 0.0;
      atom_stress[i][4] = 0.0;
      atom_stress[i][5] = 0.0;
    }

  for (i = 0; i < elem_ntotal; i++)
    for (j = 0; j < npe; j++) {
      node_stress[i][j][0] = 0.0;
      node_stress[i][j][1] = 0.0;
      node_stress[i][j][2] = 0.0;
      node_stress[i][j][3] = 0.0;
      node_stress[i][j][4] = 0.0;
      node_stress[i][j][5] = 0.0;
    }

  // add in per-atom and per-node contributions from each force

  if (pairflag && force->pair) {
    double **vatom = force->pair->vatom;
    for (i = 0; i < atom_ntotal; i++) {
      atom_stress[i][0] += vatom[i][0];
      atom_stress[i][1] += vatom[i][1];
      atom_stress[i][2] += vatom[i][2];
      atom_stress[i][3] += vatom[i][3];
      atom_stress[i][4] += vatom[i][4];
      atom_stress[i][5] += vatom[i][5];
    }
    double ***vnode = force->pair->vnode;
    for (i = 0; i < elem_ntotal; i++)
      for (j = 0; j < npe; j++) {
        node_stress[i][j][0] += vnode[i][j][0];
        node_stress[i][j][1] += vnode[i][j][1];
        node_stress[i][j][2] += vnode[i][j][2];
        node_stress[i][j][3] += vnode[i][j][3];
        node_stress[i][j][4] += vnode[i][j][4];
        node_stress[i][j][5] += vnode[i][j][5];
      }
  }

  // add in per-atom contributions from relevant fixes
  // skip if vatom = NULL
  // possible during setup phase if fix has not initialized its vatom yet
  // e.g. fix ave/spatial defined before fix shake,
  //   and fix ave/spatial uses a per-atom stress from this compute as input

  //if (fixflag) {
  //  for (int ifix = 0; ifix < modify->nfix; ifix++)
  //    if (modify->fix[ifix]->virial_flag) {
  //      double **vatom = modify->fix[ifix]->vatom;
  //      if (vatom)
  //        for (i = 0; i < nlocal; i++)
  //          for (j = 0; j < 6; j++)
  //            atom_stress[i][j] += vatom[i][j];
  //    }
  //}

  // communicate ghost virials between neighbor procs

  if (force->newton) comm->reverse_comm_compute(this);

  // zero virial of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *atom_mask = atom->mask;
  int *elem_mask = element->mask;

  for (i = 0; i < atom_nlocal; i++)
    if (!(atom_mask[i] & groupbit)) {
      atom_stress[i][0] = 0.0;
      atom_stress[i][1] = 0.0;
      atom_stress[i][2] = 0.0;
      atom_stress[i][3] = 0.0;
      atom_stress[i][4] = 0.0;
      atom_stress[i][5] = 0.0;
    }

  for (i = 0; i < elem_nlocal; i++)
    if (!(elem_mask[i] & groupbit)) 
      for (j = 0; j < npe; i++) {
        node_stress[i][j][0] = 0.0;
        node_stress[i][j][1] = 0.0;
        node_stress[i][j][2] = 0.0;
        node_stress[i][j][3] = 0.0;
        node_stress[i][j][4] = 0.0;
        node_stress[i][j][5] = 0.0;
      } 

  // CAC: kinetic energy term is not included
  // include kinetic energy term for each atom in group
  // apply temperature bias is applicable
  // mvv2e converts mv^2 to energy

  //if (keflag) {
  //  double **v = atom->v;
  //  double *mass = atom->mass;
  //  double *rmass = atom->rmass;
  //  int *type = atom->type;
  //  double mvv2e = force->mvv2e;

  //  if (biasflag == NOBIAS) {
  //    if (rmass) {
  //      for (i = 0; i < nlocal; i++)
  //        if (mask[i] & groupbit) {
  //          onemass = mvv2e * rmass[i];
  //          atom_stress[i][0] += onemass*v[i][0]*v[i][0];
  //          atom_stress[i][1] += onemass*v[i][1]*v[i][1];
  //          atom_stress[i][2] += onemass*v[i][2]*v[i][2];
  //          atom_stress[i][3] += onemass*v[i][0]*v[i][1];
  //          atom_stress[i][4] += onemass*v[i][0]*v[i][2];
  //          atom_stress[i][5] += onemass*v[i][1]*v[i][2];
  //        }

  //    } else {
  //      for (i = 0; i < nlocal; i++)
  //        if (mask[i] & groupbit) {
  //          onemass = mvv2e * mass[type[i]];
  //          atom_stress[i][0] += onemass*v[i][0]*v[i][0];
  //          atom_stress[i][1] += onemass*v[i][1]*v[i][1];
  //          atom_stress[i][2] += onemass*v[i][2]*v[i][2];
  //          atom_stress[i][3] += onemass*v[i][0]*v[i][1];
  //          atom_stress[i][4] += onemass*v[i][0]*v[i][2];
  //          atom_stress[i][5] += onemass*v[i][1]*v[i][2];
  //        }
  //    }

  //  } else {

  // invoke temperature if it hasn't been already
  // this insures bias factor is pre-computed

  //    if (keflag && temperature->invoked_scalar != update->ntimestep)
  //      temperature->compute_scalar();

  //    if (rmass) {
  //      for (i = 0; i < nlocal; i++)
  //        if (mask[i] & groupbit) {
  //          temperature->remove_bias(i,v[i]);
  //          onemass = mvv2e * rmass[i];
  //          atom_stress[i][0] += onemass*v[i][0]*v[i][0];
  //          atom_stress[i][1] += onemass*v[i][1]*v[i][1];
  //          atom_stress[i][2] += onemass*v[i][2]*v[i][2];
  //          atom_stress[i][3] += onemass*v[i][0]*v[i][1];
  //          atom_stress[i][4] += onemass*v[i][0]*v[i][2];
  //          atom_stress[i][5] += onemass*v[i][1]*v[i][2];
  //          temperature->restore_bias(i,v[i]);
  //        }

  //    } else {
  //      for (i = 0; i < nlocal; i++)
  //        if (mask[i] & groupbit) {
  //          temperature->remove_bias(i,v[i]);
  //          onemass = mvv2e * mass[type[i]];
  //          atom_stress[i][0] += onemass*v[i][0]*v[i][0];
  //          atom_stress[i][1] += onemass*v[i][1]*v[i][1];
  //          atom_stress[i][2] += onemass*v[i][2]*v[i][2];
  //          atom_stress[i][3] += onemass*v[i][0]*v[i][1];
  //          atom_stress[i][4] += onemass*v[i][0]*v[i][2];
  //          atom_stress[i][5] += onemass*v[i][1]*v[i][2];
  //          temperature->restore_bias(i,v[i]);
  //        }
  //    }
  //  }
  //}

  // convert to stress*volume units = -pressure*volume

  double nktv2p = -force->nktv2p;
  for (i = 0; i < atom_nlocal; i++)
    if (atom_mask[i] & groupbit) {
      atom_stress[i][0] *= nktv2p;
      atom_stress[i][1] *= nktv2p;
      atom_stress[i][2] *= nktv2p;
      atom_stress[i][3] *= nktv2p;
      atom_stress[i][4] *= nktv2p;
      atom_stress[i][5] *= nktv2p;
    }

  for (i = 0; i < elem_nlocal; i++)
    if (elem_mask[i] & groupbit) 
      for (j = 0; j < npe; j++) {
        node_stress[i][j][0] *= nktv2p;
        node_stress[i][j][1] *= nktv2p;
        node_stress[i][j][2] *= nktv2p;
        node_stress[i][j][3] *= nktv2p;
        node_stress[i][j][4] *= nktv2p;
        node_stress[i][j][5] *= nktv2p;
      }

}

/* ---------------------------------------------------------------------- */

int ComputeStressAtom::pack_atom_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = atom_stress[i][0];
    buf[m++] = atom_stress[i][1];
    buf[m++] = atom_stress[i][2];
    buf[m++] = atom_stress[i][3];
    buf[m++] = atom_stress[i][4];
    buf[m++] = atom_stress[i][5];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::unpack_atom_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    atom_stress[j][0] += buf[m++];
    atom_stress[j][1] += buf[m++];
    atom_stress[j][2] += buf[m++];
    atom_stress[j][3] += buf[m++];
    atom_stress[j][4] += buf[m++];
    atom_stress[j][5] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int ComputeStressAtom::pack_elem_reverse_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (j = 0; j < npe; j++) {
      buf[m++] = node_stress[i][j][0];
      buf[m++] = node_stress[i][j][1];
      buf[m++] = node_stress[i][j][2];
      buf[m++] = node_stress[i][j][3];
      buf[m++] = node_stress[i][j][4];
      buf[m++] = node_stress[i][j][5];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::unpack_elem_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < npe; k++) { 
      node_stress[j][k][0] += buf[m++];
      node_stress[j][k][1] += buf[m++];
      node_stress[j][k][2] += buf[m++];
      node_stress[j][k][3] += buf[m++];
      node_stress[j][k][4] += buf[m++];
      node_stress[j][k][5] += buf[m++];
    }
  }
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based and node-based arrays
   ------------------------------------------------------------------------- */

double ComputeStressAtom::memory_usage()
{
  double bytes = atom_nmax * 6 * sizeof(double);
  bytes += elem_nmax * 6 * npe * sizeof(double);
  return bytes;
}
