#include "compute_centroid_stress_atom.h"
#include <cstring>
#include "atom.h"
#include "element.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
//#include "bond.h"
//#include "angle.h"
//#include "dihedral.h"
//#include "improper.h"
//#include "kspace.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

enum{NOBIAS, BIAS};

/* ---------------------------------------------------------------------- */

ComputeCentroidStressAtom::ComputeCentroidStressAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg), 
  id_temp(nullptr), atom_stress(nullptr), node_stress(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal compute centroid/stress/atom command");

  peratom_flag = 1;
  size_peratom_cols = 9;
  element_data_style = NODE;
  pressatomflag = 2;
  timeflag = 1;
  comm_atom_reverse = 9;
  comm_elem_reverse = 9 * element->maxnpe * element->maxapc;

  // store temperature ID used by stress computation
  // insure it is valid for temperature computation

  if (strcmp(arg[3], "NULL") == 0) id_temp = nullptr;
  else {

    error->all(FLERR, "Temperature calculation is not available yet");
    /*
    int n = strlen(arg[3]) + 1;
    id_temp = new char[n];
    strcpy(id_temp, arg[3]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR, "Could not find compute centroid/stress/atom temperature ID");
    if (modify->compute[icompute]->tempflag == 0)
      error->all(FLERR, 
                 "Compute centroid/stress/atom temperature ID does not "
                 "compute temperature");
    */
  }

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = 1;
    fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = 0;
    fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg], "pair") == 0) pairflag = 1;
      //else if (strcmp(arg[iarg], "bond") == 0) bondflag = 1;
      //else if (strcmp(arg[iarg], "angle") == 0) angleflag = 1;
      //else if (strcmp(arg[iarg], "dihedral") == 0) dihedralflag = 1;
      //else if (strcmp(arg[iarg], "improper") == 0) improperflag = 1;
      //else if (strcmp(arg[iarg], "kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg], "fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg], "virial") == 0) {
        pairflag = 1;
        fixflag = 1;
        //bondflag = angleflag = dihedralflag = improperflag = 1;
        //kspaceflag = fixflag = 1;
      } else error->all(FLERR, "Illegal compute centroid/stress/atom command");
      iarg++;
    }
  }

  atom_nmax = 0;
  elem_nmax = 0;
  maxnpe = 0;
  maxapc = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCentroidStressAtom::~ComputeCentroidStressAtom()
{
  delete [] id_temp;
  memory->destroy(atom_stress);
  memory->destroy(node_stress);
}

/* ---------------------------------------------------------------------- */

void ComputeCentroidStressAtom::init()
{
  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (id_temp) {
    error->all(FLERR, "Temperature calculation is not available yet");
    //int icompute = modify->find_compute(id_temp);
    //if (icompute < 0)
    //  error->all(FLERR, "Could not find compute centroid/stress/atom temperature ID");
    //temperature = modify->compute[icompute];
    //if (temperature->tempbias) biasflag = BIAS;
    //else biasflag = NOBIAS;
  } else biasflag = NOBIAS;

  // check if pair styles support centroid atom stress
  if (pairflag && force->pair)
    if (force->pair->centroidstressflag == CENTROID_NOTAVAIL)
      error->all(FLERR, "Pair style does not support compute centroid/stress/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeCentroidStressAtom::compute_peratom()
{
  double onemass;

  invoked_peratom = update->ntimestep;
  if (update->vflag_atom != invoked_peratom)
    error->all(FLERR, "Per-atom virial was not tallied on needed timestep");

  // grow local stress array if necessary
  // needs to be atom->nmax and element->nmax in length

  if (atom->nmax > atom_nmax) {
    memory->destroy(atom_stress);
    atom_nmax = atom->nmax;
    array_atom = memory->create(atom_stress, atom_nmax, 9, "centroid/stress/atom:atom_stress");
  }

  if (element->nmax > elem_nmax || maxnpe < element->maxnpe || maxapc < element->maxapc) {
    memory->destroy(node_stress);
    elem_nmax = element->nmax;
    maxnpe = element->maxnpe;
    maxapc = element->maxapc;
    array_node = memory->create(node_stress, elem_nmax, maxapc, 
        maxnpe, 9, "centroid/stress/atom:atom_stress");
  }

  comm_elem_reverse = 9 * maxnpe * maxapc;
  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  /*
  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton_bond) nbond += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;
  */

  int atom_nlocal = atom->nlocal;
  int atom_ntotal = atom_nlocal;
  int elem_nlocal = element->nlocal;
  int elem_ntotal = elem_nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  if (force->newton) {
    atom_ntotal += atom->nghost;
    elem_ntotal += element->nghost;
  }


  // clear local stress array

  memset(&atom_stress[0][0], 0, atom_ntotal * 9 * sizeof(double));
  memset(&node_stress[0][0][0][0], 0, elem_ntotal * maxapc * 
      maxnpe * 9 * sizeof(double));

  // add in per-atom contributions from each force

  // per-atom virial and per-atom centroid virial are the same for two-body
  // pair styles are either CENTROID_SAME or CENTROID_AVAIL or CENTROID_NOTAVAIL
  if (pairflag && force->pair && force->pair->compute_flag) {
    if (force->pair->centroidstressflag == CENTROID_AVAIL) {
      double **cvatom = force->pair->cvatom;
      double ****cvnode = force->pair->cvnode;
      for (int i = 0; i < atom_ntotal; i++)
        for (int j = 0; j < 9; j++)
          atom_stress[i][j] += cvatom[i][j];
      for (int i = 0; i < elem_ntotal; i++) 
        for (int j = 0; j < apc[etype[i]]; j++) 
          for (int k = 0; k < npe[etype[i]]; k++) 
            for (int l = 0; l < 9; l++)
              node_stress[i][j][k][l] += cvnode[i][j][k][l];

    } else {
      double **vatom = force->pair->vatom;
      double ****vnode = force->pair->vnode;
      for (int i = 0; i < atom_ntotal; i++) {
        for (int j = 0; j < 6; j++)
          atom_stress[i][j] += vatom[i][j];
        for (int j = 6; j < 9; j++)
          atom_stress[i][j] += vatom[i][j-3];
      }
      for (int i = 0; i < elem_ntotal; i++) 
        for (int j = 0; j < apc[etype[i]]; j++) 
          for (int k = 0; k < npe[etype[i]]; k++) {
            for (int l = 0; l < 6; l++)
              node_stress[i][j][k][l] += vnode[i][j][k][l];
            for (int l = 6; l < 9; l++)
              node_stress[i][j][k][l] += vnode[i][j][k][l-3];
          }
    }
  }

  /*
  // per-atom virial and per-atom centroid virial are the same for bonds
  if (bondflag && force->bond) {
    double **vatom = force->bond->vatom;
    for (i = 0; i < nbond; i++) {
      for (j = 0; j < 6; j++)
        stress[i][j] += vatom[i][j];
      for (j = 6; j < 9; j++)
        stress[i][j] += vatom[i][j-3];
    }
  }

  if (angleflag && force->angle) {
    double **cvatom = force->angle->cvatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 9; j++)
        stress[i][j] += cvatom[i][j];
  }

  if (dihedralflag && force->dihedral) {
    double **cvatom = force->dihedral->cvatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 9; j++)
        stress[i][j] += cvatom[i][j];
  }

  if (improperflag && force->improper) {
    double **cvatom = force->improper->cvatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 9; j++)
        stress[i][j] += cvatom[i][j];
  }

  if (kspaceflag && force->kspace) {
    double **vatom = force->kspace->vatom;
    for (i = 0; i < nkspace; i++) {
      for (j = 0; j < 6; j++)
        stress[i][j] += vatom[i][j];
      for (j = 6; j < 9; j++)
        stress[i][j] += vatom[i][j-3];
    }
  }
  */
  // add in per-atom contributions from relevant fixes
  // skip if vatom = nullptr
  // possible during setup phase if fix has not initialized its vatom yet
  // e.g. fix ave/spatial defined before fix shake, 
  //   and fix ave/spatial uses a per-atom stress from this compute as input
  /*
  if (fixflag) {
    for (int ifix = 0; ifix < modify->nfix; ifix++)
      if (modify->fix[ifix]->virial_flag) {
        double **vatom = modify->fix[ifix]->vatom;
        if (vatom)
          for (i = 0; i < nlocal; i++) {
            for (j = 0; j < 6; j++)
              stress[i][j] += vatom[i][j];
            for (j = 6; j < 9; j++)
              stress[i][j] += vatom[i][j-3];
          }
      }
  }
  */
  // communicate ghost virials between neighbor procs

  //if (force->newton || (force->kspace && force->kspace->tip4pflag))
  if (force->newton)
    comm->reverse_comm_compute(this);

  // zero virial of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *atom_mask = atom->mask;
  int *elem_mask = element->mask;

  for (int i = 0; i < atom_nlocal; i++)
    if (!(atom_mask[i] & groupbit)) 
      for (int j = 0; j < 9; j++)
        atom_stress[i][j] = 0.0;


  for (int i = 0; i < elem_nlocal; i++) 
    if (!(elem_mask[i] & groupbit)) 
      for (int j = 0; j < apc[etype[i]]; j++) 
        for (int k = 0; k < npe[etype[i]]; k++) 
          for (int l = 0; l < 9; l++)
            node_stress[i][j][k][l] = 0.0;

  // include kinetic energy term for each atom in group
  // apply temperature bias is applicable
  // mvv2e converts mv^2 to energy
  /*
  if (keflag) {
    double **v = atom->v;
    double *mass = atom->mass;
    double *rmass = atom->rmass;
    int *type = atom->type;
    double mvv2e = force->mvv2e;

    if (biasflag == NOBIAS) {
      if (rmass) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            onemass = mvv2e * rmass[i];
            stress[i][0] += onemass*v[i][0]*v[i][0];
            stress[i][1] += onemass*v[i][1]*v[i][1];
            stress[i][2] += onemass*v[i][2]*v[i][2];
            stress[i][3] += onemass*v[i][0]*v[i][1];
            stress[i][4] += onemass*v[i][0]*v[i][2];
            stress[i][5] += onemass*v[i][1]*v[i][2];
            stress[i][6] += onemass*v[i][1]*v[i][0];
            stress[i][7] += onemass*v[i][2]*v[i][0];
            stress[i][8] += onemass*v[i][2]*v[i][1];
          }

      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            onemass = mvv2e * mass[type[i]];
            stress[i][0] += onemass*v[i][0]*v[i][0];
            stress[i][1] += onemass*v[i][1]*v[i][1];
            stress[i][2] += onemass*v[i][2]*v[i][2];
            stress[i][3] += onemass*v[i][0]*v[i][1];
            stress[i][4] += onemass*v[i][0]*v[i][2];
            stress[i][5] += onemass*v[i][1]*v[i][2];
            stress[i][6] += onemass*v[i][1]*v[i][0];
            stress[i][7] += onemass*v[i][2]*v[i][0];
            stress[i][8] += onemass*v[i][2]*v[i][1];
          }
      }

    } else {

      // invoke temperature if it hasn't been already
      // this insures bias factor is pre-computed

      if (keflag && temperature->invoked_scalar != update->ntimestep)
        temperature->compute_scalar();

      if (rmass) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            temperature->remove_bias(i, v[i]);
            onemass = mvv2e * rmass[i];
            stress[i][0] += onemass*v[i][0]*v[i][0];
            stress[i][1] += onemass*v[i][1]*v[i][1];
            stress[i][2] += onemass*v[i][2]*v[i][2];
            stress[i][3] += onemass*v[i][0]*v[i][1];
            stress[i][4] += onemass*v[i][0]*v[i][2];
            stress[i][5] += onemass*v[i][1]*v[i][2];
            stress[i][6] += onemass*v[i][1]*v[i][0];
            stress[i][7] += onemass*v[i][2]*v[i][0];
            stress[i][8] += onemass*v[i][2]*v[i][1];
            temperature->restore_bias(i, v[i]);
          }

      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            temperature->remove_bias(i, v[i]);
            onemass = mvv2e * mass[type[i]];
            stress[i][0] += onemass*v[i][0]*v[i][0];
            stress[i][1] += onemass*v[i][1]*v[i][1];
            stress[i][2] += onemass*v[i][2]*v[i][2];
            stress[i][3] += onemass*v[i][0]*v[i][1];
            stress[i][4] += onemass*v[i][0]*v[i][2];
            stress[i][5] += onemass*v[i][1]*v[i][2];
            stress[i][6] += onemass*v[i][1]*v[i][0];
            stress[i][7] += onemass*v[i][2]*v[i][0];
            stress[i][8] += onemass*v[i][2]*v[i][1];
            temperature->restore_bias(i, v[i]);
          }
      }
    }
  }
  */
  // convert to stress * volume units = -pressure * volume

  double nktv2p = -force->nktv2p;
  for (int i = 0; i < atom_nlocal; i++)
    if (atom_mask[i] & groupbit) 
      for (int j = 0; j < 9; j++)
        atom_stress[i][j] *= nktv2p;

  for (int i = 0; i < elem_nlocal; i++)
    if (elem_mask[i] & groupbit) 
      for (int j = 0; j < apc[etype[i]]; j++) 
        for (int k = 0; k < npe[etype[i]]; k++) 
          for (int l = 0; l < 9; l++)
            node_stress[i][j][k][l] *= nktv2p;

}

/*  ----------------------------------------------------------------------  */

int ComputeCentroidStressAtom::pack_atom_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = atom_stress[i][0];
    buf[m++] = atom_stress[i][1];
    buf[m++] = atom_stress[i][2];
    buf[m++] = atom_stress[i][3];
    buf[m++] = atom_stress[i][4];
    buf[m++] = atom_stress[i][5];
    buf[m++] = atom_stress[i][6];
    buf[m++] = atom_stress[i][7];
    buf[m++] = atom_stress[i][8];
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

void ComputeCentroidStressAtom::unpack_atom_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    atom_stress[j][0] += buf[m++];
    atom_stress[j][1] += buf[m++];
    atom_stress[j][2] += buf[m++];
    atom_stress[j][3] += buf[m++];
    atom_stress[j][4] += buf[m++];
    atom_stress[j][5] += buf[m++];
    atom_stress[j][6] += buf[m++];
    atom_stress[j][7] += buf[m++];
    atom_stress[j][8] += buf[m++];
  }
}

/*  ----------------------------------------------------------------------  */

int ComputeCentroidStressAtom::pack_elem_reverse_comm(int n, int first, double *buf)
{
  int i, j, k, m, last;
  int *etype = element->etype;
  int *npe = element->npe;
  int *apc = element->apc;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) {
        buf[m++] = node_stress[i][j][k][0];
        buf[m++] = node_stress[i][j][k][1];
        buf[m++] = node_stress[i][j][k][2];
        buf[m++] = node_stress[i][j][k][3];
        buf[m++] = node_stress[i][j][k][4];
        buf[m++] = node_stress[i][j][k][5];
        buf[m++] = node_stress[i][j][k][6];
        buf[m++] = node_stress[i][j][k][7];
        buf[m++] = node_stress[i][j][k][8];
      }
  }  
  return m;
}

/*  ----------------------------------------------------------------------  */

void ComputeCentroidStressAtom::unpack_elem_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, l, m;
  int *etype = element->etype;
  int *npe = element->npe;
  int *apc = element->apc;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < apc[etype[j]]; k++) 
      for (l = 0; l < npe[etype[j]]; l++) {
        node_stress[j][k][l][0] += buf[m++];
        node_stress[j][k][l][1] += buf[m++];
        node_stress[j][k][l][2] += buf[m++];
        node_stress[j][k][l][3] += buf[m++];
        node_stress[j][k][l][4] += buf[m++];
        node_stress[j][k][l][5] += buf[m++];
        node_stress[j][k][l][6] += buf[m++];
        node_stress[j][k][l][7] += buf[m++];
        node_stress[j][k][l][8] += buf[m++];
      }
  }
}
/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeCentroidStressAtom::memory_usage()
{
  double bytes = atom_nmax * 9 * sizeof(double);
  bytes += elem_nmax * 9 * maxapc * maxnpe * sizeof(double);
  return bytes;
}
