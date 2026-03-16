#include <cmath>
#include <cstring>
#include "compute_displace_atom.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "input.h"
//#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------------- */

ComputeDisplaceAtom::ComputeDisplaceAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  atomdisplace(NULL), nodedisplace(NULL), id_fix(NULL)
{

  if (narg != 3) error->all(FLERR,"Illegal compute displace/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;
  create_attribute = 1;

  // optional args

  //refreshflag = 0;
  //rvar = NULL;

  //int iarg = 3;
  //while (iarg < narg) {
  //  if (strcmp(arg[iarg],"refresh") == 0) {
  //    if (iarg+2 > narg)
  //      error->all(FLERR,"Illegal compute displace/atom command");
  //    refreshflag = 1;
  //    delete [] rvar;
  //    int n = strlen(arg[iarg+1]) + 1;
  //    rvar = new char[n];
  //    strcpy(rvar,arg[iarg+1]);
  //    iarg += 2;
  //  } else error->all(FLERR,"Illegal compute displace/atom command");
  //}

  // error check

  //if (refreshflag) {
  //  ivar = input->variable->find(rvar);
  //  if (ivar < 0)
  //    error->all(FLERR,"Variable name for compute displace/atom does not exist");
  //  if (input->variable->atomstyle(ivar) == 0)
  //    error->all(FLERR,"Compute displace/atom variable "
  //               "is not atom-style variable");
  //}

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "node";
  newarg[5] = (char *) "1";
  newarg[6] = (char *) "3";
  modify->add_fix(7,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;
    double **x = atom->x;
    int *amask = atom->mask;
    imageint *aimage = atom->image;
    int nalocal = atom->nlocal;

    double ***nodexoriginal = fix->a3store;
    double ***nodex = element->nodex;
    int *emask = element->mask;
    imageint *eimage = element->image;
    int nelocal = element->nlocal;
    int npe = element->npe;

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) domain->unmap(x[i],aimage[i],xoriginal[i]);
      else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
    for (int i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) 
        for (int j = 0; j < npe; j++)
          domain->unmap(nodex[i][j],eimage[i],nodexoriginal[i][j]);
      else 
        for (int j = 0; j < npe; j++)
          nodexoriginal[i][j][0] = nodexoriginal[i][j][1] = nodexoriginal[i][j][2] = 0.0;
  }

  // per-atom displacement array

  namax = nemax = 0;
  //nvmax = 0;
  //varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDisplaceAtom::~ComputeDisplaceAtom()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  memory->destroy(atomdisplace);
  memory->destroy(nodedisplace);
  //delete [] rvar;
  //memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeDisplaceAtom::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute displace/atom fix ID");
  fix = (FixStore *) modify->fix[ifix];

  //if (refreshflag) {
  //  ivar = input->variable->find(rvar);
  //  if (ivar < 0)
  //    error->all(FLERR,"Variable name for compute displace/atom does not exist");
  //}
}

/* ---------------------------------------------------------------------- */

void ComputeDisplaceAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  int npe = element->npe;

  // grow local displacement arrays if necessary

  if (atom->nmax > namax) {
    memory->destroy(atomdisplace);
    namax = atom->nmax;
    memory->create(atomdisplace,namax,4,"displace/atom:atomdisplace");
    array_atom = atomdisplace;
  }

  if (element->nmax > nemax) {
    memory->destroy(nodedisplace);
    nemax = element->nmax;
    memory->create(nodedisplace,nemax,npe,4,"displace/atom:nodedisplace");
    array_node = nodedisplace;
  }

  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;
  double ***nodexoriginal = fix->a3store;

  double **x = atom->x;
  int *amask = atom->mask;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;

  double ***nodex = element->nodex;
  int *emask = element->mask;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int xbox,ybox,zbox;
  double dx,dy,dz;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        xbox = (aimage[i] & IMGMASK) - IMGMAX;
        ybox = (aimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (aimage[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - xoriginal[i][2];
        atomdisplace[i][0] = dx;
        atomdisplace[i][1] = dy;
        atomdisplace[i][2] = dz;
        atomdisplace[i][3] = sqrt(dx*dx + dy*dy + dz*dz);
      } else atomdisplace[i][0] = atomdisplace[i][1] =
        atomdisplace[i][2] = atomdisplace[i][3] = 0.0;

    for (int i = 0; i < nelocal; i++)
      if (emask[i] & groupbit) {
        xbox = (eimage[i] & IMGMASK) - IMGMAX;
        ybox = (eimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (eimage[i] >> IMG2BITS) - IMGMAX;
        for (int j = 0; j < npe; j++) {
          dx = nodex[i][j][0] + xbox*xprd - nodexoriginal[i][j][0];
          dy = nodex[i][j][1] + ybox*yprd - nodexoriginal[i][j][1];
          dz = nodex[i][j][2] + zbox*zprd - nodexoriginal[i][j][2];
          nodedisplace[i][j][0] = dx;
          nodedisplace[i][j][1] = dy;
          nodedisplace[i][j][2] = dz;
          nodedisplace[i][j][3] = sqrt(dx*dx + dy*dy + dz*dz);
        }
      } else 
        for (int j = 0; j < npe; j++) 
          nodedisplace[i][j][0] = nodedisplace[i][j][1] =
            nodedisplace[i][j][2] = nodedisplace[i][j][3] = 0.0;

  } else {
    for (int i = 0; i < nalocal; i++)
      if (amask[i] & groupbit) {
        xbox = (aimage[i] & IMGMASK) - IMGMAX;
        ybox = (aimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (aimage[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
        dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
        atomdisplace[i][0] = dx;
        atomdisplace[i][1] = dy;
        atomdisplace[i][2] = dz;
        atomdisplace[i][3] = sqrt(dx*dx + dy*dy + dz*dz);
      } else atomdisplace[i][0] = atomdisplace[i][1] =
        atomdisplace[i][2] = atomdisplace[i][3] = 0.0;

    for (int i = 0; i < nelocal; i++) 
      if (emask[i] & groupbit) {
        xbox = (eimage[i] & IMGMASK) - IMGMAX;
        ybox = (eimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (eimage[i] >> IMG2BITS) - IMGMAX;
        for (int j = 0; j < npe; j++) {
          dx = nodex[i][j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - nodexoriginal[i][j][0];
          dy = nodex[i][j][1] + h[1]*ybox + h[3]*zbox - nodexoriginal[i][j][1];
          dz = nodex[i][j][2] + h[2]*zbox - nodexoriginal[i][j][2];
          nodedisplace[i][j][0] = dx;
          nodedisplace[i][j][1] = dy;
          nodedisplace[i][j][2] = dz;
          nodedisplace[i][j][3] = sqrt(dx*dx + dy*dy + dz*dz);
        }
      } else 
        for (int j = 0; j < npe; j++) 
          nodedisplace[i][j][0] = nodedisplace[i][j][1] =
            nodedisplace[i][j][2] = nodedisplace[i][j][3] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
   ------------------------------------------------------------------------- */

void ComputeDisplaceAtom::set_atom_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/* ----------------------------------------------------------------------
   initialize one element's storage values, called when element is created
   ------------------------------------------------------------------------- */

void ComputeDisplaceAtom::set_elem_arrays(int i)
{
  double ***nodexoriginal = fix->a3store;
  double ***nodex = element->nodex;
  int npe = element->npe;
  for (int j = 0; j < npe; j++) {
    nodexoriginal[i][j][0] = nodex[i][j][0];
    nodexoriginal[i][j][1] = nodex[i][j][1];
    nodexoriginal[i][j][2] = nodex[i][j][2];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values passed on from another element, called when atom is created from an element
   ------------------------------------------------------------------------- */

void ComputeDisplaceAtom::pass_on_atom_arrays(int i, int i_old, int iintpl)
{
  double **xoriginal = fix->astore;
  double ***nodexoriginal = fix->a3store;
  element->evec->interpolate(xoriginal[i],nodexoriginal,i_old,iintpl,3);
}

/* ----------------------------------------------------------------------
   initialize one element's storage values passed on from another element, called when element is created from another element
   ------------------------------------------------------------------------- */

void ComputeDisplaceAtom::pass_on_elem_arrays(int i, int i_old, int *iintpl_list)
{
  double ***nodexoriginal = fix->a3store;
  int npe = element->npe;
  for (int j = 0; j < npe; j++) 
    element->evec->interpolate(nodexoriginal[i][j],nodexoriginal,i_old,iintpl_list[j],3);
}

/* ----------------------------------------------------------------------
   reset per-atom storage values, based on atom-style variable evaluation
   called by dump when dump_modify refresh is set
   ------------------------------------------------------------------------- */

//void ComputeDisplaceAtom::refresh()
//{
//  if (atom->namax > nvmax) {
//    nvmax = atom->namax;
//    memory->destroy(varatom);
//    memory->create(varatom,nvmax,"displace/atom:varatom");
//  }
//
//  input->variable->compute_atom(ivar,igroup,varatom,1,0);
//
//  double **xoriginal = fix->astore;
//  double **x = atom->x;
//  imageint *image = atom->image;
//  int nlocal = atom->nlocal;
//
//  for (int i = 0; i < nlocal; i++)
//    if (varatom[i]) domain->unmap(x[i],image[i],xoriginal[i]);
//}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeDisplaceAtom::memory_usage()
{
  double bytes = (namax + nemax*element->npe) * 4 * sizeof(double);
  //bytes += nvmax * sizeof(double);
  return bytes;
}
