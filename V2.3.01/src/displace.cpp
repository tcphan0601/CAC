#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "displace.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "irregular_comm.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "math_const.h"
#include "universe.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include "variable.h"

using namespace CAC_NS;
using namespace MathConst;

enum{MOVE,RAMP,RANDOM,ROTATE,DISLOCATION};
enum{NODE,ELEMENT};

/* ---------------------------------------------------------------------- */

Displace::Displace(CAC *cac) : Pointers(cac)
{
  mvec = NULL;
}

/* ---------------------------------------------------------------------- */

Displace::~Displace()
{
  memory->destroy(mvec);
}

/* ---------------------------------------------------------------------- */

void Displace::command(int narg, char **arg)
{
  int i,j;

  if (domain->box_exist == 0)
    error->all(FLERR,"Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace command");
  /* no read/write restart yet
     if (modify->nfix_restart_peratom)
     error->all(FLERR,"Cannot displace after "
     "reading restart file with per-atom info");
     */
  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");
  if (comm->me == 0 && logfile) fprintf(logfile,"Displacing atoms ...\n");

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace group ID");
  groupbit = group->bitmask[igroup];

  /* no fix rigid yet
     if (modify->check_rigid_group_overlap(groupbit))
     error->warning(FLERR,"Attempting to displace atoms in rigid bodies");
     */
  int style = -1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  //  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  //  else if (strcmp(arg[1],"rotate") == 0) style = ROTATE;
  else if (strcmp(arg[1],"dislocation") == 0) style = DISLOCATION;
  else error->all(FLERR,"Illegal displace command");

  // set option defaults

  scaleflag = 0;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
//  else if (style == RANDOM) options(narg-6,&arg[6]);
//  else if (style == ROTATE) options(narg-9,&arg[9]);
  else if (style == DISLOCATION) options(narg-7,&arg[7]);

  // setup scaling


  if (scaleflag) {
    xscale[0] = domain->xprd;
    xscale[1] = domain->yprd;
    xscale[2] = domain->zprd;
  }
  else MathExtra::set3(xscale,1.0);

  // move atoms by 3-vector or specified variable(s)

  if (style == MOVE) 
    for (int idim = 0; idim < 3; idim++) move(idim,arg[idim+2]);

  // move atoms in ramped fashion

  if (style == RAMP) ramp(&arg[2]);

  // displace atom according to theoretical displacement field of dislocation

  if (style == DISLOCATION) dislocation(&arg[2]); 

  double **ax = atom->x;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

  double **ex = element->x;
  double ***nodex = element->nodex;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  for (i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i]);


  if (domain->triclinic) {
    domain->x2lambda(atom->nlocal,atom->x);
    domain->x2lambda(element->nlocal,element->x);
  }
  domain->reset_box();
  IrregularComm *irregular_comm = new IrregularComm(cac);
  irregular_comm->migrate(1);
  delete irregular_comm;
  if (domain->triclinic) {
    domain->lambda2x(atom->nlocal,atom->x);
    domain->lambda2x(element->nlocal,element->x);
  }

  // check if any atoms were lost

  bigint natoms,nelements,nblocal;

  nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0) {
    char str1[128];
    sprintf(str1,"Lost atoms via displace: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->warning(FLERR,str1);
  }

  nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && comm->me == 0) {
    char str2[128];
    sprintf(str2,"Lost elements via displace: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,element->nelements,nelements);
    error->warning(FLERR,str2);
  }

}

/* ----------------------------------------------------------------------
   move atoms/elements/nodes in ramped fashion   
   ------------------------------------------------------------------------- */

void Displace::ramp(char **arg){
  int d_dim,i,j;
  if (strcmp(arg[0],"x") == 0) d_dim = 0;
  else if (strcmp(arg[0],"y") == 0) d_dim = 1;
  else if (strcmp(arg[0],"z") == 0) d_dim = 2;
  else error->all(FLERR,"Illegal displace ramp command");

  double d_lo,d_hi;
  d_lo = xscale[d_dim]*universe->numeric(FLERR,arg[1]);
  d_hi = xscale[d_dim]*universe->numeric(FLERR,arg[2]);

  int coord_dim;
  if (strcmp(arg[3],"x") == 0) coord_dim = 0;
  else if (strcmp(arg[3],"y") == 0) coord_dim = 1;
  else if (strcmp(arg[3],"z") == 0) coord_dim = 2;
  else error->all(FLERR,"Illegal displace ramp command");

  double coord_lo,coord_hi;
  coord_lo = xscale[coord_dim]*universe->numeric(FLERR,arg[4]);
  coord_hi = xscale[coord_dim]*universe->numeric(FLERR,arg[5]);

  double **ax = atom->x;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int npe = element->npe;

  double fraction,dramp;

  if (nalocal){
    for (i = 0; i < nalocal; i++) 
      if (amask[i] & groupbit) ramp_coord(coord_lo,coord_hi,d_lo,d_hi,ax[i],coord_dim,d_dim);
  }

  if (nelocal){
    for (i = 0; i < nelocal; i++) {
      if (group_flag){
        if (emask[i] & groupbit)
          for (j = 0; j < npe; j++) ramp_coord(coord_lo,coord_hi,d_lo,d_hi,nodex[i][j],coord_dim,d_dim);
      } else {
        for (j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit)
            ramp_coord(coord_lo,coord_hi,d_lo,d_hi,nodex[i][j],coord_dim,d_dim);
      }
    }

    element->evec->update_center_coord();
  }

}


/* ----------------------------------------------------------------------
   change coord according to ramped fashion
   ------------------------------------------------------------------------- */

void Displace::ramp_coord(double coord_lo, double coord_hi, double d_lo, double d_hi, double *coord, int coord_dim, int d_dim){
  double fraction = (coord[coord_dim] - coord_lo) / (coord_hi - coord_lo);
  fraction = MAX(fraction,0.0);
  fraction = MIN(fraction,1.0);
  double dramp = d_lo + fraction*(d_hi - d_lo);
  coord[d_dim] += dramp;
}

/* ----------------------------------------------------------------------
   move atoms/elements/nodes either by specified numeric displacement or variable evaluation
   ------------------------------------------------------------------------- */

void Displace::move(int idim, char *arg)
{
  if (strstr(arg,"v_") != arg) {
    double delta = xscale[idim] * universe->numeric(FLERR,arg);
    move(idim,delta);
  } else {
    int ivar = input->variable->find(arg+2);
    if (ivar < 0)
      error->all(FLERR,"Variable name for displace does not exist");

    if (input->variable->equalstyle(ivar)) {
      double delta = xscale[idim] * input->variable->compute_equal(ivar);
      move(idim,delta);
    } else if (input->variable->atomstyle(ivar)) {
      error->all(FLERR,"Atom style variable not working yet");
      //if (mvec == NULL) memory->create(mvec,nlocal,"displace:mvec");
      //input->variable->compute_atom(ivar,igroup,mvec,1,0);
      //for (int i = 0; i < nlocal; i++)
      //  if (mask[i] & groupbit) x[i][idim] += scale*mvec[i];
    } else error->all(FLERR,"Variable for displace is invalid style");
  }
}
/* ----------------------------------------------------------------------
   parse optional parameters at end of displace input line
   ------------------------------------------------------------------------- */

void Displace::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix setforce command");
      if (strcmp(arg[iarg+1],"node") == 0) group_flag = NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = ELEMENT;
      else error->all(FLERR, "Illegal fix setforce command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace command");
  }
}

/* ----------------------------------------------------------------------
   move atoms/elements according to theoretical dislocation displacement field 
   ------------------------------------------------------------------------- */


void Displace::dislocation(char **arg){

  int d_dim,b_dim,i,j,dim_flag[3];
  dim_flag[0] = dim_flag[1] = dim_flag[2] = 1;

  // dislocation line direction
  if (strcmp(arg[0],"x") == 0) {
    d_dim = 0;
    dim_flag[0] = 0;
  } else if (strcmp(arg[0],"y") == 0) {
    d_dim = 1;
    dim_flag[1] = 0;
  } else if (strcmp(arg[0],"z") == 0) {
    d_dim = 2;
    dim_flag[2] = 0;
  } else error->all(FLERR,"Illegal displace dislocation command");

  // burger vector direction
  if (strcmp(arg[3],"x") == 0) {
    b_dim = 0;
    dim_flag[0] = 0;
  } else if (strcmp(arg[3],"y") == 0) {
    b_dim = 1;
    dim_flag[1] = 0;
  } else if (strcmp(arg[3],"z") == 0) {
    b_dim = 2;
    dim_flag[2] = 0;
  } else error->all(FLERR,"Illegal displace dislocation command");

  if (d_dim == b_dim) 
    error->all(FLERR,"Dislocation line direction must be different from Burger vector direction");


  int dim1 = b_dim;
  int dim2 = dim_flag[0] * 0 + dim_flag[1] * 1 + dim_flag[2] * 2;
  double center1,center2;
  double b; // Burgers vector
  double v; // Poisson's ratio
  double u1,u2,x1,x2; 

  if (dim1 == (d_dim + 1)%3){
    center1 = universe->numeric(FLERR,arg[1]);
    center2 = universe->numeric(FLERR,arg[2]);
  } else {
    center1 = universe->numeric(FLERR,arg[2]);
    center2 = universe->numeric(FLERR,arg[1]);
  }

  b = universe->numeric(FLERR,arg[4]);
  v = universe->numeric(FLERR,arg[5]);

  double **ax = atom->x;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int npe = element->npe;

  // displace atoms
  if (nalocal) 
    for (i = 0; i < nalocal; i++){
      if (amask[i] & groupbit){
        x1 = ax[i][dim1] - center1;
        x2 = ax[i][dim2] - center2;
        dislc_field(&u1,&u2,x1,x2,b,v);
        ax[i][dim1] += u1;
        ax[i][dim2] += u2;
      }
    }

  // displace elements
  if (nelocal){
    for (i = 0; i < nelocal; i++){
      if (group_flag) {
        if (emask[i] & groupbit)
          for (j = 0; j < npe; j++) {
            x1 = nodex[i][j][dim1] - center1;
            x2 = nodex[i][j][dim2] - center2;
            dislc_field(&u1,&u2,x1,x2,b,v);
            nodex[i][j][dim1] += u1;
            nodex[i][j][dim2] += u2;
          }
      } else {
        for (j = 0; j < npe; j++) 
          if (nodemask[i][j] & groupbit){
            x1 = nodex[i][j][dim1] - center1;
            x2 = nodex[i][j][dim2] - center2;
            dislc_field(&u1,&u2,x1,x2,b,v);
            nodex[i][j][dim1] += u1;
            nodex[i][j][dim2] += u2;
          }
      }
    }
    element->evec->update_center_coord();
  }
}

/* ----------------------------------------------------------------------
   calculate theoretical dislocation displacement field 
   ------------------------------------------------------------------------- */

void Displace::dislc_field(double *u1, double *u2, double x1, double x2, double b, double v){
  double rsq = x1*x1 + x2*x2;

  // for u1: need range of inverse tan to be from 0 to 2*pi
  // range of atan2 is from -pi to +pi so add a pi to offset

  *u1 = b/MY_2PI*((atan2(x2,x1) + MY_PI) + x1*x2/(2.0*(1-v)*rsq));
  *u2 = -b/MY_8PI/(1-v)*((1-2*v)*log(rsq) + (x1*x1-x2*x2)/rsq);
}

/* ----------------------------------------------------------------------
   move all atoms/elements/nodes in group by delta along dimension dim
   ------------------------------------------------------------------------- */

void Displace::move(int dim, double delta)
{
  double **x = atom->x;
  int *mask = atom->mask;
  for (int i = 0; i < atom->nlocal; i++)
    if (mask[i] & groupbit) x[i][dim] += delta;

  double ***nodex = element->nodex;
  mask = element->mask;
  int **nodemask = element->nodemask;
  for (int i = 0; i < element->nlocal; i++){
    if (group_flag) {
      if (mask[i] & groupbit)
        for (int j = 0; j < element->npe; j++) nodex[i][j][dim] += delta;
    } else {
      for (int j = 0; j < element->npe; j++)
        if (nodemask[i][j] & groupbit) nodex[i][j][dim] += delta;
    }
  }
  if (element->nlocal) element->evec->update_center_coord();
}
