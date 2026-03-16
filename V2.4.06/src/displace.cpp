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
#include "lattice.h"
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

#define   EPSILON 1.0e-4

enum{MOVE,RAMP,RANDOM,ROTATE,MIRROR,DISLOCATION,MULTIDISL};

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
  if (domain->box_exist == 0)
    error->all(FLERR,"Displace command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace command");
  /* no read/write restart yet
     if (modify->nfix_restart_peratom)
     error->all(FLERR,"Cannot displace after "
     "reading restart file with per-atom info");
     */
  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms/elements ...\n");
  if (comm->me == 0 && logfile) fprintf(logfile,"Displacing atoms/elements ...\n");

  /* no fix rigid yet
     if (modify->check_rigid_group_overlap(groupbit))
     error->warning(FLERR,"Attempting to displace atoms in rigid bodies");
     */
  style = -1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  //  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  else if (strcmp(arg[1],"rotate") == 0) style = ROTATE;
  else if (strcmp(arg[1],"mirror") == 0) style = MIRROR;
  else if (strcmp(arg[1],"dislocation") == 0) style = DISLOCATION;
  else if (strcmp(arg[1],"multidisl") == 0) style = MULTIDISL;
  else error->all(FLERR,"Illegal displace command");

  // group and style
  // ignore if dislocation/multidisl style
  
  if (style != DISLOCATION && style != MULTIDISL) {
    igroup = group->find(arg[0]);
    if (igroup == -1) error->all(FLERR,"Could not find displace group ID");
    groupbit = group->bitmask[igroup];
  }

  // set option defaults

  scaleflag = 0;
  group_flag = Group::ELEMENT;
  sequence_flag = 0;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  //  else if (style == RANDOM) options(narg-6,&arg[6]);
  else if (style == ROTATE) options(narg-9,&arg[9]);
  else if (style == MIRROR) options(narg-4,&arg[4]);
  else if (style == DISLOCATION) options(narg-9,&arg[9]);
  else if (style == MULTIDISL) options(narg-12,&arg[12]);

  // setup scaling

  if (scaleflag) {
    xscale[0] = domain->lattice->xlattice;
    xscale[1] = domain->lattice->ylattice;
    xscale[2] = domain->lattice->zlattice;
  }
  else xscale[0] = xscale[1] = xscale[2] = 1.0;

  // move atoms/elements by 3-vector or specified variable(s)

  if (style == MOVE) 
    for (int idim = 0; idim < 3; idim++) move(idim,arg[idim+2]);

  // move atoms/elements in ramped fashion

  else if (style == RAMP) ramp(&arg[2]);

  // rotate atoms/elements by right-hand rule by theta around R

  else if (style == ROTATE) rotate(&arg[2]);

  // reflect atoms/elements accross a mirror
  
  else if (style == MIRROR) mirror(&arg[2]);

  // displace atoms/elements according to theoretical displacement field of dislocation
  // create 1 dislocation

  else if (style == DISLOCATION) dislocation(&arg[2]); 

  // displace atoms/elements according to theoretical displacement field of dislocation
  // create multiple dislocation

  else if (style == MULTIDISL) multi_dislocations(&arg[2]); 

  // move atoms/elements back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms/elements moved a long distance
  // use irregular() in case atoms/elements moved a long distance

  double **ax = atom->x;
  imageint *aimage = atom->image;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

  double **ex = element->x;
  double ***nodex = element->nodex;
  imageint *eimage = element->image;
  int nelocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i],npe[etype[i]]);

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->reset_box();
  IrregularComm *irregular_comm = new IrregularComm(cac);
  irregular_comm->migrate(1);
  delete irregular_comm;
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2nodex(element->nlocal,element->nodex);
  }

  // check if any atoms/elements were lost

  bigint natoms,nelements,nblocal;

  nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0) {
    char str[128];
    sprintf(str,"Lost atoms via displace: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->warning(FLERR,str);
  }

  nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && comm->me == 0) {
    char str[128];
    sprintf(str,"Lost elements via displace: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,element->nelements,nelements);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   move atoms/elements/nodes in ramped fashion   
   ------------------------------------------------------------------------- */

void Displace::ramp(char **arg)
{
  int d_dim,i,j;

  // displace direction 

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
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int *npe = element->npe;

  double fraction,dramp;

  for (i = 0; i < nalocal; i++) 
    if (amask[i] & groupbit) ramp_coord(coord_lo,coord_hi,d_lo,d_hi,ax[i],coord_dim,d_dim);

  if (group_flag == Group::ATOM) return;

  for (i = 0; i < nelocal; i++) 
    if (group_flag == Group::ELEMENT) {
      if (emask[i] & groupbit)
        for (j = 0; j < npe[etype[i]]; j++) ramp_coord(coord_lo,coord_hi,d_lo,d_hi,nodex[i][j],coord_dim,d_dim);
    } else if (group_flag == Group::NODE) {
      for (j = 0; j < npe[etype[i]]; j++) 
        if (nodemask[i][j] & groupbit)
          ramp_coord(coord_lo,coord_hi,d_lo,d_hi,nodex[i][j],coord_dim,d_dim);
    }
  element->evec->update_center_coord();

}


/* ----------------------------------------------------------------------
   change coord according to ramped fashion
   ------------------------------------------------------------------------- */

void Displace::ramp_coord(double coord_lo, double coord_hi, 
    double d_lo, double d_hi, double *coord, int coord_dim, int d_dim)
{
  double fraction = (coord[coord_dim] - coord_lo) / (coord_hi - coord_lo);
  fraction = MAX(fraction,0.0);
  fraction = MIN(fraction,1.0);
  double dramp = d_lo + fraction*(d_hi - d_lo);
  coord[d_dim] += dramp;
}

/* ----------------------------------------------------------------------
   move atoms/elements/nodes by reflecting accross a mirror (x,y, or z plane)
   ------------------------------------------------------------------------- */

void Displace::mirror(char **arg)
{
  int dim;
  if (strcmp(arg[0],"x") == 0) dim = 0;
  else if (strcmp(arg[0],"y") == 0) dim = 1;
  else if (strcmp(arg[0],"z") == 0) dim = 2;
  double mirror_point = universe->numeric(FLERR,arg[1]);
  for (int i = 0; i < atom->nlocal; i++) 
    atom->x[i][dim] = 2.0*mirror_point - atom->x[i][dim];
  for (int i = 0; i < element->nlocal; i++) {
    element->x[i][dim] = 2.0*mirror_point - element->x[i][dim];
    for (int j = 0; j < element->npe[element->etype[i]]; j++)
      element->nodex[i][j][dim] = 2.0*mirror_point - element->nodex[i][j][dim];
  } 

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
  while (iarg < narg) 
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace command");
      if (strcmp(arg[iarg+1],"atom") == 0) group_flag = Group::ATOM;
      else if (strcmp(arg[iarg+1],"node") == 0) group_flag = Group::NODE;
      else if (strcmp(arg[iarg+1],"element") == 0) group_flag = Group::ELEMENT;
      else error->all(FLERR,"Illegal displace command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sequence") == 0) {
      if (style != MULTIDISL) error->all(FLERR,"Illegal displace command");
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace command");
      if (strcmp(arg[iarg+1],"yes") == 0) sequence_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) sequence_flag = 0;
      else error->all(FLERR,"Illegal displace command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace command");
}

/* ----------------------------------------------------------------------
   rotate atoms/elements by right-hand rule by theta around R
   P = point = vector = point of rotation
   R = vector = axis of rotation
   R0 = runit = unit vector for R
   D = X - P = vector from P to X
   C = (D dot R0) R0 = projection of atom coord onto R line
   A = D - C = vector from R line to X
   B = R0 cross A = vector perp to A in plane of rotation
   A,B define plane of circular rotation around R line
   X = P + C + A cos(theta) + B sin(theta)
   ------------------------------------------------------------------------- */

void Displace::rotate(char **arg)
{
  double theta_new;
  double axis[3],point[3],runit[3];

  int dim = domain->dimension;
  point[0] = xscale[0]*universe->numeric(FLERR,arg[0]);
  point[1] = xscale[1]*universe->numeric(FLERR,arg[1]);
  point[2] = xscale[2]*universe->numeric(FLERR,arg[2]);
  axis[0] = universe->numeric(FLERR,arg[3]);
  axis[1] = universe->numeric(FLERR,arg[4]);
  axis[2] = universe->numeric(FLERR,arg[5]);
  double theta = universe->numeric(FLERR,arg[6]);
  if (dim == 2 && (axis[0] != 0.0 || axis[1] != 0.0))
    error->all(FLERR,"Invalid displace rotate axis for 2d");

  double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
  if (len == 0.0)
    error->all(FLERR,"Zero length rotation vector with displace");
  runit[0] = axis[0]/len;
  runit[1] = axis[1]/len;
  runit[2] = axis[2]/len;

  double angle = MY_PI*theta/180.0;
  double cosine = cos(angle);
  double sine = sin(angle);

  double **ax = atom->x;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  imageint *aimage = atom->image;

  for (int i = 0; i < nalocal; i++) {
    if (amask[i] & groupbit) {
      // unwrap coordinate and reset image flags accordingly
      domain->unmap(ax[i],aimage[i]);
      aimage[i] = ((imageint) IMGMAX << IMG2BITS) |
        ((imageint) IMGMAX << IMGBITS) | IMGMAX;
      rotate_point(ax[i],point,runit,cosine,sine);
    }
  }

  if (group_flag == Group::ATOM) return;

  double **ex = element->x;
  double ***nodex = element->nodex;

  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *npe = element->npe;
  imageint *eimage = element->image;

  for (int i = 0; i < nelocal; i++) {
    if (emask[i] & groupbit) {
      // unwrap coordinate and reset image flags accordingly
      domain->unmap(ex[i],nodex[i],eimage[i],npe[etype[i]]);
      eimage[i] = ((imageint) IMGMAX << IMG2BITS) |
        ((imageint) IMGMAX << IMGBITS) | IMGMAX;

      rotate_point(ex[i],point,runit,cosine,sine);
      for (int j = 0; j < npe[etype[i]]; j++)
        rotate_point(nodex[i][j],point,runit,cosine,sine);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Displace::rotate_point(double *x, double *point, double *runit, double cosine, double sine)
{
  double a[3],b[3],c[3],d[3],disp[3];
  double ddotr;
  d[0] = x[0] - point[0];
  d[1] = x[1] - point[1];
  d[2] = x[2] - point[2];
  ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
  c[0] = ddotr*runit[0];
  c[1] = ddotr*runit[1];
  c[2] = ddotr*runit[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1]*a[2] - runit[2]*a[1];
  b[1] = runit[2]*a[0] - runit[0]*a[2];
  b[2] = runit[0]*a[1] - runit[1]*a[0];
  disp[0] = a[0]*cosine  + b[0]*sine;
  disp[1] = a[1]*cosine  + b[1]*sine;
  disp[2] = a[2]*cosine  + b[2]*sine;
  x[0] = point[0] + c[0] + disp[0];
  x[1] = point[1] + c[1] + disp[1];
  if (domain->dimension == 3) x[2] = point[2] + c[2] + disp[2];
}

/* ----------------------------------------------------------------------
   move atoms/elements according to theoretical dislocation displacement field 
   create 1 dislocation
   ------------------------------------------------------------------------- */

void Displace::dislocation(char **arg)
{
  int d_dim,slip_dim,i,j,dim_flag[3];
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
  } else error->all(FLERR,"Illegal displace dislocation command:"
      "invalid dislocation line direction");

  // slip direction

  int direction;
  if (strcmp(arg[3],"+x") == 0) {
    slip_dim = 0;
    dim_flag[0] = 0;
    direction = 1;
  } else if (strcmp(arg[3],"+y") == 0) {
    slip_dim = 1;
    dim_flag[1] = 0;
    direction = 1;
  } else if (strcmp(arg[3],"+z") == 0) {
    slip_dim = 2;
    dim_flag[2] = 0;
    direction = 1;
  } else if (strcmp(arg[3],"-x") == 0) {
    slip_dim = 0;
    dim_flag[0] = 0;
    direction = -1;
  } else if (strcmp(arg[3],"-y") == 0) {
    slip_dim = 1;
    dim_flag[1] = 0;
    direction = -1;
  } else if (strcmp(arg[3],"-z") == 0) {
    slip_dim = 2;
    dim_flag[2] = 0;
    direction = -1;
  } else error->all(FLERR,"Illegal displace dislocation command:"
      "invalid splip direction");

  if (slip_dim == d_dim) 
    error->all(FLERR,"Illegal displace dislocation command:"
        "splip direction and dislocation line direction must be different");

  int dim1 = slip_dim;
  int dim2 = dim_flag[0] * 0 + dim_flag[1] * 1 + dim_flag[2] * 2;
  double center1,center2;
  double b;           // Burgers vector
  double disl_angle;  // angle (in degrees) between Burgers vector and 
  // dislocation line direction (toward slip direction)
  double v;           // Poisson's ratio
  double u,u1,u2,x1,x2;

  if (dim1 == (d_dim + 1)%3) {
    center1 = universe->numeric(FLERR,arg[1]);
    center2 = universe->numeric(FLERR,arg[2]);
  } else {
    center1 = universe->numeric(FLERR,arg[2]);
    center2 = universe->numeric(FLERR,arg[1]);
  }

  b = universe->numeric(FLERR,arg[4]);
  disl_angle = universe->numeric(FLERR,arg[5]);
  v = universe->numeric(FLERR,arg[6]);

  // check for errors

  if (b <= 0) error->all(FLERR,"Illegal displace dislocation command:"
      "Burger vector must be positive");
  if (v >= 1 || v < 0) 
    error->all(FLERR,"Illegal displace dislocation command:"
        "Invalid poisson's ratio");

  double bedge = b*sin(disl_angle/180*MY_PI);
  double bscrew = b*cos(disl_angle/180*MY_PI);

  double **ax = atom->x;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int *npe = element->npe;

  // displace atoms

  if (nalocal) 
    for (i = 0; i < nalocal; i++) {
      x1 = direction*(ax[i][dim1] - center1);
      x2 = direction*(ax[i][dim2] - center2);
      u = u1 = u2 = 0.0;
      edge_dislocation_field(u1,u2,x1,x2,bedge,v);
      screw_dislocation_field(u,x1,x2,bscrew);
      ax[i][dim1] += u1;
      ax[i][dim2] += u2;
      ax[i][d_dim] += u;
    }

  // displace elements

  if (nelocal) {
    for (i = 0; i < nelocal; i++)
      for (j = 0; j < npe[etype[i]]; j++) {
        x1 = direction*(nodex[i][j][dim1] - center1);
        x2 = direction*(nodex[i][j][dim2] - center2);
        u = u1 = u2 = 0.0;
        edge_dislocation_field(u1,u2,x1,x2,bedge,v);
        screw_dislocation_field(u,x1,x2,bscrew);
        nodex[i][j][dim1] += u1;
        nodex[i][j][dim2] += u2;
        nodex[i][j][d_dim] += u;
      }
    element->evec->update_center_coord();
  }
}

/* ----------------------------------------------------------------------
   move atoms/elements according to theoretical dislocation displacement field 
   create multiple dislocations
   ------------------------------------------------------------------------- */

void Displace::multi_dislocations(char **arg)
{
  int d_dim,slip_dim,r_dim,i,j,dim_flag[3];
  dim_flag[0] = dim_flag[1] = dim_flag[2] = 1;

  int ndisl = universe->inumeric(FLERR,arg[0]); // number of dislocations

  // dislocation line direction

  if (strcmp(arg[1],"x") == 0) {
    d_dim = 0;
    dim_flag[0] = 0;
  } else if (strcmp(arg[1],"y") == 0) {
    d_dim = 1;
    dim_flag[1] = 0;
  } else if (strcmp(arg[1],"z") == 0) {
    d_dim = 2;
    dim_flag[2] = 0;
  } else error->all(FLERR,"Illegal displace dislocation command:"
      "invalid dislocation line direction");

  // slip direction

  int direction;
  if (strcmp(arg[4],"+x") == 0) {
    slip_dim = 0;
    dim_flag[0] = 0;
    direction = 1;
  } else if (strcmp(arg[4],"+y") == 0) {
    slip_dim = 1;
    dim_flag[1] = 0;
    direction = 1;
  } else if (strcmp(arg[4],"+z") == 0) {
    slip_dim = 2;
    dim_flag[2] = 0;
    direction = 1;
  } else if (strcmp(arg[4],"-x") == 0) {
    slip_dim = 0;
    dim_flag[0] = 0;
    direction = -1;
  } else if (strcmp(arg[4],"-y") == 0) {
    slip_dim = 1;
    dim_flag[1] = 0;
    direction = -1;
  } else if (strcmp(arg[4],"-z") == 0) {
    slip_dim = 2;
    dim_flag[2] = 0;
    direction = -1;
  } else error->all(FLERR,"Illegal displace dislocation command:"
      "invalid splip direction");

  if (slip_dim == d_dim) 
    error->all(FLERR,"Illegal displace dislocation command:"
        "splip direction and dislocation line direction must be different");

  int dim1 = slip_dim;
  int dim2 = dim_flag[0] * 0 + dim_flag[1] * 1 + dim_flag[2] * 2;
  double center1,center2;
  double offset,ratio;
  double b;           // Burgers vector
  double disl_angle;  // angle (in degrees) between Burgers vector and 
  // dislocation line direction (toward slip direction)
  double v;           // Poisson's ratio
  double u,u1,u2,x1,x2; 

  if (dim1 == (d_dim + 1)%3) {
    center1 = universe->numeric(FLERR,arg[2]);
    center2 = universe->numeric(FLERR,arg[3]);
  } else {
    center1 = universe->numeric(FLERR,arg[3]);
    center2 = universe->numeric(FLERR,arg[2]);
  }

  offset = universe->numeric(FLERR,arg[5]);     // offset distance between dislocation
  ratio = universe->numeric(FLERR,arg[6]);      // ratio between each offset
  b = universe->numeric(FLERR,arg[7]);
  disl_angle = universe->numeric(FLERR,arg[8]);
  v = universe->numeric(FLERR,arg[9]);

  // check for errors

  if (offset == 0) error->all(FLERR,"Illegal displace dislocation command:"
      "offset between dislocation can't be 0");
  if (ratio <= 0) error->all(FLERR,"Illegal displace dislocation command:"
      "ratio must be positive");
  if (b <= 0) error->all(FLERR,"Illegal displace dislocation command:"
      "Burger vector must be positive");
  if (v >= 1 || v < 0) 
    error->all(FLERR,"Illegal displace dislocation command:"
        "Invalid poisson's ratio");

  double bedge = b*sin(disl_angle/180*MY_PI);
  double bscrew = b*cos(disl_angle/180*MY_PI);

  double *list_center1 = new double[ndisl];

  list_center1[0] = center1;
  for (int idisl = 1; idisl < ndisl; idisl++) 
    list_center1[idisl] = list_center1[idisl-1] + offset*pow(ratio,idisl-1);

  double **ax = atom->x;
  int *amask = atom->mask;
  int nalocal = atom->nlocal;
  double ***nodex = element->nodex;
  int *emask = element->mask;
  int *etype = element->etype;
  int **nodemask = element->nodemask;
  int nelocal = element->nlocal;
  int *npe = element->npe;

  // displace atoms

  if (nalocal) 
    for (i = 0; i < nalocal; i++) {
      if (sequence_flag) {
        for (int idisl = 0; idisl < ndisl; idisl++) {
          u = u1 = u2 = 0.0;
          x1 = direction*(ax[i][dim1] - list_center1[idisl]);
          x2 = direction*(ax[i][dim2] - center2);
          edge_dislocation_field(u1,u2,x1,x2,bedge,v);
          screw_dislocation_field(u,x1,x2,bscrew);
          ax[i][dim1] += u1;
          ax[i][dim2] += u2;
          ax[i][d_dim] += u;
        }
      } else {
        u = u1 = u2 = 0.0;
        x2 = direction*(ax[i][dim2] - center2);
        for (int idisl = 0; idisl < ndisl; idisl++) {
          x1 = direction*(ax[i][dim1] - list_center1[idisl]);
          edge_dislocation_field(u1,u2,x1,x2,bedge,v);
          screw_dislocation_field(u,x1,x2,bscrew);
        }
        ax[i][dim1] += u1;
        ax[i][dim2] += u2;
        ax[i][d_dim] += u;
      }
    }

  // displace elements

  if (nelocal) {
    for (i = 0; i < nelocal; i++)
      for (j = 0; j < npe[etype[i]]; j++) {
        if (sequence_flag) {
          for (int idisl = 0; idisl < ndisl; idisl++) {
            u = u1 = u2 = 0.0;
            x1 = direction*(nodex[i][j][dim1] - list_center1[idisl]);
            x2 = direction*(nodex[i][j][dim2] - center2);
            edge_dislocation_field(u1,u2,x1,x2,bedge,v);
            screw_dislocation_field(u,x1,x2,bscrew);
            nodex[i][j][dim1] += u1;
            nodex[i][j][dim2] += u2;
            nodex[i][j][d_dim] += u;
          }
        } else {
          u = u1 = u2 = 0.0;
          x2 = direction*(nodex[i][j][dim2] - center2);
          for (int idisl = 0; idisl < ndisl; idisl++) {
            x1 = direction*(nodex[i][j][dim1] - list_center1[idisl]);
            edge_dislocation_field(u1,u2,x1,x2,bedge,v);
            screw_dislocation_field(u,x1,x2,bscrew);
          }
          nodex[i][j][dim1] += u1;
          nodex[i][j][dim2] += u2;
          nodex[i][j][d_dim] += u;
        }
      }
    element->evec->update_center_coord();
  }
  delete list_center1;
}


/* ----------------------------------------------------------------------
   calculate theoretical dislocation displacement field 
   ------------------------------------------------------------------------- */

void Displace::edge_dislocation_field(double &u1, double &u2, 
    double x1, double x2, double b, double v)
{

  double rsq = x1*x1 + x2*x2;

  if (rsq < EPSILON)
    error->one(FLERR,"Atom/Node too close to singularity ");

  // for u1: need range of inverse tan to be from 0 to 2*pi
  // range of atan2 is from -pi to +pi so add a pi to offset

  u1 += b / MY_2PI * ((atan2(x2,x1) + MY_PI) + x1*x2/(2.0*(1-v)*rsq));
  u2 += -b / (MY_8PI*(1-v)) * ((1-2*v)*log(rsq) + (x1*x1-x2*x2)/rsq);
}

/* ----------------------------------------------------------------------
   calculate theoretical screw dislocation displacement field 
   ------------------------------------------------------------------------- */

void Displace::screw_dislocation_field(double &u, double x1, double x2, double b)
{
  // need range of inverse tan to be from 0 to 2*pi
  // range of atan2 is from -pi to +pi so add a pi to offset

  u += b/MY_2PI*(atan2(x2,x1) + MY_PI);
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

  if (group_flag == Group::ATOM) return;

  double ***nodex = element->nodex;
  mask = element->mask;
  int **nodemask = element->nodemask;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < element->nlocal; i++)
    if (group_flag == Group::ELEMENT) {
      if (mask[i] & groupbit)
        for (int j = 0; j < npe[etype[i]]; j++) nodex[i][j][dim] += delta;
    } else if (group_flag == Group::NODE) {
      for (int j = 0; j < npe[etype[i]]; j++)
        if (nodemask[i][j] & groupbit) nodex[i][j][dim] += delta;
    }
  element->evec->update_center_coord();
}
