#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "create_elements.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "universe.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "error.h"

using namespace CAC_NS;
using namespace MathConst;

#define BIG 1.0e30
#define EPSILON 1.0e-6

enum{BOX,REGION,RANDOM};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

CreateElements::CreateElements(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void CreateElements::command(int narg, char **arg)
{
  
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_elements command before simulation box is defined");
  if (domain->lattice->nbasis == 0)
    error->all(FLERR,"Cannot create elements with undefined lattice");

  //if (modify->nfix_restart_peratom)
  //  error->all(FLERR,"Cannot create_atoms after "
  //             "reading restart file with per-atom info");

  // parse arguments

  if (narg < 2) error->all(FLERR,"Illegal create_elements command");
  netype = universe->inumeric(FLERR,arg[0]);
  nctype = universe->inumeric(FLERR,arg[1]);

  int iarg;
  if (strcmp(arg[2],"box") == 0) {
    style = BOX;
    iarg = 3;
  } else if (strcmp(arg[2],"region") == 0) {
    style = REGION;
    if (narg < 4) error->all(FLERR,"Illegal create_elements command");
    iregion = domain->find_region(arg[3]);
    if (iregion == -1) error->all(FLERR,
        "Create_elements region ID does not exist");
    domain->regions[iregion]->init();
    domain->regions[iregion]->prematch();
    iarg = 4;
  //} else if (strcmp(arg[2],"single") == 0) {
  //  style = SINGLE;
  //  if (narg < 6) error->all(FLERR,"Illegal create_elements command");
  //  xone[0] = universe->numeric(FLERR,arg[3]);
  //  xone[1] = universe->numeric(FLERR,arg[4]);
  //  xone[2] = universe->numeric(FLERR,arg[5]);
  //  iarg = 6;
  } else error->all(FLERR,"Illegal create_elements command");

  // process optional keywords

  int scaleflag = 0;
  remapflag = 0;
  int molseed;
  varflag = 0;
  vstr = xstr = ystr = zstr = NULL;
  quatone[0] = quatone[1] = quatone[2] = 0.0;

  nbasis = domain->lattice->nbasis;
  basisctype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) {
    basisctype[i] = nctype;
  }

  while (iarg < narg) {
    if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_elements command");
      int ibasis = universe->inumeric(FLERR,arg[iarg+1]);
      int ictype = universe->inumeric(FLERR,arg[iarg+2]);
      if (ibasis <= 0 || ibasis > nbasis || ictype <= 0 || ictype > atom->ntypes)
        error->all(FLERR,"Invalid basis setting in create_elements command");
      basisctype[ibasis-1] = ictype;
      iarg += 3;
    } else if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_elements command");
      if (strcmp(arg[iarg+1],"yes") == 0) remapflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) remapflag = 0;
      else error->all(FLERR,"Illegal create_elements command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_elements command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal create_elements command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"var") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_elements command");
      delete [] vstr;
      int n = strlen(arg[iarg+1]) + 1;
      vstr = new char[n];
      strcpy(vstr,arg[iarg+1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"set") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_elements command");
      if (strcmp(arg[iarg+1],"x") == 0) {
        delete [] xstr;
        int n = strlen(arg[iarg+2]) + 1;
        xstr = new char[n];
        strcpy(xstr,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"y") == 0) {
        delete [] ystr;
        int n = strlen(arg[iarg+2]) + 1;
        ystr = new char[n];
        strcpy(ystr,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"z") == 0) {
        delete [] zstr;
        int n = strlen(arg[iarg+2]) + 1;
        zstr = new char[n];
        strcpy(zstr,arg[iarg+2]);
      } else error->all(FLERR,"Illegal create_elements command");
      iarg += 3;
    //} else if (strcmp(arg[iarg],"rotate") == 0) {
    //  if (style != SINGLE)
    //    error->all(FLERR,"Cannot use create_elements rotate unless single style");
    //  if (iarg+5 > narg) error->all(FLERR,"Illegal create_elements command");
    //  double thetaone;
    //  double axisone[3];
    //  thetaone = universe->numeric(FLERR,arg[iarg+1]);
    //  axisone[0] = universe->numeric(FLERR,arg[iarg+2]);
    //  axisone[1] = universe->numeric(FLERR,arg[iarg+3]);
    //  axisone[2] = universe->numeric(FLERR,arg[iarg+4]);
    //  if (axisone[0] == 0.0 && axisone[1] == 0.0 && axisone[2] == 0.0)
    //    error->all(FLERR,"Illegal create_elements command");
    //  if (domain->dimension == 2 && (axisone[0] != 0.0 || axisone[1] != 0.0))
    //    error->all(FLERR,"Invalid create_elements rotation vector for 2d model");
    //  MathExtra::norm3(axisone);
    //  MathExtra::axisangle_to_quat(axisone,thetaone,quatone);
    //  iarg += 5;
    } else error->all(FLERR,"Illegal create_elements command");
  }

  // error checks

  if (nctype <= 0 || nctype > atom->ntypes)
    error->all(FLERR,"Invalid atom type in create_elements command");
  if (netype <= 0 || netype > element->netypes)
    error->all(FLERR,"Invalid element type in create_elements command");
  if (element->apc[netype] != nbasis) 
    error->all(FLERR,"Number of apc for this element type must match number of basis atoms in the unit cell");

  //if (style == RANDOM) {
  //  if (nrandom < 0) error->all(FLERR,"Illegal create_elements command");
  //  if (seed <= 0) error->all(FLERR,"Illegal create_elements command");
  //}

  // error check and further setup for variable test

  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR,"Incomplete use of variables in create_elements command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR,"Incomplete use of variables in create_elements command");

  if (varflag) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for create_elements does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for create_elements is invalid style");

    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0)
        error->all(FLERR,"Variable name for create_elements does not exist");
      if (!input->variable->internalstyle(xvar))
        error->all(FLERR,"Variable for create_elements is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0)
        error->all(FLERR,"Variable name for create_elements does not exist");
      if (!input->variable->internalstyle(yvar))
        error->all(FLERR,"Variable for create_elements is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0)
        error->all(FLERR,"Variable name for create_elements does not exist");
      if (!input->variable->internalstyle(zvar))
        error->all(FLERR,"Variable for create_elements is invalid style");
    }
  }

  // demand non-none lattice be defined for BOX and REGION
  // else setup scaling for SINGLE and RANDOM
  // could use domain->lattice->lattice2box() to do conversion of
  //   lattice to box, but not consistent with other uses of units=lattice
  // triclinic remapping occurs in add_single()


  //if (style == SINGLE && scaleflag == 1) {
  //  xone[0] *= domain->lattice->xlattice;
  //  xone[1] *= domain->lattice->ylattice;
  //  xone[2] *= domain->lattice->zlattice;
  //}

  // set bounds for my proc in sublo[3] & subhi[3]
  // if periodic and style = BOX or REGION, i.e. using lattice:
  //   should create exactly 1 element when 2 images are both "on" the boundary
  //   either image may be slightly inside/outside true box due to round-off
  //   if I am lo proc, decrement lower bound by EPSILON
  //     this will insure lo image is created
  //   if I am hi proc, decrement upper bound by 2.0*EPSILON
  //     this will insure hi image is not created
  //   thus insertion box is EPSILON smaller than true box
  //     and is shifted away from true boundary
  //     which is where elements are likely to be generated

  triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (style == BOX || style == REGION) {
    if (comm->layout != LAYOUT_TILED) {
      if (domain->xperiodic) {
        if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
        if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] -= 2.0*epsilon[0];
      }
      if (domain->yperiodic) {
        if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
        if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] -= 2.0*epsilon[1];
      }
      if (domain->zperiodic) {
        if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
        if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] -= 2.0*epsilon[2];
      }
    } else {
      if (domain->xperiodic) {
        if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
        if (comm->mysplit[0][1] == 1.0) subhi[0] -= 2.0*epsilon[0];
      }
      if (domain->yperiodic) {
        if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
        if (comm->mysplit[1][1] == 1.0) subhi[1] -= 2.0*epsilon[1];
      }
      if (domain->zperiodic) {
        if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
        if (comm->mysplit[2][1] == 1.0) subhi[2] -= 2.0*epsilon[2];
      }
    }
  }

  // Record wall time for element creation

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // clear ghost count 
  // do it now b/c creating elements will overwrite ghost elements

  element->nghost = 0;
  //atom->evec->clear_bonus();

  // add elements in one of 3 ways

  bigint nelements_previous = element->nelements;
  int nlocal_previous = element->nlocal;

  //if (style == SINGLE) add_single();
  //else if (style == RANDOM) add_random();
  //else add_lattice();
  add_lattice();

  // init per-element fix/compute/variable values for created atoms

  element->data_fix_compute_variable(nlocal_previous,element->nlocal);

  // set new total # of elements and error check

  bigint nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&element->nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (element->nelements < 0 || element->nelements >= MAXBIGINT)
    error->all(FLERR,"Too many total elements");

  bigint nelements_added = element->nelements - nelements_previous;

  element->nintpls += nelements_added * element->nintpl[netype];

  // add IDs for newly created elements
  // check that element IDs are valid

  if (element->tag_enable) element->tag_extend();
  element->tag_check();

  // assign element ids to clusters
  
  if (element->element_cluster_flag) {
    tagint **element_clusters = element->element_clusters;
    tagint *tag = element->tag;
    int i;
    for (int i = nlocal_previous; i < element->nlocal; i += nbasis) {
      for (int j = 0; j < nbasis; j++) {
        for (int k = 0; k < nbasis; k++)
          element_clusters[i+j][k] = tag[i+k];
        for (int k = nbasis; k < element->max_apc; k++)
          element_clusters[i+j][k] = 0;
      }
    }
  }

  // if global map exists, reset it
  // invoke map_init() b/c element count has grown
  // set nghost to 0 so old ghosts won't be mapped

  if (element->map_style) {
    element->nghost = 0;
    element->map_init();
    element->map_set();
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // clean up

  delete [] basisctype;
  delete [] vstr;
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;

  // print status

  FILE *out;
  if (comm->me == 0) {
    for (int i = 0; i < 2; i++) {
      if (i == 0) out = screen;
      else out = logfile;
      if (out) {
        fprintf(out,"Created " BIGINT_FORMAT " elements\n",nelements_added);
        fprintf(out,"  Total " BIGINT_FORMAT " elements\n",element->nelements);
        fprintf(out,"  Time spent = %g secs\n",time2-time1);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add single element with coords at xone if it's in my sub-box
   if triclinic, xone is in lamda coords
   ------------------------------------------------------------------------- */

void CreateElements::add_single()
{
  // remap element if requested

//  if (remapflag) {
//    imageint imagetmp = ((imageint) IMGMAX << IMG2BITS) |
//      ((imageint) IMGMAX << IMGBITS) | IMGMAX;
//    domain->remap(xone,imagetmp);
//  }
//
//  // if triclinic, convert to lamda coords (0-1)
//
//  double lamda[3],*coord;
//  if (triclinic) {
//    domain->x2lamda(xone,lamda);
//    coord = lamda;
//  } else coord = xone;
//
//  // if element is in my subbox, create it
//
//  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
//      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
//      coord[2] >= sublo[2] && coord[2] < subhi[2]) {
//    element->evec->create_element(xone,netype,nctype,0);
//    for
//  }
}

/* ----------------------------------------------------------------------
   add Nrandom elements at random locations
   ------------------------------------------------------------------------- */
/*
   void CreateElements::add_random()
   {
   double xlo,ylo,zlo,xhi,yhi,zhi,zmid;
   double lamda[3],*coord;
   double *boxlo,*boxhi;

// random number generator, same for all procs

RanPark *random = new RanPark(cac,seed);

// bounding box for element creation
// in real units, even if triclinic
// only limit bbox by region if its bboxflag is set (interior region)

if (triclinic == 0) {
xlo = domain->boxlo[0]; xhi = domain->boxhi[0];
ylo = domain->boxlo[1]; yhi = domain->boxhi[1];
zlo = domain->boxlo[2]; zhi = domain->boxhi[2];
zmid = zlo + 0.5*(zhi-zlo);
} else {
xlo = domain->boxlo_bound[0]; xhi = domain->boxhi_bound[0];
ylo = domain->boxlo_bound[1]; yhi = domain->boxhi_bound[1];
zlo = domain->boxlo_bound[2]; zhi = domain->boxhi_bound[2];
zmid = zlo + 0.5*(zhi-zlo);
boxlo = domain->boxlo_lamda;
boxhi = domain->boxhi_lamda;
}

if (iregion >= 0 && domain->regions[iregion]->bboxflag) {
xlo = MAX(xlo,domain->regions[iregion]->extent_xlo);
xhi = MIN(xhi,domain->regions[iregion]->extent_xhi);
ylo = MAX(ylo,domain->regions[iregion]->extent_ylo);
yhi = MIN(yhi,domain->regions[iregion]->extent_yhi);
zlo = MAX(zlo,domain->regions[iregion]->extent_zlo);
zhi = MIN(zhi,domain->regions[iregion]->extent_zhi);
}

// generate random positions for each new element within bounding box
// iterate until element is within region, variable, and triclinic simulation box
// if final element position is in my subbox, create it

if (xlo > xhi || ylo > yhi || zlo > zhi)
error->all(FLERR,"No overlap of box and region for create_elements");

int valid;
for (int i = 0; i < nrandom; i++) {
while (1) {
xone[0] = xlo + random->uniform() * (xhi-xlo);
xone[1] = ylo + random->uniform() * (yhi-ylo);
xone[2] = zlo + random->uniform() * (zhi-zlo);
if (domain->dimension == 2) xone[2] = zmid;

valid = 1;
if (iregion >= 0 &&
domain->regions[iregion]->match(xone) == 0)
valid = 0;
if (varflag && vartest(xone) == 0) valid = 0;
if (triclinic) {
domain->x2lamda(xone,lamda);
coord = lamda;
if (coord[0] < boxlo[0] || coord[0] >= boxhi[0] ||
coord[1] < boxlo[1] || coord[1] >= boxhi[1] ||
coord[2] < boxlo[2] || coord[2] >= boxhi[2]) valid = 0;
} else coord = xone;

if (valid) break;
}

// if triclinic, coord is now in lamda units

if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
    coord[1] >= sublo[1] && coord[1] < subhi[1] &&
    coord[2] >= sublo[2] && coord[2] < subhi[2]) {
  element->evec->create_element(xone,ntype);
}
}

// clean-up

delete random;
}
*/
/* ----------------------------------------------------------------------
   add many elements by looping over lattice
   ------------------------------------------------------------------------- */

void CreateElements::add_lattice()
{

  // convert 8 corners of my box domain from box coords to lattice coords
  // for orthogonal, use corner pts of my subbox
  // for triclinic, use bounding box of my subbox
  // xyz min to max = bounding box around the domain corners in lattice space

  double bboxlo[3],bboxhi[3];

  if (triclinic == 0) {
    bboxlo[0] = domain->sublo[0]; bboxhi[0] = domain->subhi[0];
    bboxlo[1] = domain->sublo[1]; bboxhi[1] = domain->subhi[1];
    bboxlo[2] = domain->sublo[2]; bboxhi[2] = domain->subhi[2];
  } else domain->bbox(domain->sublo_lamda,domain->subhi_lamda,bboxlo,bboxhi);

  double xmin,ymin,zmin,xmax,ymax,zmax;
  xmin = ymin = zmin = BIG;
  xmax = ymax = zmax = -BIG;

  domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxlo[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxlo[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxlo[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxlo[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxhi[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxhi[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxhi[2],
      xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxhi[2],
      xmin,ymin,zmin,xmax,ymax,zmax);

  // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my subbox
  // overlap = any part of a unit cell (face,edge,pt) in common with my subbox
  // in lattice space, subbox is a tilted box
  // but bbox of subbox is aligned with lattice axes
  // so ilo:khi unit cells should completely tile bounding box
  // decrement lo, increment hi to avoid round-off issues in lattice->bbox(),
  //   which can lead to missing elements in rare cases
  // extra decrement of lo if min < 0, since static_cast(-1.5) = -1

  int ilo,ihi,jlo,jhi,klo,khi;
  ilo = static_cast<int> (xmin) - 1;
  jlo = static_cast<int> (ymin) - 1;
  klo = static_cast<int> (zmin) - 1;
  ihi = static_cast<int> (xmax) + 1;
  jhi = static_cast<int> (ymax) + 1;
  khi = static_cast<int> (zmax) + 1;

  if (xmin < 0.0) ilo--;
  if (ymin < 0.0) jlo--;
  if (zmin < 0.0) klo--;

  
  int nx,ny,nz;
  element->evec->get_element_size(netype,nx,ny,nz);

  // reset ilo,jlo,klo lattice sites to be multiple of element sizes so that element positions are in sync with other procs
  
  ilo = (ilo/nx - 1)*nx;
  jlo = (jlo/ny - 1)*ny;
  klo = (klo/nz - 1)*nz;

  // iterate through lattice sites
  // add element if it meets all criteria

  int i,j,k,m,dim;
  double x[3],lamda[3];
  double *coord;
  int nbasis = domain->lattice->nbasis;
  double **basis = domain->lattice->basis;
  double basis_center[3];
  
  for (dim = 0; dim < 3; dim++) {
    basis_center[dim] = 0;
    for (m = 0; m < nbasis; m++) 
      basis_center[dim] += basis[m][dim];
    basis_center[dim] /= nbasis;
  }

  for (k = klo; k <= khi; k += nz) {
    for (j = jlo; j <= jhi; j += ny) {
      for (i = ilo; i <= ihi; i += nx) {
        x[0] = i + basis_center[0] + (nx - 1)/2.0;
        x[1] = j + basis_center[1] + (ny - 1)/2.0;
        x[2] = k + basis_center[2] + (nz - 1)/2.0;

        domain->lattice->lattice2box(x[0],x[1],x[2]);

        // if a region was specified, test if element is in it

        if (style == REGION)
          if (!domain->regions[iregion]->match(x)) continue;

        // if variable test specified, eval variable

        if (varflag && vartest(x) == 0) continue;

        // test if element position is in my subbox

        if (triclinic) {
          domain->x2lamda(x,lamda);
          coord = lamda;
        } else coord = x;

        if (coord[0] < sublo[0] || coord[0] >= subhi[0] ||
            coord[1] < sublo[1] || coord[1] >= subhi[1] ||
            coord[2] < sublo[2] || coord[2] >= subhi[2]) continue;

        // add the element(s) to my list of elements

        if (nbasis == 1) {
          element->evec->create_element(x,netype,basisctype[0],0);
        } else {
          for (m = 0; m < nbasis; j++) {
            x[0] = i + basis[m][0] + (nx - 1)/2.0;
            x[1] = j + basis[m][1] + (ny - 1)/2.0;
            x[2] = k + basis[m][2] + (nz - 1)/2.0;

            domain->lattice->lattice2box(x[0],x[1],x[2]);

            element->evec->create_element(x,netype,basisctype[m],0);
          }
        }
      }
    }
  }

  //if (nbasis > 1) domain->reset_box();

}

/* ----------------------------------------------------------------------
   test a generated element position against variable evaluation
   first set x,y,z values in internal variables
   ------------------------------------------------------------------------- */

int CreateElements::vartest(double *x)
{
  if (xstr) input->variable->internal_set(xvar,x[0]);
  if (ystr) input->variable->internal_set(yvar,x[1]);
  if (zstr) input->variable->internal_set(zvar,x[2]);

  double value = input->variable->compute_equal(vvar);

  if (value == 0.0) return 0;
  return 1;
}
