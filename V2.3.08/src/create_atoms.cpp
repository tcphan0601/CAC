#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "create_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
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

enum{BOX,REGION,SINGLE,RANDOM};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

CreateAtoms::CreateAtoms(CAC *cac) : Pointers(cac) {}

/* ---------------------------------------------------------------------- */

void CreateAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_atoms command before simulation box is defined");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot create_atoms after "
               "reading restart file with per-atom info");

  // parse arguments

  if (narg < 2) error->all(FLERR,"Illegal create_atoms command");
  ntype = universe->inumeric(FLERR,arg[0]);

  int iarg;
  if (strcmp(arg[1],"box") == 0) {
    style = BOX;
    iarg = 2;
  } else if (strcmp(arg[1],"region") == 0) {
    style = REGION;
    if (narg < 3) error->all(FLERR,"Illegal create_atoms command");
    nregion = domain->find_region(arg[2]);
    if (nregion == -1) error->all(FLERR,
                                  "Create_atoms region ID does not exist");
    domain->regions[nregion]->init();
    domain->regions[nregion]->prematch();
    iarg = 3;;
  } else if (strcmp(arg[1],"single") == 0) {
    style = SINGLE;
    if (narg < 5) error->all(FLERR,"Illegal create_atoms command");
    xone[0] = universe->numeric(FLERR,arg[2]);
    xone[1] = universe->numeric(FLERR,arg[3]);
    xone[2] = universe->numeric(FLERR,arg[4]);
    iarg = 5;
  } else if (strcmp(arg[1],"random") == 0) {
    style = RANDOM;
    if (narg < 5) error->all(FLERR,"Illegal create_atoms command");
    nrandom = universe->inumeric(FLERR,arg[2]);
    seed = universe->inumeric(FLERR,arg[3]);
    if (strcmp(arg[4],"NULL") == 0) nregion = -1;
    else {
      nregion = domain->find_region(arg[4]);
      if (nregion == -1) error->all(FLERR,
                                    "Create_atoms region ID does not exist");
      domain->regions[nregion]->init();
      domain->regions[nregion]->prematch();
    }
    iarg = 5;
  } else error->all(FLERR,"Illegal create_atoms command");

  // process optional keywords

  int scaleflag = 0;
  remapflag = 0;
  int molseed;
  varflag = 0;
  vstr = xstr = ystr = zstr = NULL;
  quatone[0] = quatone[1] = quatone[2] = 0.0;

  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) basistype[i] = ntype;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_atoms command");
      int ibasis = universe->inumeric(FLERR,arg[iarg+1]);
      int itype = universe->inumeric(FLERR,arg[iarg+2]);
      if (ibasis <= 0 || ibasis > nbasis || itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Invalid basis setting in create_atoms command");
      basistype[ibasis-1] = itype;
      iarg += 3;
    } else if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_atoms command");
      if (strcmp(arg[iarg+1],"yes") == 0) remapflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) remapflag = 0;
      else error->all(FLERR,"Illegal create_atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal create_atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"var") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_atoms command");
      delete [] vstr;
      int n = strlen(arg[iarg+1]) + 1;
      vstr = new char[n];
      strcpy(vstr,arg[iarg+1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"set") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_atoms command");
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
      } else error->all(FLERR,"Illegal create_atoms command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (style != SINGLE)
        error->all(FLERR,"Cannot use create_atoms rotate unless single style");
      if (iarg+5 > narg) error->all(FLERR,"Illegal create_atoms command");
      double thetaone;
      double axisone[3];
      thetaone = universe->numeric(FLERR,arg[iarg+1]);
      axisone[0] = universe->numeric(FLERR,arg[iarg+2]);
      axisone[1] = universe->numeric(FLERR,arg[iarg+3]);
      axisone[2] = universe->numeric(FLERR,arg[iarg+4]);
      if (axisone[0] == 0.0 && axisone[1] == 0.0 && axisone[2] == 0.0)
        error->all(FLERR,"Illegal create_atoms command");
      if (domain->dimension == 2 && (axisone[0] != 0.0 || axisone[1] != 0.0))
        error->all(FLERR,"Invalid create_atoms rotation vector for 2d model");
      MathExtra::norm3(axisone);
      MathExtra::axisangle_to_quat(axisone,thetaone,quatone);
      iarg += 5;
    } else error->all(FLERR,"Illegal create_atoms command");
  }

  // error checks

  if (ntype <= 0 || ntype > atom->ntypes)
    error->all(FLERR,"Invalid atom type in create_atoms command");

  if (style == RANDOM) {
    if (nrandom < 0) error->all(FLERR,"Illegal create_atoms command");
    if (seed <= 0) error->all(FLERR,"Illegal create_atoms command");
  }

  // error check and further setup for variable test

  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR,"Incomplete use of variables in create_atoms command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR,"Incomplete use of variables in create_atoms command");

  if (varflag) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for create_atoms does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for create_atoms is invalid style");

    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(xvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(yvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0)
        error->all(FLERR,"Variable name for create_atoms does not exist");
      if (!input->variable->internalstyle(zvar))
        error->all(FLERR,"Variable for create_atoms is invalid style");
    }
  }

  // demand non-none lattice be defined for BOX and REGION
  // else setup scaling for SINGLE and RANDOM
  // could use domain->lattice->lattice2box() to do conversion of
  //   lattice to box, but not consistent with other uses of units=lattice
  // triclinic remapping occurs in add_single()

  if (style == BOX || style == REGION) {
    if (nbasis == 0)
      error->all(FLERR,"Cannot create atoms with undefined lattice");
  } else if (scaleflag == 1) {
    xone[0] *= domain->lattice->xlattice;
    xone[1] *= domain->lattice->ylattice;
    xone[2] *= domain->lattice->zlattice;
  }

  // set bounds for my proc in sublo[3] & subhi[3]
  // if periodic and style = BOX or REGION, i.e. using lattice:
  //   should create exactly 1 atom when 2 images are both "on" the boundary
  //   either image may be slightly inside/outside true box due to round-off
  //   if I am lo proc, decrement lower bound by EPSILON
  //     this will insure lo image is created
  //   if I am hi proc, decrement upper bound by 2.0*EPSILON
  //     this will insure hi image is not created
  //   thus insertion box is EPSILON smaller than true box
  //     and is shifted away from true boundary
  //     which is where atoms are likely to be generated

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

  // Record wall time for atom creation

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // clear ghost count and any ghost bonus data internal to AtomVec
  // same logic as beginning of Comm::exchange()
  // do it now b/c creating atoms will overwrite ghost atoms

  atom->nghost = 0;
  //atom->avec->clear_bonus();

  // add atoms in one of 3 ways

  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  if (style == SINGLE) add_single();
  else if (style == RANDOM) add_random();
  else add_lattice();

  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_variable(nlocal_previous,atom->nlocal);

  // set new total # of atoms and error check

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR,"Too many total atoms");

  // add IDs for newly created atoms
  // check that atom IDs are valid

  if (atom->tag_enable) atom->tag_extend();
  atom->tag_check();

  // if global map exists, reset it
  // invoke map_init() b/c atom count has grown
  // set nghost to 0 so old ghosts won't be mapped

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // clean up

  if (domain->lattice) delete [] basistype;
  delete [] vstr;
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;

  // print status

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Created " BIGINT_FORMAT " atoms\n",
              atom->natoms-natoms_previous);
      fprintf(screen,"  Time spent = %g secs\n",time2-time1);
    }
    if (logfile) {
      fprintf(logfile,"Created " BIGINT_FORMAT " atoms\n",
              atom->natoms-natoms_previous);
      fprintf(logfile,"  Time spent = %g secs\n",time2-time1);
    }
  }
}

/* ----------------------------------------------------------------------
   add single atom with coords at xone if it's in my sub-box
   if triclinic, xone is in lamda coords
------------------------------------------------------------------------- */

void CreateAtoms::add_single()
{
  // remap atom if requested

  if (remapflag) {
    imageint imagetmp = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(xone,imagetmp);
  }

  // if triclinic, convert to lamda coords (0-1)

  double lamda[3],*coord;
  if (triclinic) {
    domain->x2lamda(xone,lamda);
    coord = lamda;
  } else coord = xone;

  // if atom is in my subbox, create it

  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2]) 
    atom->avec->create_atom(xone,ntype,0);
}

/* ----------------------------------------------------------------------
   add Nrandom atoms at random locations
------------------------------------------------------------------------- */

void CreateAtoms::add_random()
{
  double xlo,ylo,zlo,xhi,yhi,zhi,zmid;
  double lamda[3],*coord;
  double *boxlo,*boxhi;

  // random number generator, same for all procs

  RanPark *random = new RanPark(cac,seed);

  // bounding box for atom creation
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

  if (nregion >= 0 && domain->regions[nregion]->bboxflag) {
    xlo = MAX(xlo,domain->regions[nregion]->extent_xlo);
    xhi = MIN(xhi,domain->regions[nregion]->extent_xhi);
    ylo = MAX(ylo,domain->regions[nregion]->extent_ylo);
    yhi = MIN(yhi,domain->regions[nregion]->extent_yhi);
    zlo = MAX(zlo,domain->regions[nregion]->extent_zlo);
    zhi = MIN(zhi,domain->regions[nregion]->extent_zhi);
  }

  // generate random positions for each new atom within bounding box
  // iterate until atom is within region, variable, and triclinic simulation box
  // if final atom position is in my subbox, create it

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"No overlap of box and region for create_atoms");

  int valid;
  for (int i = 0; i < nrandom; i++) {
    while (1) {
      xone[0] = xlo + random->uniform() * (xhi-xlo);
      xone[1] = ylo + random->uniform() * (yhi-ylo);
      xone[2] = zlo + random->uniform() * (zhi-zlo);
      if (domain->dimension == 2) xone[2] = zmid;

      valid = 1;
      if (nregion >= 0 &&
          domain->regions[nregion]->match(xone) == 0)
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
      atom->avec->create_atom(xone,ntype,0);
    }
  }

  // clean-up

  delete random;
}

/* ----------------------------------------------------------------------
   add many atoms by looping over lattice
------------------------------------------------------------------------- */

void CreateAtoms::add_lattice()
{
  // convert 8 corners of my subdomain from box coords to lattice coords
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
  //   which can lead to missing atoms in rare cases
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

  // iterate on 3d periodic lattice of unit cells using loop bounds
  // iterate on nbasis atoms in each unit cell
  // convert lattice coords to box coords
  // add atom if it meets all criteria

  double **basis = domain->lattice->basis;
  double x[3],lamda[3];
  double *coord;

  int i,j,k,m;
  for (k = klo; k <= khi; k++)
    for (j = jlo; j <= jhi; j++)
      for (i = ilo; i <= ihi; i++)
        for (m = 0; m < nbasis; m++) {

          x[0] = i + basis[m][0];
          x[1] = j + basis[m][1];
          x[2] = k + basis[m][2];

          // convert from lattice coords to box coords

          domain->lattice->lattice2box(x[0],x[1],x[2]);

          // if a region was specified, test if atom is in it

          if (style == REGION)
            if (!domain->regions[nregion]->match(x)) continue;

          // if variable test specified, eval variable

          if (varflag && vartest(x) == 0) continue;

          // test if atom position is in my subbox

          if (triclinic) {
            domain->x2lamda(x,lamda);
            coord = lamda;
          } else coord = x;

          if (coord[0] < sublo[0] || coord[0] >= subhi[0] ||
              coord[1] < sublo[1] || coord[1] >= subhi[1] ||
              coord[2] < sublo[2] || coord[2] >= subhi[2]) continue;

          // add the atom to my list of atoms

          atom->avec->create_atom(x,basistype[m],0);
        }
}

/* ----------------------------------------------------------------------
   test a generated atom position against variable evaluation
   first set x,y,z values in internal variables
------------------------------------------------------------------------- */

int CreateAtoms::vartest(double *x)
{
  if (xstr) input->variable->internal_set(xvar,x[0]);
  if (ystr) input->variable->internal_set(yvar,x[1]);
  if (zstr) input->variable->internal_set(zvar,x[2]);

  double value = input->variable->compute_equal(vvar);

  if (value == 0.0) return 0;
  return 1;
}
