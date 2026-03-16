//#define BALANCE_DEBUG 1

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "balance.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "rcb.h"
#include "irregular_comm.h"
#include "domain.h"
#include "universe.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

enum{XYZ,SHIFT,BISECTION};
enum{NONE,UNIFORM,USER};
enum{X,Y,Z};

/* ---------------------------------------------------------------------- */

Balance::Balance(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  user_xsplit = user_ysplit = user_zsplit = NULL;
  shift_allocate = 0;
  proccost = allproccost = NULL;

  rcb = NULL;

  fp = NULL;
  firststep = 1;

  user_element_weight = 0.0;
}

/* ---------------------------------------------------------------------- */

Balance::~Balance()
{
  memory->destroy(proccost);
  memory->destroy(allproccost);

  delete [] user_xsplit;
  delete [] user_ysplit;
  delete [] user_zsplit;

  if (shift_allocate) {
    delete [] bdim;
    delete [] onecost;
    delete [] allcost;
    delete [] sum;
    delete [] target;
    delete [] lo;
    delete [] hi;
    delete [] losum;
    delete [] hisum;
  }

  delete rcb;
 
  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   called as balance command in input script
------------------------------------------------------------------------- */

void Balance::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Balance command before simulation box is defined");

  // insure particles are in current box & update box via shrink-wrap
  // init entire system since comm->setup is done
  
  // must reset atom map after exchange() since it clears it

  MPI_Barrier(world);
  double start_time = MPI_Wtime();

  if (me == 0 && screen) fprintf(screen,"System init for balancing work load ...\n");
  if (me == 0 && screen) fprintf(screen,"System init for balancing work load ...\n");
  cac->init();

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();
  comm->exchange();
  if (atom->map_style) atom->map_set();
  if (element->map_style) element->map_set();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2nodex(element->nlocal,element->nodex);
  }

  if (me == 0 && screen) fprintf(screen,"Balancing work load ...\n");
  if (me == 0 && logfile) fprintf(logfile,"Balancing work load ...\n");

  // parse required arguments

  if (narg < 2) error->all(FLERR,"Illegal balance command");

  thresh = universe->numeric(FLERR,arg[0]);

  int dimension = domain->dimension;
  int *procgrid = comm->procgrid;
  style = -1;
  xflag = yflag = zflag = NONE;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        xflag = UNIFORM;
        iarg += 2;
      } else {
        if (procgrid[0] > narg)
          error->all(FLERR,"Illegal balance command");
        xflag = USER;
        delete [] user_xsplit;
        user_xsplit = new double[procgrid[0]+1];
        user_xsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[0]; i++)
          user_xsplit[i] = universe->numeric(FLERR,arg[iarg++]);
        user_xsplit[procgrid[0]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        yflag = UNIFORM;
        iarg += 2;
      } else {
        if (procgrid[1] > narg)
          error->all(FLERR,"Illegal balance command");
        yflag = USER;
        delete [] user_ysplit;
        user_ysplit = new double[procgrid[1]+1];
        user_ysplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[1]; i++)
          user_ysplit[i] = universe->numeric(FLERR,arg[iarg++]);
        user_ysplit[procgrid[1]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        zflag = UNIFORM;
        iarg += 2;
      } else {
        if (procgrid[2] > narg)
          error->all(FLERR,"Illegal balance command");
        zflag = USER;
        delete [] user_zsplit;
        user_zsplit = new double[procgrid[2]+1];
        user_zsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[2]; i++)
          user_zsplit[i] = universe->numeric(FLERR,arg[iarg++]);
        user_zsplit[procgrid[2]] = 1.0;
      }

    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (style != -1) error->all(FLERR,"Illegal balance command");
      if (iarg+4 > narg) error->all(FLERR,"Illegal balance command");
      style = SHIFT;
      if (strlen(arg[iarg+1]) > 3) error->all(FLERR,"Illegal balance command");
      strcpy(bstr,arg[iarg+1]);
      nitermax = universe->inumeric(FLERR,arg[iarg+2]);
      if (nitermax <= 0) error->all(FLERR,"Illegal balance command");
      stopthresh = universe->numeric(FLERR,arg[iarg+3]);
      if (stopthresh < 1.0) error->all(FLERR,"Illegal balance command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"rcb") == 0) {
      if (style != -1) error->all(FLERR,"Illegal balance command");
      style = BISECTION;
      iarg++;

    } else break;
  }

  // error checks

  if (style == XYZ) {
    if (zflag != NONE  && dimension == 2)
      error->all(FLERR,"Cannot balance in z dimension for 2d simulation");

    if (xflag == USER)
      for (int i = 1; i <= procgrid[0]; i++)
        if (user_xsplit[i-1] >= user_xsplit[i])
          error->all(FLERR,"Illegal balance command");
    if (yflag == USER)
      for (int i = 1; i <= procgrid[1]; i++)
        if (user_ysplit[i-1] >= user_ysplit[i])
          error->all(FLERR,"Illegal balance command");
    if (zflag == USER)
      for (int i = 1; i <= procgrid[2]; i++)
        if (user_zsplit[i-1] >= user_zsplit[i])
          error->all(FLERR,"Illegal balance command");
  }

  if (style == SHIFT) {
    const int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Balance shift string is invalid");
    }
  }

  if (style == BISECTION && comm->style == 0)
    error->all(FLERR,"Balance rcb cannot be used with comm_style brick");

  // process remaining optional args

  options(iarg,narg,arg);

  // set weight for elements

  double maxinit,mininit;
  double init_imbalance = imbalance_factor(maxinit,mininit);

  // no load-balance if imbalance doesn't exceed threshold
  // unless switching from tiled to non tiled layout, then force rebalance

  if (comm->layout == Comm::LAYOUT_TILED && style != BISECTION) {
  } else if (init_imbalance < thresh) return;

  // debug output of initial state

#ifdef BALANCE_DEBUG
  if (outflag) dumpout(update->ntimestep);
#endif

  int niter = 0;

  // perform load-balance
  // style XYZ = explicit setting of cutting planes of logical 3d grid

  if (style == XYZ) {
    if (comm->layout == Comm::LAYOUT_UNIFORM) {
      if (xflag == USER || yflag == USER || zflag == USER)
        comm->layout = Comm::LAYOUT_NONUNIFORM;
    } else if (comm->style == Comm::LAYOUT_NONUNIFORM) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
        comm->layout = Comm::LAYOUT_UNIFORM;
    } else if (comm->style == Comm::LAYOUT_TILED) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
        comm->layout = Comm::LAYOUT_UNIFORM;
      else comm->layout = Comm::LAYOUT_NONUNIFORM;
    }

    if (xflag == UNIFORM) {
      for (int i = 0; i < procgrid[0]; i++)
        comm->xsplit[i] = i * 1.0/procgrid[0];
      comm->xsplit[procgrid[0]] = 1.0;
    } else if (xflag == USER)
      for (int i = 0; i <= procgrid[0]; i++) comm->xsplit[i] = user_xsplit[i];

    if (yflag == UNIFORM) {
      for (int i = 0; i < procgrid[1]; i++)
        comm->ysplit[i] = i * 1.0/procgrid[1];
      comm->ysplit[procgrid[1]] = 1.0;
    } else if (yflag == USER)
      for (int i = 0; i <= procgrid[1]; i++) comm->ysplit[i] = user_ysplit[i];

    if (zflag == UNIFORM) {
      for (int i = 0; i < procgrid[2]; i++)
        comm->zsplit[i] = i * 1.0/procgrid[2];
      comm->zsplit[procgrid[2]] = 1.0;
    } else if (zflag == USER)
      for (int i = 0; i <= procgrid[2]; i++) comm->zsplit[i] = user_zsplit[i];
  }

  // style SHIFT = adjust cutting planes of logical 3d grid

  if (style == SHIFT) {
    comm->layout = Comm::LAYOUT_NONUNIFORM;
    shift_setup_static(bstr);
    niter = shift();
  }

  // style BISECTION = recursive coordinate bisectioning

  if (style == BISECTION) {
    comm->layout = Comm::LAYOUT_TILED;
    bisection(1);
  }

  // reset proc sub-domains
  // for either brick or tiled comm style

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();

  // move particles to new processors via IrregularComm

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
  }

  IrregularComm *irr_comm = new IrregularComm(cac);

  if (style == BISECTION) irr_comm->migrate(1,1,rcb->sendproc);
  else irr_comm->migrate(1);

  delete irr_comm;

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
  }

  // output of final result

  if (outflag) dumpout(update->ntimestep);

  // check if any atoms/elements were lost

  bigint n;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&n,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (n != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via balance: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,atom->natoms,n);
    error->all(FLERR,str);
  }


  nblocal = element->nlocal;
  MPI_Allreduce(&nblocal,&n,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (n != element->nelements) {
    char str[128];
    sprintf(str,"Lost elements via balance: original " BIGINT_FORMAT
        " current " BIGINT_FORMAT,element->nelements,n);
    error->all(FLERR,str);
  }

  // final_imbalance = final imbalance

  double maxfinal,minfinal;
  double final_imbalance = imbalance_factor(maxfinal,minfinal);

  // stats output

  double stop_time = MPI_Wtime();

  FILE *out;
  if (me == 0) {
    for (int i = 0; i < 2; i++) {
      if (i == 0) out = screen;
      else out = logfile;
      if (out) {
        fprintf(out,"  rebalancing time: %g seconds\n",stop_time-start_time);
        fprintf(out,"  iteration count = %d\n",niter);
        fprintf(out,"  initial/final max load/proc = %g %g\n",
            maxinit,maxfinal);
        fprintf(out,"  initial/final min load/proc = %g %g\n",
            mininit,minfinal);
        fprintf(out,"  initial/final imbalance factor = %g %g\n",
            init_imbalance,final_imbalance);

        if (style != BISECTION) {
          fprintf(out,"  x cuts:");
          for (int i = 0; i <= comm->procgrid[0]; i++)
            fprintf(out," %g",comm->xsplit[i]);
          fprintf(out,"\n");
          fprintf(out,"  y cuts:");
          for (int i = 0; i <= comm->procgrid[1]; i++)
            fprintf(out," %g",comm->ysplit[i]);
          fprintf(out,"\n");
          fprintf(out,"  z cuts:");
          for (int i = 0; i <= comm->procgrid[2]; i++)
            fprintf(out," %g",comm->zsplit[i]);
          fprintf(out,"\n");
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   process optional command args for Balance and FixBalance
   ------------------------------------------------------------------------- */

int Balance::options(int iarg, int narg, char **arg)
{
  int screenflag = 0;
  outflag = 0;
  int outarg = 0;
  fp = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal (fix) balance command");
      if (strcmp(arg[iarg+1],"pizza") == 0) outflag = 1;
      else if (strcmp(arg[iarg+1],"tecplot") == 0) outflag = 2;
      else error->all(FLERR,"Illegal (fix) balance command");
      outarg = iarg+2;
      iarg += 3;
    } else if (strcmp(arg[iarg],"eweight") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal (fix) balance command");
      user_element_weight = universe->numeric(FLERR,arg[iarg+1]);
      if (user_element_weight < 0) error->all(FLERR,"Illegal (fix) balance command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) screenflag = 0;
      else error->all(FLERR,"Illegal (fix) balance command");
      iarg += 2;
    } else error->all(FLERR,"Illegal (fix) balance command");
  }

  // output file

  if (outflag && me == 0) {
    fp = fopen(arg[outarg],"w");
    if (fp == NULL) error->one(FLERR,"Cannot open (fix) balance output file");
    if (outflag == 2) {
      fprintf(fp,"title=\"Tiled Proc Map\"\n");
      fprintf(fp,"variables=\"x\", \"y\", \"z\"\n");
    }
  }

  return screenflag;
}


/* ----------------------------------------------------------------------
   calculate imbalance factor based on particle count or particle weights
   return max = max load per proc
   return imbalance = max load per proc / ave load per proc
   ------------------------------------------------------------------------- */

double Balance::imbalance_factor(double &maxcost, double &mincost)
{
  double mycost,totalcost;
  mycost = atom->nlocal;

  if (user_element_weight) 
    mycost += user_element_weight * element->nlocal;
  else {
    int *nintg = element->nintg;
    int *etype = element->etype;
    for (int i = 0; i < element->nlocal; i++)
      mycost += nintg[etype[i]];
  }


  MPI_Allreduce(&mycost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&mycost,&mincost,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&mycost,&totalcost,1,MPI_DOUBLE,MPI_SUM,world);


  double imbalance = 1.0;

  if (maxcost > 0.0) imbalance = maxcost / (totalcost/nprocs);

  return imbalance;

}

/* ----------------------------------------------------------------------
   perform balancing via RCB class
   sortflag = flag for sorting order of received messages by proc ID
   return list of procs to send my atoms to
   ------------------------------------------------------------------------- */

int *Balance::bisection(int sortflag)
{
  if (!rcb) rcb = new RCB(cac);

  // NOTE: this logic is specific to orthogonal boxes, not triclinic

  int dim = domain->dimension;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;

  // shrink-wrap simulation box around atoms/nodes for input to RCB
  // leads to better-shaped sub-boxes when atoms are far from box boundaries

  double shrink[6],shrinkall[6];
  shrink[0] = boxhi[0]; shrink[1] = boxhi[1]; shrink[2] = boxhi[2];
  shrink[3] = boxlo[0]; shrink[4] = boxlo[1]; shrink[5] = boxlo[2];

  double **x = atom->x;
  for (int i = 0; i < atom->nlocal; i++) {
    shrink[0] = MIN(shrink[0],x[i][0]);
    shrink[1] = MIN(shrink[1],x[i][1]);
    shrink[2] = MIN(shrink[2],x[i][2]);
    shrink[3] = MAX(shrink[3],x[i][0]);
    shrink[4] = MAX(shrink[4],x[i][1]);
    shrink[5] = MAX(shrink[5],x[i][2]);
  }

  double ***nodex = element->nodex;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < element->nlocal; i++) 
    for (int j = 0; j < npe[etype[i]]; j++) {
      shrink[0] = MIN(shrink[0],nodex[i][j][0]);
      shrink[1] = MIN(shrink[1],nodex[i][j][1]);
      shrink[2] = MIN(shrink[2],nodex[i][j][2]);
      shrink[3] = MAX(shrink[3],nodex[i][j][0]);
      shrink[4] = MAX(shrink[4],nodex[i][j][1]);
      shrink[5] = MAX(shrink[5],nodex[i][j][2]);
    }

  shrink[3] = -shrink[3]; shrink[4] = -shrink[4]; shrink[5] = -shrink[5];
  MPI_Allreduce(shrink,shrinkall,6,MPI_DOUBLE,MPI_MIN,world);
  shrinkall[3] = -shrinkall[3];
  shrinkall[4] = -shrinkall[4];
  shrinkall[5] = -shrinkall[5];

  double *shrinklo = &shrinkall[0];
  double *shrinkhi = &shrinkall[3];

  // merge atom coords and element coords into dots array
  // and assign appropriate weights (element coords go first)
  // if no element, ignore per-particle weights (wt = NULL)
  // otherwise, grow it to at least size of 1

  int ntotal = atom->nlocal + element->nlocal;
  double **dots;
  double *wt;

  memory->create(dots,ntotal,3,"balance:dots");
  if (element->nelements) {
    if (ntotal) memory->create(wt,ntotal,"balance:wt");
    else memory->create(wt,1,"balance:wt");
  } else wt = NULL;

  x = element->x;
  for (int i = 0; i < element->nlocal; i++) {
    dots[i][0] = x[i][0]; 
    dots[i][1] = x[i][1]; 
    dots[i][2] = x[i][2]; 
    if (wt) {
      if (user_element_weight) wt[i] = user_element_weight;
      else wt[i] = element->nintg[element->etype[i]];
    }
  }

  x = atom->x;
  int j = 0;
  for (int i = element->nlocal; i < ntotal; i++) {
    dots[i][0] = x[j][0]; 
    dots[i][1] = x[j][1]; 
    dots[i][2] = x[j][2]; 
    j++;
    if (wt) wt[i] = 1.0;
  }

  // invoke RCB
  // then invert() to create list of proc assignments for my atoms

  rcb->compute(dim,ntotal,dots,wt,shrinklo,shrinkhi);

  rcb->invert(sortflag);

  // reset RCB lo/hi bounding box to full simulation box as needed

  double *lo = rcb->lo;
  double *hi = rcb->hi;

  if (lo[0] == shrinklo[0]) lo[0] = boxlo[0];
  if (lo[1] == shrinklo[1]) lo[1] = boxlo[1];
  if (lo[2] == shrinklo[2]) lo[2] = boxlo[2];
  if (hi[0] == shrinkhi[0]) hi[0] = boxhi[0];
  if (hi[1] == shrinkhi[1]) hi[1] = boxhi[1];
  if (hi[2] == shrinkhi[2]) hi[2] = boxhi[2];

  // store RCB cut, dim, lo/hi box in CommTiled
  // cut and lo/hi need to be in fractional form so can
  // OK if changes by epsilon from what RCB used since atoms
  //   will subsequently migrate to new owning procs by exchange() anyway
  // ditto for atoms exactly on lo/hi RCB box boundaries due to ties

  comm->rcbnew = 1;

  int idim = rcb->cutdim;
  if (idim >= 0) comm->rcbcutfrac = (rcb->cut - boxlo[idim]) / prd[idim];
  else comm->rcbcutfrac = 0.0;
  comm->rcbcutdim = idim;

  double (*mysplit)[2] = comm->mysplit;

  mysplit[0][0] = (lo[0] - boxlo[0]) / prd[0];
  if (hi[0] == boxhi[0]) mysplit[0][1] = 1.0;
  else mysplit[0][1] = (hi[0] - boxlo[0]) / prd[0];

  mysplit[1][0] = (lo[1] - boxlo[1]) / prd[1];
  if (hi[1] == boxhi[1]) mysplit[1][1] = 1.0;
  else mysplit[1][1] = (hi[1] - boxlo[1]) / prd[1];

  mysplit[2][0] = (lo[2] - boxlo[2]) / prd[2];
  if (hi[2] == boxhi[2]) mysplit[2][1] = 1.0;
  else mysplit[2][1] = (hi[2] - boxlo[2]) / prd[2];

  // clean up 

  memory->destroy(dots);
  memory->destroy(wt); 

  // return list of procs to send my atoms to

  return rcb->sendproc;
}

/* ----------------------------------------------------------------------
   setup static load balance operations
   called from command and indirectly initially from fix balance
   set rho = 0 for static balancing
   ------------------------------------------------------------------------- */

void Balance::shift_setup_static(char *str)
{
  shift_allocate = 1;

  memory->create(proccost,nprocs,"balance:proccost");
  memory->create(allproccost,nprocs,"balance:allproccost");

  ndim = strlen(str);
  bdim = new int[ndim];

  for (int i = 0; i < ndim; i++) {
    if (str[i] == 'x') bdim[i] = X;
    if (str[i] == 'y') bdim[i] = Y;
    if (str[i] == 'z') bdim[i] = Z;
  }

  // max = max number of processors in each direction

  int max = MAX(comm->procgrid[0],comm->procgrid[1]);
  max = MAX(max,comm->procgrid[2]);

  onecost = new double[max];
  allcost = new double[max];
  sum = new double[max+1];
  target = new double[max+1];
  lo = new double[max+1];
  hi = new double[max+1];
  losum = new double[max+1];
  hisum = new double[max+1];

  // if current layout is TILED, set initial uniform splits in Comm
  // this gives starting point to subsequent shift balancing

  if (comm->layout == Comm::LAYOUT_TILED) {
    int *procgrid = comm->procgrid;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;

    for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
    for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
    for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];
    xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;
  }

  rho = 0;
}

/* ----------------------------------------------------------------------
   setup shift load balance operations
   called from fix balance
   set rho = 1 to do dynamic balancing after call to shift_setup_static()
   ------------------------------------------------------------------------- */

void Balance::shift_setup(char *str, int nitermax_in, double thresh_in)
{
  shift_setup_static(str);
  nitermax = nitermax_in;
  stopthresh = thresh_in;
  rho = 1;
}

/* ----------------------------------------------------------------------
   load balance by changing xyz split proc boundaries in Comm
   called one time from input script command or many times from fix balance
   return niter = iteration count
   ------------------------------------------------------------------------- */

int Balance::shift()
{
  int i,j,k,m,np;
  double mycost,totalcost;
  double *split;

  // no balancing if no atom and element

  bigint natoms = atom->natoms;
  bigint nelements = element->nelements;
  if (natoms+nelements == 0) return 0;

  // set delta for 1d balancing = root of threshold
  // root = # of dimensions being balanced on

  double delta = pow(stopthresh,1.0/ndim) - 1.0;
  int *procgrid = comm->procgrid;

  // all balancing done in lamda coords
  // only use element center coords so no need to convert
  // node coords to lamda coords

  domain->x2lamda(atom->nlocal,atom->x);
  domain->x2lamda(element->nlocal,element->x);

  // loop over dimensions in balance string

  int niter = 0;
  for (int idim = 0; idim < ndim; idim++) {

    // split = ptr to xyz split in Comm

    if (bdim[idim] == X) split = comm->xsplit;
    else if (bdim[idim] == Y) split = comm->ysplit;
    else if (bdim[idim] == Z) split = comm->zsplit;
    else continue;

    // initial count and sum

    np = procgrid[bdim[idim]];
    tally(bdim[idim],np,split);

    // target[i] = desired sum at split I

    mycost = atom->nlocal;
    if (user_element_weight) 
      mycost += user_element_weight * element->nlocal;
    else {
      int *nintg = element->nintg;
      int *etype = element->etype;
      for (i = 0; i < element->nlocal; i++)
        mycost += nintg[etype[i]];
    }
    MPI_Allreduce(&mycost,&totalcost,1,MPI_DOUBLE,MPI_SUM,world);

    for (i = 0; i < np; i++) target[i] = totalcost/np * i;
    target[np] = totalcost;

    // lo[i] = closest split <= split[i] with a sum <= target
    // hi[i] = closest split >= split[i] with a sum >= target

    lo[0] = hi[0] = 0.0;
    lo[np] = hi[np] = 1.0;
    losum[0] = hisum[0] = 0.0;
    losum[np] = hisum[np] = totalcost;

    for (i = 1; i < np; i++) {
      for (j = i; j >= 0; j--)
        if (sum[j] <= target[i]) {
          lo[i] = split[j];
          losum[i] = sum[j];
          break;
        }
      for (j = i; j <= np; j++)
        if (sum[j] >= target[i]) {
          hi[i] = split[j];
          hisum[i] = sum[j];
          break;
        }
    }

    // iterate until balanced

#ifdef BALANCE_DEBUG
    if (me == 0) debug_shift_output(idim,0,np,split);
#endif

    int doneflag;
    int change = 1;
    for (m = 0; m < nitermax; m++) {
      change = adjust(np,split);
      tally(bdim[idim],np,split);
      niter++;

#ifdef BALANCE_DEBUG
      if (me == 0) debug_shift_output(idim,m+1,np,split);
      if (outflag) dumpout(update->ntimestep);
#endif

      // stop if no change in splits, b/c all targets are met exactly

      if (!change) break;

      // stop if all split sums are within delta of targets
      // this is a 1d test of particle count per slice
      // assumption is that this is sufficient accuracy
      //   for 3d imbalance factor to reach threshold

      doneflag = 1;
      for (i = 1; i < np; i++)
        if (fabs(1.0*(sum[i]-target[i]))/target[i] > delta) doneflag = 0;
      if (doneflag) break;
    }

    // eliminate final adjacent splits that are duplicates
    // can happen if particle distribution is narrow and Nitermax is small
    // set lo = midpt between splits
    // spread duplicates out evenly between bounding midpts with non-duplicates
    // i,j = lo/hi indices of set of duplicate splits
    // delta = new spacing between duplicates
    // bounding midpts = lo[i-1] and lo[j]

    int duplicate = 0;
    for (i = 1; i < np-1; i++)
      if (split[i] == split[i+1]) duplicate = 1;
    if (duplicate) {
      for (i = 0; i < np; i++)
        lo[i] = 0.5 * (split[i] + split[i+1]);
      i = 1;
      while (i < np-1) {
        j = i + 1;
        while (split[j] == split[i]) j++;
        j--;
        if (j > i) {
          delta = (lo[j] - lo[i-1]) / (j-i+2);
          for (k = i; k <= j; k++)
            split[k] = lo[i-1] + (k-i+1)*delta;
        }
        i = j + 1;
      }
    }

    // sanity check on bad duplicate or inverted splits
    // zero or negative width sub-domains will break Comm class
    // should never happen if recursive multisection algorithm is correct

    int bad = 0;
    for (i = 0; i < np; i++)
      if (split[i] >= split[i+1]) bad = 1;
    if (bad) error->all(FLERR,"Balance produced bad splits");
    /*
       if (me == 0) {
       printf("BAD SPLITS %d %d %d\n",np+1,niter,delta);
       for (i = 0; i < np+1; i++)
       printf(" %g",split[i]);
       printf("\n");
       }
       */

    // stop at this point in bstr if imbalance factor < threshold
    // this is a true 3d test of particle count per processor

    if (imbalance_splits() <= stopthresh) break;
  }

  // restore real coords

  domain->lamda2x(atom->nlocal,atom->x);
  domain->lamda2x(element->nlocal,element->x);

  return niter;
}

/* ----------------------------------------------------------------------
   count atoms in each slice, based on their dim coordinate
   N = # of slices
   split = N+1 cuts between N slices
   return updated count = particles per slice
   return updated sum = cumulative count below each of N+1 splits
   use binary search to find which slice each atom is in
   ------------------------------------------------------------------------- */

void Balance::tally(int dim, int n, double *split)
{
  for (int i = 0; i < n; i++) onecost[i] = 0.0;
  int index;

  double **x;

  x = atom->x;
  for (int i = 0; i < atom->nlocal; i++) {
    index = binary_search(x[i][dim],n,split);
    onecost[index] += 1.0;
  }

  x = element->x;
  if (user_element_weight) {
    for (int i = 0; i < element->nlocal; i++) {
      index = binary_search(x[i][dim],n,split);
      onecost[index] += user_element_weight;
    }
  } else {
    int *nintg = element->nintg;
    int *etype = element->etype;
    for (int i = 0; i < element->nlocal; i++) {
      index = binary_search(x[i][dim],n,split);
      onecost[index] += nintg[etype[i]];
    }
  }

  MPI_Allreduce(onecost,allcost,n,MPI_DOUBLE,MPI_SUM,world);

  sum[0] = 0.0;
  for (int i = 1; i < n+1; i++)
    sum[i] = sum[i-1] + allcost[i-1];
}

/* ----------------------------------------------------------------------
   adjust cuts between N slices in a dim via recursive multisectioning method
   split = current N+1 cuts, with 0.0 and 1.0 at end points
   sum = cumulative count up to each split
   target = desired cumulative count up to each split
   lo/hi = split values that bound current split
   update lo/hi to reflect sums at current split values
   overwrite split with new cuts
   guaranteed that splits will remain in ascending order,
   though adjacent values may be identical
   recursive bisectioning zooms in on each cut by halving lo/hi
   return 0 if no changes in any splits, b/c they are all perfect
   ------------------------------------------------------------------------- */

int Balance::adjust(int n, double *split)
{
  int i;
  double fraction;

  // reset lo/hi based on current sum and splits
  // insure lo is monotonically increasing, ties are OK
  // insure hi is monotonically decreasing, ties are OK
  // this effectively uses info from nearby splits
  // to possibly tighten bounds on lo/hi

  for (i = 1; i < n; i++) {
    if (sum[i] <= target[i]) {
      lo[i] = split[i];
      losum[i] = sum[i];
    }
    if (sum[i] >= target[i]) {
      hi[i] = split[i];
      hisum[i] = sum[i];
    }
  }
  for (i = 1; i < n; i++)
    if (lo[i] < lo[i-1]) {
      lo[i] = lo[i-1];
      losum[i] = losum[i-1];
    }
  for (i = n-1; i > 0; i--)
    if (hi[i] > hi[i+1]) {
      hi[i] = hi[i+1];
      hisum[i] = hisum[i+1];
    }

  int change = 0;
  for (int i = 1; i < n; i++)
    if (sum[i] != target[i]) {
      change = 1;
      if (rho == 0) split[i] = 0.5 * (lo[i]+hi[i]);
      else {
        fraction = 1.0*(target[i]-losum[i]) / (hisum[i]-losum[i]);
        split[i] = lo[i] + fraction * (hi[i]-lo[i]);
      }
    }
  return change;
}

/* ----------------------------------------------------------------------
   calculate imbalance based on processor splits in 3 dims
   atoms must be in lamda coords (0-1) before called
   map particles to 3d grid of procs
   return imbalance factor = max load per proc / ave load per proc
   ------------------------------------------------------------------------- */

double Balance::imbalance_splits()
{
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  int nx = comm->procgrid[0];
  int ny = comm->procgrid[1];
  int nz = comm->procgrid[2];

  for (int i = 0; i < nprocs; i++) proccost[i] = 0.0;

  double **x;
  int ix,iy,iz;

  x = atom->x;
  for (int i = 0; i < atom->nlocal; i++) {
    ix = binary_search(x[i][0],nx,xsplit);
    iy = binary_search(x[i][1],ny,ysplit);
    iz = binary_search(x[i][2],nz,zsplit);
    proccost[iz*nx*ny + iy*nx + ix] += 1.0;
  }

  x = element->x;
  if (user_element_weight) {
    for (int i = 0; i < element->nlocal; i++) {
      ix = binary_search(x[i][0],nx,xsplit);
      iy = binary_search(x[i][1],ny,ysplit);
      iz = binary_search(x[i][2],nz,zsplit);
      proccost[iz*nx*ny + iy*nx + ix] += user_element_weight;
    }
  } else {
    int *nintg = element->nintg;
    int *etype = element->etype;
    for (int i = 0; i < element->nlocal; i++) {
      ix = binary_search(x[i][0],nx,xsplit);
      iy = binary_search(x[i][1],ny,ysplit);
      iz = binary_search(x[i][2],nz,zsplit);
      proccost[iz*nx*ny + iy*nx + ix] += nintg[etype[i]];
    }
  }

  // one proc's particles may map to many partitions, so must Allreduce

  MPI_Allreduce(proccost,allproccost,nprocs,MPI_DOUBLE,MPI_SUM,world);

  double maxcost = 0.0;
  double totalcost = 0.0;
  for (int i = 0; i < nprocs; i++) {
    maxcost = MAX(maxcost,allproccost[i]);
    totalcost += allproccost[i];
  }

  double imbalance = 1.0;
  if (maxcost > 0.0) imbalance = maxcost / (totalcost/nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   binary search for where value falls in N-length vec
   note that vec actually has N+1 values, but ignore last one
   values in vec are monotonically increasing, but adjacent values can be ties
   value may be outside range of vec limits
   always return index from 0 to N-1 inclusive
   return 0 if value < vec[0]
   reutrn N-1 if value >= vec[N-1]
   return index = 1 to N-2 inclusive if vec[index] <= value < vec[index+1]
   note that for adjacent tie values, index of lower tie is not returned
   since never satisfies 2nd condition that value < vec[index+1]
   ------------------------------------------------------------------------- */

int Balance::binary_search(double value, int n, double *vec)
{
  int lo = 0;
  int hi = n-1;

  if (value < vec[lo]) return lo;
  if (value >= vec[hi]) return hi;

  // insure vec[lo] <= value < vec[hi] at every iteration
  // done when lo,hi are adjacent

  int index = (lo+hi)/2;
  while (lo < hi-1) {
    if (value < vec[index]) hi = index;
    else if (value >= vec[index]) lo = index;
    index = (lo+hi)/2;
  }

  return index;
}

/* ----------------------------------------------------------------------
   write dump snapshot of line segments in Pizza.py mdump mesh format
   write xy lines around each proc's sub-domain for 2d
   write xyz cubes around each proc's sub-domain for 3d
   only called by proc 0
NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */

void Balance::dumpout(bigint tstep)
{
  int dimension = domain->dimension;
  int triclinic = domain->triclinic;

  // Allgather each proc's sub-box
  // could use Gather, but that requires MPI to alloc memory

  double *lo,*hi;
  if (triclinic == 0) {
    lo = domain->sublo;
    hi = domain->subhi;
  } else {
    lo = domain->sublo_lamda;
    hi = domain->subhi_lamda;
  }

  double box[6];
  box[0] = lo[0]; box[1] = lo[1]; box[2] = lo[2];
  box[3] = hi[0]; box[4] = hi[1]; box[5] = hi[2];

  double **boxall;
  memory->create(boxall,nprocs,6,"balance:dumpout");
  MPI_Allgather(box,6,MPI_DOUBLE,&boxall[0][0],6,MPI_DOUBLE,world);

  if (me) {
    memory->destroy(boxall);
    return;
  }

  if (outflag == 1) {
    // proc 0 writes out nodal coords
    // some will be duplicates

    double *boxlo = domain->boxlo;
    double *boxhi = domain->boxhi;

    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",tstep);
    fprintf(fp,"ITEM: NUMBER OF NODES\n");
    if (dimension == 2) fprintf(fp,"%d\n",4*nprocs);
    else fprintf(fp,"%d\n",8*nprocs);
    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxlo[0],boxhi[0]);
    fprintf(fp,"%g %g\n",boxlo[1],boxhi[1]);
    fprintf(fp,"%g %g\n",boxlo[2],boxhi[2]);
    fprintf(fp,"ITEM: NODES\n");

    if (triclinic == 0) {
      if (dimension == 2) {
        int m = 0;
        for (int i = 0; i < nprocs; i++) {
          fprintf(fp,"%d %d %g %g %g\n",m+1,1,boxall[i][0],boxall[i][1],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+2,1,boxall[i][3],boxall[i][1],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+3,1,boxall[i][3],boxall[i][4],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+4,1,boxall[i][0],boxall[i][4],0.0);
          m += 4;
        }
      } else {
        int m = 0;
        for (int i = 0; i < nprocs; i++) {
          fprintf(fp,"%d %d %g %g %g\n",m+1,1,
              boxall[i][0],boxall[i][1],boxall[i][2]);
          fprintf(fp,"%d %d %g %g %g\n",m+2,1,
              boxall[i][3],boxall[i][1],boxall[i][2]);
          fprintf(fp,"%d %d %g %g %g\n",m+3,1,
              boxall[i][3],boxall[i][4],boxall[i][2]);
          fprintf(fp,"%d %d %g %g %g\n",m+4,1,
              boxall[i][0],boxall[i][4],boxall[i][2]);
          fprintf(fp,"%d %d %g %g %g\n",m+5,1,
              boxall[i][0],boxall[i][1],boxall[i][5]);
          fprintf(fp,"%d %d %g %g %g\n",m+6,1,
              boxall[i][3],boxall[i][1],boxall[i][5]);
          fprintf(fp,"%d %d %g %g %g\n",m+7,1,
              boxall[i][3],boxall[i][4],boxall[i][5]);
          fprintf(fp,"%d %d %g %g %g\n",m+8,1,
              boxall[i][0],boxall[i][4],boxall[i][5]);
          m += 8;
        }
      }

    } else {
      double (*bc)[3] = domain->corners;

      if (dimension == 2) {
        int m = 0;
        for (int i = 0; i < nprocs; i++) {
          domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
          fprintf(fp,"%d %d %g %g %g\n",m+1,1,bc[0][0],bc[0][1],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+2,1,bc[1][0],bc[1][1],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+3,1,bc[2][0],bc[2][1],0.0);
          fprintf(fp,"%d %d %g %g %g\n",m+4,1,bc[3][0],bc[3][1],0.0);
          m += 4;
        }
      } else {
        int m = 0;
        for (int i = 0; i < nprocs; i++) {
          domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
          fprintf(fp,"%d %d %g %g %g\n",m+1,1,bc[0][0],bc[0][1],bc[0][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+2,1,bc[1][0],bc[1][1],bc[1][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+3,1,bc[2][0],bc[2][1],bc[2][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+4,1,bc[3][0],bc[3][1],bc[3][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+5,1,bc[4][0],bc[4][1],bc[4][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+6,1,bc[5][0],bc[5][1],bc[5][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+7,1,bc[6][0],bc[6][1],bc[6][1]);
          fprintf(fp,"%d %d %g %g %g\n",m+8,1,bc[7][0],bc[7][1],bc[7][1]);
          m += 8;
        }
      }
    }

    // write out one square/cube per processor for 2d/3d

    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",tstep);
    if (dimension == 2) fprintf(fp,"ITEM: NUMBER OF SQUARES\n");
    else fprintf(fp,"ITEM: NUMBER OF CUBES\n");
    fprintf(fp,"%d\n",nprocs);
    if (dimension == 2) fprintf(fp,"ITEM: SQUARES\n");
    else fprintf(fp,"ITEM: CUBES\n");

    if (dimension == 2) {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %d %d %d %d\n",i+1,1,m+1,m+2,m+3,m+4);
        m += 4;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %d %d %d %d %d %d %d %d\n",
            i+1,1,m+1,m+2,m+3,m+4,m+5,m+6,m+7,m+8);
        m += 8;
      }
    }
  } else {
    fprintf(fp,"zone t=\"TIMESTEP " BIGINT_FORMAT "\", n = %d e = %d datapacking = point, zonetype = febrick\n",tstep,
        static_cast<int> (nprocs*pow(2,dimension)),nprocs);
    if (triclinic == 0) {
      if (dimension == 2) 
        for (int i = 0; i < nprocs; i++) {
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][1],0.0);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][1],0.0);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][4],0.0);
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][4],0.0);
        }
      else 
        for (int i = 0; i < nprocs; i++) {
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][1],boxall[i][2]);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][1],boxall[i][2]);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][4],boxall[i][2]);
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][4],boxall[i][2]);
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][1],boxall[i][5]);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][1],boxall[i][5]);
          fprintf(fp,"%g %g %g\n",boxall[i][3],boxall[i][4],boxall[i][5]);
          fprintf(fp,"%g %g %g\n",boxall[i][0],boxall[i][4],boxall[i][5]);
        }

    } else {
      double (*bc)[3] = domain->corners;

      if (dimension == 2) 
        for (int i = 0; i < nprocs; i++) {
          domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
          fprintf(fp,"%g %g %g\n",bc[0][0],bc[0][1],0.0);
          fprintf(fp,"%g %g %g\n",bc[1][0],bc[1][1],0.0);
          fprintf(fp,"%g %g %g\n",bc[2][0],bc[2][1],0.0);
          fprintf(fp,"%g %g %g\n",bc[3][0],bc[3][1],0.0);
        }
      else 
        for (int i = 0; i < nprocs; i++) {
          domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
          fprintf(fp,"%g %g %g\n",bc[0][0],bc[0][1],bc[0][1]);
          fprintf(fp,"%g %g %g\n",bc[1][0],bc[1][1],bc[1][1]);
          fprintf(fp,"%g %g %g\n",bc[2][0],bc[2][1],bc[2][1]);
          fprintf(fp,"%g %g %g\n",bc[3][0],bc[3][1],bc[3][1]);
          fprintf(fp,"%g %g %g\n",bc[4][0],bc[4][1],bc[4][1]);
          fprintf(fp,"%g %g %g\n",bc[5][0],bc[5][1],bc[5][1]);
          fprintf(fp,"%g %g %g\n",bc[6][0],bc[6][1],bc[6][1]);
          fprintf(fp,"%g %g %g\n",bc[7][0],bc[7][1],bc[7][1]);
        }
    }

    if (dimension == 2) {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %d %d\n",m+1,m+2,m+3,m+4);
        m += 4;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
            m+1,m+2,m+3,m+4,m+5,m+6,m+7,m+8);
        m += 8;
      }
    }

  }
  memory->destroy(boxall);
}

/* ----------------------------------------------------------------------
   debug output for Idim and count
   only called by proc 0
   ------------------------------------------------------------------------- */

#ifdef BALANCE_DEBUG
void Balance::debug_shift_output(int idim, int m, int np, double *split)
{
  int i;
  const char *dim = NULL;

  double *boxlo = domain->boxlo;
  double *prd = domain->prd;

  if (bdim[idim] == X) dim = "X";
  else if (bdim[idim] == Y) dim = "Y";
  else if (bdim[idim] == Z) dim = "Z";
  fprintf(stderr,"Dimension %s, Iteration %d\n",dim,m);

  fprintf(stderr,"  Count:");
  for (i = 0; i < np; i++) fprintf(stderr," " BIGINT_FORMAT,count[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Sum:");
  for (i = 0; i <= np; i++) fprintf(stderr," " BIGINT_FORMAT,sum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Target:");
  for (i = 0; i <= np; i++) fprintf(stderr," " BIGINT_FORMAT,target[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Actual cut:");
  for (i = 0; i <= np; i++)
    fprintf(stderr," %g",boxlo[bdim[idim]] + split[i]*prd[bdim[idim]]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Split:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",split[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Low:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",lo[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Low-sum:");
  for (i = 0; i <= np; i++) fprintf(stderr," " BIGINT_FORMAT,losum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Hi:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",hi[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Hi-sum:");
  for (i = 0; i <= np; i++) fprintf(stderr," " BIGINT_FORMAT,hisum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Delta:");
  for (i = 0; i < np; i++) fprintf(stderr," %g",split[i+1]-split[i]);
  fprintf(stderr,"\n");

  bigint max = 0;
  for (i = 0; i < np; i++) max = MAX(max,count[i]);
  fprintf(stderr,"  Imbalance factor: %g\n",1.0*max*np/target[np]);
}
#endif
