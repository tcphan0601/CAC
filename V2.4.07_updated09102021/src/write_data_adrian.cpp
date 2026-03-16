#include <mpi.h>
#include <string.h>
#include "write_data_adrian.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "element_vec.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

//enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteDataAdrian::WriteDataAdrian(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data_adrian command in input script
------------------------------------------------------------------------- */

void WriteDataAdrian::command(int narg, char **arg)
{
  if (nprocs>1) error->all(FLERR,"Cannot use more than 1 proc ro write data");
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_adrian command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_data_adrian command");

  // if filename contains a "*", replace with current timestep


  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  if (narg > 1) error->all(FLERR,"Illegal write_data_adrian command");
  if (element->max_apc > 1 && comm->nprocs != 1) error->all(FLERR,"Can only use write_data_adrian command for 1 processor for multiple atoms per node");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch
  /*
     pairflag = II;
     coeffflag = 1;
     int noinit = 0;

     int iarg = 1;
     while (iarg < narg) {
     if (strcmp(arg[iarg],"pair") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal write_data_adrian command");
     if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
     else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
     else error->all(FLERR,"Illegal write_data_adrian command");
     iarg += 2;
     } else if (strcmp(arg[iarg],"noinit") == 0) {
     noinit = 1;
     iarg++;
     } else if (strcmp(arg[iarg],"nocoeff") == 0) {
     coeffflag = 0;
     iarg++;
     } else error->all(FLERR,"Illegal write_data_adrian command");
     }
     */
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data_adrian immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {
  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_data_adrian ...\n");
  cac->init();

  // move atoms to new processors before writing file
  // do setup_pre_exchange to force update of per-atom info if needed
  // enforce PBC in case atoms are outside box
  // call borders() to rebuild atom map since exchange() destroys map

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
	  domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();
  comm->exchange();
  comm->setup_borders();
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
	  domain->lamda2nodex(element->nlocal,element->nodex);
  }

  write(file);

  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteDataAdrian::write(char *file)
{
  // special case where reneighboring is not done in integrator
  //   on timestep data file is written (due to build_once being set)
  // if box is changing, must be reset, else data file will have
  //   wrong box size and atoms will be lost when data file is read
  // other calls to pbc and domain and comm are not made,
  //   b/c they only make sense if reneighboring is actually performed

  //if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write data file

  bigint nalocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nalocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && output->thermo->elemlostflag == Thermo::ERROR)
    error->all(FLERR,"Element count is inconsistent, cannot write data file");


  // sum up bond,angle counts
  // may be different than atom->nbonds,nangles if broken/turned-off
  /*
     if (atom->molecular == 1 && (atom->nbonds || atom->nbondtypes)) {
     nbonds_local = atom->avec->pack_bond(NULL);
     MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_CAC_BIGINT,MPI_SUM,world);
     }
     if (atom->molecular == 1 && (atom->nangles || atom->nangletypes)) {
     nangles_local = atom->avec->pack_angle(NULL);
     MPI_Allreduce(&nangles_local,&nangles,1,MPI_CAC_BIGINT,MPI_SUM,world);
     }
     */
  // open data file

  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open data file %s",file);
      error->one(FLERR,str);
    }
  }

  // proc 0 writes header, ntype-length arrays, force fields

  if (me == 0) {
    header();

    fprintf(fp,"\nCAC Elements\n\n");
  }


  // per element info

  if (natoms)
    atoms();
  if (nelements)
    elements();


  /*

  // extra sections managed by fixes

  for (int i = 0; i < modify->nfix; i++)
  if (modify->fix[i]->wd_section)
  for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i,m);
  */
  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
   ------------------------------------------------------------------------- */

void WriteDataAdrian::header()
{
  fprintf(fp,"CAC data file via write_data_adrian, "
      "timestep = " BIGINT_FORMAT "\n",
      update->ntimestep);

  fprintf(fp,"\n");

  fprintf(fp,BIGINT_FORMAT " cac elements\n",atom->natoms+element->nelements/element->apc[1]);
  fprintf(fp,"%d atom types\n",atom->ntypes);

  fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",domain->boxlo[0],domain->boxhi[0]);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",domain->boxlo[1],domain->boxhi[1]);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",domain->boxlo[2],domain->boxhi[2]);

  if (domain->triclinic)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e xy xz yz\n",
        domain->xy,domain->xz,domain->yz);

}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteDataAdrian::atoms()
{
  tagint maxtag,tagoffset;
  maxtag = 0;
  for (int i = 0; i < element->nelements; i++)
    maxtag = MAX(maxtag,element->tag[i]);
  MPI_Allreduce(&maxtag,&tagoffset,1,MPI_CAC_TAGINT,MPI_MAX,world);
  for (int i = 0; i < atom->natoms; i++)
    fprintf(fp, TAGINT_FORMAT " Atom 1 1 1 1\n1 1 %d %g %g %g\n",atom->tag[i]+tagoffset,atom->type[i],
        atom->x[i][0],
        atom->x[i][1],
        atom->x[i][2]);

}

/* ----------------------------------------------------------------------
   write out Elements section of data file
   ------------------------------------------------------------------------- */

void WriteDataAdrian::elements()
{

  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc

  double **buf;

  int ncol = 5 + 8*3;

  int sendrow = element->nlocal;
  int maxrow;

  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);


  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_adrian:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_adrian:buf");


  // pack my element data into buf

  element->evec->pack_element_data_adrian(buf);
  // write one chunk of elements per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      element->evec->write_element_data_adrian(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);

}


