#include <mpi.h>
#include <string.h>
#include "write_data_ghost_debug.h"
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

enum{IGNORE,WARN,ERROR};                    // same as thermo.cpp
//enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteDataGhostDebug::WriteDataGhostDebug(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data_ghost_debug_ghost_debug command in input script
   all processor write out its own file contain local+ghost atoms/elements/nodes 
------------------------------------------------------------------------- */

void WriteDataGhostDebug::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_ghost_debug command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_data_ghost_debug command");

  // if filename contains a "*", replace with current timestep


  char *ptr;
  int n = strlen(arg[0]) + 20;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s_proc_%d",arg[0],update->ntimestep,ptr+1,me);
  } else sprintf(file,"%s_proc_%d",arg[0],me);

  if (narg > 1) error->all(FLERR,"Illegal write_data_ghost_debug command");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch
  /*
     pairflag = II;
     coeffflag = 1;
     int noinit = 0;

     int iarg = 1;
     while (iarg < narg) {
     if (strcmp(arg[iarg],"pair") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal write_data_ghost_debug command");
     if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
     else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
     else error->all(FLERR,"Illegal write_data_ghost_debug command");
     iarg += 2;
     } else if (strcmp(arg[iarg],"noinit") == 0) {
     noinit = 1;
     iarg++;
     } else if (strcmp(arg[iarg],"nocoeff") == 0) {
     coeffflag = 0;
     iarg++;
     } else error->all(FLERR,"Illegal write_data_ghost_debug command");
     }
     */
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data_ghost_debug immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {
  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_data_ghost_debug ...\n");
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

void WriteDataGhostDebug::write(char *file)
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
  if (natoms != atom->natoms && output->thermo->atomlostflag == ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && output->thermo->elemlostflag == ERROR)
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

  fp = fopen(file,"w");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open data file %s",file);
    error->one(FLERR,str);
  }


  header();

  // per atom info
  
  fprintf(fp,"\nAtoms\n\n");
  if (natoms) 
    atoms();

  // per element info

  if (nelements) {
    elements();

    // reset all node tag before writing node to file

    //element->reset_node_tags(1,0);
    //nodes(nnode_locals);
  }


  fclose(fp);
}

/* ----------------------------------------------------------------------
   proc me writes out data file header
   ------------------------------------------------------------------------- */

void WriteDataGhostDebug::header()
{
  fprintf(fp,"CAC data file via write_data_ghost_debug, "
      "timestep = " BIGINT_FORMAT " from proc %d\n",
      update->ntimestep,me);

  fprintf(fp,"\n");

  int n = atom->nlocal+element->nlocal+atom->nghost+element->nghost;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < element->nlocal+element->nghost; i++)
    n += npe[etype[i]];
  fprintf(fp,BIGINT_FORMAT " atoms\n",n);
  fprintf(fp,"%d atom types\n",3);

  fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",domain->sublo[0],domain->subhi[0]);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",domain->sublo[1],domain->subhi[1]);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",domain->sublo[2],domain->subhi[2]);

  if (domain->triclinic)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e xy xz yz\n",
        domain->xy,domain->xz,domain->yz);

}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteDataGhostDebug::atoms()
{
  // write local+ghost atoms
  
  int *tag = atom->tag;
  double **x = atom->x;
  
  for (int i = 0; i < atom->nlocal+atom->nghost; i++) {
    fprintf(fp,"%d %d %g %g %g\n",tag[i],1,x[i][0],x[i][1],x[i][2]);
  }
}

/* ----------------------------------------------------------------------
   write out Elements section of data file
   ------------------------------------------------------------------------- */

void WriteDataGhostDebug::elements()
{

  // write local+ghost atoms
  
  int *tag = element->tag;
  double **x = element->x;
  double ***nodex = element->nodex;
  int *npe = element->npe;
  int *etype = element->etype;
  
  for (int i = 0; i < element->nlocal+element->nghost; i++) {
    fprintf(fp,"%d %d %g %g %g\n",tag[i],2,x[i][0],x[i][1],x[i][2]);
    for (int j = 0; j < npe[etype[i]]; j++)
      fprintf(fp,"%d %d %g %g %g\n",tag[i],3,nodex[i][j][0],nodex[i][j][1],nodex[i][j][2]);
  }

}


/* ----------------------------------------------------------------------
   write out Nodes section of data file
   ------------------------------------------------------------------------- */

void WriteDataGhostDebug::nodes()
{

}




