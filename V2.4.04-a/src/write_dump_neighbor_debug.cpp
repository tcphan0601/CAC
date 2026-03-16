#include <mpi.h>
#include <string.h>
#include "write_dump_neighbor_debug.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "element_vec.h"
#include "neigh_list.h"
#include "neighbor.h"
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

WriteDumpNeighborDebug::WriteDumpNeighborDebug(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_dump_neighbor_debug_neighbor_debug command in input script
   all processor write out its own file contain local+ghost atoms/elements/nodes 
------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_neighbor_debug command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_dump_neighbor_debug command");

  // if filename contains a "*", replace with current timestep


  char *ptr;
  int n = strlen(arg[0]) + 20;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s_proc_%d",arg[0],update->ntimestep,ptr+1,me);
  } else sprintf(file,"%s_proc_%d",arg[0],me);

  if (narg > 1) error->all(FLERR,"Illegal write_dump_neighbor_debug command");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch
  /*
     pairflag = II;
     coeffflag = 1;
     int noinit = 0;

     int iarg = 1;
     while (iarg < narg) {
     if (strcmp(arg[iarg],"pair") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal write_dump_neighbor_debug command");
     if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
     else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
     else error->all(FLERR,"Illegal write_dump_neighbor_debug command");
     iarg += 2;
     } else if (strcmp(arg[iarg],"noinit") == 0) {
     noinit = 1;
     iarg++;
     } else if (strcmp(arg[iarg],"nocoeff") == 0) {
     coeffflag = 0;
     iarg++;
     } else error->all(FLERR,"Illegal write_dump_neighbor_debug command");
     }
     */
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_dump_neighbor_debug immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {

 
  
  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_dump_neighbor_debug ...\n");
  cac->init();

  // move atoms to new processors before writing file
  // do setup_pre_exchange to force update of per-atom info if needed
  // enforce PBC in case atoms are outside box
  // call borders() to rebuild atom map since exchange() destroys map
  
  atom->setup();
  element->setup();
  modify->setup_pre_exchange();

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
  if (neighbor->style) neighbor->setup_bins(); 
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost,atom->x);
    domain->lamda2x(element->nlocal+element->nghost,element->x);
    domain->lamda2nodex(element->nlocal+element->nghost,element->nodex);
  }

  modify->setup_pre_neighbor();
  neighbor->build();

  if (force->pair)
    list = force->pair->list;
  else error->all(FLERR,"pair has not been defined");
 
  if (list == NULL) error->all(FLERR,"neighbor list not build yet");


  write(file);

  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::write(char *file)
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

  fp = fopen(file,"w");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open data file %s",file);
    error->one(FLERR,str);
  }

  header();
  int count = 1;

  // per atom info

  if (atom->nlocal) {
    int *tag = atom->tag;
    double **x = atom->x;
    int *numneigha2a = list->numneigha2a;
    int *numneigha2ia = list->numneigha2ia;

    if (list->ainum < atom->nlocal) error->one(FLERR,"number of neighbors for atoms not defined yet");
    for (int i = 0; i < atom->nlocal; i++) 
      fprintf(fp,"%d 1 %g %g %g %d %d %d %d -1\n",
          count++,x[i][0],x[i][1],x[i][2],tag[i],numneigha2a[i]+numneigha2ia[i],numneigha2a[i],numneigha2ia[i]);

  }

  // per element info

  if (element->nlocal) {

    int *tag = element->tag;
    double ***nodex = element->nodex;
    int *nintg = element->nintg;
    int *etype = element->etype;
    int **i2ia = element->i2ia;
    double coord[3];
    int *numneighi2a = list->numneighi2a;
    int *numneighi2ia = list->numneighi2ia;
    int *e2ilist = list->e2ilist;

    ElementVec *evec = element->evec;

    if (list->einum < element->nlocal) error->one(FLERR,"number of neighbors for integration points not defined yet");
    for (int i = 0; i < element->nlocal; i++) {
      int itype = etype[i];
      for (int j = 0; j < nintg[itype]; j++) {
        int iintpl = i2ia[itype][j];
        int iintg_local = e2ilist[i] + j;
        evec->interpolate(coord,nodex,i,iintpl,3);
        fprintf(fp,"%d 3 %g %g %g %d %d %d %d %d\n"
            ,count++,coord[0],coord[1],coord[2],tag[i],numneighi2a[iintg_local]+numneighi2ia[iintg_local],numneighi2a[iintg_local],numneighi2ia[iintg_local],j);
      }
    }
  }


  fclose(fp);
}

/* ----------------------------------------------------------------------
   proc me writes out data file header
   ------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::header()
{
  int n = atom->nlocal;
  int *nintg = element->nintg;
  int *etype = element->etype;
  for (int i = 0; i < element->nlocal; i++)
    n += nintg[etype[i]];

  fprintf(fp,"ITEM: TIMESTEP\n");

  fprintf(fp,"%d\n",me);

  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",n);
  char boundstr[9];
  domain->boundary_string(boundstr);
  if (domain->triclinic == 0) {
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
  } else {
    fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
  }
  fprintf(fp,"ITEM: ATOMS id type x y z tag numneigh numneigha numneighia iintg\n");

}


