#include <mpi.h>
#include <string.h>
#include "write_dump_ghost_debug.h"
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

WriteDumpGhostDebug::WriteDumpGhostDebug(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_dump_ghost_debug_neighbor_debug command in input script
   all processor write out its own file contain local+ghost atoms/elements/nodes 
------------------------------------------------------------------------- */

void WriteDumpGhostDebug::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_neighbor_debug command before simulation box is defined");

  if (narg != 1) error->all(FLERR,"Illegal write_dump_ghost_debug command");

  char *ptr;
  int n = strlen(arg[0]) + 20;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s_proc_%d",arg[0],update->ntimestep,ptr+1,me);
  } else sprintf(file,"%s_proc_%d",arg[0],me);

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_dump_ghost_debug immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {



  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_dump_ghost_debug ...\n");
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

  write(file);
  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteDumpGhostDebug::write(char *file)
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



  fp = fopen(file,"w");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open data file %s",file);
    error->one(FLERR,str);
  }

  header();

  int i,j;
  double ***nodex = element->nodex;

  int nlocal = element->nlocal;
  int nghost = element->nghost;
  tagint *tag = element->tag;
  int *element_shape_ids = element->element_shape_ids;
  int *etype = element->etype;
  int *npe = element->npe;
  int n = nlocal+nghost;
  int nnode = 0;
  for (int i = 0; i < n; i++) {
    nnode += npe[etype[i]];
  }

  if (n) {
    fprintf(fp,"zone t = \"Coarse Element\"");
    fprintf(fp," n = %d e = %d",nnode,n);
    if (domain->dimension == 3) 
      fprintf(fp," datapacking = point, zonetype = febrick\n");
    else 
      fprintf(fp," datapacking = point, zonetype = fequadrilateral\n");

    for (i = 0; i < n; i++) 
      for (j = 0; j < npe[etype[i]]; j++)
        fprintf(fp,"%g %g %g %d %d %d\n",nodex[i][j][0],nodex[i][j][1],nodex[i][j][2],tag[i],etype[i],i>=nlocal);

    int count = 1;

    for (i = 0; i < n; i++) 
      if (element_shape_ids[etype[i]] == Element::QUADRILATERAL) {
        fprintf(fp,"%d %d %d %d\n",count,count+1,count+2,count+3);
        count += 4;
      } else if (element_shape_ids[etype[i]] == Element::TRIANGLE) {
        fprintf(fp,"%d %d %d %d\n",count,count+1,count+2,count+2);
        count += 3;
      } else if (element_shape_ids[etype[i]] == Element::HEXAHEDRON) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+3,count+4,count+5,count+6,count+7);
        count += 8;
      } else if (element_shape_ids[etype[i]] == Element::PYRAMID) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+3,count+4,count+4,count+4,count+4);
        count += 5;
      } else if (element_shape_ids[etype[i]] == Element::TETRAHEDRON) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+2,count+3,count+3,count+3,count+3);
        count += 4;
      } else if (element_shape_ids[etype[i]] == Element::WEDGE) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+2,count+3,count+4,count+5,count+5);
        count += 6;
      }

  }
  nlocal = atom->nlocal;
  nghost = atom->nghost;
  double **x = atom->x;
  tag = atom->tag;
  n = nlocal+nghost;
  if (n) {
    fprintf(fp, "zone t=\"Discrete Atoms\", f = point\n");
    for (i = 0; i < n; i++) 
      fprintf(fp,"%g %g %g %d %d %d\n",x[i][0],x[i][1],x[i][2],tag[i],0,i>=nlocal);
  }
  fclose(fp);
}

/* ----------------------------------------------------------------------
   proc me writes out data file header
   ------------------------------------------------------------------------- */

void WriteDumpGhostDebug::header()
{

  fprintf(fp,"title = \"CAC Simulation\"");
  fprintf(fp,", variables = \"x\", \"y\", \"z\", \"tag\", \"etype\", \"ghost\"\n");

}


