#include <mpi.h>
#include <string.h>
#include "write_data_elem.h"
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

WriteDataElem::WriteDataElem(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data_elem command in input script
------------------------------------------------------------------------- */

void WriteDataElem::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_elem command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_data_elem command");

  // if filename contains a "*", replace with current timestep


  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  //if (narg > 1) error->all(FLERR,"Illegal write_data_elem command");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch
  /*
     pairflag = II;
     coeffflag = 1;
     int noinit = 0;
     */
  nodeflag = 1;
  ghostflag = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"node") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_data_elem command");
      if (strcmp(arg[iarg+1],"yes") == 0) nodeflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) nodeflag = 0;
      else error->all(FLERR,"Illegal write_data_elem command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ghost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_data_elem command");
      if (comm->nprocs != 1) error->all(FLERR,"Ghost option for write_data_elem doesn't allow multiple procs");
      if (strcmp(arg[iarg+1],"yes") == 0) ghostflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) ghostflag = 0;
      else error->all(FLERR,"Illegal write_data_elem command");
      iarg += 2;
    } else error->all(FLERR,"Illegal write_data_elem command");

  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {
  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_data_elem ...\n");
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
    domain->lamda2nodex(element->nlocal, element->nodex);
  }

  write(file);

  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteDataElem::write(char *file)
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
    type_arrays();
    //    if (coeffflag) force_fields();
  }

  tagint max_atom_tag = 0;
  tagint max_atom_tag_all = 0;
  tagint max_elem_tag = 0;
  tagint max_elem_tag_all = 0;

  int nnode_locals = 0;
  // per atom/element info

  if (comm->me == 0)
    fprintf(fp,"\nAtoms\n\n");

  if (natoms) {
    atoms();
    if (ghostflag) {
      max_atom_tag_all = atom->nlocal + atom->nghost;
    } else {
      for (int i = 0; i < atom->nlocal; i++) 
        max_atom_tag = MAX(max_atom_tag,atom->tag[i]);
      MPI_Allreduce(&max_atom_tag,&max_atom_tag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
    }
  }
  if (nelements) {
    elements(max_atom_tag_all);
    if (nodeflag) {
      if (ghostflag) {
        max_elem_tag_all = element->nlocal + element->nghost;
        nnode_locals = 0;
        for (int i = 0; i < element->nlocal + element->nghost; i++)
          nnode_locals += element->npe[element->etype[i]];
      } else {
        for (int i = 0; i < element->nlocal; i++) 
          max_elem_tag = MAX(max_elem_tag,element->tag[i]);
        MPI_Allreduce(&max_elem_tag,&max_elem_tag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
        element->reset_node_tags(1,0);
        nnode_locals = element->count_nodes(0);
      }
      nodes(nnode_locals,max_atom_tag_all + max_elem_tag_all);
    }
  }

  // per atom/element velocity info

  //if (comm->me == 0)
  //  fprintf(fp,"\nVelocities\n\n");
  //if (natoms) atom_velocities();
  //if (nelements) {
  //  elem_velocities(max_atom_tag_all);
  //  if (nodeflag) 
  //    node_velocities(nnode_locals,max_atom_tag_all+max_elem_tag_all);
  //}

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

void WriteDataElem::header()
{
  fprintf(fp,"CAC data file via write_data_elem and should not be used as an input data, "
      "timestep = " BIGINT_FORMAT "\n",
      update->ntimestep);

  fprintf(fp,"\n");

  bigint n;
  if (ghostflag) {
    n = atom->nlocal + atom->nghost + element->nlocal + element->nghost;
    if (nodeflag)
      for (int i = 0; i < element->nlocal + element->nghost; i++)
        n += element->npe[element->etype[i]];
  } else {
    n = atom->natoms + element->nelements;
    if (nodeflag) {
      element->count_nodes(1);
      n += element->nnodes;
    }
  }
  fprintf(fp,BIGINT_FORMAT " atoms\n",n);
  fprintf(fp,"%d atom types\n",atom->ntypes+nodeflag);

  fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",domain->boxlo[0],domain->boxhi[0]);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",domain->boxlo[1],domain->boxhi[1]);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",domain->boxlo[2],domain->boxhi[2]);

  if (domain->triclinic)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e xy xz yz\n",
        domain->xy,domain->xz,domain->yz);

}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
   ------------------------------------------------------------------------- */

void WriteDataElem::type_arrays()
{
  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp,"\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp,"%d %g\n",i,mass[i]);
    if (nodeflag) 
      fprintf(fp,"%d %g\n",atom->ntypes+1,0);
  }
}

/* ----------------------------------------------------------------------
   write out atom info in Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::atoms()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_data_atom;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_elem:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_elem:buf");

  // pack my atom data into buf

  atom->avec->pack_data_elem(buf,ghostflag);

  // write one chunk of atoms per proc to file
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
      atom->avec->write_data(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out element info in Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::elements(tagint tag_offset)
{

  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc

  double **buf;

  int ncol = element->evec->size_data_elem;

  int sendrow = element->nlocal;
  if (ghostflag) sendrow += element->nghost;
  int maxrow;

  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_elem:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_elem:buf");


  // pack my element data into buf

  element->evec->pack_element_data_elem(buf,tag_offset,ghostflag);

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

      element->evec->write_data_elem(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);

}

/* ----------------------------------------------------------------------
   write out node info in Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::nodes(int n, tagint tag_offset)
{

  // communication buffer for all my Node info
  // max_size = largest buffer needed by any proc

  double **buf;

  int ncol = element->evec->size_data_elem;

  int sendrow = n;
  int maxrow;

  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_elem:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_elem:buf");


  // pack my element data into buf

  element->evec->pack_node_data_elem(buf,tag_offset,ghostflag);

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
      element->evec->write_data_elem(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);

}

/* ----------------------------------------------------------------------
   write out atom velocities into Velocities section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::atom_velocities()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_data_vel;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my velocity data into buf

  atom->avec->pack_vel(buf);

  // write one chunk of velocities per proc to file
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

      atom->avec->write_vel(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out element velocities into Velocities section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::elem_velocities(tagint tag_offset)
{
  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc

  int ncol = element->evec->size_data_elem_vel;

  int sendrow = element->nlocal*1;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_elem:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_elem:buf");

  // pack my velocity data into buf

  element->evec->pack_element_vel_data_elem(buf,tag_offset);

  // write one chunk of velocities per proc to file
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

      element->evec->write_vel_data_elem(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out node velocities into Velocities section of data file
   ------------------------------------------------------------------------- */

void WriteDataElem::node_velocities(int n, tagint tag_offset)
{
  // communication buffer for all my Node info
  // max_size = largest buffer needed by any proc

  int ncol = element->evec->size_data_elem_vel;

  int sendrow = n;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data_elem:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data_elem:buf");

  // pack my velocity data into buf

  element->evec->pack_node_vel_data_elem(buf,tag_offset);

  // write one chunk of velocities per proc to file
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

      element->evec->write_vel_data_elem(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);
}

