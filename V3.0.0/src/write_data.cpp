#include <mpi.h>
#include <string.h>
#include "write_data.h"
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

//enum{II, IJ};

/* ---------------------------------------------------------------------- */

WriteData::WriteData(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data command in input script
------------------------------------------------------------------------- */

void WriteData::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Write_data command before simulation box is defined");

  if (narg < 1) error->all(FLERR, "Illegal write_data command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0], '*'))) {
    *ptr = '\0';
    sprintf(file, "%s" BIGINT_FORMAT "%s", arg[0], update->ntimestep, ptr+1);
  } else strcpy(file, arg[0]);

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch

  //pairflag = II;
  //coeffflag = 1;
  velocityflag = elementflag = atomflag = 1;
  int noinit = 0;

  int iarg = 1;
  while (iarg < narg) {
    //if (strcmp(arg[iarg], "pair") == 0) {
    //  if (iarg+2 > narg) error->all(FLERR, "Illegal write_data command");
    //  if (strcmp(arg[iarg+1], "ii") == 0) pairflag = II;
    //  else if (strcmp(arg[iarg+1], "ij") == 0) pairflag = IJ;
    //  else error->all(FLERR, "Illegal write_data command");
    //  iarg += 2;
    if (strcmp(arg[iarg], "noinit") == 0) {
      noinit = 1;
      iarg++;
    //} else if (strcmp(arg[iarg], "nocoeff") == 0) {
    //  coeffflag = 0;
    //  iarg++;
    } else if (strcmp(arg[iarg], "novelocity") == 0) {
      velocityflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg], "noelement") == 0) {
      elementflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg], "noatom") == 0) {
      atomflag = 0;
      iarg++;
    } else error->all(FLERR, "Illegal write_data command");
  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  if (noinit == 0) {
    if (comm->me == 0) {
      if (screen) fprintf(screen, "System init for write_data ...\n");
      if (logfile) fprintf(logfile, "System init for write_data ...\n");
    }
    cac->init();

    // move atoms to new processors before writing file
    // do setup_pre_exchange to force update of per-atom info if needed
    // enforce PBC in case atoms are outside box
    // call borders() to rebuild atom map since exchange() destroys map

    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal, atom->x);
      domain->x2lamda(element->nlocal, element->x);
      domain->nodex2lamda(element->nlocal, element->nodex);
    }
    domain->pbc();
    domain->reset_box();
    comm->setup_exchange();
    comm->exchange();
    comm->setup_borders();
    comm->borders();
    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal, atom->x);
      domain->lamda2x(element->nlocal, element->x);
      domain->lamda2nodex(element->nlocal, element->nodex);
    }
  }

  write(file);

  delete [] file;

}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteData::write(char *file)
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
  MPI_Allreduce(&nalocal, &natoms, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == Thermo::ERROR) {
    if (comm->me == 0) printf("natoms = %d atom->natoms = %d\n", natoms, atom->natoms);
    error->all(FLERR, "Atom count is inconsistent, cannot write data file");
  }
  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal, &nelements, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (nelements != element->nelements && output->thermo->elemlostflag == Thermo::ERROR) {
    printf("nelement_count = %d nelements = %d \n", nelements, element->nelements);
    error->all(FLERR, "Element count is inconsistent, cannot write data file");
    }

  // sum up bond, angle counts
  // may be different than atom->nbonds, nangles if broken/turned-off
  /*
     if (atom->molecular == 1 && (atom->nbonds || atom->nbondtypes)) {
     nbonds_local = atom->avec->pack_bond(NULL);
     MPI_Allreduce(&nbonds_local, &nbonds, 1, MPI_CAC_BIGINT, MPI_SUM, world);
     }
     if (atom->molecular == 1 && (atom->nangles || atom->nangletypes)) {
     nangles_local = atom->avec->pack_angle(NULL);
     MPI_Allreduce(&nangles_local, &nangles, 1, MPI_CAC_BIGINT, MPI_SUM, world);
     }
     */
  // open data file

  if (me == 0) {
    fp = fopen(file, "w");
    if (fp == NULL) {
      char str[128];
      sprintf(str, "Cannot open data file %s", file);
      error->one(FLERR, str);
    }
  }

  // proc 0 writes header, ntype-length arrays, force fields

  if (me == 0) {
    header();
    type_arrays();
    //    if (coeffflag) force_fields();
  }

  // per atom info

  if (natoms && atomflag) {
    atoms();
    if (velocityflag) atom_velocities();
  }

  // per element info

  if (nelements && elementflag) {
    elements();
    int nnode_locals = element->count_nodes(0);

    // reset all node tag before writing node to file

    //element->reset_node_tags(1, 0);
    nodes(nnode_locals);
    if (velocityflag) node_velocities(nnode_locals);

  }

  /*

  // extra sections managed by fixes

  for (int i = 0; i < modify->nfix; i++)
  if (modify->fix[i]->wd_section)
  for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i, m);
  */
  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
   ------------------------------------------------------------------------- */

void WriteData::header()
{
  fprintf(fp, "CAC data file via write_data, "
      "timestep = " BIGINT_FORMAT "\n", 
      update->ntimestep);

  fprintf(fp, "\n");

  fprintf(fp, BIGINT_FORMAT " atoms\n", atom->natoms);
  if (element->nelements && elementflag) fprintf(fp, BIGINT_FORMAT " elements\n", element->nelements);
  fprintf(fp, "%d atom types\n", atom->ntypes);
  if (element->nelements && elementflag) fprintf(fp, "%d element types\n", element->netypes);
  fprintf(fp, "\n");
  fprintf(fp, "%-1.16e %-1.16e xlo xhi\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e ylo yhi\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e zlo zhi\n", domain->boxlo[2], domain->boxhi[2]);
  if (domain->triclinic)
    fprintf(fp, "%-1.16e %-1.16e %-1.16e xy xz yz\n", 
        domain->xy, domain->xz, domain->yz);
}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
   ------------------------------------------------------------------------- */

void WriteData::type_arrays()
{
  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp, "\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g\n", i, mass[i]);
  }

  // write out Element Types section if there exists elements in model

  if (element->nelements && elementflag) {
    fprintf(fp, "\nElement Types\n\n");
    element->evec->write_element_types(fp);
  }

}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
   ------------------------------------------------------------------------- */

void WriteData::atoms()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_data_atom;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);

  double **buf;
  if (me == 0) memory->create(buf, MAX(1, maxrow), ncol, "write_data:buf");
  else memory->create(buf, MAX(1, sendrow), ncol, "write_data:buf");

  // pack my atom data into buf

  atom->avec->pack_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp, "\nAtoms\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_data(fp, recvrow, buf);
    }

  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0], sendrow * ncol, MPI_DOUBLE, 0, 0, world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Elements section of data file
   ------------------------------------------------------------------------- */

void WriteData::elements()
{

  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc

  double **buf;

  int ncol = element->evec->size_data_element;

  int sendrow = element->nlocal;
  int maxrow;

  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);


  if (me == 0) memory->create(buf, MAX(1, maxrow), ncol, "write_data:buf");
  else memory->create(buf, MAX(1, sendrow), ncol, "write_data:buf");


  // pack my element data into buf

  element->evec->pack_element_data(buf);

  // write one chunk of elements per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;
    fprintf(fp, "\nElements\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      element->evec->write_element_data(fp, recvrow, buf);
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0], sendrow * ncol, MPI_DOUBLE, 0, 0, world);
  }
  memory->destroy(buf);

}


/* ----------------------------------------------------------------------
   write out Nodes section of data file
   ------------------------------------------------------------------------- */

void WriteData::nodes(int n)
{

  // communication buffer for all my Node info
  // max_size = largest buffer needed by any proc

  double **buf;

  int ncol = element->evec->size_data_node;

  bigint sendrow = n;
  bigint maxrow;

  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);


  if (me == 0) memory->create(buf, (1, maxrow), ncol, "write_data:buf");
  else memory->create(buf, MAX(1, sendrow), ncol, "write_data:buf");


  // pack my Node data into buf

  element->evec->pack_node_data(buf);

  // write one chunk of nodes per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;


  if (me == 0) {
    MPI_Status status;
    MPI_Request request;


    fprintf(fp, "\nNodes\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;
      element->evec->write_node_data(fp, recvrow, buf);
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0], sendrow * ncol, MPI_DOUBLE, 0, 0, world);
  }
  memory->destroy(buf);

}



/* ----------------------------------------------------------------------
   write out Atom Velocities section of data file
   ------------------------------------------------------------------------- */

void WriteData::atom_velocities()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_data_vel;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);

  double **buf;
  if (me == 0) memory->create(buf, MAX(1, maxrow), ncol, "write_data:buf");
  else memory->create(buf, MAX(1, sendrow), ncol, "write_data:buf");

  // pack my velocity data into buf

  atom->avec->pack_vel(buf);

  // write one chunk of velocities per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp, "\nVelocities\n\n");

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_vel(fp, recvrow, buf);
    }

  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0], sendrow * ncol, MPI_DOUBLE, 0, 0, world);
  }

  memory->destroy(buf);
}
/* ----------------------------------------------------------------------
   write out Node Velocities section of data file
   ------------------------------------------------------------------------- */

void WriteData::node_velocities(int n)
{
  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc

  int ncol = element->evec->size_data_vel;

  int sendrow = n;
  int maxrow;
  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);

  double **buf;
  if (me == 0) memory->create(buf, MAX(1, maxrow), ncol, "write_data:buf");
  else memory->create(buf, MAX(1, sendrow), ncol, "write_data:buf");

  // pack my velocity data into buf

  element->evec->pack_vel_data(buf);

  // write one chunk of velocities per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp, "\nNode Velocities\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      element->evec->write_vel_data(fp, recvrow, buf);
    }

  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0], sendrow * ncol, MPI_DOUBLE, 0, 0, world);
  }

  memory->destroy(buf);
}

