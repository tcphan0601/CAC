#include "TECIO.h"
#include <mpi.h>
#include <string.h>
#include "write_tecplot.h"
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

enum{PLT = 0, SZPLT = 1, ASCII = 2};        // output file format ASCII (.dat), binary (.plt), or SZL (.szplt)
//enum{II, IJ};

/* ---------------------------------------------------------------------- */

WriteTecplot::WriteTecplot(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
}

/* ----------------------------------------------------------------------
   called as write_tecplot command in input script
------------------------------------------------------------------------- */

void WriteTecplot::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Write_tecplot command before simulation box is defined");

  if (narg < 1) error->all(FLERR, "Illegal write_tecplot command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0], '*'))) {
    *ptr = '\0';
    sprintf(file, "%s" BIGINT_FORMAT "%s", arg[0], update->ntimestep, ptr+1);
  } else strcpy(file, arg[0]);

//  if (narg > 1) error->all(FLERR, "Illegal write_tecplot command");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch

  //pairflag = II;
  //coeffflag = 1;
  int noinit = 0;
  
  fileformat = PLT;
  aveflag = 1;
  debug = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "format") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal write_tecplot command");
      if (strcmp(arg[iarg+1], "ascii") == 0) fileformat = ASCII;
      else if (strcmp(arg[iarg+1], "plt") == 0) fileformat = PLT;
      else if (strcmp(arg[iarg+1], "szplt") == 0) fileformat = SZPLT;
      else error->all(FLERR, "Illegal write_tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "debug") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal write_tecplot command");
      if (strcmp(arg[iarg+1], "yes") == 0) debug = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) debug = 0;
      else error->all(FLERR, "Illegal write_tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal write_tecplot command");
      if (strcmp(arg[iarg+1], "yes") == 0) aveflag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) aveflag = 0;
      else error->all(FLERR, "Illegal write_tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "noinit") == 0) {
      noinit = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "nocoeff") == 0) {
      //  coeffflag = 0;
      //  iarg++;
    } else error->all(FLERR, "Illegal write_tecplot command");
  }


  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_tecplot immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  if (noinit == 0) {
    if (comm->me == 0) {
      if (screen) fprintf(screen, "System init for write_tecplot ...\n");
      if (logfile) fprintf(logfile, "System init for write_tecplot ...\n");
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

void WriteTecplot::write(char *file)
{
  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write dat file

  bigint nalocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nalocal, &natoms, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == Thermo::ERROR)
    error->all(FLERR, "Atom count is inconsistent, cannot write dat file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal, &nelements, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (nelements != element->nelements && output->thermo->elemlostflag == Thermo::ERROR)
    error->all(FLERR, "Element count is inconsistent, cannot write dat file");

  // open data file and write headers

  if (me == 0) {
    if (fileformat == ASCII) {
      fp = fopen(file, "w");
      if (fp == NULL) {
        char str[128];
        sprintf(str, "Cannot open data file %s", file);
        error->one(FLERR, str);
      }

      fprintf(fp, "title=\"TIMESTEP " BIGINT_FORMAT "\"\n", update->ntimestep);
      fprintf(fp, "variables = \"x\", \"y\", \"z\", \"ID\", \"etype\", \"ctype\"\n");
    } else {
      soltime = update->ntimestep;
      dummy = 0;
      dummy1 = 1;
      int success = TECINI142((char *) "OUTPUT FILE FROM WRITE_TECPLOT COMMAND IN CAC", 
          (char *) "x y z ID etype ctype", 
          file, 
          (char *) ".", 
          &fileformat, 
          &dummy, 
          &debug, 
          &dummy1);
      if (success == -1) {
        char str[128];
        sprintf(str, "Cannot open data file %s", file);
        error->one(FLERR, str);
      }
    }
  }

  // proc 0 writes header, ntype-length arrays, force fields

  // per atom info

  if (natoms) {
    atoms();
  }

  // per element info

  if (nelements) {
    nodes();

    // reset all node tag before writing node connectivity to file

    if (aveflag) {
      memory->create(nodecell_ids, element->nlocal, element->maxnpe, "write_tecplot:nodecell_ids");
      element->set_node_connectivities(1, nodecell_ids);
    } else {
      memory->create(node_ids, element->nlocal, element->maxapc, element->maxnpe, "write_tecplot:nodecell_ids");
      element->set_node_connectivities(1, node_ids);
    }
    node_connect();
    if (aveflag) memory->destroy(nodecell_ids);
    else memory->destroy(node_ids);

  }

  /*

  // extra sections managed by fixes

  for (int i = 0; i < modify->nfix; i++)
  if (modify->fix[i]->wd_section)
  for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i, m);
  */

  // close data file

  if (me == 0) {
    if (fileformat == ASCII)
      fclose(fp);
    else {
      int success = TECEND142();
      if (success == -1) 
        error->one(FLERR, "Cannot close data file");
    }
  }
}

/* ----------------------------------------------------------------------
   write out Discrete Atoms section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::atoms()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_tecplot_atom;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);


  int tmp, recvrow;

  if (fileformat == ASCII) {
    double **buf;
    if (me == 0) memory->create(buf, MAX(1, maxrow), ncol, "write_tecplot:buf");
    else memory->create(buf, MAX(1, sendrow), ncol, "write_tecplot:buf");

    // pack my atom data into buf

    atom->avec->pack_tecplot(buf);

    // write one chunk of atoms per proc to file
    // proc 0 pings each proc, receives its chunk, writes to file
    // all other procs wait for ping, send their chunk to proc 0


    if (me == 0) {
      MPI_Status status;
      MPI_Request request;

      fprintf(fp, "zone t=\"Discrete Atoms\", f = point\n");
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(&buf[0][0], maxrow*ncol, MPI_DOUBLE, iproc, 0, world, &request);
          MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
          MPI_Wait(&request, &status);
          MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
          recvrow /= ncol;
        } else recvrow = sendrow;

        atom->avec->write_tecplot(fp, recvrow, buf);
      }

    } else {
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(&buf[0][0], sendrow*ncol, MPI_DOUBLE, 0, 0, world);
    }


    memory->destroy(buf);
  } else {

    double *buf;
    if (me == 0) memory->create(buf, MAX(1, maxrow), "write_tecplot:bbuf");
    else memory->create(buf, MAX(1, sendrow), "write_tecplot:bbuf");
    if (me == 0) {
      zonetype = 0;       // Ordered zone
      int imax = atom->natoms;
      int jmax = 1;
      int kmax = 1;

      // write zone header information

      int success = TECZNE142((char*)"Discrete Atoms", 
          &zonetype, 
          &imax, 
          &jmax, 
          &kmax, 
          &dummy, 
          &dummy, 
          &dummy, 
          &soltime, 
          &dummy, 
          &dummy, 
          &dummy1, 
          &dummy, 
          &dummy, 
          &dummy1, 
          &dummy1, 
          &dummy1, 
          NULL, 
          NULL, 
          NULL, 
          &dummy);

      if (success == -1) error->one(FLERR, "Cannot write Discrete Atoms Zone");

      MPI_Status status;
      MPI_Request request;

      for (int ival = 0; ival < ncol; ival++) {

        atom->avec->pack_tecplot_binary(buf, ival);
        for (int iproc = 0; iproc < nprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(&buf[0], maxrow, MPI_DOUBLE, iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
          } else recvrow = sendrow;

          success = TECDAT142(&recvrow, buf, &dummy1);
          if (success == -1) error->one(FLERR, "Cannot write Discrete Atoms Zone");
        }
      }
    } else {
      for (int ival = 0; ival < ncol; ival++) {
        atom->avec->pack_tecplot_binary(buf, ival);
        MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(&buf[0], sendrow, MPI_DOUBLE, 0, 0, world);
      }
    }

    memory->destroy(buf);
  }

}

/* ----------------------------------------------------------------------
   write out node connectivity section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::node_connect()
{

  // communication buffer for all my Element info
  // max_size = largest buffer needed by any proc


  int ncol = element->evec->size_tecplot_node_connect;

  int sendrow;
  if (aveflag) 
    sendrow = element->nlocal;
  else {
    sendrow = element->count_polygons(0);
  }

  int maxrow;

  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);

  int *buf;
  if (me == 0) memory->create(buf, MAX(1, maxrow) * ncol, "write_tecplot:buf");
  else memory->create(buf, MAX(1, sendrow)*ncol, "write_tecplot:buf");

  // pack my element data into buf

  if (aveflag) 
    element->evec->pack_node_connect_tecplot(buf, nodecell_ids);
  else 
    element->evec->pack_node_connect_tecplot(buf, node_ids);

  // write one chunk of elements per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp, recvrow;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0], maxrow * ncol, MPI_INT, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_INT, &recvrow);
      } else recvrow = sendrow*ncol;

      if (fileformat == ASCII) {
        element->evec->write_node_connect_tecplot(fp, recvrow / ncol, buf);
      } else {
        int success = tecnode142(&recvrow, buf);
        if (success == -1) error->one(FLERR, "Cannot write node connectivity section");
      }
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0], sendrow * ncol, MPI_INT, 0, 0, world);
  }

  memory->destroy(buf);
}


/* ----------------------------------------------------------------------
   write out Nodes section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::nodes()
{

  // communication buffer for all my Node info
  // max_size = largest buffer needed by any proc


  int ncol = element->evec->size_tecplot_node;

  bigint sendrow;
  if (aveflag) sendrow = element->count_node_cells(0);
  else sendrow = element->count_nodes(0); 


  bigint maxrow;

  MPI_Allreduce(&sendrow, &maxrow, 1, MPI_INT, MPI_MAX, world);


  int tmp, recvrow;

  if (fileformat == ASCII) {

    double **buf;
    if (me == 0) memory->create(buf, (1, maxrow), ncol, "write_tecplot:buf");
    else memory->create(buf, MAX(1, sendrow), ncol, "write_tecplot:buf");

    // pack my Node data into buf

    element->evec->pack_node_tecplot(buf, aveflag);

    // write one chunk of nodes per proc to file
    // proc 0 pings each proc, receives its chunk, writes to file
    // all other procs wait for ping, send their chunk to proc 0
    bigint nelements;
    bigint nnodes;

    if (aveflag) {
      nelements = element->nelements;
      nnodes = element->nncells;
    } else {
      nelements = element->count_polygons(1);
      nnodes = element->nnodes;
    }

    if (me == 0) {
      MPI_Status status;
      MPI_Request request;
      fprintf(fp, "zone t=\"Coarse Element\", ");
      fprintf(fp, " n = " BIGINT_FORMAT " e = " BIGINT_FORMAT, nnodes, nelements);
      fprintf(fp, " datapacking=point");
      if (domain->dimension == 3) 
        fprintf(fp, ", zonetype=febrick\n");
      else 
        fprintf(fp, ", zonetype=fequadrilateral\n");

      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(&buf[0][0], maxrow * ncol, MPI_DOUBLE, iproc, 0, world, &request);
          MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
          MPI_Wait(&request, &status);
          MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
          recvrow /= ncol;
        } else recvrow = sendrow;
        element->evec->write_node_tecplot(fp, recvrow, buf);
      }
    } else {
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(&buf[0][0], sendrow*ncol, MPI_DOUBLE, 0, 0, world);
    }
    memory->destroy(buf);

  } else {

    double *buf;
    int nelements;
    int nnodes;
    if (aveflag) {
      nelements = element->nelements;
      element->count_node_cells(1);
      nnodes = element->nncells;
    } else {
      nelements = element->count_polygons(1);
      element->count_nodes(1);
      nnodes = element->nnodes;
    }

    if (me == 0) memory->create(buf, (1, maxrow), "write_tecplot:buf");
    else memory->create(buf, MAX(1, sendrow), "write_tecplot:buf");

    if (me == 0) {
      MPI_Status status;
      MPI_Request request;
      if (domain->dimension == 3) zonetype = 5;    // Brick zone
      else zonetype = 3;                // Quadrilateral zone

      int valuelocation[] = {1, 1, 1, 1, 1, 1};

      // write header section for Coarse Elements Zone

      int success = TECZNE142((char *) "Coarse Elements", 
          &zonetype, 
          &nnodes, 
          &nelements, 
          &dummy, 
          &dummy, 
          &dummy, 
          &dummy, 
          &soltime, 
          &dummy, 
          &dummy, 
          &dummy1, 
          &dummy, 
          &dummy, 
          0,              /* totalnumfacenodes */
          0,              /* numconnectedboundaryfaces */
          0,              /* totalnumboundaryconnections */
          NULL,           /* passivevarlist */
          valuelocation,  /* valuelocation = nodal */
          NULL,           /* sharvarfromzone */
          &dummy);

      if (success == -1) error->one(FLERR, "Cannot write Coarse Element zone");
      for (int ival = 0; ival < ncol; ival++) {

        // pack my Node data into buf

        element->evec->pack_node_tecplot_binary(buf, ival, aveflag);

        for (int iproc = 0; iproc < nprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(&buf[0], maxrow, MPI_DOUBLE, iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &recvrow);
          } else recvrow = sendrow;
          success = TECDAT142(&recvrow, buf, &dummy1);
          if (success == -1) error->one(FLERR, "Cannot write Coarse Elements Zones");
        }
      }
    } else {
      for (int ival = 0; ival < ncol; ival++) {

        // pack my Node data into buf

        element->evec->pack_node_tecplot_binary(buf, ival, aveflag);

        MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(&buf[0], sendrow, MPI_DOUBLE, 0, 0, world);
      }
    }

    memory->destroy(buf);
  }
}



