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

enum {IGNORE,WARN,ERROR};                    // same as thermo.cpp
//enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteTecplot::WriteTecplot(CAC *cac) : Pointers(cac)
{
	MPI_Comm_rank(world,&me);
	MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
   called as write_tecplot command in input script
------------------------------------------------------------------------- */
=======
   called as write_dat command in input script
   ------------------------------------------------------------------------- */
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

void WriteTecplot::command(int narg, char **arg)
{
<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_tecplot command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_tecplot command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  if (narg > 1) error->all(FLERR,"Illegal write_tecplot command");

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch
  /*
     pairflag = II;
     coeffflag = 1;
     int noinit = 0;

     int iarg = 1;
     while (iarg < narg) {
     if (strcmp(arg[iarg],"pair") == 0) {
     if (iarg+2 > narg) error->all(FLERR,"Illegal write_tecplot command");
     if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
     else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
     else error->all(FLERR,"Illegal write_tecplot command");
     iarg += 2;
     } else if (strcmp(arg[iarg],"noinit") == 0) {
     noinit = 1;
     iarg++;
     } else if (strcmp(arg[iarg],"nocoeff") == 0) {
     coeffflag = 0;
     iarg++;
     } else error->all(FLERR,"Illegal write_tecplot command");
     }
     */
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_tecplot immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {
  if (comm->me == 0 && screen)
    fprintf(screen,"System init for write_tecplot ...\n");

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
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
	  domain->lamda2nodex(element->nlocal,element->nodex);
  }

  write(file);

  delete [] file;
=======
	if (domain->box_exist == 0)
		error->all(FLERR,"Write_dat command before simulation box is defined");

	if (narg < 1) error->all(FLERR,"Illegal write_dat command");

	// if filename contains a "*", replace with current timestep

	char *ptr;
	int n = strlen(arg[0]) + 16;
	char *file = new char[n];

	if ((ptr = strchr(arg[0],'*'))) {
		*ptr = '\0';
		sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
	} else strcpy(file,arg[0]);

	if (narg > 1) error->all(FLERR,"Illegal write_dat command");

	// read optional args
	// noinit is a hidden arg, only used by -r command-line switch
	/*
	   pairflag = II;
	   coeffflag = 1;
	   int noinit = 0;

	   int iarg = 1;
	   while (iarg < narg) {
	   if (strcmp(arg[iarg],"pair") == 0) {
	   if (iarg+2 > narg) error->all(FLERR,"Illegal write_dat command");
	   if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
	   else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
	   else error->all(FLERR,"Illegal write_dat command");
	   iarg += 2;
	   } else if (strcmp(arg[iarg],"noinit") == 0) {
	   noinit = 1;
	   iarg++;
	   } else if (strcmp(arg[iarg],"nocoeff") == 0) {
	   coeffflag = 0;
	   iarg++;
	   } else error->all(FLERR,"Illegal write_dat command");
	   }
	 */
	// init entire system since comm->exchange is done
	// comm::init needs neighbor::init needs pair::init needs kspace::init, etc
	// exception is when called by -r command-line switch
	//   then write_dat immediately follows reading of restart file
	//   assume that read_restart initialized necessary values
	//   if don't make exception:
	//     pair->init() can fail due to various unset values:
	//     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

	//  if (noinit == 0) {
	if (comm->me == 0 && screen)
		fprintf(screen,"System init for write_dat ...\n");

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
	comm->setup();
	comm->exchange();
	comm->borders();
	if (domain->triclinic) {
		domain->lamda2x(atom->nlocal,atom->x);
		domain->lamda2x(element->nlocal,element->x);
		domain->lamda2nodex(element->nlocal,element->nodex);
	}

	write(file);

	delete [] file;
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteTecplot::write(char *file)
{
<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write dat file

  bigint nalocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nalocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write dat file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && output->thermo->elemlostflag == ERROR)
    error->all(FLERR,"Element count is inconsistent, cannot write dat file");

  // open data file

  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open data file %s",file);
      error->one(FLERR,str);
    }

    fprintf(fp,"title=\"TIMESTEP "BIGINT_FORMAT"\"\n",update->ntimestep);
    fprintf(fp,"variables = \"x\", \"y\", \"z\"\n");
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
    
    element->reset_node_tags(1,0);
    node_connect();
  }

  /*

    // extra sections managed by fixes

  for (int i = 0; i < modify->nfix; i++)
  if (modify->fix[i]->wd_section)
  for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i,m);
  */
  // close data file

  if (me == 0) fclose(fp);
=======
	// natoms = sum of nlocal = value to write into data file
	// if unequal and thermo lostflag is "error", don't write dat file

	bigint nalocal = atom->nlocal;
	bigint natoms;
	MPI_Allreduce(&nalocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
	if (natoms != atom->natoms && output->thermo->atomlostflag == ERROR)
		error->all(FLERR,"Atom count is inconsistent, cannot write dat file");

	bigint nelocal = element->nlocal;
	bigint nelements;
	MPI_Allreduce(&nelocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
	if (nelements != element->nelements && output->thermo->elemlostflag == ERROR)
		error->all(FLERR,"Element count is inconsistent, cannot write dat file");

	// open data file

	if (me == 0) {
		fp = fopen(file,"w");
		if (fp == NULL) {
			char str[128];
			sprintf(str,"Cannot open data file %s",file);
			error->one(FLERR,str);
		}

		fprintf(fp,"title=\"TIMESTEP " BIGINT_FORMAT "\"\n",update->ntimestep);
		fprintf(fp,"variables = \"x\", \"y\", \"z\"\n");
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

		element->reset_node_tags(1);
		node_connect();
	}

	/*

	   // extra sections managed by fixes

	   for (int i = 0; i < modify->nfix; i++)
	   if (modify->fix[i]->wd_section)
	   for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i,m);
	 */
	// close data file

	if (me == 0) fclose(fp);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp
}

/* ----------------------------------------------------------------------
   write out Discrete Atoms section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::atoms()
{
	// communication buffer for all my Atom info
	// max_size = largest buffer needed by any proc

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  int ncol = atom->avec->size_tecplot_atom;
=======
	int ncol = atom->avec->size_dat_atom;
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	int sendrow = atom->nlocal;
	int maxrow;
	MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_tecplot:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_tecplot:buf");
=======
	double **buf;
	if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_dat:buf");
	else memory->create(buf,MAX(1,sendrow),ncol,"write_dat:buf");
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	// pack my atom data into buf

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  atom->avec->pack_tecplot(buf);
=======
	atom->avec->pack_dat(buf);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	// write one chunk of atoms per proc to file
	// proc 0 pings each proc, receives its chunk, writes to file
	// all other procs wait for ping, send their chunk to proc 0

	int tmp,recvrow;

	if (me == 0) {
		MPI_Status status;
		MPI_Request request;

		fprintf(fp,"zone t=\"Discrete Atoms\", f = point\n");
		for (int iproc = 0; iproc < nprocs; iproc++) {
			if (iproc) {
				MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
				MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
				MPI_Wait(&request,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
				recvrow /= ncol;
			} else recvrow = sendrow;

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
      atom->avec->write_tecplot(fp,recvrow,buf);
    }
=======
			atom->avec->write_dat(fp,recvrow,buf);
		}
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
		MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
	}

	memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out node connectivity section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::node_connect()
{

	// communication buffer for all my Element info
	// max_size = largest buffer needed by any proc

	double **buf;

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  int ncol = element->evec->size_tecplot_node_connect;
=======
	int ncol = element->evec->size_dat_node_connect;
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	int sendrow = element->nlocal;
	int maxrow;

	MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);


<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_tecplot:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_tecplot:buf");
=======
	if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_dat:buf");
	else memory->create(buf,MAX(1,sendrow),ncol,"write_dat:buf");
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	// pack my element data into buf

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  element->evec->pack_node_connect_tecplot(buf);
=======
	element->evec->pack_node_connect_dat(buf);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

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

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
      element->evec->write_node_connect_tecplot(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);
=======
			element->evec->write_node_connect_dat(fp,recvrow,buf);
		}
	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
		MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
	}
	memory->destroy(buf);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

}


/* ----------------------------------------------------------------------
   write out Nodes section of dat file
   ------------------------------------------------------------------------- */

void WriteTecplot::nodes()
{

	// communication buffer for all my Node info
	// max_size = largest buffer needed by any proc

	double **buf;

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  int ncol = element->evec->size_tecplot_node;
=======
	int ncol = element->evec->size_dat_node;
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	bigint sendrow = element->count_nodes(1);
	bigint maxrow;

	MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);


<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  if (me == 0) memory->create(buf,(1,maxrow),ncol,"write_tecplot:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_tecplot:buf");
=======
	if (me == 0) memory->create(buf,(1,maxrow),ncol,"write_dat:buf");
	else memory->create(buf,MAX(1,sendrow),ncol,"write_dat:buf");
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp


	// pack my Node data into buf

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  element->evec->pack_node_tecplot(buf);
=======
	element->evec->pack_node_dat(buf);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

	// write one chunk of nodes per proc to file
	// proc 0 pings each proc, receives its chunk, writes to file
	// all other procs wait for ping, send their chunk to proc 0

	int tmp,recvrow;

<<<<<<< HEAD:V2.4.0-a1/src/write_tecplot.cpp
  double xtmp,ytmp,ztmp;
  for (int i = 0; i < element->nsubelem[1]; i++) {
    xtmp = ytmp = ztmp = 0.0;
    for (int j = 0; j < 5; j++) {
      xtmp += element->nodex[0][j][0]*element->shape_array_center_subelem[1][i][j];
      ytmp += element->nodex[0][j][1]*element->shape_array_center_subelem[1][i][j];
      ztmp += element->nodex[0][j][2]*element->shape_array_center_subelem[1][i][j];
    }
    fprintf(screen,"%d 1 %g %g %g\n",i+1,xtmp,ytmp,ztmp);
  }
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;
    bigint nelements = element->nelements;
    bigint nnodes = element->nnodes;

    fprintf(fp,"zone t=\"Coarse Element\",");
    fprintf(fp," n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,nnodes,nelements);
    fprintf(fp," datapacking=point");
    if (domain->dimension == 3) 
      fprintf(fp,", zonetype=febrick\n");
    else 
      fprintf(fp,", zonetype=fequadrilateral\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;
      element->evec->write_node_tecplot(fp,recvrow,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }
  memory->destroy(buf);
=======

	if (me == 0) {
		MPI_Status status;
		MPI_Request request;
		bigint nelements = element->nelements;
		bigint nnodes = element->nnodes;

		fprintf(fp,"zone t=\"Coarse Element\",");
		fprintf(fp," n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,nnodes,nelements);
		fprintf(fp," datapacking=point, zonetype=febrick\n");
		for (int iproc = 0; iproc < nprocs; iproc++) {
			if (iproc) {
				MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
				MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
				MPI_Wait(&request,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
				recvrow /= ncol;
			} else recvrow = sendrow;
			element->evec->write_node_dat(fp,recvrow,buf);
		}
	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
		MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
	}
	memory->destroy(buf);
>>>>>>> 35213ed516bb7dded2274611730a1df3fc3c02a5:V2.4.0-a1/src/write_dat.cpp

}
