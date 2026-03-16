#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "dump.h"
#include "comm.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "universe.h"

using namespace CAC_NS;

Dump *Dump::dumpptr;

#define BIG 1.0e20
#define EPSILON 1.0e-6

enum{ASCEND,DESCEND};

/* ---------------------------------------------------------------------- */

Dump::Dump(CAC *cac, int narg, char **arg) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  comm_atom_forward = comm_atom_reverse = 0;
  comm_elem_forward = comm_elem_reverse = 0;

  format = NULL;
  first_flag = 0;
  flush_flag = 1;
  format_elem = NULL;
  format_atom = NULL;
  format_node = NULL;
  format_user = NULL;
  format_default = NULL;
  clearstep = 0;
  sort_flag = 0;
  append_flag = 0;
  buffer_allow = 0;
  buffer_flag = 0;
  padflag = 0;
  natom_me = nelem_me = nnode_me = nintpl_me = 0;
  natom_total = nelem_total = nnode_total = nintpl_total = 0;

  maxabuf = maxebuf = maxnbuf =  maxids = maxsort = maxproc = 0;
  abuf = ebuf = nbuf = bufsort = NULL;
  ids = idsort = NULL;
  index = proclist = NULL;

  maxsabuf = 0;
  sabuf = NULL;
  maxsebuf = 0;
  sebuf = NULL;
  maxsnbuf = 0;
  snbuf = NULL;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //  if ends in .bin = binary file
  //  else if ends in .gz = gzipped text file
  //  else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;

  multiproc = 0;
  nclusterprocs = nprocs;
  filewriter = 0;
  if (me == 0) filewriter = 1;
  fileproc = 0;
  multiname = NULL;

  // NOTE: multiproc dump not tested yet
  
  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    if (strstr(style,"mpiio"))
      error->all(FLERR,
          "Dump file MPI-IO output not allowed with \% in filename");
    multiproc = 1;
    nclusterprocs = 1;
    filewriter = 1;
    fileproc = me;
    MPI_Comm_split(world,me,0,&clustercomm);
    multiname = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(multiname,"%s%d%s",filename,me,ptr+1);
    *ptr = '%';
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;

}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;
  delete [] multiname;

  delete [] format;
  delete [] format_atom;
  delete [] format_node;
  delete [] format_elem;
  delete [] format_default;
  delete [] format_user;

  memory->destroy(abuf);
  memory->destroy(ebuf);
  memory->destroy(nbuf);
  memory->destroy(bufsort);
  memory->destroy(ids);
  memory->destroy(idsort);
  memory->destroy(index);
  memory->destroy(proclist);
  memory->destroy(sabuf);
  memory->destroy(sebuf);
  memory->destroy(snbuf);

  if (multiproc) MPI_Comm_free(&clustercomm);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (filewriter) pclose(fp);
    } else {
      if (filewriter) fclose(fp);
    }
    fp = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();

  if (!sort_flag) {
    memory->destroy(bufsort);
    memory->destroy(ids);
    memory->destroy(idsort);
    memory->destroy(index);
    memory->destroy(proclist);

    maxids = maxsort = maxproc = 0;
    bufsort = NULL;
    ids = idsort = NULL;
    index = proclist = NULL; 
  } 
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
   ------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  char *filecurrent = filename;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    sprintf(filecurrent, "%s" BIGINT_FORMAT "%s",
        filestar,update->ntimestep, ptr+1);
    *ptr = '*';

  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    fp = fopen(filecurrent,"w");

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{

  int nmax;
  int nsmin, nsmax;
  int tmp,nlines,nchars;
  bigint bnatom_me,bnelem_me,bnnode_me,bnintpl_me;
  MPI_Status status;
  MPI_Request request;
  bigint natom_header,nelem_header,nnode_header;

  
  setup_pre_write();

  // if file per timestep, open new file
  
  if (multifile) openfile();

  // natom_me = # of dump lines in discrete atom section this proc contributes to dump 
  
  natom_me = count_atoms();
  bnatom_me = natom_me;
  MPI_Allreduce(&bnatom_me,&natom_total,1,MPI_CAC_BIGINT,MPI_SUM,world);

  natom_header = natom_total;
  if (multiproc)
    MPI_Allreduce(&bnatom_me,&natom_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);

  int npe = element->npe;

  // nelem_me = # of dump lines in element connection section this proc contributes to dump
  // nnode_me = # of dump lines in node section this proc contributes to dump
  // nintpl_me = # of dump lines in addition to atom section this proc contributes to dump (only used for counting total atom for dump atom style)
  
  nelem_me = count_elements();
  nintpl_me = count_intpl_atoms();
  nnode_me = count_nodes();

  // define bigint variables to use in MPI_Allreduce
  
  bnintpl_me = nintpl_me;
  bnelem_me = nelem_me;
  bnnode_me = nnode_me;

  MPI_Allreduce(&bnintpl_me,&nintpl_total,1,MPI_CAC_BIGINT,MPI_SUM,world); 
  MPI_Allreduce(&bnelem_me,&nelem_total,1,MPI_CAC_BIGINT,MPI_SUM,world); 
  MPI_Allreduce(&bnnode_me,&nnode_total,1,MPI_CAC_BIGINT,MPI_SUM,world); 


  nelem_header = nelem_total;
  nnode_header = nnode_total;
  if (multiproc) {
    MPI_Allreduce(&bnelem_me,&nelem_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);
    MPI_Allreduce(&bnnode_me,&nnode_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);
  }

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // write common headers
  
  write_header(natom_total + nintpl_total);



  // dump atoms
  
  if (natom_total) {

    // write header specific for atom section
    
    if (filewriter) write_atom_header(natom_header);

    if (multiproc != nprocs) // not multiproc
      MPI_Allreduce(&natom_me,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = natom_me; // multiproc

    // insure buf is sized for packing and communicating
    // use nmax to insure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls

    if (nmax > maxabuf){
      if ((bigint) nmax * atom_size_one > MAXSMALLINT) {
        char errstr[100];  
        sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
        error->all(FLERR,errstr);
      }
      maxabuf = nmax;
      memory->destroy(abuf);
      memory->create(abuf,maxabuf*atom_size_one,"dump:abuf");
    }

    // pack atoms into buf

    pack_atom(NULL);

    // if buffering, convert doubles into strings
    // insure sbuf is sized for communicating
    // cannot buffer if output is to binary file

    if (buffer_flag && !binary) {
      nsatom_me = convert_atom_string(natom_me,atom_size_one,abuf);

      MPI_Allreduce(&nsatom_me,&nsmin,1,MPI_INT,MPI_MIN,world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs)
        MPI_Allreduce(&nsatom_me,&nsmax,1,MPI_INT,MPI_MAX,world);
      else nsmax = nsatom_me;

      if (nsmax>maxsabuf) {
        maxsabuf = nsmax;
        memory->grow(sabuf,maxsabuf,"dump:sabuf");
      }
    }

    // comm and output buf of doubles

    if (buffer_flag == 0 || binary) {
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(abuf,maxabuf*atom_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_DOUBLE,&nlines);
            nlines /= atom_size_one;
          } else nlines = natom_me;

          write_atom_data(nlines,abuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(abuf,natom_me*atom_size_one,MPI_DOUBLE,fileproc,0,world);
      }

      // comm and output sbuf = one big string of formatted values per proc

    } else {

      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(sabuf,maxsabuf,MPI_CHAR,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_CHAR,&nchars);
          } else nchars = nsatom_me;

          write_atom_data(nchars,(double *) sabuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(sabuf,nsatom_me,MPI_CHAR,fileproc,0,world);

      }
    }
  }

  // dump elements

  if (element->nelements) {

    // dump node section
    // should dump node before dumping element
    // write header specific for element section

    if (filewriter) write_node_header(nnode_header);
    if (multiproc != nprocs) 
      MPI_Allreduce(&nnode_me,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = nnode_me;

    // insure buf is sized for packing and communicating
    // use nmax to insure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls

    if (nmax > maxnbuf){
      if ((bigint) nmax * node_size_one > MAXSMALLINT) {
        char errstr[100];  
        sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
        error->all(FLERR,errstr);
      }

      maxnbuf = nmax;
      memory->destroy(nbuf);
      memory->create(nbuf,maxnbuf*node_size_one,"dump:nbuf");
    }

    // pack nodes into buf

    pack_node(NULL);

    // if buffering, convert doubles into strings
    // insure sbuf is sized for communicating
    // cannot buffer if output is to binary file

    if (buffer_flag && !binary) {
      nsnode_me = convert_node_string(nnode_me,node_size_one,nbuf);

      MPI_Allreduce(&nsnode_me,&nsmin,1,MPI_INT,MPI_MIN,world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs)
        MPI_Allreduce(&nsnode_me,&nsmax,1,MPI_INT,MPI_MAX,world);
      else nsmax = nsnode_me;

      if (nsmax>maxsnbuf) {
        maxsnbuf = nsmax;
        memory->grow(snbuf,maxsnbuf,"dump:snbuf");
      }
    }

    // comm and output buf of doubles

    if (buffer_flag == 0 || binary) {
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(nbuf,maxnbuf*node_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_DOUBLE,&nlines);
            nlines /= node_size_one;
          } else nlines = nnode_me;

          write_node_data(nlines,nbuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(nbuf,nnode_me*node_size_one,MPI_DOUBLE,fileproc,0,world);
      }

      // comm and output sbuf = one big string of formatted values per proc

    } else {
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(snbuf,maxsnbuf,MPI_CHAR,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_CHAR,&nchars);
          } else nchars = nsnode_me;

          write_node_data(nchars,(double *) snbuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(snbuf,nsnode_me,MPI_CHAR,fileproc,0,world);
      }
    }

    // dump element section

    if (filewriter) write_elem_header(nelem_header);
    if (multiproc != nprocs) 
      MPI_Allreduce(&nelem_me,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = nelem_me;

    // insure buf is sized for packing and communicating
    // use nmax to insure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls

    if (nmax > maxebuf) {
      if ((bigint) nmax * elem_size_one > MAXSMALLINT) {
        char errstr[100];  
        sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
        error->all(FLERR,errstr);
      }

      maxebuf = nmax;
      memory->destroy(ebuf);
      memory->create(ebuf,maxebuf*elem_size_one,"dump:ebuf");
    }

    // pack elements into buf

    pack_elem(NULL);

    // if buffering, convert doubles into strings
    // insure sbuf is sized for communicating
    // cannot buffer if output is to binary file

    if (buffer_flag && !binary) {
      nselem_me = convert_elem_string(nelem_me,elem_size_one,ebuf);
      MPI_Allreduce(&nselem_me,&nsmin,1,MPI_INT,MPI_MIN,world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs)
        MPI_Allreduce(&nselem_me,&nsmax,1,MPI_INT,MPI_MAX,world);
      else nsmax = nselem_me;

      if (nsmax>maxsebuf) {
        maxsebuf = nsmax;
        memory->grow(sebuf,maxsebuf,"dump:sebuf");
      }
    }

    // comm and output buf of doubles

    if (buffer_flag == 0 || binary) {
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(ebuf,maxebuf*elem_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_DOUBLE,&nlines);
            nlines /= elem_size_one;
          } else nlines = nelem_me;

          write_elem_data(nlines,ebuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(ebuf,nelem_me*elem_size_one,MPI_DOUBLE,fileproc,0,world);
      }

      // comm and output sbuf = one big string of formatted values per proc

    } else {
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(sebuf,maxsebuf,MPI_CHAR,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_CHAR,&nchars);
          } else nchars = nselem_me;

          write_elem_data(nchars,(double *) sebuf);
        }
        if (flush_flag && fp) fflush(fp);

      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(sebuf,nselem_me,MPI_CHAR,fileproc,0,world);
      }
    }
  }
  // if file per timestep, close file if I am filewriter

  if (multifile) {
    if (compressed) {
      if (filewriter && fp != NULL) pclose(fp);
    } else {
      if (filewriter && fp != NULL) fclose(fp);
    }
    fp = NULL;
  }

}
/* ---------------------------------------------------------------------- */

int Dump::count_atoms() 
{
  if (igroup == 0 || igroup == 2) 
    return atom->nlocal;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

int Dump::count_elements() 
{
  if (igroup == 0 || igroup == 1)
    return element->nlocal;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;

  return m;
}

/* ---------------------------------------------------------------------- */

int Dump::count_nodes()
{
  return count_elements()*element->npe;
}
/* ---------------------------------------------------------------------- */

int Dump::count_intpl_atoms() 
{
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *nintpl = element->nintpl;
  int *etype = element->etype;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m += nintpl[etype[i]];

  return m;
}

/* ----------------------------------------------------------------------
   process params common to all dumps here
   if unknown param, call modify_param specific to the dump
   ------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) append_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) append_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"buffer") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) buffer_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) buffer_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      if (buffer_flag && buffer_allow == 0)
        error->all(FLERR,"Dump_modify buffer yes not allowed for this style");
      iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
        if (strcmp(id,output->dump[idump]->id) == 0) break;
      int n;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        delete [] output->var_dump[idump];
        n = strlen(&arg[iarg+1][2]) + 1;
        output->var_dump[idump] = new char[n];
        strcpy(output->var_dump[idump],&arg[iarg+1][2]);
        n = 0;
      } else {
        n = universe->inumeric(FLERR,arg[iarg+1]);
        if (n <= 0) error->all(FLERR,"Illegal dump_modify command");
      }
      output->every_dump[idump] = n;
      iarg += 2;

    } else if (strcmp(arg[iarg],"fileper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify fileper "
            "without % in dump file name");
      int nper = universe->inumeric(FLERR,arg[iarg+1]);
      if (nper <= 0) error->all(FLERR,"Illegal dump_modify command");

      multiproc = nprocs/nper;
      if (nprocs % nper) multiproc++;
      fileproc = me/nper * nper;
      int fileprocnext = MIN(fileproc+nper,nprocs);
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;
      int icluster = fileproc/nper;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) first_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) first_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");

      if (strcmp(arg[iarg+1],"none") == 0) {
        delete [] format_line_user;
        delete [] format_int_user;
        delete [] format_bigint_user;
        delete [] format_float_user;
        format_line_user = NULL;
        format_int_user = NULL;
        format_bigint_user = NULL;
        format_float_user = NULL;
        // pass format none to child classes which may use it
        // not an error if they don't
        modify_param(narg-iarg,&arg[iarg]);
        iarg += 2;
        continue;
      }

      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");

      if (strcmp(arg[iarg+1],"line") == 0) {
        delete [] format_line_user;
        int n = strlen(arg[iarg+2]) + 1;
        format_line_user = new char[n];
        strcpy(format_line_user,arg[iarg+2]);
        iarg += 3;
      } else {   // pass other format options to child classes
        int n = modify_param(narg-iarg,&arg[iarg]);
        if (n == 0) error->all(FLERR,"Illegal dump_modify command");
        iarg += n;
      }

    } else if (strcmp(arg[iarg],"nfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (!multiproc)
        error->all(FLERR,"Cannot use dump_modify nfile "
            "without % in dump file name");
      int nfile = universe->inumeric(FLERR,arg[iarg+1]);
      if (nfile <= 0) error->all(FLERR,"Illegal dump_modify command");
      nfile = MIN(nfile,nprocs);

      multiproc = nfile;
      int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
      fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
      int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
      if (fcluster < icluster) fileproc++;
      int fileprocnext =
        static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
      fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
      if (fcluster < icluster+1) fileprocnext++;
      nclusterprocs = fileprocnext - fileproc;
      if (me == fileproc) filewriter = 1;
      else filewriter = 0;

      MPI_Comm_free(&clustercomm);
      MPI_Comm_split(world,icluster,0,&clustercomm);

      delete [] multiname;
      multiname = new char[strlen(filename) + 16];
      char *ptr = strchr(filename,'%');
      *ptr = '\0';
      sprintf(multiname,"%s%d%s",filename,icluster,ptr+1);
      *ptr = '%';
      iarg += 2;

    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      padflag = universe->inumeric(FLERR,arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"off") == 0) sort_flag = 0;
      else if (strcmp(arg[iarg+1],"id") == 0) {
        sort_flag = 1;
        sortcol = 0;
        sortorder = ASCEND;
      } else {
        sort_flag = 1;
        sortcol = universe->inumeric(FLERR,arg[iarg+1]);
        sortorder = ASCEND;
        if (sortcol == 0) error->all(FLERR,"Illegal dump_modify command");
        if (sortcol < 0) {
          sortorder = DESCEND;
          sortcol = -sortcol;
        }
        sortcolm1 = sortcol - 1;
      }
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += n;
    }
  }
}


/* ---------------------------------------------------------------------- */
bigint Dump::memory_usage()
{
  bigint bytes = memory->usage(ebuf,elem_size_one*maxebuf);
  bytes += memory->usage(sebuf,maxsebuf);
  bytes += memory->usage(abuf,atom_size_one*maxabuf);
  bytes += memory->usage(sabuf,maxsabuf);
  bytes += memory->usage(nbuf,node_size_one*maxnbuf);
  bytes += memory->usage(snbuf,maxsnbuf);
  return bytes;
}
