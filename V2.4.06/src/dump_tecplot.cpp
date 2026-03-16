#include "TECIO.h"
#include <string.h>
#include "dump_tecplot.h"
#include "compute.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

enum{PLT = 0, SZPLT = 1, ASCII = 2};        // output file format ASCII (.dat), binary (.plt), or SZL (.szplt)


enum{COORDINATE,VELOCITY,FORCE,STRESS,DISPLACEMENT,PE};

/*--------------------------------------------------------------------------------*/

DumpTecplot::DumpTecplot(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump tecplot command");
  nevery = universe->inumeric(FLERR,arg[3]);

  if (domain->dimension == 3) 
    npe_connect = 8;
  else 
    npe_connect = 4;

  title = NULL;
  values = NULL;
  valuelocation = NULL;

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;

  int iarg = 5;
  stress_flag = 0;
  velocity_flag = 0;
  force_flag = 0;
  displace_flag = 0;
  pe_flag = 0;
  average_flag = 0;
  fileformat = PLT;
  debug = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) stress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) stress_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"velocity") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) velocity_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) velocity_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) force_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) force_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"displace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) displace_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) displace_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pe") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) pe_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pe_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) average_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) average_flag = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      if (average_flag && !element->element_cluster_flag) {
        average_flag = 0;
        if (me == 0) error->warning(FLERR,"No elements with multiple atoms per unit cell, ave keyword is ignored");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"ascii") == 0) fileformat = ASCII;
      else if (strcmp(arg[iarg+1],"plt") == 0) fileformat = PLT;
      else if (strcmp(arg[iarg+1],"szplt") == 0) fileformat = SZPLT;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"debug") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"yes") == 0) debug = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) debug = 0;
      else error->all(FLERR,"Illegal dump tecplot command");
      iarg += 2;
 
    } else error->all(FLERR,"Illegal dump tecplot command");
  }

  int n;

  if (stress_flag) {
    n = strlen(id) + strlen("_DUMP_STRESS") + 2;
    id_stress = new char[n];
    strcpy(id_stress,id);
    strcat(id_stress,"_DUMP_STRESS");
    char **newarg = new char*[4];
    newarg[0] = id_stress;
    newarg[1] = arg[1];
    newarg[2] = (char *) "stress/atom";
    newarg[3] = (char *) "NULL";
    modify->add_compute(4,newarg);
    delete [] newarg;
  }

  if (displace_flag) {
    n = strlen(id) + strlen("_DUMP_DISPLACE") + 2;
    id_displace = new char[n];
    strcpy(id_displace,id);
    strcat(id_displace,"_DUMP_DISPLACE");
    char **newarg = new char*[3];
    newarg[0] = id_displace;
    newarg[1] = arg[1];
    newarg[2] = (char *) "displace/atom";
    modify->add_compute(3,newarg);
    delete [] newarg;
  }

  if (pe_flag) {
    n = strlen(id) + strlen("_DUMP_PE") + 2;
    id_pe = new char[n];
    strcpy(id_pe,id);
    strcat(id_pe,"_DUMP_PE");
    char **newarg = new char*[3];
    newarg[0] = id_pe;
    newarg[1] = arg[1];
    newarg[2] = (char *) "pe/atom";
    modify->add_compute(3,newarg);
    delete [] newarg;
  }

  if (average_flag) {
    if (stress_flag) comm_elem_forward += 6*element->max_npe;
    if (velocity_flag) comm_elem_forward += 3*element->max_npe;
    if (force_flag) comm_elem_forward += 3*element->max_npe;
    if (displace_flag) comm_elem_forward += 4*element->max_npe;
  }

  char *suffix = filename + strlen(filename) - strlen(".plt");
  if (suffix > filename && strcmp(suffix,".plt") == 0) fileformat = PLT;
  suffix = filename + strlen(filename) - strlen(".szplt");
  if (suffix > filename && strcmp(suffix,".szplt") == 0) fileformat = SZPLT;
  suffix = filename + strlen(filename) - strlen(".dat");
  if (suffix > filename && strcmp(suffix,".dat") == 0) fileformat = ASCII;

}

/*--------------------------------------------------------------------------------*/

DumpTecplot::~DumpTecplot()
{
  if (modify->ncompute && stress_flag) modify->delete_compute(id_stress);
  if (modify->ncompute && displace_flag) modify->delete_compute(id_displace);
  if (modify->ncompute && pe_flag) modify->delete_compute(id_pe);
  delete [] title;
  delete [] values;
  delete [] valuelocation;

  // fp is closed in dump.cpp destructor
  
  if (multifile == 0 && filewriter && fileformat != ASCII) {
    int success = TECEND142();
    if (success == -1)
      error->all(FLERR,"Cannot close dump file");
  }
}

/*--------------------------------------------------------------------------------*/

void DumpTecplot::init_style() 
{
  atom_size_one = node_size_one = 3 + 
    3*(velocity_flag + force_flag) + 6*stress_flag + 4*displace_flag + pe_flag;
  if (domain->dimension == 3) elem_size_one = 8;
  else elem_size_one = 4;

  max_size_one = MAX(atom_size_one,elem_size_one);

  int n;

  if (fileformat == ASCII) {

    // format for writing atom and node section 

    delete [] format;
    char *str;
    str = (char *) "%g %g %g";
    n = strlen(str) + 2;
    format = new char[n];
    strcpy(format,str);

    if (velocity_flag) {
      delete [] format_velocity;
      str = (char *) " %g %g %g";
      n = strlen(str) + 2;
      format_velocity = new char[n];
      strcpy(format_velocity,str);
    }

    if (force_flag) {
      delete [] format_force;
      str = (char *) " %g %g %g";
      n = strlen(str) + 2;
      format_force = new char[n];
      strcpy(format_force,str);
    }

    if (stress_flag) {
      delete [] format_stress;
      str = (char *) " %g %g %g %g %g %g";
      n = strlen(str) + 2;
      format_stress = new char[n];
      strcpy(format_stress,str);
    }

    if (displace_flag) {
      delete [] format_displace;
      str = (char *) " %g %g %g %g";
      n = strlen(str) + 2;
      format_displace = new char[n];
      strcpy(format_displace,str);
    }

    if (pe_flag) {
      delete [] format_pe;
      str = (char *) " %g";
      n = strlen(str) + 2;
      format_pe = new char[n];
      strcpy(format_pe,str);
    }

    // format for writing element section (node connectivity)

    delete [] format_elem;
    str = (char *) TAGINT_FORMAT " ";
    n = strlen(str) + 2;
    format_elem = new char[n];
    strcpy(format_elem, str);

    // setup column string 

    delete [] columns; 
    columns = new char[500];
    columns = strcpy(columns,"\"x\", \"y\", \"z\"");
    if (velocity_flag) 
      strcat(columns,", \"vx\", \"vy\", \"vz\"");
    if (force_flag) 
      strcat(columns,", \"fx\", \"fy\", \"fz\"");
    if (stress_flag) 
      strcat(columns,", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"");
    if (displace_flag) 
      strcat(columns,", \"dx\", \"dy\", \"dz\", \"d\"");
    if (pe_flag) 
      strcat(columns,", \"pe\"");

    header_choice = &DumpTecplot::header_item;
    atom_header_choice = &DumpTecplot::atom_header_item;
    elem_header_choice = &DumpTecplot::elem_header_item;
    node_header_choice = &DumpTecplot::node_header_item;

    pack_elem_choice = &DumpTecplot::pack_elem_noscale_noimage;
    pack_atom_choice = &DumpTecplot::pack_atom_noscale_noimage;
    pack_node_choice = &DumpTecplot::pack_node_noscale_noimage;

    convert_atom_choice = &DumpTecplot::convert_atom_noimage;
    convert_node_choice = &DumpTecplot::convert_node_noimage;
    convert_elem_choice = &DumpTecplot::convert_elem_noimage;

    if (buffer_flag == 1) {
      write_atom_choice = &DumpTecplot::write_string;
      write_node_choice = &DumpTecplot::write_string;
      write_elem_choice = &DumpTecplot::write_string;
    } else {
      write_atom_choice = &DumpTecplot::write_lines_noimage;
      write_elem_choice = &DumpTecplot::write_elem_lines_noimage;
      write_node_choice = &DumpTecplot::write_lines_noimage;
    }
  } else {
    pack_atom_choice = &DumpTecplot::pack_atom_noscale_noimage_binary;
    pack_node_choice = &DumpTecplot::pack_node_noscale_noimage_binary;
    pack_elem_choice = &DumpTecplot::pack_elem_noscale_noimage_binary;

    soltime = 0;
    dummy = 0;
    dummy1 = 1;

    title = new char[100];
    sprintf(title,"Output file from CAC");
    values = new char[500];
    values = strcpy(values,"x y z");
    if (velocity_flag)
      strcat(values," vx vy vz");
    if (force_flag)
      strcat(values," fx fy fz");
    if (stress_flag)
      strcat(values," sxx syy szz sxy sxz syz");
    if (displace_flag)
      strcat(values," dx dy dz d");
    if (pe_flag)
      strcat(values," pe");
    valuelocation = new int[node_size_one];
    for (int i = 0; i < node_size_one; i++) 
      valuelocation[i] = 1;
  }

  // set stress/displace compute pointer

  if (stress_flag) {
    int icompute = modify->find_compute(id_stress);
    if (icompute < 0) error->all(FLERR,"stress/atom ID for dump tecplot does not exist");
    stress = modify->compute[icompute];
  } else stress = NULL;

  if (displace_flag) {
    int icompute = modify->find_compute(id_displace);
    if (icompute < 0) error->all(FLERR,"displace/atom ID for dump tecplot does not exist");
    displace = modify->compute[icompute];
  } else displace = NULL;

  if (pe_flag) {
    int icompute = modify->find_compute(id_pe);
    if (icompute < 0) error->all(FLERR,"pe/atom ID for dump tecplot does not exist");
    pe = modify->compute[icompute];
  } else pe = NULL;


  // open single file, one time only

  if (multifile == 0) {
    openfile();
    if (fileformat == ASCII) write_header(0);
  }
}

/* ----------------------------------------------------------------------
   opening of a dump tecplot file
   ASCII or PLT or SZPLT
   ------------------------------------------------------------------------- */

void DumpTecplot::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  char *filecurrent = filename;

  // if one file per timestep, replace '*' with current timestep

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
          filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent) + 2];
        strcpy(nameslist[numfiles],filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent) + 1];
        strcpy(nameslist[fileidx],filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (fileformat == ASCII) {
      fp = fopen(filecurrent,"w");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open dump file %s",filecurrent);
        error->one(FLERR,str);
      }
    } else {
      int success = TECINI142(title,values,filecurrent,(char *) ".",&fileformat,&dummy,&debug,&dummy1);
      if (success == -1) {
        char str[128];
        sprintf(str,"Cannot open dump file %s",filecurrent);
        error->one(FLERR,str);
      }
    }
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}


/* -----------------------------------------
   overwrite write() function from dump class 
   ------------------------------------------ */

void DumpTecplot::write()
{
  // if file per timestep, open new file

  if (multifile) openfile();

  bigint bnatom_me,bnelem_me,bnnode_me;
  bigint natom_header,nelem_header,nnode_header;

  // natom_me = # of dump lines in discrete atom section this proc contributes to dump 

  natom_me = count_atoms();
  bnatom_me = natom_me;
  MPI_Allreduce(&bnatom_me,&natom_total,1,MPI_CAC_BIGINT,MPI_SUM,world);

  natom_header = natom_total;
  if (multiproc)
    MPI_Allreduce(&bnatom_me,&natom_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);

  // nelem_me = # of dump lines in element connection section this proc contributes to dump
  // nnode_me = # of dump lines in node section this proc contributes to dump

  nelem_me = count_elements();
  nnode_me = count_nodes();

  // define bigint variables to use in MPI_Allreduce

  bnelem_me = nelem_me;
  bnnode_me = nnode_me;

  MPI_Allreduce(&bnelem_me,&nelem_total,1,MPI_CAC_BIGINT,MPI_SUM,world); 
  MPI_Allreduce(&bnnode_me,&nnode_total,1,MPI_CAC_BIGINT,MPI_SUM,world); 

  nelem_header = nelem_total;
  nnode_header = nnode_total;
  if (multiproc) {
    MPI_Allreduce(&bnelem_me,&nelem_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);
    MPI_Allreduce(&bnnode_me,&nnode_header,1,MPI_CAC_BIGINT,MPI_SUM,clustercomm);
  }

  int nmax;
  MPI_Status status;
  MPI_Request request;
  int nsmin,nsmax;
  int tmp,nlines,nchars;

  if (fileformat == ASCII) {

    // if file per timestep, write common header

    if (multifile) write_header(0);

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

      if (nmax > maxbuf) {
        if ((bigint) nmax * atom_size_one > MAXSMALLINT) {
          char errstr[100];  
          sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
          error->all(FLERR,errstr);
        }
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*max_size_one,"dump:buf");
      }

      // pack atoms into buf

      pack_atom(NULL);

      // if buffering, convert doubles into strings
      // insure sbuf is sized for communicating
      // cannot buffer if output is to binary file

      if (buffer_flag && !binary) {
        nsatom_me = convert_atom_string(natom_me,atom_size_one,buf);

        MPI_Allreduce(&nsatom_me,&nsmin,1,MPI_INT,MPI_MIN,world);
        if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
        if (multiproc != nprocs)
          MPI_Allreduce(&nsatom_me,&nsmax,1,MPI_INT,MPI_MAX,world);
        else nsmax = nsatom_me;

        if (nsmax > maxsbuf) {
          maxsbuf = nsmax;
          memory->grow(sbuf,maxsbuf,"dump:sbuf");
        }
      }

      // comm and output buf of doubles

      if (buffer_flag == 0 || binary) {
        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(buf,maxbuf*max_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_DOUBLE,&nlines);
              nlines /= atom_size_one;
            } else nlines = natom_me;

            write_atom_data(nlines,buf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(buf,natom_me*atom_size_one,MPI_DOUBLE,fileproc,0,world);
        }

        // comm and output sbuf = one big string of formatted values per proc

      } else {

        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_CHAR,&nchars);
            } else nchars = nsatom_me;

            write_atom_data(nchars,(double *) sbuf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(sbuf,nsatom_me,MPI_CHAR,fileproc,0,world);

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

      if (nmax > maxbuf) {
        if ((bigint) nmax * node_size_one > MAXSMALLINT) {
          char errstr[100];  
          sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
          error->all(FLERR,errstr);
        }

        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*max_size_one,"dump:buf");
      }

      // pack nodes into buf

      pack_node(NULL);

      // if buffering, convert doubles into strings
      // insure sbuf is sized for communicating
      // cannot buffer if output is to binary file

      if (buffer_flag && !binary) {
        nsnode_me = convert_node_string(nnode_me,node_size_one,buf);

        MPI_Allreduce(&nsnode_me,&nsmin,1,MPI_INT,MPI_MIN,world);
        if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
        if (multiproc != nprocs)
          MPI_Allreduce(&nsnode_me,&nsmax,1,MPI_INT,MPI_MAX,world);
        else nsmax = nsnode_me;

        if (nsmax > maxsbuf) {
          maxsbuf = nsmax;
          memory->grow(sbuf,maxsbuf,"dump:sbuf");
        }
      }

      // comm and output buf of doubles

      if (buffer_flag == 0 || binary) {
        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(buf,maxbuf*max_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_DOUBLE,&nlines);
              nlines /= node_size_one;
            } else nlines = nnode_me;

            write_node_data(nlines,buf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(buf,nnode_me*node_size_one,MPI_DOUBLE,fileproc,0,world);
        }

        // comm and output sbuf = one big string of formatted values per proc

      } else {
        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_CHAR,&nchars);
            } else nchars = nsnode_me;

            write_node_data(nchars,(double *) sbuf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(sbuf,nsnode_me,MPI_CHAR,fileproc,0,world);
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

      if (nmax > maxbuf) {
        if ((bigint) nmax * max_size_one > MAXSMALLINT) {
          char errstr[100];  
          sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
          error->all(FLERR,errstr);
        }

        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*elem_size_one,"dump:buf");
      }

      // pack elements into buf

      pack_elem(NULL);

      // if buffering, convert doubles into strings
      // insure sbuf is sized for communicating
      // cannot buffer if output is to binary file

      if (buffer_flag && !binary) {
        nselem_me = convert_elem_string(nelem_me,elem_size_one,buf);
        MPI_Allreduce(&nselem_me,&nsmin,1,MPI_INT,MPI_MIN,world);
        if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
        if (multiproc != nprocs)
          MPI_Allreduce(&nselem_me,&nsmax,1,MPI_INT,MPI_MAX,world);
        else nsmax = nselem_me;

        if (nsmax > maxsbuf) {
          maxsbuf = nsmax;
          memory->grow(sbuf,maxsbuf,"dump:sbuf");
        }
      }

      // comm and output buf of doubles

      if (buffer_flag == 0 || binary) {
        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(buf,maxbuf*max_size_one,MPI_DOUBLE,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_DOUBLE,&nlines);
              nlines /= elem_size_one;
            } else nlines = nelem_me;

            write_elem_data(nlines,buf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(buf,nelem_me*elem_size_one,MPI_DOUBLE,fileproc,0,world);
        }

        // comm and output sbuf = one big string of formatted values per proc

      } else {
        if (filewriter) {
          for (int iproc = 0; iproc < nclusterprocs; iproc++) {
            if (iproc) {
              MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
              MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
              MPI_Wait(&request,&status);
              MPI_Get_count(&status,MPI_CHAR,&nchars);
            } else nchars = nselem_me;

            write_elem_data(nchars,(double *) sbuf);
          }
          if (flush_flag && fp) fflush(fp);

        } else {
          MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(sbuf,nselem_me,MPI_CHAR,fileproc,0,world);
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
  } else {

    int success;

    // dump atoms

    if (natom_total) {

      // write header specific for atom section

      if (filewriter) {
        zonetype = 0;       // Ordered zone
        int imax = natom_header;
        int jmax = 1;
        int kmax = 1;
        success = TECZNE142((char*)"Discrete Atoms",&zonetype,&imax,&jmax,&kmax,
            &dummy,&dummy,&dummy,
            &soltime,
            &dummy,&dummy,
            &dummy1,
            &dummy,&dummy,
            &dummy1,&dummy1,&dummy1,
            NULL,NULL,NULL,
            &dummy);
        if (success == -1) error->one(FLERR,"Cannot write Discrete Atoms Zone");
      }

      grow_buf(natom_me);

      // communicate buf and dump out each value one by one

      valueflag = COORDINATE;
      for (ival = 0; ival < 3; ival++) {
        pack_atom(NULL);
        comm_buf_tecio(natom_me);
      }
      if (velocity_flag) {
        valueflag = VELOCITY;
        for (ival = 0; ival < 3; ival++) {
          pack_atom(NULL);
          comm_buf_tecio(natom_me);
        }
      }
      if (force_flag) {
        valueflag = FORCE;
        for (ival = 0; ival < 3; ival++) {
          pack_atom(NULL);
          comm_buf_tecio(natom_me);
        }
      }
      if (stress_flag) {
        valueflag = STRESS;
        for (ival = 0; ival < 6; ival++) {
          pack_atom(NULL);
          comm_buf_tecio(natom_me);
        }
      }
      if (displace_flag) {
        valueflag = DISPLACEMENT;
        for (ival = 0; ival < 4; ival++) {
          pack_atom(NULL);
          comm_buf_tecio(natom_me);
        }
      }
      if (pe_flag) {
        valueflag = PE;
        pack_atom(NULL);
        comm_buf_tecio(natom_me);
      }
    }

    // dump elements

    if (element->nelements) {

      // dump node info

      if (filewriter) {
        if (domain->dimension == 3) zonetype = 5;    // Brick zone
        else zonetype = 3;                           // Quadrilateral zone
        int nnodes = nnode_total;
        int nelements = nelem_total;
        success = TECZNE142((char *) "Coarse Elements",&zonetype,&nnodes,&nelements,
            &dummy,&dummy,&dummy,&dummy,
            &soltime,
            &dummy,&dummy,
            &dummy1,
            &dummy,&dummy,
            0,0,0,
            NULL,valuelocation,NULL,
            &dummy);
        if (success == -1) error->one(FLERR,"Cannot write Coarse Element zone");
      }

      grow_buf(nnode_me);

      // communicate buf and dump out each value one by one

      valueflag = COORDINATE;
      for (ival = 0; ival < 3; ival++) {
        pack_node(NULL);
        comm_buf_tecio(nnode_me);
      }
      if (velocity_flag) {
        valueflag = VELOCITY;
        for (ival = 0; ival < 3; ival++) {
          pack_node(NULL);
          comm_buf_tecio(nnode_me);
        }
      }
      if (force_flag) {
        valueflag = FORCE;
        for (ival = 0; ival < 3; ival++) {
          pack_node(NULL);
          comm_buf_tecio(nnode_me);
        }
      }
      if (stress_flag) {
        valueflag = STRESS;
        for (ival = 0; ival < 6; ival++) {
          pack_node(NULL);
          comm_buf_tecio(nnode_me);
        }
      }
      if (displace_flag) {
        valueflag = DISPLACEMENT;
        for (ival = 0; ival < 4; ival++) {
          pack_node(NULL);
          comm_buf_tecio(nnode_me);
        }
      }
      if (pe_flag) {
        valueflag = PE;
        pack_node(NULL);
        comm_buf_tecio(nnode_me);
      }

      // dump node connectivity
      // use int buf for node connectivity 
      // due to tecio function requirement
      // nlines for TECNODE142 fuctions = # of values to write

      nelem_me *= elem_size_one;
      if (multiproc != nprocs) 
        MPI_Allreduce(&nelem_me,&nmax,1,MPI_INT,MPI_MAX,world);
      else nmax = nelem_me;

      if (nmax > maxibuf) {
        if ((bigint) nmax > MAXSMALLINT) {
          char errstr[100];  
          sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
          error->all(FLERR,errstr);
        }

        maxibuf = nmax;
        memory->destroy(ibuf);
        memory->create(ibuf,maxibuf,"dump:ibuf");
      }
      pack_elem(NULL);
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(ibuf,maxibuf,MPI_INT,me+iproc,0,world,&request);
            MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
            MPI_Wait(&request,&status);
            MPI_Get_count(&status,MPI_INT,&nlines);
          } else nlines = nelem_me;
          success = TECNODE142(&nlines,ibuf);
          if (success == -1) error->one(FLERR,"Cannot write to dump file");
        }
      } else {
        MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        MPI_Rsend(ibuf,nelem_me,MPI_INT,fileproc,0,world);
      }
    }

    // if file per timestep, close file if I am filewriter

    if (multifile && filewriter) {
      success = TECEND142();
      if (success == -1)
        error->all(FLERR,"Cannot close dump file");
    }

  }
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::header_item(bigint ndump)
{
  fprintf(fp,"title = \"CAC Simulation\"");
  fprintf(fp,"variables = %s\n",columns);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::node_header_item(bigint ndump)
{
  fprintf(fp,"zone t = \"Coarse Element\"");
  fprintf(fp," n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,ndump,nelem_total);
  if (domain->dimension == 3) 
    fprintf(fp," datapacking = point, zonetype = febrick\n");
  else 
    fprintf(fp," datapacking = point, zonetype = fequadrilateral\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::atom_header_item(bigint ndump)
{
  fprintf(fp,"zone t = \"Discrete Atoms\", datapacking = point\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::elem_header_item(bigint ndump)
{
  // no header for element section (node connectivity)
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_header(bigint ndump)
{
  if (multiproc) (this->*elem_header_choice)(ndump);
  else if (filewriter) (this->*elem_header_choice)(ndump);
}
/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_atom_header(bigint ndump)
{
  if (multiproc) (this->*atom_header_choice)(ndump);
  else if (filewriter) (this->*atom_header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecplot::write_node_header(bigint ndump)
{
  if (multiproc) (this->*node_header_choice)(ndump);
  else if (filewriter) (this->*node_header_choice)(ndump);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_elem(tagint *ids)
{
  (this->*pack_elem_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_node(tagint *ids)
{
  (this->*pack_node_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecplot::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_elem_noscale_noimage(tagint *ids)
{
  int *mask = element->mask;
  int nlocal = element->nlocal;
  tagint *tag = element->tag;
  tagint **nodetag = element->nodetag;

  int *npe = element->npe;
  int *etype = element->etype;
  int **element_clusters = element->element_clusters;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;
  int inpe,itype;

  element->reset_node_tags(groupbit,average_flag);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (average_flag && tag[i] != element_clusters[i][0]) 
        continue;

      itype = etype[i];
      if (element_shape_ids[itype] == Element::QUADRILATERAL) {
        for (int j = 0; j < 4; j++) 
          buf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == Element::TRIANGLE) {

        // node connectivity scheme: 1 2 3 3

        for (int j = 0; j < 3; j++) 
          buf[m++] = nodetag[i][j];
        buf[m++] = nodetag[i][2];

      } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

        for (int j = 0; j < 8; j++) 
          buf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == Element::PYRAMID) {

        // node connectivity scheme: 1 2 3 4 5 5 5 5

        for (int j = 0; j < 4; j++) 
          buf[m++] = nodetag[i][j];
        for (int j = 4; j < 8; j++)
          buf[m++] = nodetag[i][4];

      } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

        // node connectivity scheme: 1 2 3 3 4 4 4 4

        for (int j = 0; j < 3; j++) 
          buf[m++] = nodetag[i][j];
        buf[m++] = nodetag[i][2];
        for (int j = 4; j < 8; j++) 
          buf[m++] = nodetag[i][3];
      } else if (element_shape_ids[itype] == Element::WEDGE) {

        // node connectivity scheme: 1 2 3 3 4 5 6 6

        for (int j = 0; j < 3; j++) 
          buf[m++] = nodetag[i][j];
        for (int j = 3; j < 7; j++) 
          buf[m++] = nodetag[i][j-1];
        buf[m++] = nodetag[i][5];

      }
    }
} 

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_elem_noscale_noimage_binary(tagint *ids)
{
  int *mask = element->mask;
  int nlocal = element->nlocal;
  tagint *tag = element->tag;
  tagint **nodetag = element->nodetag;

  int *npe = element->npe;
  int *etype = element->etype;
  int **element_clusters = element->element_clusters;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;
  int inpe,itype;

  element->reset_node_tags(groupbit,average_flag);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (average_flag && tag[i] != element_clusters[i][0]) 
        continue;

      itype = etype[i];
      if (element_shape_ids[itype] == Element::QUADRILATERAL) {
        for (int j = 0; j < 4; j++) 
          ibuf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == Element::TRIANGLE) {

        // node connectivity scheme: 1 2 3 3

        for (int j = 0; j < 3; j++) 
          ibuf[m++] = nodetag[i][j];
        ibuf[m++] = nodetag[i][2];

      } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

        for (int j = 0; j < 8; j++) 
          ibuf[m++] = nodetag[i][j];

      } else if (element_shape_ids[itype] == Element::PYRAMID) {

        // node connectivity scheme: 1 2 3 4 5 5 5 5

        for (int j = 0; j < 4; j++) 
          ibuf[m++] = nodetag[i][j];
        for (int j = 4; j < 8; j++)
          ibuf[m++] = nodetag[i][4];

      } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

        // node connectivity scheme: 1 2 3 3 4 4 4 4

        for (int j = 0; j < 3; j++) 
          ibuf[m++] = nodetag[i][j];
        ibuf[m++] = nodetag[i][2];
        for (int j = 4; j < 8; j++) 
          ibuf[m++] = nodetag[i][3];

      } else if (element_shape_ids[itype] == Element::WEDGE) {

        // node connectivity scheme: 1 2 3 3 4 5 6 6

        for (int j = 0; j < 3; j++) 
          ibuf[m++] = nodetag[i][j];
        for (int j = 3; j < 7; j++) 
          ibuf[m++] = nodetag[i][j-1];
        ibuf[m++] = nodetag[i][5];

      }
    }
} 

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_atom_noscale_noimage(tagint *ids)
{
  double **s,**d,*e;
  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    s = stress->array_atom;
  } 
  if (displace_flag) {
    if (displace->invoked_peratom != update->ntimestep) {
      displace->compute_peratom();
      displace->addstep(update->ntimestep+nevery);
    }
    d = displace->array_atom;
  }
  if (pe_flag) {
    if (pe->invoked_peratom != update->ntimestep) {
      pe->compute_peratom();
      pe->addstep(update->ntimestep+nevery);
    }
    e = pe->vector_atom;
  }

  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = x[i][0]; 
      buf[m++] = x[i][1]; 
      buf[m++] = x[i][2];
      if (velocity_flag) {
        buf[m++] = v[i][0]; 
        buf[m++] = v[i][1]; 
        buf[m++] = v[i][2];
      }
      if (force_flag) {
        buf[m++] = f[i][0]; 
        buf[m++] = f[i][1]; 
        buf[m++] = f[i][2];
      }
      if (stress_flag) {
        buf[m++] = s[i][0]; 
        buf[m++] = s[i][1]; 
        buf[m++] = s[i][2];
        buf[m++] = s[i][3]; 
        buf[m++] = s[i][4]; 
        buf[m++] = s[i][5];
      }
      if (displace_flag) {
        buf[m++] = d[i][0]; 
        buf[m++] = d[i][1]; 
        buf[m++] = d[i][2];
        buf[m++] = d[i][3]; 
      }
      if (pe_flag) {
        buf[m++] = e[i]; 
      }

      if (ids) ids[n++] = i;
    }
  }
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_atom_noscale_noimage_binary(tagint *ids)
{
  double **s,**d,*e;
  if (valueflag == STRESS) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    s = stress->array_atom;
  } 
  if (valueflag == DISPLACEMENT) {
    if (displace->invoked_peratom != update->ntimestep) {
      displace->compute_peratom();
      displace->addstep(update->ntimestep+nevery);
    }
    d = displace->array_atom;
  }
  if (valueflag == PE) {
    if (pe->invoked_peratom != update->ntimestep) {
      pe->compute_peratom();
      pe->addstep(update->ntimestep+nevery);
    }
    e = pe->vector_atom;
  }

  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (valueflag == COORDINATE) 
        buf[m++] = x[i][ival];
      else if (valueflag == VELOCITY) 
        buf[m++] = v[i][ival];
      else if (valueflag == FORCE)
        buf[m++] = f[i][ival];
      else if (valueflag == STRESS)
        buf[m++] = s[i][ival];
      else if (valueflag == DISPLACEMENT)
        buf[m++] = d[i][ival];
      else if (valueflag == PE)
        buf[m++] = e[i];
    }
    if (ids) ids[n++] = i;
  }
}

/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_node_noscale_noimage(tagint *ids)
{
  double ***nodes,***noded,**nodee;
  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    nodes = stress->array_node;
  }
  if (displace_flag) {
    if (displace->invoked_peratom != update->ntimestep) {
      displace->compute_peratom();
      displace->addstep(update->ntimestep+nevery);
    }
    noded = displace->array_node;
  }
  if (pe_flag) {
    if (pe->invoked_peratom != update->ntimestep) {
      pe->compute_peratom();
      pe->addstep(update->ntimestep+nevery);
    }
    nodee = pe->vector_node;
  }

  int m,n,ii,i,j,k;
  int *mask = element->mask;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  m = n = 0;

  if (average_flag) {
    int **element_clusters = element->element_clusters;
    tagint *tag = element->tag;
    int *apc = element->apc;

    // forward comm of ghost stress or force/velocity 

    comm->forward_comm_dump(this);

    // build element map to find elements in clusters

    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }

    double ix,iy,iz,s0,s1,s2,s3,s4,s5,e0;
    double vx,vy,vz,fx,fy,fz,d0,d1,d2,d3;
    for (i = 0; i < nlocal; i++) {
      if (tag[i] != element_clusters[i][0]) continue;
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (j = 0; j < inpe; j++) {
        ix = nodex[i][j][0];
        iy = nodex[i][j][1];
        iz = nodex[i][j][2];
        if (velocity_flag) {
          vx = nodev[i][j][0]; 
          vy = nodev[i][j][1]; 
          vz = nodev[i][j][2];
        }
        if (force_flag) {
          fx = nodef[i][j][0]; 
          fy = nodef[i][j][1]; 
          fz = nodef[i][j][2];
        }
        if (stress_flag) {
          s0 = nodes[i][j][0]; 
          s1 = nodes[i][j][1];
          s2 = nodes[i][j][2]; 
          s3 = nodes[i][j][3];
          s4 = nodes[i][j][4]; 
          s5 = nodes[i][j][5];
        }
        if (displace_flag) {
          d0 = noded[i][j][0]; 
          d1 = noded[i][j][1];
          d2 = noded[i][j][2]; 
          d3 = noded[i][j][3];
        } 
        if (pe_flag) {
          e0 = nodee[i][j]; 
        } 

        for (k = 1; k < iapc; k++) {
          ii = element->map(element_clusters[i][k]);
          ix += nodex[ii][j][0]; 
          iy += nodex[ii][j][1]; 
          iz += nodex[ii][j][2];
          if (velocity_flag) {
            vx += nodev[ii][j][0]; 
            vy += nodev[ii][j][1]; 
            vz += nodev[ii][j][2];
          }
          if (force_flag) {
            fx += nodef[ii][j][0]; 
            fy += nodef[ii][j][1]; 
            fz += nodef[ii][j][2];
          }
          if (stress_flag) {
            s0 += nodes[ii][j][0]; 
            s1 += nodes[ii][j][1];
            s2 += nodes[ii][j][2]; 
            s3 += nodes[ii][j][3];
            s4 += nodes[ii][j][4]; 
            s5 += nodes[ii][j][5];
          }
          if (displace_flag) {
            d0 += noded[ii][j][0]; 
            d1 += noded[ii][j][1];
            d2 += noded[ii][j][2]; 
            d3 += noded[ii][j][3];
          } 
          if (pe_flag) {
            e0 += nodee[ii][j]; 
          } 

        }
        buf[m++] = ix/iapc; 
        buf[m++] = iy/iapc; 
        buf[m++] = iz/iapc;
        if (velocity_flag) {
          buf[m++] = vx/iapc; 
          buf[m++] = vy/iapc; 
          buf[m++] = vz/iapc;
        }
        if (force_flag) {
          buf[m++] = fx/iapc; 
          buf[m++] = fy/iapc; 
          buf[m++] = fz/iapc;
        }
        if (stress_flag) {
          buf[m++] = s0/iapc; 
          buf[m++] = s1/iapc; 
          buf[m++] = s2/iapc; 
          buf[m++] = s3/iapc; 
          buf[m++] = s4/iapc; 
          buf[m++] = s5/iapc;
        }
        if (displace_flag) {
          buf[m++] = d0/iapc; 
          buf[m++] = d1/iapc; 
          buf[m++] = d2/iapc; 
          buf[m++] = d3/iapc; 
        }
        if (pe_flag) {
          buf[m++] = e0/iapc; 
        }
      }
    }

    // clear map if needed

    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }

  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (j = 0; j < npe[etype[i]]; j++) {
          buf[m++] = nodex[i][j][0]; 
          buf[m++] = nodex[i][j][1]; 
          buf[m++] = nodex[i][j][2];
          if (velocity_flag) {
            buf[m++] = nodev[i][j][0]; 
            buf[m++] = nodev[i][j][1]; 
            buf[m++] = nodev[i][j][2]; 
          }
          if (force_flag) {
            buf[m++] = nodef[i][j][0]; 
            buf[m++] = nodef[i][j][1]; 
            buf[m++] = nodef[i][j][2];
          }
          if (stress_flag) {
            buf[m++] = nodes[i][j][0]; 
            buf[m++] = nodes[i][j][1]; 
            buf[m++] = nodes[i][j][2]; 
            buf[m++] = nodes[i][j][3]; 
            buf[m++] = nodes[i][j][4]; 
            buf[m++] = nodes[i][j][5];
          }
          if (displace_flag) {
            buf[m++] = noded[i][j][0]; 
            buf[m++] = noded[i][j][1]; 
            buf[m++] = noded[i][j][2]; 
            buf[m++] = noded[i][j][3]; 
          }
          if (pe_flag) {
            buf[m++] = nodee[i][j]; 
          }
        }
        if (ids) ids[n++] = i;
      }
    }
  }
}


/*------------------------------------------------------------------------------------------*/

void DumpTecplot::pack_node_noscale_noimage_binary(tagint *ids)
{
  double ***nodes,***noded,**nodee;
  if (valueflag == STRESS) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    nodes = stress->array_node;
  }
  if (valueflag == DISPLACEMENT) {
    if (displace->invoked_peratom != update->ntimestep) {
      displace->compute_peratom();
      displace->addstep(update->ntimestep+nevery);
    }
    noded = displace->array_node;
  }
  if (valueflag == PE) {
    if (pe->invoked_peratom != update->ntimestep) {
      pe->compute_peratom();
      pe->addstep(update->ntimestep+nevery);
    }
    nodee = pe->vector_node;
  }

  int m,n,ii,i,j,k;
  int *mask = element->mask;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  m = n = 0;

  if (average_flag) {
    int **element_clusters = element->element_clusters;
    tagint *tag = element->tag;
    int *apc = element->apc;

    // forward comm of ghost stress or force/velocity 

    comm->forward_comm_dump(this);

    // build element map to find elements in clusters

    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }

    double ix,iy,iz,s0,s1,s2,s3,s4,s5;
    double vx,vy,vz,fx,fy,fz,d0,d1,d2,d3;
    for (i = 0; i < nlocal; i++) {
      if (tag[i] != element_clusters[i][0]) continue;
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (j = 0; j < inpe; j++) {
        if (valueflag == COORDINATE) 
          buf[m] = nodex[i][j][ival];
        else if (valueflag == VELOCITY)
          buf[m] = nodev[i][j][ival];
        else if (valueflag == FORCE)
          buf[m] = nodef[i][j][ival];
        else if (valueflag == STRESS)
          buf[m] = nodes[i][j][ival];
        else if (valueflag == DISPLACEMENT)
          buf[m] = noded[i][j][ival];
        else if (valueflag == PE)
          buf[m] = nodee[i][j];
        for (k = 1; k < iapc; k++) {
          ii = element->map(element_clusters[i][k]);
          if (valueflag == COORDINATE) 
            buf[m] += nodex[ii][j][ival];
          else if (valueflag == VELOCITY)
            buf[m] += nodev[ii][j][ival];
          else if (valueflag == FORCE)
            buf[m] += nodef[ii][j][ival];
          else if (valueflag == STRESS)
            buf[m] += nodes[ii][j][ival];
          else if (valueflag == DISPLACEMENT)
            buf[m] += noded[ii][j][ival];
          else if (valueflag == PE)
            buf[m] += nodee[ii][j];
        }
        buf[m] /= iapc; 
        m++;
      }
    }

    // clear map if needed

    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }

  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (j = 0; j < npe[etype[i]]; j++) {
          if (valueflag == COORDINATE) 
            buf[m++] = nodex[i][j][ival];
          else if (valueflag == VELOCITY)
            buf[m++] = nodev[i][j][ival];
          else if (valueflag == FORCE)
            buf[m++] = nodef[i][j][ival];
          else if (valueflag == STRESS)
            buf[m++] = nodes[i][j][ival];
          else if (valueflag == DISPLACEMENT)
            buf[m++] = noded[i][j][ival];
          else if (valueflag == PE)
            buf[m++] = nodee[i][j];
        }
      }
      if (ids) ids[n++] = i;
    }
  }
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_node_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_node_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_elem_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_elem_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_node_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
    offset += sprintf(&sbuf[offset],format,
        mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (velocity_flag) {
      offset += sprintf(&sbuf[offset],format_velocity,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
    if (force_flag) {
      offset += sprintf(&sbuf[offset],format_force,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
    if (stress_flag) {
      offset += sprintf(&sbuf[offset],format_stress,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      offset += sprintf(&sbuf[offset],format_displace,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    if (pe_flag) {
      offset += sprintf(&sbuf[offset],format_pe,mybuf[m]);
      m += 1;
    }

    offset += sprintf(&sbuf[offset],"\n");
  }
  return offset;
}
/*--------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_atom_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
    offset += sprintf(&sbuf[offset],format,
        mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (velocity_flag) {
      offset += sprintf(&sbuf[offset],format_velocity,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
    if (force_flag) {
      offset += sprintf(&sbuf[offset],format_force,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
    if (stress_flag) {
      offset += sprintf(&sbuf[offset],format_stress,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      offset += sprintf(&sbuf[offset],format_displace,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    if (pe_flag) {
      offset += sprintf(&sbuf[offset],format_pe,mybuf[m]);
      m += 1;
    }
    offset += sprintf(&sbuf[offset],"\n");
  }
  return offset;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::convert_elem_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
    for (int j = 0; j < npe_connect; j++)
      offset += sprintf(&sbuf[offset],format_elem,
          static_cast<tagint>(mybuf[m++]));

    offset += sprintf(&sbuf[offset],"\n");
  }
  return offset;
}	

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_node_data(int n, double *mybuf)
{
  (this->*write_node_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_data(int n, double *mybuf)
{
  (this->*write_elem_choice) (n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (velocity_flag) {
      fprintf(fp,format_velocity,mybuf[m],mybuf[m+1],
          mybuf[m+2]);
      m += 3;
    }
    if (force_flag) {
      fprintf(fp,format_force,mybuf[m],mybuf[m+1],
          mybuf[m+2]);
      m += 3;
    }
    if (stress_flag) {
      fprintf(fp,format_stress,mybuf[m],mybuf[m+1],
          mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      fprintf(fp,format_displace,mybuf[m],mybuf[m+1],
          mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    if (pe_flag) {
      fprintf(fp,format_pe,mybuf[m]);
      m += 1;
    }

    fprintf(fp,"\n");
  }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::write_elem_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  int k;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < npe_connect; j++)
      fprintf(fp,format_elem, static_cast<tagint>(mybuf[m++]));
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpTecplot::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int ii,i,j,m;
  int npe = element->max_npe;
  m = 0;

  if (velocity_flag) {
    double ***nodev = element->nodev;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = nodev[i][j][0];
        buf[m++] = nodev[i][j][1];
        buf[m++] = nodev[i][j][2];
      }
    }
  }
  if (force_flag) {
    double ***nodef = element->nodef;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = nodef[i][j][0];
        buf[m++] = nodef[i][j][1];
        buf[m++] = nodef[i][j][2];
      }
    }
  }
  if (stress_flag) {
    double ***nodes = stress->array_node;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = nodes[i][j][0];
        buf[m++] = nodes[i][j][1];
        buf[m++] = nodes[i][j][2];
        buf[m++] = nodes[i][j][3];
        buf[m++] = nodes[i][j][4];
        buf[m++] = nodes[i][j][5];
      }
    }
  } 
  if (displace_flag) {
    double ***noded = displace->array_node;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = noded[i][j][0];
        buf[m++] = noded[i][j][1];
        buf[m++] = noded[i][j][2];
        buf[m++] = noded[i][j][3];
      }
    }
  }
  if (pe_flag) {
    double **nodee = pe->vector_node;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = nodee[i][j];
      }
    }
  }


  return m;
}

/* ---------------------------------------------------------------------- */

void DumpTecplot::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i,j,m,last;
  last = first + n;
  int npe = element->max_npe;
  m = 0;
  if (velocity_flag) {
    double ***nodev = element->nodev;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        nodev[i][j][0] = buf[m++];
        nodev[i][j][1] = buf[m++];
        nodev[i][j][2] = buf[m++];
      }
  }
  if (force_flag) {
    double ***nodef = element->nodef;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        nodef[i][j][0] = buf[m++];
        nodef[i][j][1] = buf[m++];
        nodef[i][j][2] = buf[m++];
      }
  }
  if (stress_flag) {
    double ***nodes = stress->array_node;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        nodes[i][j][0] = buf[m++];
        nodes[i][j][1] = buf[m++];
        nodes[i][j][2] = buf[m++];
        nodes[i][j][3] = buf[m++];
        nodes[i][j][4] = buf[m++];
        nodes[i][j][5] = buf[m++];
      }
  } 
  if (displace_flag) {
    double ***noded = displace->array_node;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        noded[i][j][0] = buf[m++];
        noded[i][j][1] = buf[m++];
        noded[i][j][2] = buf[m++];
        noded[i][j][3] = buf[m++];
      }
  } 
  if (pe_flag) {
    double **nodee = pe->vector_node;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        nodee[i][j] = buf[m++];
      }
  } 

}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::count_elements()
{
  if (average_flag) return element->count_element_clusters();
  else if (igroup == 0 || igroup == 1) return element->nlocal;
  int m = 0;
  int *mask = element->mask;
  for (int i = 0; i < element->nlocal; i++) 
    if (mask[i] & groupbit) m++;
  return m;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecplot::count_nodes()
{
  int n = 0;
  int *etype = element->etype;
  int *npe = element->npe;
  for (int i = 0; i < element->nlocal; i++)
    n += npe[etype[i]];
  return n;
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::comm_buf_tecio(int num_me)
{
  MPI_Status status;
  MPI_Request request;
  int nlines,tmp,success;
  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
      } else nlines = num_me;
      success = TECDAT142(&nlines,buf,&dummy1);
      if (success == -1) error->one(FLERR,"Cannot write to dump file");
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,num_me,MPI_DOUBLE,fileproc,0,world);
  }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecplot::grow_buf(int num_me)
{
  int nmax;
  if (multiproc != nprocs) 
    MPI_Allreduce(&num_me,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = num_me;

  if (nmax > maxbuf) {
    if ((bigint) nmax > MAXSMALLINT) {
      char errstr[100];  
      sprintf(errstr,"Too much per-proc info for dump ID %s, nmax = %d",id,nmax);
      error->all(FLERR,errstr);
    }

    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf,"dump:buf");
  }
}
