#include <string.h>
#include "dump_atom.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "compute.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

enum{INTERPOLATE,INTEGRATION,NODE,CENTER};

/*--------------------------------------------------------------------------------*/

DumpAtom::DumpAtom(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump atom command");
  nevery = universe->inumeric(FLERR,arg[3]);

  int iarg = 5;
  stress_flag = 0;
  force_flag = 0;
  velocity_flag = 0;
  displace_flag = 0;
  wrap_flag = 0;
  dump_element_style = INTERPOLATE;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) stress_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) stress_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"velocity") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) velocity_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) velocity_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) force_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) force_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"displace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) displace_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) displace_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pe") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) pe_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pe_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"style") == 0) {
      if (strcmp(arg[iarg+1],"interpolate") == 0)
        dump_element_style = INTERPOLATE;
      else if (strcmp(arg[iarg+1],"integration") == 0)
        dump_element_style = INTEGRATION;
      else if (strcmp(arg[iarg+1],"node") == 0)
        dump_element_style = NODE;
      else if (strcmp(arg[iarg+1],"center") == 0)
        dump_element_style = CENTER;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"wrap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) wrap_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wrap_flag = 0;
      else error->all(FLERR,"Illegal dump atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal dump atom command");
  }

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;

  int n;

  if (stress_flag) {
    n = strlen(id) + strlen("_DUMP_STRESS") + 2;
    id_stress = new char[n];
    strcpy(id_stress,id);
    strcat(id_stress,"_DUMP_STRESS");
    char **newarg = new char*[4];
    newarg[0] = id_stress;
    newarg[1] = arg[1] ;
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

}

/*--------------------------------------------------------------------------------*/

DumpAtom::~DumpAtom()
{
  // check ncompute in case all fixes have already been deleted

  if (modify->ncompute && stress_flag) modify->delete_compute(id_stress);
  if (modify->ncompute && displace_flag) modify->delete_compute(id_displace);
  if (modify->ncompute && pe_flag) modify->delete_compute(id_pe);
}
/*--------------------------------------------------------------------------------*/

void DumpAtom::init_style()
{
  atom_size_one = 5 + 6*stress_flag + 3*(velocity_flag + force_flag) + 4*displace_flag + pe_flag;
  if (image_flag) atom_size_one += 3;
  elem_size_one = node_size_one = 0;
  max_size_one = atom_size_one;

  int n;

  delete [] format;
  char *str;
  str = (char *) TAGINT_FORMAT " %d %g %g %g";
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


  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string

  delete [] columns;
  columns = new char[100];
  columns = strcpy(columns,"id type x y z");
  if (velocity_flag) 
    strcat(columns," vx vy vz");
  if (force_flag) 
    strcat(columns," fx fy fz");
  if (stress_flag) 
    strcat(columns," sxx syy szz sxy sxz syz");
  if (displace_flag) 
    strcat(columns," dx dy dz d");
  if (pe_flag) 
    strcat(columns," pe");

  // setup function ptr
  
  if (domain->triclinic == 0)
    header_choice = &DumpAtom::header_item;
  else
    header_choice = &DumpAtom::header_item_triclinic;

  pack_atom_choice = &DumpAtom::pack_atom_noscale_noimage;

  convert_atom_choice = &DumpAtom::convert_atom_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpAtom::write_string;
  }
  else {
    write_atom_choice = &DumpAtom::write_atom_lines_noimage;
  }

  // set stress/displace/pe compute ptr

  if (stress_flag) {
    int icompute = modify->find_compute(id_stress);
    if (icompute < 0) error->all(FLERR,"stress/atom ID for dump atom does not exist");
    stress = modify->compute[icompute];
  } else stress = NULL;

  if (displace_flag) {
    int icompute = modify->find_compute(id_displace);
    if (icompute < 0) error->all(FLERR,"displace/atom ID for dump atom does not exist");
    displace = modify->compute[icompute];
  } else displace = NULL;

  if (pe_flag) {
    int icompute = modify->find_compute(id_pe);
    if (icompute < 0) error->all(FLERR,"pe/atom ID for dump atom does not exist");
    pe = modify->compute[icompute];
  } else pe = NULL;


  // open single file, one time only

  if (multifile == 0) openfile();

}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%-1.16e %-1.16e\n",boxxlo,boxxhi);
  fprintf(fp,"%-1.16e %-1.16e\n",boxylo,boxyhi);
  fprintf(fp,"%-1.16e %-1.16e\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);

}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::header_item_triclinic(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxxlo,boxxhi,boxxy);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxylo,boxyhi,boxxz);
  fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxzlo,boxzhi,boxyz);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}


/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_header(bigint ndump)
{
}
/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_header(bigint ndump)
{
}

/*-------------------------------------------------------------------------------------------*/

void DumpAtom::write_node_header(bigint ndump)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_elem(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_node(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpAtom::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpAtom::pack_atom_noscale_noimage(tagint *ids)
{
  double **s,**d,*e;
  double ***nodes,***noded,**nodee;

  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    s = stress->array_atom;
    nodes = stress->array_node;
  } 
  if (displace_flag) {
    if (displace->invoked_peratom != update->ntimestep) {
      displace->compute_peratom();
      displace->addstep(update->ntimestep+nevery);
    }
    d = displace->array_atom;
    noded = displace->array_node;
  }
  if (pe_flag) {
    if (pe->invoked_peratom != update->ntimestep) {
      pe->compute_peratom();
      pe->addstep(update->ntimestep+nevery);
    }
    e = pe->vector_atom;
    nodee = pe->vector_node;
  }

  int m,n;
  int i,j,iintpl,iintg,inode;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  tagint maxtag,maxtag_all;
  maxtag = 0;

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    maxtag = MAX(maxtag,tag[i]);
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
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

  if (atom->natoms) {
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
  } else maxtag_all = 0;

  // pack atoms from elements 
  
  nlocal = element->nlocal;
  x = element->x;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int *nintpl = element->nintpl;
  int *nintg = element->nintg;
  mask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tag = element->tag;
  int *npe = element->npe;
  double ***shape_array = element->shape_array;
  int **i2ia = element->i2ia;
  double tmp;
  double *coord = new double[3];
  int ietype;

  // set IDs for interpolated atoms from element tag IDs to have consistent tag IDs across different simulations of the same system

  tagint itag = maxtag_all + 1;
  int maxintpl = element->maxintpl;
  int maxintg = element->maxintg;
  int maxnpe = element->max_npe;
  for (i = 0; i < nlocal; i++) {
    ietype = etype[i];
    if (mask[i] & groupbit) {
      if (dump_element_style == INTERPOLATE) {
        for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
          buf[m++] = itag + (tag[i]-1)*maxintpl + iintpl;
          buf[m++] = ctype[i];
          for (j = 0; j < 3; j++) {
            coord[j] = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              coord[j] += shape_array[ietype][iintpl][inode]*nodex[i][inode][j];
          }
          if (wrap_flag) domain->remap(coord);
          buf[m++] = coord[0];
          buf[m++] = coord[1];
          buf[m++] = coord[2];
          if (velocity_flag) {
            for (j = 0; j < 3; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodev[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (force_flag) {
            for (j = 0; j < 3; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodef[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (stress_flag) {
            for (j = 0; j < 6; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodes[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (displace_flag) {
            for (j = 0; j < 4; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*noded[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (pe_flag) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += shape_array[ietype][iintpl][inode]*nodee[i][inode];
            buf[m++] = tmp;
          }

        }
      } else if (dump_element_style == INTEGRATION) {
        for (iintg = 0; iintg < nintg[ietype]; iintg++) {
          iintpl = i2ia[ietype][iintg];
          buf[m++] = itag + (tag[i]-1)*maxintg + iintg;
          buf[m++] = ctype[i];
          for (j = 0; j < 3; j++) {
            coord[j] = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              coord[j] += shape_array[ietype][iintpl][inode]*nodex[i][inode][j];
          }
          if (wrap_flag) domain->remap(coord);
          buf[m++] = coord[0];
          buf[m++] = coord[1];
          buf[m++] = coord[2];

          if (velocity_flag) {
            for (j = 0; j < 3; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodev[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (force_flag) {
            for (j = 0; j < 3; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodef[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (stress_flag) {
            for (j = 0; j < 6; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*nodes[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (displace_flag) {
            for (j = 0; j < 4; j++) {
              tmp = 0.0;
              for (inode = 0; inode < npe[ietype]; inode++) 
                tmp += shape_array[ietype][iintpl][inode]*noded[i][inode][j];
              buf[m++] = tmp;
            }
          }
          if (pe_flag) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += shape_array[ietype][iintpl][inode]*nodee[i][inode];
            buf[m++] = tmp;
          }

        }
      } else if (dump_element_style == NODE) {
        for (inode = 0; inode < npe[ietype]; inode++) {
          buf[m++] = itag + (tag[i]-1)*maxnpe + inode;
          buf[m++] = ctype[i];
          coord[0] = nodex[i][inode][0];
          coord[1] = nodex[i][inode][1];
          coord[2] = nodex[i][inode][2];
          if (wrap_flag) domain->remap(coord);
          buf[m++] = coord[0];
          buf[m++] = coord[1];
          buf[m++] = coord[2];

          if (velocity_flag) {
            buf[m++] = nodev[i][inode][0];
            buf[m++] = nodev[i][inode][1];
            buf[m++] = nodev[i][inode][2];
          }
          if (force_flag) {
            buf[m++] = nodef[i][inode][0];
            buf[m++] = nodef[i][inode][1];
            buf[m++] = nodef[i][inode][2];
          }
          if (stress_flag) {
            buf[m++] = nodes[i][inode][0];
            buf[m++] = nodes[i][inode][1];
            buf[m++] = nodes[i][inode][2];
            buf[m++] = nodes[i][inode][3];
            buf[m++] = nodes[i][inode][4];
            buf[m++] = nodes[i][inode][5];
          }
          if (displace_flag) {
            buf[m++] = noded[i][inode][0];
            buf[m++] = noded[i][inode][1];
            buf[m++] = noded[i][inode][2];
            buf[m++] = noded[i][inode][3];
          }
          if (pe_flag) {
            buf[m++] = nodee[i][inode];
          }

        }
      } else {
        buf[m++] = itag + tag[i] - 1;
        buf[m++] = ctype[i];
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2];
        if (velocity_flag) {
          for (j = 0; j < 3; j++) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += nodev[i][inode][j];
            buf[m++] = tmp/npe[ietype];
          }
        }
        if (force_flag) {
          for (j = 0; j < 3; j++) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += nodef[i][inode][j];
            buf[m++] = tmp/npe[ietype];
          }
        }
        if (stress_flag) {
          for (j = 0; j < 6; j++) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += nodes[i][inode][j];
            buf[m++] = tmp/npe[ietype];
          }
        }
        if (displace_flag) {
          for (j = 0; j < 4; j++) {
            tmp = 0.0;
            for (inode = 0; inode < npe[ietype]; inode++) 
              tmp += noded[i][inode][j];
            buf[m++] = tmp/npe[ietype];
          }
        }
        if (pe_flag) {
          tmp = 0.0;
          for (inode = 0; inode < npe[ietype]; inode++) 
            tmp += nodee[i][inode];
          buf[m++] = tmp/npe[ietype];
        }
      }
    }
  }  

  delete [] coord;
}

/*--------------------------------------------------------------------------------------------*/
int DumpAtom::convert_node_string(int n, int size_one, double *mybuf)
{
  return 0;
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_elem_string(int n, int size_one, double *mybuf)
{
  return 0;
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_atom_noimage(int n, int size_one, double *mybuf)
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
        static_cast<tagint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += 5;
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

void DumpAtom::write_node_data(int n, double *mybuf)
{
}

/*----------------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_data(int n, double *mybuf)
{
}

/*---------------------------------------------------------------------------------------------------*/

void DumpAtom::write_string(int n, double *mybuf)
{
  fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
        static_cast<tagint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += 5;
    if (velocity_flag) {
      fprintf(fp,format_velocity,mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
    if (force_flag) {
      fprintf(fp,format_force,mybuf[m],mybuf[m+1],mybuf[m+2]);
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

/* ---------------------------------------------------------------------- */

int DumpAtom::count_atoms()
{
  int m = 0;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == 0 || igroup == 2) {
    m = atom->nlocal;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) m++;
  }

  mask = element->mask;
  nlocal = element->nlocal;
  int *nintpl = element->nintpl;
  int *nintg = element->nintg;
  int *npe = element->npe;
  int *etype = element->etype;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (dump_element_style == INTERPOLATE) 
        m += nintpl[etype[i]];
      else if (dump_element_style == INTEGRATION)
        m += nintg[etype[i]];
      else if (dump_element_style == NODE)
        m += npe[etype[i]];
      else m++;
    }
  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::count_elements()
{
  return 0;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::count_nodes()
{
  return 0;
}
