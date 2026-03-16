#include <string.h>
#include "dump_tec_plot.h"
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

/*--------------------------------------------------------------------------------*/

DumpTecPlot::DumpTecPlot(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump cac command");
  nevery = universe->inumeric(FLERR,arg[3]);
  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;

  nclusters = 0;
  maxclusters = 0;
  clusters = NULL;
  cutlink = 0.0;

  int iarg = 5;
  stress_flag = 0;
  force_flag = 0;
  displace_flag = 0;
  link_flag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) stress_flag = 1;
    else if (strcmp(arg[iarg],"force") == 0) force_flag = 1;
    else if (strcmp(arg[iarg],"displace") == 0) displace_flag = 1;
    else if (strcmp(arg[iarg],"link") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"cut") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal dump tecplot command");
        cutlink = universe->numeric(FLERR,arg[iarg+2]);
        link_flag = 1;
        iarg++;
      } else if (strcmp(arg[iarg+1],"id") == 0) {
        link_flag = 2;
      } else error->all(FLERR,"Illegal dump tecplot command");
      if (igroup == 2) link_flag = 0;
      else if (igroup != 0 && igroup != 1) error->all(FLERR,"Group all or element required for element link");
      iarg++;
    } else error->all(FLERR,"Illegal dump tecplot command");
    iarg++;
  }

  apc = element->apc;

  if (apc == 1) link_flag = 0;

  int n;

  if (stress_flag) {
    n = strlen(id) + strlen("_DUMP_STRESS") + 1;
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
    n = strlen(id) + strlen("_DUMP_DISPLACE") + 1;
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

  if (link_flag) {
    if (stress_flag) comm_elem_forward += 6*element->npe;
    if (force_flag) comm_elem_forward += 6*element->npe;
    if (displace_flag) comm_elem_forward += 4*element->npe;
  }
}

/*--------------------------------------------------------------------------------*/
DumpTecPlot::~DumpTecPlot()
{
  if (link_flag) memory->destroy(clusters);
  if (modify->ncompute && stress_flag) modify->delete_compute(id_stress);
  if (modify->ncompute && displace_flag) modify->delete_compute(id_displace);

}
/*--------------------------------------------------------------------------------*/

void DumpTecPlot::setup_pre_write() 
{
  nclusters = 0;
  int a = comm->me;
  int b;
  //MPI_Allreduce(&bnclusters,&totalclusters,1,MPI_CAC_BIGINT,MPI_SUM,world);
  // setup link flag and create cluster list
  
  if (element->nelements == 0 || !link_flag) return;
  if (element->nmax > maxclusters*apc) {
    memory->destroy(clusters);
    maxclusters = element->nmax/apc+1;
    memory->create(clusters,maxclusters,apc,"dump:clusters");
  }
  int nlocal = element->nlocal;
  int nghost = element->nghost;
  int nall = nlocal + nghost;
  int i,j,k;
  tagint *tag = element->tag;
  tagint itag;
  double **x = element->x;
  int *mask = element->mask;
  double ix,iy,iz,delx,dely,delz,rsq;

  int own_flag;

  if (link_flag == 1) {
    if (nlocal) {
      int *set_flag = new int[nlocal+nghost];
      for (i = 0; i < nall; i++)
        set_flag[i] = 0;

      double cutsq = cutlink*cutlink;
      for (i = 0; i < nlocal; i++) {
        if (set_flag[i]) continue;
        itag = tag[i];
        set_flag[i] = 1;
        ix = x[i][0]; 
        iy = x[i][1]; 
        iz = x[i][2]; 
        clusters[nclusters][0] = i;
        int n = 1;
        for (j = 0; j < nall; j++) {
          if (set_flag[j]) continue;
          delx = ix - x[j][0]; 
          dely = iy - x[j][1]; 
          delz = iz - x[j][2]; 
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq >= cutsq) continue;
          if (n >= apc) error->one(FLERR,"Too many elements in a cluster, reduce cut");
          clusters[nclusters][n++] = j;
          set_flag[j] = 1;
        }
        if (n != apc) error->one(FLERR,"Not enough elements in a cluster, increase cut");
        own_flag = 1;

        // own this cluster if itag is the smallest

        for (j = 1; j < apc; j++) {
          if (tag[clusters[nclusters][j]] <= itag) 
            own_flag = 0;
        }
        if (own_flag) {
          nclusters++;
        } else {
          for (j = 0; j < apc; j++)
            set_flag[clusters[nclusters][j]] = 0;
        }
      }
      delete [] set_flag;
    }
  } else {
    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }

    for (i = 0; i < nlocal; i++) 
      if ((mask[i] & groupbit)) {
        itag = tag[i];
        if ((itag-1) % apc) continue;
        clusters[nclusters][0] = i;
        for (k = 1; k < apc; k++) {
          j = element->map(itag+k);
          if (j < 0) error->one(FLERR,"ghost element not correct for link");
          clusters[nclusters][k] = j;
        }
        nclusters++;
      }
    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }
  }

  bigint bnclusters = nclusters;
  bigint totalclusters;

  MPI_Allreduce(&bnclusters,&totalclusters,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (totalclusters != element->nelements/apc) {
    error->all(FLERR,"Linking element incorrecly");
  }
}

/*--------------------------------------------------------------------------------*/

void DumpTecPlot::init_style() 
{

  atom_size_one = node_size_one = 
    3 + 6*(stress_flag + force_flag) + 4*displace_flag;
  elem_size_one = element->npe;

  if (link_flag)
    if (element->nelements%apc != 0)
      error->all(FLERR,"Invalid number of elements for linking");
  int n;

  // format for writing atom and node section 

  delete [] format;
  char *str;
  str = (char *) "%g %g %g";
  n = strlen(str);
  format = new char[n];
  strcpy(format,str);

  if (stress_flag) {
    delete [] format_stress;
    str = (char *) " %g %g %g %g %g %g";
    n = strlen(str);
    format_stress = new char[n];
    strcpy(format_stress,str);
  }

  if (force_flag) {
    delete [] format_force;
    str = (char *) " %g %g %g %g %g %g";
    n = strlen(str);
    format_force = new char[n];
    strcpy(format_force,str);
  }

  if (displace_flag) {
    delete [] format_displace;
    str = (char *) " %g %g %g %g";
    n = strlen(str);
    format_displace = new char[n];
    strcpy(format_displace,str);
  }

  // format for writing element section (node connectivity)

  delete [] format_elem;
  str = (char *) TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT;
  n = strlen(str) + 2;
  format_elem = new char[n];
  strcpy(format_elem, str);
  strcat(format_elem,"\n");

  // setup column string 

  delete [] columns; 
  columns = new char[500];
  columns = strcpy(columns,"\"x\", \"y\", \"z\"");
  if (stress_flag) 
    strcat(columns,", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"");
  if (force_flag) 
    strcat(columns,", \"vx\", \"vy\", \"vz\", \"fx\", \"fy\", \"fz\"");
  if (displace_flag) 
    strcat(columns,", \"dx\", \"dy\", \"dz\", \"d\"");

  header_choice = &DumpTecPlot::header_item;
  atom_header_choice = &DumpTecPlot::atom_header_item;
  elem_header_choice = &DumpTecPlot::elem_header_item;
  node_header_choice = &DumpTecPlot::node_header_item;

  pack_atom_choice = &DumpTecPlot::pack_atom_noscale_noimage;
  pack_node_choice = &DumpTecPlot::pack_node_noscale_noimage;
  pack_elem_choice = &DumpTecPlot::pack_elem_noscale_noimage;

  convert_atom_choice = &DumpTecPlot::convert_atom_noimage;
  convert_node_choice = &DumpTecPlot::convert_node_noimage;
  convert_elem_choice = &DumpTecPlot::convert_elem_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpTecPlot::write_string;
    write_node_choice = &DumpTecPlot::write_string;
    write_elem_choice = &DumpTecPlot::write_string;
  } else {
    write_atom_choice = &DumpTecPlot::write_lines_noimage;
    write_elem_choice = &DumpTecPlot::write_elem_lines_noimage;
    write_node_choice = &DumpTecPlot::write_lines_noimage;
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
}

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::header_item(bigint ndump)
{
  fprintf(fp, "title = \"TIMESTEP "BIGINT_FORMAT"\"\n",update->ntimestep);
  fprintf(fp, "variables = %s\n", columns);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::node_header_item(bigint ndump)
{
  int npe = element->npe;
  bigint nelem = ndump/npe;
  fprintf(fp,"zone t = \"Coarse Element\",");
  fprintf(fp," n = "BIGINT_FORMAT" e = "BIGINT_FORMAT, ndump, nelem);
  fprintf(fp," datapacking = point, zonetype = febrick\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::atom_header_item(bigint ndump)
{
  fprintf(fp,"zone t = \"Discrete Atoms\", f = point\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::elem_header_item(bigint ndump)
{
  // no header for element section (node connectivity)
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_elem_header(bigint ndump)
{
  if (multiproc) (this->*elem_header_choice)(ndump);
  else if (filewriter) (this->*elem_header_choice)(ndump);
}
/*-------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_atom_header(bigint ndump)
{
  if (multiproc) (this->*atom_header_choice)(ndump);
  else if (filewriter) (this->*atom_header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_node_header(bigint ndump)
{
  if (multiproc) (this->*node_header_choice)(ndump);
  else if (filewriter) (this->*node_header_choice)(ndump);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_elem(tagint *ids)
{
  (this->*pack_elem_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_node(tagint *ids)
{
  (this->*pack_node_choice)(ids);
}

/*-----------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_elem_noscale_noimage(tagint *ids)
{
  int *mask = element->mask;
  int nlocal = element->nlocal;
  tagint **nodetag = element->nodetag;

  int npe = element->npe;
  int m = 0;

  if (link_flag) {
    tagint maxtag = nclusters * npe;
    tagint maxtag_sum;
    MPI_Scan(&maxtag,&maxtag_sum,1,MPI_CAC_TAGINT,MPI_SUM,world);
    tagint itag = maxtag_sum - maxtag + 1;

    for (int i = 0; i < nclusters; i++)
      for (int j = 0; j < npe; j++)
        ebuf[m++] = itag++;

  } else {

    // reset node tag to for node connectivity

    element->reset_node_tags(groupbit);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        for (int j = 0; j < npe; j++) 
          ebuf[m++] = nodetag[i][j];

  }
} 

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_atom_noscale_noimage(tagint *ids)
{
  double **s,**d;
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
  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      abuf[m++] = x[i][0]; 
      abuf[m++] = x[i][1]; 
      abuf[m++] = x[i][2];
      if (stress_flag) {
        abuf[m++] = s[i][0]; 
        abuf[m++] = s[i][1]; 
        abuf[m++] = s[i][2];
        abuf[m++] = s[i][3]; 
        abuf[m++] = s[i][4]; 
        abuf[m++] = s[i][5];
      }
      if (force_flag) {
        abuf[m++] = v[i][0]; 
        abuf[m++] = v[i][1]; 
        abuf[m++] = v[i][2];
        abuf[m++] = f[i][0]; 
        abuf[m++] = f[i][1]; 
        abuf[m++] = f[i][2];
      }
      if (displace_flag) {
        abuf[m++] = d[i][0]; 
        abuf[m++] = d[i][1]; 
        abuf[m++] = d[i][2];
        abuf[m++] = d[i][3]; 
      }
      if (ids) ids[n++] = i;
    }
  }
}


/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_node_noscale_noimage(tagint *ids)
{
  double ***nodes,***noded;
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

  int m,n,ii,i,j,k;
  int *mask = element->mask;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int nlocal = element->nlocal;
  int npe = element->npe;

  m = n = 0;

  if (link_flag) {

    // forward comm of ghost stress or force/velocity 

    comm->forward_comm_dump(this);

    double ix,iy,iz,s0,s1,s2,s3,s4,s5;
    double vx,vy,vz,fx,fy,fz,d0,d1,d2,d3;
    for (ii = 0; ii < nclusters; ii++) 
      for (j = 0; j < npe; j++) {
        ix = iy = iz = 0.0;
        vx = vy = vz = 0.0;
        fx = fy = fz = 0.0;
        s0 = s1 = s2 = 0.0;
        s3 = s4 = s5 = 0.0;
        d0 = d1 = d2 = d3 = 0.0;
        for (k = 0; k < apc; k++) {
          i = clusters[ii][k];
          ix += nodex[i][j][0]; 
          iy += nodex[i][j][1]; 
          iz += nodex[i][j][2];
          if (stress_flag) {
            s0 += nodes[i][j][0]; 
            s1 += nodes[i][j][1];
            s2 += nodes[i][j][2]; 
            s3 += nodes[i][j][3];
            s4 += nodes[i][j][4]; 
            s5 += nodes[i][j][5];
          }
          if (force_flag) {
            vx += nodev[i][j][0]; 
            vy += nodev[i][j][1]; 
            vz += nodev[i][j][2];
            fx += nodef[i][j][0]; 
            fy += nodef[i][j][1]; 
            fz += nodef[i][j][2];
          }
          if (displace_flag) {
            d0 += noded[i][j][0]; 
            d1 += noded[i][j][1];
            d2 += noded[i][j][2]; 
            d3 += noded[i][j][3];
          } 
        }
        nbuf[m++] = ix/apc; 
        nbuf[m++] = iy/apc; 
        nbuf[m++] = iz/apc;
        if (stress_flag) {
          nbuf[m++] = s0/apc; 
          nbuf[m++] = s1/apc; 
          nbuf[m++] = s2/apc; 
          nbuf[m++] = s3/apc; 
          nbuf[m++] = s4/apc; 
          nbuf[m++] = s5/apc;
        }
        if (force_flag) {
          nbuf[m++] = vx/apc; 
          nbuf[m++] = vy/apc; 
          nbuf[m++] = vz/apc;
          nbuf[m++] = fx/apc; 
          nbuf[m++] = fy/apc; 
          nbuf[m++] = fz/apc;
        }
        if (displace_flag) {
          nbuf[m++] = d0/apc; 
          nbuf[m++] = d1/apc; 
          nbuf[m++] = d2/apc; 
          nbuf[m++] = d3/apc; 
        }
      }
  } else 
    for (i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
        for (j = 0; j < npe; j++) {
          nbuf[m++] = nodex[i][j][0]; 
          nbuf[m++] = nodex[i][j][1]; 
          nbuf[m++] = nodex[i][j][2];
          if (stress_flag) {
            nbuf[m++] = nodes[i][j][0]; 
            nbuf[m++] = nodes[i][j][1]; 
            nbuf[m++] = nodes[i][j][2]; 
            nbuf[m++] = nodes[i][j][3]; 
            nbuf[m++] = nodes[i][j][4]; 
            nbuf[m++] = nodes[i][j][5];
          }
          if (force_flag) {
            nbuf[m++] = nodev[i][j][0]; 
            nbuf[m++] = nodev[i][j][1]; 
            nbuf[m++] = nodev[i][j][2]; 
            nbuf[m++] = nodef[i][j][0]; 
            nbuf[m++] = nodef[i][j][1]; 
            nbuf[m++] = nodef[i][j][2];
          }
          if (displace_flag) {
            nbuf[m++] = noded[i][j][0]; 
            nbuf[m++] = noded[i][j][1]; 
            nbuf[m++] = noded[i][j][2]; 
            nbuf[m++] = noded[i][j][3]; 
          }
        }
        if (ids) ids[n++] = i;
      }
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_node_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_node_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_elem_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_elem_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_node_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsnbuf) {
      if ((bigint) maxsnbuf + DELTA > MAXSMALLINT) return -1;
      maxsnbuf += DELTA;
      memory->grow(snbuf,maxsnbuf,"dump:snbuf");
    }
    offset += sprintf(&snbuf[offset],format,
        mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (stress_flag) {
      offset += sprintf(&snbuf[offset],format_stress,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (force_flag) {
      offset += sprintf(&snbuf[offset],format_force,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      offset += sprintf(&snbuf[offset],format_displace,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    offset += sprintf(&snbuf[offset],"\n");
  }
  return offset;
}
/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_atom_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsabuf) {
      if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
      maxsabuf += DELTA;
      memory->grow(sabuf,maxsabuf,"dump:sabuf");
    }
    offset += sprintf(&sabuf[offset],format,
        mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (stress_flag) {
      offset += sprintf(&sabuf[offset],format_stress,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (force_flag) {
      offset += sprintf(&sabuf[offset],format_force,
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      offset += sprintf(&sabuf[offset],format_displace,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    offset += sprintf(&sabuf[offset],"\n");
  }
  return offset;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_elem_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsebuf) {
      if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
      maxsebuf += DELTA;
      memory->grow(sebuf,maxsebuf,"dump:sebuf");
    }
    offset += sprintf(&sebuf[offset],format_elem,
        static_cast<tagint>(mybuf[m]),
        static_cast<tagint>(mybuf[m+1]),
        static_cast<tagint>(mybuf[m+2]),
        static_cast<tagint>(mybuf[m+3]),
        static_cast<tagint>(mybuf[m+4]),
        static_cast<tagint>(mybuf[m+5]),
        static_cast<tagint>(mybuf[m+6]),
        static_cast<tagint>(mybuf[m+7]));
    m += 8;
  }
  return offset;
}	

/*----------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_node_data(int n, double *mybuf)
{
  (this->*write_node_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_elem_data(int n, double *mybuf)
{
  (this->*write_elem_choice) (n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (stress_flag) {
      fprintf(fp,format_stress,mybuf[m],mybuf[m+1],
          mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (force_flag) {
      fprintf(fp,format_force,mybuf[m],mybuf[m+1],
          mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (displace_flag) {
      fprintf(fp,format_displace,mybuf[m],mybuf[m+1],
          mybuf[m+2],mybuf[m+3]);
      m += 4;
    }
    fprintf(fp,"\n");
  }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpTecPlot::write_elem_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_elem, static_cast<tagint>(mybuf[m]),
        static_cast<tagint>(mybuf[m+1]),
        static_cast<tagint>(mybuf[m+2]),
        static_cast<tagint>(mybuf[m+3]),
        static_cast<tagint>(mybuf[m+4]),
        static_cast<tagint>(mybuf[m+5]),
        static_cast<tagint>(mybuf[m+6]),
        static_cast<tagint>(mybuf[m+7]));
    m += elem_size_one;
  }
}

/* ---------------------------------------------------------------------- */

int DumpTecPlot::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int ii,i,j,m;
  int npe = element->npe;
  m = 0;

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
  if (force_flag) {
    double ***nodef = element->nodef;
    double ***nodev = element->nodev;
    for (ii = 0; ii < n; ii++) {
      i = list[ii];
      for (j = 0; j < npe; j++) {
        buf[m++] = nodev[i][j][0];
        buf[m++] = nodev[i][j][1];
        buf[m++] = nodev[i][j][2];
        buf[m++] = nodef[i][j][0];
        buf[m++] = nodef[i][j][1];
        buf[m++] = nodef[i][j][2];
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

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpTecPlot::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i,j,m,last;
  last = first + n;
  int npe = element->npe;
  m = 0;

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
  if (force_flag) {
    double ***nodef = element->nodef;
    double ***nodev = element->nodev;
    for (i = first; i < last; i++) 
      for (j = 0; j < npe; j++) { 
        nodev[i][j][0] = buf[m++];
        nodev[i][j][1] = buf[m++];
        nodev[i][j][2] = buf[m++];
        nodef[i][j][0] = buf[m++];
        nodef[i][j][1] = buf[m++];
        nodef[i][j][2] = buf[m++];
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
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::count_elements()
{
  if (link_flag) return nclusters;
  else if (igroup == 0 || igroup == 1) return element->nlocal;
  int m = 0;
  int *mask = element->mask;
  for (int i = 0; i < element->nlocal; i++) 
    if (mask[i] & groupbit) m++;
  return m;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::count_nodes()
{
  return count_elements()*element->npe;
}
