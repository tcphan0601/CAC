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
  format_default = NULL;
  format = NULL;

  nclusters = 0;
  maxclusters = 0;
  clusters = NULL;
  cutlink = 0.0;

  int iarg = 5;
  stress_flag = 0;
  force_flag = 0;
  link_flag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) {
      stress_flag = 1;
      force_flag = 0;
    } else if (strcmp(arg[iarg],"force") == 0) {
      stress_flag = 0;
      force_flag = 1;
    } else if (strcmp(arg[iarg],"link") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump tecplot command");
      if (strcmp(arg[iarg+1],"cut") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal dump tecplot command");
        cutlink = universe->numeric(FLERR,arg[iarg+2]);
        link_flag = 1;
        iarg++;
      } else if (strcmp(arg[iarg+1],"tag") == 0) {
        link_flag = 2;
      } else error->all(FLERR,"Illegal dump tecplot command");
      if (igroup != 0 && igroup != 1) error->all(FLERR,"Group all or element required for element link");
      iarg++;
    } else error->all(FLERR,"Illegal dump tecplot command");
    iarg++;
  }

  apc = element->apc;

  if (apc == 1) link_flag = 0;
  int n = strlen(id) + 8;
  id_stress = new char[n];
  strcpy(id_stress,id);
  strcat(id_stress,"_stress");

  if (stress_flag) {
    char **newarg = new char*[4];
    newarg[0] = id_stress;
    newarg[1] = (char *) "all";
    newarg[2] = (char *) "stress/atom";
    newarg[3] = (char *) "NULL";
    modify->add_compute(4,newarg);
    delete [] newarg;
  }


}

/*--------------------------------------------------------------------------------*/
DumpTecPlot::~DumpTecPlot()
{
  if (stress_flag) delete [] id_stress;
  if (link_flag) memory->destroy(clusters);
}
/*--------------------------------------------------------------------------------*/

void DumpTecPlot::setup_pre_write() 
{
  if (element->nelements == 0 || !link_flag) return;
  if (element->nlocal == 0) return;
  if (element->nmax/apc > maxclusters) {
    memory->destroy(clusters);
    maxclusters = element->nmax/apc;
    memory->create(clusters,maxclusters,apc,"dump:clusters");
  }
  int nlocal = element->nlocal;
  int nghost = element->nghost;
  int nall = nlocal + nghost;
  int *set_flag = new int[nlocal+nghost];
  int i,j,k;
  tagint *tag = element->tag;
  tagint itag;
  double **x = element->x;
  double ix,iy,iz,delx,dely,delz,rsq;
  for (i = 0; i < nall; i++)
    set_flag[i] = 0;
  nclusters = 0;

  int own_flag;
  if (link_flag == 1) {
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
  } else {
    error->all(FLERR,"ID Link not finished yet");
  }
  bigint bnclusters = nclusters;
  bigint totalclusters;
  MPI_Allreduce(&bnclusters,&totalclusters,1,MPI_CAC_BIGINT,MPI_SUM,world);

  if (totalclusters != element->nelements/apc)
    error->all(FLERR,"Linking element incorrecly");
  delete [] set_flag;
}

/*--------------------------------------------------------------------------------*/

void DumpTecPlot::init_style() 
{
  atom_size_one = node_size_one = 3;
  if (stress_flag || force_flag) {
    atom_size_one += 6;
    node_size_one += 6;
  }
  elem_size_one = element->npe;

  if (link_flag)
    if (element->nelements%apc != 0)
      error->all(FLERR,"Invalid number of elements for linking");
  int n;

  // format for writing atom and node section 

  delete [] format;
  char *str;
  if (stress_flag || force_flag)
    str = (char *) "%g %g %g %g %g %g %g %g %g";
  else
    str = (char *) "%g %g %g";
  n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // format for writing element section (node connectivity)

  delete [] format_elem;
  str = (char *) TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT;
  n = strlen(str) + 2;
  format_elem = new char[n];
  strcpy(format_elem, str);
  strcat(format_elem,"\n");

  if (force_flag)
    columns = (char *) "\"x\", \"y\", \"z\", \"vx\", \"vy\", \"vz\", \"fx\", \"fy\", \"fz\"";
  else if (stress_flag)
    columns = (char *) "\"x\", \"y\", \"z\", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"";
  else
    columns = (char *) "\"x\", \"y\", \"z\"";

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

  // set stress compute pointer

  if (stress_flag) {
    int icompute = modify->find_compute(id_stress);
    if (icompute < 0) error->all(FLERR,"Stress/atom ID for dump tecplot do not exist");
    stress = modify->compute[icompute];
  }
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

void DumpTecPlot::write_header(bigint ndump){
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

    int *mask = element->mask;
    int nlocal = element->nlocal;
    tagint **nodetag = element->nodetag;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        for (int j = 0; j < npe; j++)
          ebuf[m++] = nodetag[i][j];
  }
} 

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_atom_noscale_noimage(tagint *ids)
{
  double **s;
  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    s = stress->array_atom;
  }

  //double **s = stress->array_atom;
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
      if (ids) ids[n++] = i;

    }
  }

}


/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_node_noscale_noimage(tagint *ids)
{
  double ***nodes;
  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    nodes = stress->array_node;
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
    double ix,iy,iz,s0,s1,s2,s3,s4,s5,vx,vy,vz,fx,fy,fz;
    for (ii = 0; ii < nclusters; ii++) 
      for (j = 0; j < npe; j++) {
        ix = iy = iz = 0.0;
        vx = vy = vz = 0.0;
        fx = fy = fz = 0.0;
        s0 = s1 = s2 = 0.0;
        s3 = s4 = s5 = 0.0;
        for (k = 0; k < apc; k++) {
          i = clusters[ii][k];
          ix += nodex[i][j][0];
          iy += nodex[i][j][1];
          iz += nodex[i][j][2];
          if (stress_flag) {
            s0 += nodes[i][j][0]; s1 += nodes[i][j][1];
            s2 += nodes[i][j][2]; s3 += nodes[i][j][3];
            s4 += nodes[i][j][4]; s5 += nodes[i][j][5];
          }
          if (force_flag) {
            vx += nodev[i][j][0]; fx += nodev[i][j][0];
            vy += nodev[i][j][1]; fy += nodev[i][j][1];
            vz += nodev[i][j][2]; fz += nodev[i][j][2];
          }
        }
        nbuf[m++] = ix/apc; 
        nbuf[m++] = iy/apc; 
        nbuf[m++] = iz/apc;
        if (stress_flag) {
          nbuf[m++] = s0/apc; nbuf[m++] = s1/apc;
          nbuf[m++] = s2/apc; nbuf[m++] = s3/apc;
          nbuf[m++] = s4/apc; nbuf[m++] = s5/apc;
        }
        if (force_flag) {
          nbuf[m++] = vx/apc; nbuf[m++] = fx/apc;
          nbuf[m++] = vy/apc; nbuf[m++] = fy/apc;
          nbuf[m++] = vz/apc; nbuf[m++] = fz/apc;
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
            nbuf[m++] = nodes[i][j][0]; nbuf[m++] = nodes[i][j][1];
            nbuf[m++] = nodes[i][j][2]; nbuf[m++] = nodes[i][j][3];
            nbuf[m++] = nodes[i][j][4]; nbuf[m++] = nodes[i][j][5];
          }
          if (force_flag) {
            nbuf[m++] = nodev[i][j][0]; nbuf[m++] = nodev[i][j][1];
            nbuf[m++] = nodev[i][j][2]; nbuf[m++] = nodef[i][j][0];
            nbuf[m++] = nodef[i][j][1]; nbuf[m++] = nodef[i][j][2];
          }
        }
        if (ids) ids[n++] = i;
      }
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_atom_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_node_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_node_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_elem_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_elem_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_node_noimage(int n,int size_one,double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsnbuf) {
      if ((bigint) maxsnbuf + DELTA > MAXSMALLINT) return -1;
      maxsnbuf += DELTA;
      memory->grow(snbuf, maxsnbuf,"dump:snbuf");
    }
    if (stress_flag || force_flag) {
      offset += sprintf(&snbuf[offset],format,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],
          mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],
          mybuf[m+8]);
      m += 9;
    } else {
      offset += sprintf(&snbuf[offset],format,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
  }
  return offset;
}
/*--------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_atom_noimage(int n,int size_one,double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsabuf) {
      if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
      maxsabuf += DELTA;
      memory->grow(sabuf, maxsabuf,"dump:sabuf");
    }
    if (stress_flag || force_flag) {
      offset += sprintf(&sabuf[offset],format,
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],
          mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],
          mybuf[m+8]);
      m += 9;
    } else {
      offset += sprintf(&sabuf[offset],format,
          mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
  }
  return offset;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::convert_elem_noimage(int n,int size_one,double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i <n;i++) {
    if (offset + ONELINE > maxsebuf) {
      if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
      maxsebuf += DELTA;
      memory->grow(sebuf, maxsebuf,"dump:sabuf");
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
    m += size_one;
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
    if (stress_flag || force_flag) {
      fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],
          mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8]);
      m += 9;
    } else {
      fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += 3;
    }
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

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::count_elements()
{
  if (link_flag) return nclusters;
  return element->nlocal;
}

/*----------------------------------------------------------------------------------------------------*/

int DumpTecPlot::count_nodes()
{
  if (link_flag) return nclusters*element->npe;
  return element->nlocal*element->npe;
}
