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
  if (narg != 5) error->all(FLERR,"Illegal dump cac command");
  nevery = universe->inumeric(FLERR,arg[3]);
  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;
  format_default = NULL;
  format = NULL;

  int n = strlen(id) + 8;
  id_stress = new char[n];
  strcpy(id_stress,id);
  strcat(id_stress,"_stress");

  char **newarg = new char*[4];
  newarg[0] = id_stress;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "stress/atom";
  newarg[3] = (char *) "NULL";
  modify->add_compute(4,newarg);
  delete [] newarg;

}

/*--------------------------------------------------------------------------------*/
DumpTecPlot::~DumpTecPlot()
{
  delete [] id_stress;
}
/*--------------------------------------------------------------------------------*/

void DumpTecPlot::init_style(){
  if (image_flag == 0) atom_size_one = node_size_one = 9;
  else atom_size_one = node_size_one = 12;

  elem_size_one = element->npe;

  int n;

  // format for writing atom and node section 
  
  delete [] format;
  char *str;
  str = (char *) "%g %g %g %g %g %g %g %g %g";
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

  columns = (char *) "\"x\", \"y\", \"z\", \"sxx\", \"syy\", \"szz\", \"sxy\", \"sxz\", \"syz\"";

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
  
  int icompute = modify->find_compute(id_stress);
  if (icompute < 0) error->all(FLERR,"Stress/atom ID for dump tecplot do not exist");
  stress = modify->compute[icompute];
}
/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::header_item(bigint ndump){
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
  // reset node tag to for node connectivity
  
  element->reset_node_tags(groupbit);

  int *mask = element->mask;
  int nlocal = element->nlocal;
  int npe = element->npe;
  tagint **nodetag = element->nodetag;
  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int j = 0; j < npe; j++)
        ebuf[m++] = nodetag[i][j];
} 

/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_atom_noscale_noimage(tagint *ids)
{
  if (stress->invoked_peratom != update->ntimestep) {
    stress->compute_peratom();
    stress->addstep(update->ntimestep+nevery);
  }

  double **s = stress->array_atom;
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
      abuf[m++] = s[i][0];
      abuf[m++] = s[i][1];
      abuf[m++] = s[i][2];
      abuf[m++] = s[i][3];
      abuf[m++] = s[i][4];
      abuf[m++] = s[i][5];
      if (ids) ids[n++] = i;

    }
  }

}


/*------------------------------------------------------------------------------------------*/

void DumpTecPlot::pack_node_noscale_noimage(tagint *ids)
{
  if (stress->invoked_peratom != update->ntimestep) {
    stress->compute_peratom();
    stress->addstep(update->ntimestep+nevery);
  }

  double ***nodes = stress->array_node;
  int m,n;
  int *mask = element->mask;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int nlocal = element->nlocal;
  int npe = element->npe;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < npe; j++) {
        nbuf[m++] = nodex[i][j][0];
        nbuf[m++] = nodex[i][j][1];
        nbuf[m++] = nodex[i][j][2];
        nbuf[m++] = nodes[i][j][0];
        nbuf[m++] = nodes[i][j][1];
        nbuf[m++] = nodes[i][j][2];
        nbuf[m++] = nodes[i][j][3];
        nbuf[m++] = nodes[i][j][4];
        nbuf[m++] = nodes[i][j][5];
      }
      if (ids) ids[n++] = i;
    }
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
    offset += sprintf(&snbuf[offset], format,
        mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],
        mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],
        mybuf[m+8]);
    m += size_one;
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
    offset += sprintf(&sabuf[offset], format,
        mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],
        mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],
        mybuf[m+8]);
    m += size_one;
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
    fprintf(fp,format,mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],
        mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8]);
    m += 9;
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



