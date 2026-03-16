#include <string.h>
#include "dump_ids_atom.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

// same as compute_ida_atom.cpp

enum{OTHER,CUBIC2,HEX2,CUBIC1,HEX1,CUBIC,HEX,UNKNOWN};

/*--------------------------------------------------------------------------------*/

DumpIDSAtom::DumpIDSAtom(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
  
  if (narg < 6 || narg > 11) error->all(FLERR,"Illegal dump ids command");

  int n = strlen(id) + 8;
  id_ids = new char[n];
  strcpy(id_ids,id);
  strcat(id_ids,"_ids_atom");
  char **newarg = new char*[3];
  newarg[0] = id_ids;
  newarg[1] = arg[1];
  newarg[2] = (char *) "ids/atom";
  modify->add_compute(3,newarg);

  int iarg = 5;

  allflag = 0;
  nkeep = 0;
  while (iarg < narg ) {
    if (strcmp(arg[iarg],"c2nd") == 0) keep[nkeep++] = CUBIC2;
    else if (strcmp(arg[iarg],"c1st") == 0) keep[nkeep++] = CUBIC1;
    else if (strcmp(arg[iarg],"cubic") == 0) keep[nkeep++] = CUBIC;
    else if (strcmp(arg[iarg],"h2nd") == 0) keep[nkeep++] = HEX2;
    else if (strcmp(arg[iarg],"h1st") == 0) keep[nkeep++] = HEX1;
    else if (strcmp(arg[iarg],"hex") == 0) keep[nkeep++] = HEX;
    else if (strcmp(arg[iarg],"other") == 0) keep[nkeep++] = OTHER;
    else if (strcmp(arg[iarg],"all") == 0) allflag = 1;
    else error->all(FLERR,"Illegal structure type to keep");
    iarg++;
  }
  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;
}

/*--------------------------------------------------------------------------------*/

DumpIDSAtom::~DumpIDSAtom()
{
  delete [] id_ids;
}


/*--------------------------------------------------------------------------------*/

void DumpIDSAtom::init_style() 
{
  if (image_flag == 0) atom_size_one = 6;
  else atom_size_one = 9;

  node_size_one = elem_size_one =  0;

  int n;

  delete [] format;
  char *str;
  if (image_flag == 0) str = (char *) TAGINT_FORMAT " %d %g %g %g %d";
  n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string
  
  delete [] columns;
  columns = new char[500];
  if (image_flag == 0) 
    columns = strcpy(columns,"id type x y z structure_id");
  else
    columns = strcpy(columns,"id type x y z structure_id ix iy iz");

  header_choice = &DumpIDSAtom::header_item;

  pack_atom_choice = &DumpIDSAtom::pack_atom_noscale_noimage;

  convert_atom_choice = &DumpIDSAtom::convert_atom_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpIDSAtom::write_string;
  }
  else {
    write_atom_choice = &DumpIDSAtom::write_atom_lines_noimage;
  }

  int icompute = modify->find_compute(id_ids);
  if (icompute < 0) error->all(FLERR,"compute ids/atom ID for dump ids/atom does not exist");
  compute_ids = modify->compute[icompute];

}
/*------------------------------------------------------------------------------------------*/

void DumpIDSAtom::header_item(bigint ndump) 
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

/*-------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_elem_header(bigint ndump)
{
}
/*-------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_atom_header(bigint ndump)
{
}

/*-------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_node_header(bigint ndump)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpIDSAtom::pack_elem(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpIDSAtom::pack_node(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpIDSAtom::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpIDSAtom::pack_atom_noscale_noimage(tagint *ids)
{

  // pack atom

  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *pattern_atom = compute_ids->vector_atom;
  int nlocal = atom->nlocal;
  int keep_flag;
  m = n = 0;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) 
      if (allflag) {
        abuf[m++] = tag[i];
        abuf[m++] = type[i];
        abuf[m++] = x[i][0];
        abuf[m++] = x[i][1];
        abuf[m++] = x[i][2];
        abuf[m++] = pattern_atom[i];
      } else {
        for (int j = 0; j < nkeep; j++)
          if (pattern_atom[i] == keep[j]) {
            abuf[m++] = tag[i];
            abuf[m++] = type[i];
            abuf[m++] = x[i][0];
            abuf[m++] = x[i][1];
            abuf[m++] = x[i][2];
            abuf[m++] = pattern_atom[i];
            if (ids) ids[n++] = i;
            break;
          }
      }

  // pack interpolated atoms

  tagint maxtag_all = atom->maxtag();

  nlocal = element->nlocal;
  double ***nodex = element->nodex;
  int *nintpl = element->nintpl;
  mask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tag = element->tag;
  double **pattern_intpl = compute_ids->vector_intpl_atom;
  int npe = element->npe;
  double ***shape_array = element->shape_array;
  double xtmp,ytmp,ztmp;
  int ietype;

  bigint needtag,needtag_sum;
  needtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ietype = etype[i];
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++) 
        if (allflag) needtag++;
        else 
          for (int k = 0; k < nkeep; k++)
            if (pattern_intpl[i][iintpl] == keep[k]) {
              needtag++;
              break;
            }

    }

  // compress all tags for interpolated atoms

  MPI_Scan(&needtag,&needtag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my interpolated atoms should use

  bigint itag = maxtag_all + needtag_sum - needtag + 1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ietype = etype[i];
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
        if (allflag) {
          xtmp = ytmp = ztmp = 0.0;
          for (int inode = 0; inode < npe; inode++) {
            xtmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
            ytmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
            ztmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
          }
          abuf[m++] = itag;
          abuf[m++] = ctype[i];
          abuf[m++] = xtmp;
          abuf[m++] = ytmp;
          abuf[m++] = ztmp;
          abuf[m++] = pattern_intpl[i][iintpl];
        } else {
          for (int k = 0; k < nkeep; k++)
            if (pattern_intpl[i][iintpl] == keep[k]) {
              xtmp = ytmp = ztmp = 0.0;
              for (int inode = 0; inode < npe; inode++) {
                xtmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
                ytmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
                ztmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
              }
              abuf[m++] = itag;
              abuf[m++] = ctype[i];
              abuf[m++] = xtmp;
              abuf[m++] = ytmp;
              abuf[m++] = ztmp;
              abuf[m++] = pattern_intpl[i][iintpl];
              break;
            }
        }
        itag++;
      }
    }

}

/*--------------------------------------------------------------------------------------------*/

int DumpIDSAtom::convert_node_string(int n,int size_one,double *mybuf) 
{
}

/*--------------------------------------------------------------------------------------------*/
int DumpIDSAtom::convert_atom_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpIDSAtom::convert_elem_string(int n,int size_one,double *mybuf)
{
}

/*--------------------------------------------------------------------------------------------*/

int DumpIDSAtom::convert_atom_noimage(int n,int size_one,double *mybuf)
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
        static_cast<tagint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4],
        static_cast<int> (mybuf[m+5]));
    m += 6;
  }
  return offset;
}

/*----------------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_node_data(int n, double *mybuf)
{

}

/*----------------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_elem_data(int n, double *mybuf)
{
}

/*---------------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpIDSAtom::write_atom_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
        static_cast<tagint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4],
        static_cast<int> (mybuf[m+5]));
    m += 6;
  }
}


/*-------------------------------------------------------------------------------------------------*/

int DumpIDSAtom::count_atoms()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int n = 0;

  if (compute_ids->invoked_peratom != update->ntimestep)
    compute_ids->compute_peratom();

  // count atoms

  double *pattern_atom = compute_ids->vector_atom;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit)
      if (allflag) n++;
      else 
        for (int j = 0; j < nkeep; j++) 
          if (pattern_atom[i] == keep[j]) {
            n++;
            break;
          }

  // count interpolated atoms

  nlocal = element->nlocal;
  double **pattern_intpl = compute_ids->vector_intpl_atom;

  int *nintpl = element->nintpl;
  int *etype = element->etype;
  mask = element->mask;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (allflag) n += nintpl[etype[i]];
      else
        for (int k = 0; k < nintpl[etype[i]]; k++) 
          for (int j = 0; j < nkeep; j++) 
            if (pattern_intpl[i][k] == keep[j]) {
              n++;
              break;
            }

  return n;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpIDSAtom::count_elements()
{
  return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpIDSAtom::count_nodes()
{
  return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpIDSAtom::count_intpl_atoms()
{
  return 0;
}
