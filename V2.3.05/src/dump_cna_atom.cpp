#include <string.h>
#include "dump_cna_atom.h"
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

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};

/*--------------------------------------------------------------------------------*/

DumpCNAAtom::DumpCNAAtom(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{

  if (narg < 6 || narg > 11) error->all(FLERR,"Illegal dump cna command");
  int compute_not_exist = 1;

  for (int i = 0; i < modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"cna/atom") == 0)
      compute_cna = modify->compute[i];
    compute_not_exist = 0;
  }

  if (compute_not_exist) error->all(FLERR,"dump cna/atom requires compute cna/atom command to be defined");

  int iarg = 5;
  nkeep = 0;
  while (iarg < narg ) {
    if (strcmp(arg[iarg],"unknown") == 0) keep[nkeep++] = UNKNOWN;
    else if (strcmp(arg[iarg],"fcc") == 0) keep[nkeep++] = FCC;
    else if (strcmp(arg[iarg],"hcp") == 0) keep[nkeep++] = HCP;
    else if (strcmp(arg[iarg],"bcc") == 0) keep[nkeep++] = BCC;
    else if (strcmp(arg[iarg],"icos") == 0) keep[nkeep++] = ICOS;
    else if (strcmp(arg[iarg],"other") == 0) keep[nkeep++] = OTHER;
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

void DumpCNAAtom::init_style() 
{
  if (image_flag == 0) {
    atom_size_one = 6;
    elem_size_one = 6;
  } else { 
    atom_size_one = 9;
    elem_size_one = 9;
  }
  node_size_one = 0;

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
  if (image_flag == 0) 
    columns = (char *) "id type x y z structure_id";
  else
    columns = (char *) "id type x y z structure_id ix iy iz";

  header_choice = &DumpCNAAtom::header_item;

  pack_atom_choice = &DumpCNAAtom::pack_atom_noscale_noimage;

  convert_atom_choice = &DumpCNAAtom::convert_atom_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpCNAAtom::write_string;
  }
  else {
    write_atom_choice = &DumpCNAAtom::write_atom_lines_noimage;
  }

}
/*------------------------------------------------------------------------------------------*/

void DumpCNAAtom::header_item(bigint ndump) 
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

void DumpCNAAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_elem_header(bigint ndump)
{
}
/*-------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_atom_header(bigint ndump)
{
}

/*-------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_node_header(bigint ndump)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpCNAAtom::pack_elem(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpCNAAtom::pack_node(tagint *ids)
{
}

/*-----------------------------------------------------------------------------------------*/

void DumpCNAAtom::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/*------------------------------------------------------------------------------------------*/

void DumpCNAAtom::pack_atom_noscale_noimage(tagint *ids)
{

  // pack atom

  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *pattern_atom = compute_cna->vector_atom;
  int nlocal = atom->nlocal;
  m = n = 0;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) 
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

  // pack interpolated atoms

  tagint maxtag_all = atom->maxtag();

  nlocal = element->nlocal;
  double ***nodex = element->nodex;
  int *nintpl = element->nintpl;
  mask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tag = element->tag;
  double **pattern_intpl = compute_cna->vector_intpl_atom;
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
        for (int k = 0; k < nkeep; k++)
          if (pattern_intpl[i][iintpl] == keep[k]) {
            needtag++;
            break;
          }

    }

  MPI_Scan(&needtag,&needtag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my interpolated atoms should use

  bigint itag = maxtag_all + needtag_sum - needtag + 1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ietype = etype[i];
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
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
        itag++;
      }
    }
}

/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_node_string(int n,int size_one,double *mybuf) {
}

/*--------------------------------------------------------------------------------------------*/
int DumpCNAAtom::convert_atom_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_elem_string(int n,int size_one,double *mybuf)
{
}


/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_atom_noimage(int n,int size_one,double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsabuf) {
      if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
      maxsabuf += DELTA;
      memory->grow(sabuf, maxsabuf,"dump:sabuf");
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

void DumpCNAAtom::write_node_data(int n, double *mybuf)
{

}

/*----------------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n,mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_elem_data(int n, double *mybuf)
{
}

/*---------------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_atom_lines_noimage(int n, double *mybuf)
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

int DumpCNAAtom::count_atoms()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int n = 0;

  if (compute_cna->invoked_peratom != update->ntimestep)
    compute_cna->compute_peratom();

  // count atoms

  double *pattern_atom = compute_cna->vector_atom;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit)
      for (int j = 0; j < nkeep; j++)
        if (pattern_atom[i] == keep[j]) {
          n++;
          break;
        }

  // count interpolated atoms

  nlocal = element->nlocal;
  double **pattern_intpl = compute_cna->vector_intpl_atom;

  int *nintpl = element->nintpl;
  int *etype = element->etype;
  mask = element->mask;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int k = 0; k < nintpl[etype[i]]; k++) 
        for (int j = 0; j < nkeep; j++) 
          if (pattern_intpl[i][k] == keep[j]) {
            n++;
            break;
          }

  return n;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpCNAAtom::count_elements()
{
  return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpCNAAtom::count_nodes()
{
  return 0;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpCNAAtom::count_intpl_atoms()
{
  return 0;
}
