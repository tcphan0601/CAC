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
  format_default = NULL;
  format = NULL;
  clearstep = 1;
}

/*--------------------------------------------------------------------------------*/

void DumpCNAAtom::init_style(){
  if (image_flag == 0) {
    atom_size_one = 5;
    elem_size_one = 5;
  } else { 
    atom_size_one = 8;
    elem_size_one = 8;
  }
  node_size_one = 0;

  int n;

  delete [] format;
  char *str;
  if (image_flag == 0) str = (char *) TAGINT_FORMAT " %d %g %g %g";
  n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string
  if (image_flag == 0) 
    columns = (char *) "id type x y z";
  else
    columns = (char *) "id type x y z ix iy iz";

  header_choice = &DumpCNAAtom::header_item;

  pack_atom_choice = &DumpCNAAtom::pack_atom_noscale_noimage;
  pack_elem_choice = &DumpCNAAtom::pack_elem_noscale_noimage;

  convert_atom_choice = &DumpCNAAtom::convert_atom_noimage;
  convert_elem_choice = &DumpCNAAtom::convert_elem_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpCNAAtom::write_string;
    write_elem_choice = &DumpCNAAtom::write_string;
  }
  else {
    write_atom_choice = &DumpCNAAtom::write_atom_lines_noimage;
    write_elem_choice = &DumpCNAAtom::write_elem_lines_noimage;
  }

}
/*------------------------------------------------------------------------------------------*/

void DumpCNAAtom::header_item(bigint ndump){

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
  (this->*pack_elem_choice)(ids);
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

void DumpCNAAtom::pack_elem_noscale_noimage(tagint *ids)
{
  tagint maxtag,maxtag_all;
  if (atom->natoms) {
    int *tag = atom->tag;
    for (int i = 0; i < atom->nlocal; i++) 
      maxtag = MAX(maxtag,tag[i]);
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
  } else maxtag_all = 0;
 
  int nlocal = element->nlocal;
  double ***nodex = element->nodex;
  int *nintpl = element->nintpl;
  int *mask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *tag = element->tag;
  double **pattern = compute_cna->vector_intpl_atom;
  int npe = element->npe;
  double ***shape_array = element->shape_array;
  double x,y,z;
  int m = 0;
  int keep_flag,ietype;

  bigint needtag,needtag_sum;
  needtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
    ietype = etype[i];
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++){
        keep_flag = 0;
        for (int k = 0; k < nkeep; k++)
          if (pattern[i][iintpl] == keep[k]) keep_flag = 1;
        if (keep_flag) needtag++;
      }
    }


   MPI_Scan(&needtag,&needtag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my interpolated atoms should use
 
  bigint itag = maxtag_all + needtag_sum - needtag + 1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
    ietype = etype[i];
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++){
        keep_flag = 0;
        for (int k = 0; k < nkeep; k++)
          if (pattern[i][iintpl] == keep[k]) keep_flag = 1;
        if (keep_flag) {
          x = y = z = 0.0;
          for (int inode = 0; inode < npe; inode++) {
            x += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
            y += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
            z += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
          }
          ebuf[m++] = itag;
          ebuf[m++] = ctype[i];
          ebuf[m++] = x;
          ebuf[m++] = y;
          ebuf[m++] = z;
        }
        itag++;
      }
    }
} 

/*------------------------------------------------------------------------------------------*/

void DumpCNAAtom::pack_atom_noscale_noimage(tagint *ids){
  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *pattern = compute_cna->vector_atom;
  int nlocal = atom->nlocal;
  int keep_flag;
  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    keep_flag = 0;
    if (mask[i] & groupbit) {
      for (int j = 0; j < nkeep; j++)
        if (pattern[i] == keep[j]) keep_flag = 1;
      if (keep_flag){
        abuf[m++] = tag[i];
        abuf[m++] = type[i];
        abuf[m++] = x[i][0];
        abuf[m++] = x[i][1];
        abuf[m++] = x[i][2];
        if (ids) ids[n++] = i;
      }

    }
  }
}

/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_node_string(int n,int size_one,double *mybuf){}

/*--------------------------------------------------------------------------------------------*/
int DumpCNAAtom::convert_atom_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_elem_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_elem_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_atom_noimage(int n,int size_one,double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i<n; i++) {
    if (offset + ONELINE > maxsabuf) {
      if ((bigint) maxsabuf + DELTA > MAXSMALLINT) return -1;
      maxsabuf += DELTA;
      memory->grow(sabuf, maxsabuf,"dump:sabuf");
    }
    offset += sprintf(&sabuf[offset],format,
        static_cast<tagint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += 5;
  }
  return offset;
}

/*--------------------------------------------------------------------------------------------*/

int DumpCNAAtom::convert_elem_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsebuf) {
      if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
      maxsebuf += DELTA;
      memory->grow(sebuf, maxsebuf,"dump:sebuf");
    }
    offset += sprintf(&sebuf[offset],format,
        static_cast<bigint> (mybuf[m]),
        static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += 5;

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
  (this->*write_elem_choice) (n,mybuf);
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
  int keep_flag;
  for (int i = 0; i < n; i++){
    keep_flag = 0;
    for (int j = 0; j < nkeep; j++) 
      if (static_cast<int> (mybuf[m]) == keep[j]) keep_flag = 1;

    if (keep_flag) fprintf(fp, format, mybuf[m+1], mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
    m += 6;

  }
}

/*-------------------------------------------------------------------------------------------------*/

void DumpCNAAtom::write_elem_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  int keep_flag;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < element->nintpl[1]; j++){
      keep_flag = 0;
      for (int k = 0; k < nkeep; k++) 
        if (static_cast<int> (mybuf[m]) == keep[k]) keep_flag = 1;
      if (keep_flag) fprintf(fp,format,mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
  }
}

/*-------------------------------------------------------------------------------------------------*/

int DumpCNAAtom::count_cna_atoms(){

  int nlocal = atom->nlocal;
  if (nlocal == 0) return 0;
  if (compute_cna->invoked_peratom != update->ntimestep)
    compute_cna->compute_peratom();
  int n = 0;
  double *pattern = compute_cna->vector_atom;

  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < nkeep; j++)
      if (pattern[i] == keep[j]) n++;
  return n;
}

/*-------------------------------------------------------------------------------------------------*/

int DumpCNAAtom::count_cna_intpl(){
  int nlocal = element->nlocal;
  if (nlocal == 0) return 0;
  if (compute_cna->invoked_peratom != update->ntimestep) 
    compute_cna->compute_peratom();
  int n = 0;
  double **pattern = compute_cna->vector_intpl_atom;

  int *nintpl = element->nintpl;
  int *etype = element->etype;
  int keep_flag;

  for (int i = 0; i < nlocal; i++)
    for (int k = 0; k < nintpl[etype[i]]; k++) {
      keep_flag = 0;
      for (int j = 0; j < nkeep; j++) {
        if (pattern[i][k] == keep[j]) keep_flag = 1;
        if (keep_flag) n++;
      }
    }
  return n;
}

