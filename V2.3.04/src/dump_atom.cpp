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

/*--------------------------------------------------------------------------------*/

DumpAtom::DumpAtom(CAC *cac, int narg, char **arg) : Dump(cac,narg,arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump atom command");

  nevery = universe->inumeric(FLERR,arg[3]);

  // will remove optional args after added dump custom style
  
  int iarg = 5;
  stress_flag = 0;
  force_flag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stress") == 0) {
      stress_flag = 1;
      force_flag = 0;
    } else if (strcmp(arg[iarg],"force") == 0) {
      stress_flag = 0;
      force_flag = 1;
    } else error->all(FLERR,"Illegal dump atom command");
    iarg++;
  }

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  format_default = NULL;
  format = NULL;

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

DumpAtom::~DumpAtom()
{
  delete [] id_stress;
}

/*--------------------------------------------------------------------------------*/

void DumpAtom::init_style()
{

  atom_size_one = 5;
  if (image_flag) atom_size_one += 3;
  if (stress_flag) atom_size_one += 6;
  if (force_flag) atom_size_one += 6;
  elem_size_one = atom_size_one*element->max_nintpl;
  node_size_one = 0;

  int n;

  delete [] format;
  char *str;
  if (image_flag == 0) str = (char *) TAGINT_FORMAT " %d %g %g %g";
  n = strlen(str) + 20;
  format = new char[n];
  strcpy(format,str);
  if (stress_flag) strcat(format," %g %g %g %g %g %g");
  if (force_flag) strcat(format," %g %g %g %g %g %g");
  strcat(format,"\n");

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string

  if (image_flag == 0)  {
    if (stress_flag && force_flag) 
      columns = (char *) "id type x y z vx vy vz fx fy fz sxx syy szz sxy sxz syz";
    else if (stress_flag) 
      columns = (char *) "id type x y z sxx syy szz sxy sxz syz";
    else if (force_flag)
      columns = (char *) "id type x y z vx vy vz fx fy fz";
    else
      columns = (char *) "id type x y z";
  }


  header_choice = &DumpAtom::header_item;

  pack_atom_choice = &DumpAtom::pack_atom_noscale_noimage;
  pack_elem_choice = &DumpAtom::pack_elem_noscale_noimage;

  convert_atom_choice = &DumpAtom::convert_atom_noimage;
  convert_elem_choice = &DumpAtom::convert_elem_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpAtom::write_string;
    write_elem_choice = &DumpAtom::write_string;
  }
  else {
    write_atom_choice = &DumpAtom::write_atom_lines_noimage;
    write_elem_choice = &DumpAtom::write_elem_lines_noimage;
  }

  // set stress compute pointer

  if (stress_flag) {
    int icompute = modify->find_compute(id_stress);
    if (icompute < 0) error->all(FLERR,"Stress/atom ID for dump atom do not exist");
    stress = modify->compute[icompute];
  }

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
  (this->*pack_elem_choice)(ids);
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

void DumpAtom::pack_elem_noscale_noimage(tagint *ids)
{
  double ***nodes;
  
  if (stress_flag) {
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }
    nodes = stress->array_node;
  }
  tagint maxtag,maxtag_all,*tag;
  if (atom->natoms) {
    tag = atom->tag;
    maxtag = 0;
    for (int i = 0; i < atom->nlocal; i++) 
      maxtag = MAX(maxtag,tag[i]);
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
  } else maxtag_all = 0;

  int nlocal = element->nlocal;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  double ***nodef = element->nodef;
  int *nintpl = element->nintpl;
  int *mask = element->mask;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tag = element->tag;
  int npe = element->npe;
  double ***shape_array = element->shape_array;
  double x,y,z,tmp;
  double vx,vy,vz,fx,fy,fz;
  int m = 0;
  int ietype;

  //bigint needtag,needtag_sum;
  //needtag = 0;
  //for (int i = 0; i < nlocal; i++)
  //  if (mask[i] & groupbit ) needtag += nintpl[etype[i]];

  //MPI_Scan(&needtag,&needtag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my interpolated atoms should use

  //bigint itag = maxtag_all + needtag_sum - needtag + 1;
  tagint  itag = maxtag_all + 1;
  int maxintpl = element->maxintpl;
  for (int i = 0; i < nlocal; i++) {
    ietype = etype[i];
    if (mask[i] & groupbit)
      for (int iintpl = 0; iintpl < nintpl[ietype]; iintpl++){
        ebuf[m++] = itag + (tag[i]-1)*maxintpl + iintpl;
        ebuf[m++] = ctype[i];
        x = y = z = 0.0;
        for (int inode = 0; inode < npe; inode++) {
          x += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
          y += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
          z += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
        }
        ebuf[m++] = x;
        ebuf[m++] = y;
        ebuf[m++] = z;

        if (stress_flag) {
          for (int j = 0; j < 6; j++) {
            tmp = 0.0;
            for (int inode = 0; inode < npe; inode++) 
              tmp += shape_array[ietype][iintpl][inode]*nodes[i][inode][j];
            ebuf[m++] = tmp;
          }
        }
        if (force_flag) {
          vx = vy = vz = fx = fy = fz = 0;
          for (int inode = 0; inode < npe; inode++) {
            vx += shape_array[ietype][iintpl][inode]*nodev[i][inode][0];
            vy += shape_array[ietype][iintpl][inode]*nodev[i][inode][1];
            vz += shape_array[ietype][iintpl][inode]*nodev[i][inode][2];
            fx += shape_array[ietype][iintpl][inode]*nodef[i][inode][0];
            fy += shape_array[ietype][iintpl][inode]*nodef[i][inode][1];
            fz += shape_array[ietype][iintpl][inode]*nodef[i][inode][2];
          }
          ebuf[m++] = vx;
          ebuf[m++] = vy;
          ebuf[m++] = vz;
          ebuf[m++] = fx;
          ebuf[m++] = fy;
          ebuf[m++] = fz;
        }
      }
  }
} 

/*------------------------------------------------------------------------------------------*/

void DumpAtom::pack_atom_noscale_noimage(tagint *ids)
{
  if (stress_flag) 
    if (stress->invoked_peratom != update->ntimestep) {
      stress->compute_peratom();
      stress->addstep(update->ntimestep+nevery);
    }

  int m,n;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **s;
  if (stress_flag) s = stress->array_atom;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      abuf[m++] = tag[i];
      abuf[m++] = type[i];
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

/*--------------------------------------------------------------------------------------------*/
int DumpAtom::convert_node_string(int n,int size_one,double *mybuf)
{

}

/*--------------------------------------------------------------------------------------------*/
int DumpAtom::convert_atom_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_atom_choice)(n,size_one,mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_elem_string(int n,int size_one,double *mybuf)
{
  return (this->*convert_elem_choice)(n,size_one,mybuf);
}


/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_atom_noimage(int n,int size_one,double *mybuf)
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
          static_cast<tagint> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          mybuf[m+2],mybuf[m+3],mybuf[m+4],
          mybuf[m+5],mybuf[m+6],mybuf[m+7],
          mybuf[m+8],mybuf[m+9],mybuf[m+10]);
      m += 11;
    } else {
      offset += sprintf(&sabuf[offset],format,
          static_cast<tagint> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          mybuf[m+2],mybuf[m+3],mybuf[m+4]);
      m += 5;
    }

  }
  return offset;
}

/*--------------------------------------------------------------------------------------------*/

int DumpAtom::convert_elem_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  int *etype = element->etype;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < element->nintpl[etype[i]]; j++){ // need to fix
      if (offset + ONELINE > maxsebuf) {
        if ((bigint) maxsebuf + DELTA > MAXSMALLINT) return -1;
        maxsebuf += DELTA;
        memory->grow(sebuf, maxsebuf,"dump:sebuf");
      }

      if (stress_flag || force_flag) {
        offset += sprintf(&sebuf[offset],format,
            static_cast<tagint> (mybuf[m]),
            static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4],
            mybuf[m+5],mybuf[m+6],mybuf[m+7],
            mybuf[m+8],mybuf[m+9],mybuf[m+10]);
        m += 11;
      } else { 
        offset += sprintf(&sebuf[offset],format,
            static_cast<tagint> (mybuf[m]),
            static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);

        m += 5;
      }
    }
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
  (this->*write_elem_choice) (n,mybuf);
}

/*---------------------------------------------------------------------------------------------------*/

void DumpAtom::write_string(int n, double *mybuf)
{
  int m;
  fwrite(mybuf, sizeof(char),n,fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpAtom::write_atom_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (stress_flag || force_flag) {
      fprintf(fp, format,
          static_cast<tagint> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          mybuf[m+2],mybuf[m+3],mybuf[m+4],
          mybuf[m+5],mybuf[m+6],mybuf[m+7],
          mybuf[m+8],mybuf[m+9],mybuf[m+10]);
      m += 11;
    } else {
      fprintf(fp, format,
          static_cast<tagint> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          mybuf[m+2],mybuf[m+3],mybuf[m+4]);
      m += 5;
    }
  }
}

/*-------------------------------------------------------------------------------------------------*/

void DumpAtom::write_elem_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  int *etype = element->etype;
  for (int i = 0; i < n; i++) 
    for (int j = 0; j < element->nintpl[etype[i]]; j++){
      if (stress_flag || force_flag) {
        fprintf(fp, format,
            static_cast<tagint> (mybuf[m]),
            static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4],
            mybuf[m+5],mybuf[m+6],mybuf[m+7],
            mybuf[m+8],mybuf[m+9],mybuf[m+10]);
        m += 11;
      } else {
        fprintf(fp, format,
            static_cast<tagint> (mybuf[m]),
            static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
        m += 5;
      }
    }

}

