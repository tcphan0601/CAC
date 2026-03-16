#include <string.h>
#include "dump_xyz.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
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

enum{FULLMAP, GAUSSIAN, NODE, CENTER};

/* -------------------------------------------------------------------------------- */

DumpXYZ::DumpXYZ(CAC *cac, int narg, char **arg) : Dump(cac, narg, arg)
{
  if (narg < 5) error->all(FLERR, "Illegal dump atom command");
  if (binary || multiproc) error->all(FLERR, "Invalid dump xyz filename");

  int iarg = 5;
  dump_element_style = FULLMAP;
  wrap_flag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "style") == 0) {
      if (strcmp(arg[iarg+1], "full/map") == 0)
        dump_element_style = FULLMAP;
      else if (strcmp(arg[iarg+1], "gaussian") == 0)
        dump_element_style = GAUSSIAN;
      else if (strcmp(arg[iarg+1], "node") == 0)
        dump_element_style = NODE;
      else if (strcmp(arg[iarg+1], "center") == 0)
        dump_element_style = CENTER;
      else error->all(FLERR, "Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "wrap") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "yes") == 0) wrap_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) wrap_flag = 0;
      else error->all(FLERR, "Illegal dump atom command");
      iarg += 2;
    } else error->all(FLERR, "Illegal dump atom command");
  }

  buffer_allow = 1;
  buffer_flag = 1;

  if (format_default) delete [] format_default;

  char *str = (char *) "%s %g %g %g";
  int n = strlen(str) + 2;
  format_default = new char[n];
  strcpy(format_default, str);

  ntypes = atom->ntypes;
  typenames = NULL;

}

/* -------------------------------------------------------------------------------- */

DumpXYZ::~DumpXYZ()
{
  delete[] format_default;
  format_default = NULL;

  if (typenames) {
    for (int i = 1; i <= ntypes; i++)
      delete [] typenames[i];
    delete [] typenames;
    typenames = NULL;
  }

}
/* -------------------------------------------------------------------------------- */

void DumpXYZ::init_style()
{
  atom_size_one = 4;

  if (element->nelements) {
    if (dump_element_style == FULLMAP) 
      elem_size_one = atom_size_one * element->max_nucell;
    else if (dump_element_style == GAUSSIAN) 
      elem_size_one = atom_size_one * element->max_ngcell;
    else if (dump_element_style == NODE)
      elem_size_one = atom_size_one * element->maxnpe;
  } else elem_size_one = 0;

  node_size_one = 0;
  max_size_one = MAX(atom_size_one, elem_size_one);

  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format, str);
  strcat(format, "\n");

  // initialize typenames array to be backward compatible by default
  // a 32-bit int can be maximally 10 digits plus sign

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[12];
      sprintf(typenames[itype], "%d", itype);
    }
  }

  // setup function ptr
  
  if (buffer_flag == 1) 
    write_atom_choice = &DumpXYZ::write_string;
  else 
    write_atom_choice = &DumpXYZ::write_atom_lines;

  // open single file, one time only

  if (multifile == 0) openfile();

}

/* ------------------------------------------------------------------------------------------- */

void DumpXYZ::write_header(bigint ndump)
{
  if (filewriter) {
    fprintf(fp, BIGINT_FORMAT "\n", ndump);
    fprintf(fp, "Atoms. Timestep: " BIGINT_FORMAT "\n", update->ntimestep);
  }
}

/* ------------------------------------------------------------------------------------------- */

void DumpXYZ::write_elem_header(bigint ndump)
{
}
/* ------------------------------------------------------------------------------------------- */

void DumpXYZ::write_atom_header(bigint ndump)
{
}

/* ------------------------------------------------------------------------------------------- */

void DumpXYZ::write_node_header(bigint ndump)
{
}

/* ----------------------------------------------------------------------------------------- */

void DumpXYZ::pack_node(tagint *ids)
{
}

/* ------------------------------------------------------------------------------------------ */

void DumpXYZ::pack_elem(tagint *ids)
{
} 

/* ------------------------------------------------------------------------------------------ */

void DumpXYZ::pack_atom(tagint *ids)
{
  int m, n;
  int *mask = atom->mask;
  double **x = atom->x;

  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = i;
    }
  }

  nlocal = element->nlocal;
  x = element->x;
  double ****nodex = element->nodex;
  int *nucell = element->nucell;
  int *ngcell = element->ngcell;
  mask = element->mask;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *npe = element->npe;
  int *apc = element->apc;
  double ***shape_array = element->shape_array;
  int **g2u = element->g2u;
  double tmp;
  double *coord = new double[3];
  int ietype, ibasis, iapc;

  int i, j, iucell, igcell, inode;
  for (i = 0; i < nlocal; i++) {
    ietype = etype[i];
    iapc = apc[ietype];
    if (mask[i] & groupbit) {
      if (dump_element_style == CENTER) {
        buf[m++] = 0;
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2];
      } else {
        for (ibasis = 0; ibasis < iapc; ibasis++) {
          if (dump_element_style == FULLMAP) {
            for (iucell = 0; iucell < nucell[ietype]; iucell++) {
              buf[m++] = ctype[i][ibasis];
              element->evec->interpolate(coord, nodex, i, ibasis, iucell, 3);
              if (wrap_flag) domain->remap(coord);
              buf[m++] = coord[0];
              buf[m++] = coord[1];
              buf[m++] = coord[2];
            }
          } else if (dump_element_style == GAUSSIAN) {
            for (igcell = 0; igcell < ngcell[ietype]; igcell++) {
              iucell = g2u[ietype][igcell];
              buf[m++] = ctype[i][ibasis];
              element->evec->interpolate(coord, nodex, i, ibasis, iucell, 3);
              if (wrap_flag) domain->remap(coord);
              buf[m++] = coord[0];
              buf[m++] = coord[1];
              buf[m++] = coord[2];
            }
          } else if (dump_element_style == NODE) {
            for (inode = 0; inode < npe[ietype]; inode++) {
              buf[m++] = ctype[i][ibasis];
              coord[0] = nodex[i][ibasis][inode][0];
              coord[1] = nodex[i][ibasis][inode][1];
              coord[2] = nodex[i][ibasis][inode][2];
              if (wrap_flag) domain->remap(coord);
              buf[m++] = coord[0];
              buf[m++] = coord[1];
              buf[m++] = coord[2];
            }
          }
        }
      }
    }
  }
  delete [] coord;
}

/* -------------------------------------------------------------------------------------------- */
int DumpXYZ::convert_node_string(int n, int size_one, double *mybuf)
{
  return 0;
}

/* -------------------------------------------------------------------------------------------- */

int DumpXYZ::convert_atom_string(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format, 
        typenames[static_cast<int> (mybuf[m])], 
        mybuf[m+1], mybuf[m+2], mybuf[m+3]);
    m += 4;
  }
  return offset;
}

/* -------------------------------------------------------------------------------------------- */

int DumpXYZ::convert_elem_string(int n, int size_one, double *mybuf)
{
  return 0;
}


/* ---------------------------------------------------------------------------------------------------- */

void DumpXYZ::write_node_data(int n, double *mybuf)
{

}

/* ---------------------------------------------------------------------------------------------------- */

void DumpXYZ::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n, mybuf);
}

/* ---------------------------------------------------------------------------------------------------- */

void DumpXYZ::write_elem_data(int n, double *mybuf)
{
}

/* --------------------------------------------------------------------------------------------------- */

void DumpXYZ::write_string(int n, double *mybuf)
{
  fwrite(mybuf, sizeof(char), n, fp);
}

/* ------------------------------------------------------------------------------------------------- */

void DumpXYZ::write_atom_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format, 
        typenames[static_cast<int> (mybuf[m])], 
        mybuf[m+1], mybuf[m+2], mybuf[m+3]);
    m += 4;
  }
}

/*  ----------------------------------------------------------------------  */

int DumpXYZ::count_atoms()
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
  int *nucell = element->nucell;
  int *ngcell = element->ngcell;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (dump_element_style == FULLMAP) 
        m += nucell[etype[i]] * apc[etype[i]];
      else if (dump_element_style == GAUSSIAN)
        m += ngcell[etype[i]] * apc[etype[i]];
      else if (dump_element_style == NODE)
        m += npe[etype[i]] * apc[etype[i]];
      else m++;
    }
  return m;
}

/*  ----------------------------------------------------------------------  */

int DumpXYZ::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR, "Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++)
        delete [] typenames[i];

      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 2;
      typenames[itype] = new char[n];
      strcpy(typenames[itype], arg[itype]);
    }

    return ntypes+1;
  }

  return 0;
}


