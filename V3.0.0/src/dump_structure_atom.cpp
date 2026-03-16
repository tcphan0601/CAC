#include <string.h>
#include "dump_structure_atom.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "neighbor.h"
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

enum{UNKNOWN_CNA, OTHER_CNA, FCC_CNA, HCP_CNA, BCC_CNA, ICOS_CNA};
enum{UNKNOWN_IDS, OTHER_IDS, CUBIC_IDS, HEX_IDS};
enum{IDS, CNA};

/* -------------------------------------------------------------------------------- */

DumpStructureAtom::DumpStructureAtom(CAC *cac, int narg, char **arg) : Dump(cac, narg, arg)
{

  if (narg < 7) error->all(FLERR, "Illegal dump structure/atom command");
  id_compute = arg[5];
  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR, "Could not find dump structure/atom compute ID");
  struct_compute = modify->compute[icompute];
  if (strcmp(struct_compute->style, "ids/atom") == 0) structure_style = IDS;
  else if (strcmp(struct_compute->style, "cna/atom") == 0) structure_style = CNA;
  else error->all(FLERR, "Dump structure/atom compute is not supported");

  int iarg = 6;
  nkeep = 0;
  while (iarg < narg ) {
    if (structure_style == CNA) {
      if (strcmp(arg[iarg], "unknown") == 0) keep[nkeep++] = UNKNOWN_CNA;
      else if (strcmp(arg[iarg], "fcc") == 0) keep[nkeep++] = FCC_CNA;
      else if (strcmp(arg[iarg], "hcp") == 0) keep[nkeep++] = HCP_CNA;
      else if (strcmp(arg[iarg], "bcc") == 0) keep[nkeep++] = BCC_CNA;
      else if (strcmp(arg[iarg], "icos") == 0) keep[nkeep++] = ICOS_CNA;
      else if (strcmp(arg[iarg], "other") == 0) keep[nkeep++] = OTHER_CNA;
      else if (strcmp(arg[iarg], "all") == 0) allflag = 1;
      else error->all(FLERR, "Illegal structure type to keep");
    } else {
      if (strcmp(arg[iarg], "unknown") == 0) keep[nkeep++] = UNKNOWN_IDS;
      else if (strcmp(arg[iarg], "cubic") == 0) keep[nkeep++] = CUBIC_IDS;
      else if (strcmp(arg[iarg], "hex") == 0) keep[nkeep++] = HEX_IDS;
      else if (strcmp(arg[iarg], "other") == 0) keep[nkeep++] = OTHER_IDS;
      else if (strcmp(arg[iarg], "all") == 0) allflag = 1;
      else error->all(FLERR, "Illegal structure type to keep");
    }
    iarg++;
  }
  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  clearstep = 1;
}

/* -------------------------------------------------------------------------------- */

void DumpStructureAtom::init_style() 
{
  if (image_flag == 0) atom_size_one = 9;
  else atom_size_one = 12;

  max_size_one = atom_size_one;
  node_size_one = elem_size_one =  0;

  int n;

  delete [] format;
  char *str;
  if (image_flag == 0) str = (char *) TAGINT_FORMAT " %d %g %g %g %d %d %d %d";
  else str = (char *) TAGINT_FORMAT " %d %g %g %g %d %d %d %d %d %d %d";
  n = strlen(str) + 2;
  format = new char[n];
  strcpy(format, str);
  strcat(format, "\n");

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string

  delete [] columns;
  columns = new char[500];
  if (image_flag == 0) 
    columns = strcpy(columns, "id type x y z structure_id eid basis ucell");
  else
    columns = strcpy(columns, "id type x y z structure_id eid basis ucell ix iy iz");

  if (domain->triclinic == 0)
    header_choice = &DumpStructureAtom::header_item;
  else
    header_choice = &DumpStructureAtom::header_item_triclinic;

  pack_atom_choice = &DumpStructureAtom::pack_atom_noscale_noimage;

  convert_atom_choice = &DumpStructureAtom::convert_atom_noimage;

  if (buffer_flag == 1) {
    write_atom_choice = &DumpStructureAtom::write_string;
  }
  else {
    write_atom_choice = &DumpStructureAtom::write_atom_lines_noimage;
  }

  // open single file, one time only

  if (multifile == 0) openfile();

}
/* ------------------------------------------------------------------------------------------ */

void DumpStructureAtom::header_item(bigint ndump) 
{
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, BIGINT_FORMAT "\n", ndump);
  fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
  fprintf(fp, "%-1.16e %-1.16e\n", boxxlo, boxxhi);
  fprintf(fp, "%-1.16e %-1.16e\n", boxylo, boxyhi);
  fprintf(fp, "%-1.16e %-1.16e\n", boxzlo, boxzhi);
  fprintf(fp, "ITEM: ATOMS %s\n", columns);
}

/* ------------------------------------------------------------------------------------------ */

void DumpStructureAtom::header_item_triclinic(bigint ndump)
{
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, BIGINT_FORMAT "\n", ndump);
  fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxxlo, boxxhi, boxxy);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxylo, boxyhi, boxxz);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxzlo, boxzhi, boxyz);
  fprintf(fp, "ITEM: ATOMS %s\n", columns);
}


/* ------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/* ------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_elem_header(bigint ndump)
{
}
/* ------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_atom_header(bigint ndump)
{
}

/* ------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_node_header(bigint ndump)
{
}

/* ----------------------------------------------------------------------------------------- */

void DumpStructureAtom::pack_elem(tagint *ids)
{
}

/* ----------------------------------------------------------------------------------------- */

void DumpStructureAtom::pack_node(tagint *ids)
{
}

/* ----------------------------------------------------------------------------------------- */

void DumpStructureAtom::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/* ------------------------------------------------------------------------------------------ */

void DumpStructureAtom::pack_atom_noscale_noimage(tagint *ids)
{

  // pack atom

  int m, n;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *pattern_atom = struct_compute->vector_atom;
  int nlocal = atom->nlocal;
  m = n = 0;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) {
      if (allflag) {
        buf[m++] = tag[i];
        buf[m++] = type[i];
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2];
        buf[m++] = pattern_atom[i];
        buf[m++] = -1;
        buf[m++] = -1;
        buf[m++] = -1;
        if (ids) ids[n++] = i;
      } else 
        for (int j = 0; j < nkeep; j++)
          if (pattern_atom[i] == keep[j]) {
            buf[m++] = tag[i];
            buf[m++] = type[i];
            buf[m++] = x[i][0];
            buf[m++] = x[i][1];
            buf[m++] = x[i][2];
            buf[m++] = pattern_atom[i];
            buf[m++] = -1;
            buf[m++] = -1;
            buf[m++] = -1;
            if (ids) ids[n++] = i;
            break;
          }
    }
  // pack virtual atoms

  tagint maxtag_all = atom->maxtag();

  nlocal = element->nlocal;
  double ****nodex = element->nodex;
  int *nucell = element->nucell;
  int **is_outer = element->is_outer;
  mask = element->mask;
  int *etype = element->etype;
  int **ctype = element->ctype;
  tag = element->tag;
  double ***pattern_vatom = struct_compute->vector_vatom;
  int *npe = element->npe;
  int *apc = element->apc;
  double ***shape_array = element->shape_array;
  double xtmp, ytmp, ztmp;
  int ietype, iouter;

  // first loop: scan for number of new tags needed

  bigint needtag, needtag_sum;
  needtag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ietype = etype[i];
      for (int iucell = 0; iucell < nucell[ietype]; iucell++) {
        iouter = is_outer[ietype][iucell]-1;
        if (iouter < 0) continue;
        for (int ibasis = 0; ibasis < apc[ietype]; ibasis++) {
          if (allflag) needtag++;
          else 
            for (int k = 0; k < nkeep; k++) {
              if (pattern_vatom[i][ibasis][iucell] == keep[k]) {
                needtag++;
                break;
              }
            }
        }
      }
    }

  MPI_Scan(&needtag, &needtag_sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);

  // itag = 1st new tag that my virtual atoms should use

  bigint itag = maxtag_all + needtag_sum - needtag + 1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ietype = etype[i];
      for (int iucell = 0; iucell < nucell[ietype]; iucell++) {
        if (is_outer[ietype][iucell])
          for (int ibasis = 0; ibasis < apc[ietype]; ibasis++) {
            if (allflag) {
              buf[m++] = itag++;
              buf[m++] = ctype[i][ibasis];
              buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 0);
              buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 1);
              buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 2);
              buf[m++] = pattern_vatom[i][ibasis][iucell];
              buf[m++] = tag[i];
              buf[m++] = ibasis;
              buf[m++] = iucell;
            } else 
              for (int k = 0; k < nkeep; k++) {
                if (pattern_vatom[i][ibasis][iucell] == keep[k]) {
                  buf[m++] = itag++;
                  buf[m++] = ctype[i][ibasis];
                  buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 0);
                  buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 1);
                  buf[m++] = element->evec->interpolate(nodex, i, ibasis, iucell, 2);
                  buf[m++] = pattern_vatom[i][ibasis][iucell];
                  buf[m++] = tag[i];
                  buf[m++] = ibasis;
                  buf[m++] = iucell;
                  break;
                }
              }
          }
      }
    }
}

/* -------------------------------------------------------------------------------------------- */

int DumpStructureAtom::convert_node_string(int n, int size_one, double *mybuf) 
{
  return 0;
}

/* -------------------------------------------------------------------------------------------- */
int DumpStructureAtom::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n, size_one, mybuf);
}

/* -------------------------------------------------------------------------------------------- */

int DumpStructureAtom::convert_elem_string(int n, int size_one, double *mybuf)
{
  return 0;
}

/* -------------------------------------------------------------------------------------------- */

int DumpStructureAtom::convert_atom_noimage(int n, int size_one, double *mybuf)
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
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4], 
        static_cast<int> (mybuf[m+5]), 
        static_cast<int> (mybuf[m+6]), 
        static_cast<int> (mybuf[m+7]), 
        static_cast<int> (mybuf[m+8]));
    m += 9;
  }
  return offset;
}

/* ---------------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_node_data(int n, double *mybuf)
{

}

/* ---------------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n, mybuf);
}

/* ---------------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_elem_data(int n, double *mybuf)
{
}

/* --------------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_string(int n, double *mybuf)
{
  fwrite(mybuf, sizeof(char), n, fp);
}
/* ------------------------------------------------------------------------------------------------- */

void DumpStructureAtom::write_atom_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4], 
        static_cast<int> (mybuf[m+5]), 
        static_cast<int> (mybuf[m+6]), 
        static_cast<int> (mybuf[m+7]), 
        static_cast<int> (mybuf[m+8]));
    m += 9;
  }
}


/* ------------------------------------------------------------------------------------------------- */

int DumpStructureAtom::count_atoms()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int n = 0;

  // if already invoked on this time step, 
  // check if it's invoked pre neighbor and 
  // neigh list is rebuild on this time step
  // set preneighflag = 0 before invoke

  if (struct_compute->invoked_peratom != update->ntimestep ||
      (struct_compute->invoked_peratom == update->ntimestep && struct_compute->preneighflag 
       && neighbor->lastcall == update->ntimestep)) {
    struct_compute->compute_peratom();
    struct_compute->preneighflag = 0;
  }

  // count atoms

  double *pattern_atom = struct_compute->vector_atom;
  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) {
      if (allflag) n++;
      else
        for (int j = 0; j < nkeep; j++)
          if (pattern_atom[i] == keep[j]) {
            n++;
            break;
          }
    }

  // count virtual atoms

  nlocal = element->nlocal;
  double ***pattern_vatom = struct_compute->vector_vatom;

  int *nucell = element->nucell;
  int **is_outer = element->is_outer;
  int *apc = element->apc;
  int *etype = element->etype;
  mask = element->mask;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      for (int k = 0; k < nucell[etype[i]]; k++) {
        if (is_outer[etype[i]][k])
          for (int ibasis = 0; ibasis < apc[etype[i]]; ibasis++) {
            if (allflag) n++;
            else 
              for (int j = 0; j < nkeep; j++) {
                if (pattern_vatom[i][ibasis][k] == keep[j]) {
                  n++;
                  break;
                }
              }
          }
      }

  return n;
}

/* ------------------------------------------------------------------------------------------------- */

int DumpStructureAtom::count_elements()
{
  return 0;
}

/* ------------------------------------------------------------------------------------------------- */

int DumpStructureAtom::count_nodes()
{
  return 0;
}


