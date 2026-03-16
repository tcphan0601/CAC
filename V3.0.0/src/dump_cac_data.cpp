#include <string.h>
#include "dump_cac_data.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "fix.h"
#include "fix_store.h"
#include "force.h"
#include "pair.h"
#include "error.h"
#include "comm.h"
#include "compute.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA 1048576

enum{INTERPOLATE, INTEGRATION, NODE, CENTER};

/* -------------------------------------------------------------------------------- */

DumpCACData::DumpCACData(CAC *cac, int narg, char **arg) : Dump(cac, narg, arg)
{
  if (narg < 5) error->all(FLERR, "Illegal dump cac/data command");
  nevery = universe->inumeric(FLERR, arg[3]);

  create_attribute = 1;

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  
  reference_flag = 0;
  precision_flag = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "reference") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "yes") == 0) reference_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) reference_flag = 0;
      else error->all(FLERR, "Illegal dump atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "precision") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "high") == 0) precision_flag = 1;
      else if (strcmp(arg[iarg+1], "normal") == 0) precision_flag = 0;
      else error->all(FLERR, "Illegal dump atom command");
      iarg += 2;
    } else error->all(FLERR, "Illegal dump atom command");
  }

  if (reference_flag && precision_flag)
    error->warning(FLERR, "Precision flag is disabled when reference flag is on");


  if (reference_flag) {
    int n = strlen(id) + strlen("_DUMP_STORE") + 1;
    id_fix = new char[n];
    strcpy(id_fix, id);
    strcat(id_fix, "_DUMP_STORE");

    char **newarg = new char*[6];
    newarg[0] = id_fix;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "node";
    newarg[5] = (char *) "1";
    newarg[6] = (char *) "3";
    modify->add_fix(7, newarg);
    fix = (FixStore *) modify->fix[modify->nfix-1];

    // calculate xu, yu, zu for fix store array

    double **xoriginal = fix->astore;
    double **x = atom->x;
    int *amask = atom->mask;
    imageint *aimage = atom->image;
    int nalocal = atom->nlocal;

    double ****nodexoriginal = fix->a4store;
    double ****nodex = element->nodex;
    int *emask = element->mask;
    int *etype = element->etype;
    imageint *eimage = element->image;
    int nelocal = element->nlocal;
    int *npe = element->npe;
    int *apc = element->apc;

    for (int i = 0; i < nalocal; i++) {
      if (amask[i] & groupbit) { 
        domain->unmap(x[i], aimage[i], xoriginal[i]);
      } else { 
        xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
      }
    }

    for (int i = 0; i < nelocal; i++) {
      if (emask[i] & groupbit) {
        for (int j = 0; j < apc[etype[i]]; j++)
          for (int k = 0; k < npe[etype[i]]; k++)
            domain->unmap(nodex[i][j][k], eimage[i], nodexoriginal[i][j][k]);
      } else {
        for (int j = 0; j < apc[etype[i]]; j++)
          for (int k = 0; k < npe[etype[i]]; k++) {
            nodexoriginal[i][j][k][0] = 0.0;
            nodexoriginal[i][j][k][1] = 0.0;
            nodexoriginal[i][j][k][2] = 0.0;
          }
      }
    }
  }

  // store simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo_ref = domain->boxlo[0];
    boxxhi_ref = domain->boxhi[0];
    boxylo_ref = domain->boxlo[1];
    boxyhi_ref = domain->boxhi[1];
    boxzlo_ref = domain->boxlo[2];
    boxzhi_ref = domain->boxhi[2];
  } else {
    boxxlo_ref = domain->boxlo_bound[0];
    boxxhi_ref = domain->boxhi_bound[0];
    boxylo_ref = domain->boxlo_bound[1];
    boxyhi_ref = domain->boxhi_bound[1];
    boxzlo_ref = domain->boxlo_bound[2];
    boxzhi_ref = domain->boxhi_bound[2];
    boxxy_ref = domain->xy;
    boxxz_ref = domain->xz;
    boxyz_ref = domain->yz;
  }
}

/* -------------------------------------------------------------------------------- */

DumpCACData::~DumpCACData()
{
  // check nfix in case all fixes have already been deleted

  if (reference_flag) {
    if (modify->nfix) 
      modify->delete_fix(id_fix);
    delete [] id_fix;
  }
}

/* -------------------------------------------------------------------------------- */

void DumpCACData::init_style()
{
  atom_size_one = 5;
  //if (image_flag) atom_size_one += 3;
  elem_size_one = 5;
  //if (image_flag) elem_size_one += 3;
  node_size_one = 7;
  max_size_one = 6;
  if (reference_flag) {
    atom_size_one += 3;
    node_size_one += 3;
    max_size_one += 3;
    if (image_flag) error->all(FLERR, "Cannot dump image and reference at the same time for dump cac/data");
  }

  int n;

  if (reference_flag) {
    delete [] format_atom;
    char *str;
    str = (char *) TAGINT_FORMAT " %d %g %g %g %g %g %g\n";
    n = strlen(str) + 2;
    format_atom = new char[n];
    strcpy(format_atom, str);

    delete [] format_elem;
    str = (char *) TAGINT_FORMAT " %d %g %g %g\n";
    n = strlen(str) + 2;
    format_elem = new char[n];
    strcpy(format_elem, str);

    delete [] format_node;
    str = (char *) TAGINT_FORMAT " %d %d %d %g %g %g %g %g %g\n";
    n = strlen(str) + 2;
    format_node = new char[n];
    strcpy(format_node, str);
  } else {
    if (precision_flag) {
      delete [] format_atom;
      char *str;
      str = (char *) TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_atom = new char[n];
      strcpy(format_atom, str);

      delete [] format_elem;
      str = (char *) TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_elem = new char[n];
      strcpy(format_elem, str);

      delete [] format_node;
      str = (char *) TAGINT_FORMAT " %d %d %d %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_node = new char[n];
      strcpy(format_node, str);
    } else {
      delete [] format_atom;
      char *str;
      str = (char *) TAGINT_FORMAT " %d %g %g %g\n";
      n = strlen(str) + 2;
      format_atom = new char[n];
      strcpy(format_atom, str);

      delete [] format_elem;
      str = (char *) TAGINT_FORMAT " %d %g %g %g\n";
      n = strlen(str) + 2;
      format_elem = new char[n];
      strcpy(format_elem, str);

      delete [] format_node;
      str = (char *) TAGINT_FORMAT " %d %d %d %g %g %g\n";
      n = strlen(str) + 2;
      format_node = new char[n];
      strcpy(format_node, str);
    }
  }
  // setup boundary string

  domain->boundary_string(boundstr);

  // switch node and element section so that element section is written first

  if (reference_flag) {
    if (domain->triclinic == 0)
      header_choice = &DumpCACData::header_item_reference;
    else
      header_choice = &DumpCACData::header_item_reference_triclinic;

    pack_atom_choice = &DumpCACData::pack_atom_noscale_reference;
    pack_elem_choice = &DumpCACData::pack_node_noscale_reference;
    pack_node_choice = &DumpCACData::pack_elem_noscale_reference;
    convert_atom_choice = &DumpCACData::convert_atom_reference;
    convert_elem_choice = &DumpCACData::convert_node_reference;
    convert_node_choice = &DumpCACData::convert_elem_reference;
  } else {
    if (domain->triclinic == 0)
      header_choice = &DumpCACData::header_item;
    else
      header_choice = &DumpCACData::header_item_triclinic;

    pack_atom_choice = &DumpCACData::pack_atom_noscale_noimage;
    pack_elem_choice = &DumpCACData::pack_node_noscale_noimage;
    pack_node_choice = &DumpCACData::pack_elem_noscale_noimage;
    convert_atom_choice = &DumpCACData::convert_atom_noimage;
    convert_elem_choice = &DumpCACData::convert_node_noimage;
    convert_node_choice = &DumpCACData::convert_elem_noimage;
  }

  if (buffer_flag == 1) {
    write_atom_choice = &DumpCACData::write_string;
    write_elem_choice = &DumpCACData::write_string;
    write_node_choice = &DumpCACData::write_string;
  } else {
    if (reference_flag) {
      write_atom_choice = &DumpCACData::write_atom_lines_reference;
      write_elem_choice = &DumpCACData::write_node_lines_reference;
      write_node_choice = &DumpCACData::write_elem_lines_reference;
    } else {
      write_atom_choice = &DumpCACData::write_atom_lines_noimage;
      write_elem_choice = &DumpCACData::write_node_lines_noimage;
      write_node_choice = &DumpCACData::write_elem_lines_noimage;
    }
  }

  if (!multifile) 
    error->all(FLERR, "Dump cac/data requires separate dump files");
}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::header_item(bigint ndump)
{
  fprintf(fp, "CAC data file from dump cac/data command\n");
  fprintf(fp, "# Timestep " BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "# Boundary Conditions: %s\n\n", boundstr);

  if (natom_total) 
    fprintf(fp, BIGINT_FORMAT " atoms\n", natom_total);
  if (nnode_total) 
    fprintf(fp, BIGINT_FORMAT " elements\n", nnode_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xlo xhi\n", boxxlo, boxxhi);
  fprintf(fp, "%-1.16e %-1.16e ylo yhi\n", boxylo, boxyhi);
  fprintf(fp, "%-1.16e %-1.16e zlo zhi\n", boxzlo, boxzhi);

  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp, "\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g\n", i, mass[i]);
  }

  // write out Element Types section if there exists elements in dump

  if (nelem_total) {
    fprintf(fp, "\nElement Types\n\n");
    element->evec->write_element_types(fp);
  }

}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::header_item_triclinic(bigint ndump)
{
  fprintf(fp, "CAC data file from dump cac/data command\n");
  fprintf(fp, "# Timestep " BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "# Boundary Conditions: %s\n\n", boundstr);

  if (natom_total) 
    fprintf(fp, BIGINT_FORMAT " atoms\n", natom_total);
  if (nnode_total) 
    fprintf(fp, BIGINT_FORMAT " elements\n", nnode_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xlo xhi\n", boxxlo, boxxhi);
  fprintf(fp, "%-1.16e %-1.16e ylo yhi\n", boxylo, boxyhi);
  fprintf(fp, "%-1.16e %-1.16e zlo zhi\n", boxzlo, boxzhi);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e xy xz yz\n", 
      boxxy, boxxz, boxyz);

  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp, "\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g\n", i, mass[i]);
  }

  // write out Element Types section if there exists elements in dump

  if (nelem_total) {
    fprintf(fp, "\nElement Types\n\n");
    element->evec->write_element_types(fp);
  }

}


/* ------------------------------------------------------------------------------------------ */

void DumpCACData::header_item_reference(bigint ndump)
{
  fprintf(fp, "CAC data file from dump cac/data command\n");
  fprintf(fp, "# Atom/Element positions are in reference configuration with current configuration for strain calculation\n");
  fprintf(fp, "# Timestep " BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "# Boundary Conditions: %s\n\n", boundstr);

  if (natom_total) 
    fprintf(fp, BIGINT_FORMAT " atoms\n", natom_total);
  if (nelem_total) 
    fprintf(fp, BIGINT_FORMAT " elements\n", nelem_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xloref xhiref\n", boxxlo_ref, boxxhi_ref);
  fprintf(fp, "%-2.16e %-1.16e yloref yhiref\n", boxylo_ref, boxyhi_ref);
  fprintf(fp, "%-1.16e %-1.16e zloref zhiref\n", boxzlo_ref, boxzhi_ref);
  fprintf(fp, "%-1.16e %-1.16e xlocur xhicur\n", boxxlo, boxxhi);
  fprintf(fp, "%-2.16e %-1.16e ylocur yhicur\n", boxylo, boxyhi);
  fprintf(fp, "%-1.16e %-1.16e zlocur zhicur\n", boxzlo, boxzhi);


  // write out Element Types section if there exists elements in dump

  if (nelem_total) {
    fprintf(fp, "\nElement Types\n\n");
    element->evec->write_element_types(fp);
  }

}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::header_item_reference_triclinic(bigint ndump)
{
  fprintf(fp, "CAC data file from dump cac/data command\n");
  fprintf(fp, "# Atom/Element positions are in reference configuration with current configuration for strain calculation\n");
  fprintf(fp, "# Timestep " BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "# Boundary Conditions: %s\n\n", boundstr);

  if (natom_total) 
    fprintf(fp, BIGINT_FORMAT " atoms\n", natom_total);
  if (nelem_total) 
    fprintf(fp, BIGINT_FORMAT " elements\n", nelem_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xloref xhiref\n", boxxlo_ref, boxxhi_ref);
  fprintf(fp, "%-1.16e %-1.16e yloref yhiref\n", boxylo_ref, boxyhi_ref);
  fprintf(fp, "%-1.16e %-1.16e zloref zhiref\n", boxzlo_ref, boxzhi_ref);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e xyref xzref yzref\n", 
      boxxy_ref, boxxz_ref, boxyz_ref);
  fprintf(fp, "%-1.16e %-1.16e xlocur xhicur\n", boxxlo, boxxhi);
  fprintf(fp, "%-1.16e %-1.16e ylocur yhicur\n", boxylo, boxyhi);
  fprintf(fp, "%-1.16e %-1.16e zlocur zhicur\n", boxzlo, boxzhi);
  fprintf(fp, "%-1.16e %-1.16e %-1.16e xycur xzcur yzcur\n", 
      boxxy, boxxz, boxyz);

  // write out Element Types section if there exists elements in dump

  if (nelem_total) {
    fprintf(fp, "\nElement Types\n\n");
    element->evec->write_element_types(fp);
  }

}

/* ------------------------------------------------------------------------------------------- */

void DumpCACData::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/* ------------------------------------------------------------------------------------------- */

void DumpCACData::write_elem_header(bigint ndump)
{
  fprintf(fp, "\nNodes\n\n");
}
/* ------------------------------------------------------------------------------------------- */

void DumpCACData::write_atom_header(bigint ndump)
{
  fprintf(fp, "\nAtoms\n\n");
}

/* ------------------------------------------------------------------------------------------- */

void DumpCACData::write_node_header(bigint ndump)
{
  fprintf(fp, "\nElements\n\n");
}

/* ----------------------------------------------------------------------------------------- */

void DumpCACData::pack_elem(tagint *ids)
{
  (this->*pack_elem_choice)(ids);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCACData::pack_node(tagint *ids)
{
  (this->*pack_node_choice)(ids);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCACData::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_atom_noscale_noimage(tagint *ids)
{

  int m, n;
  int i;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = i;
    }
  }

}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_elem_noscale_noimage(tagint *ids)
{

  int m, n;
  int i;
  int *mask = element->mask;
  double **x = element->x;
  int **ctype = element->ctype;
  int *etype = element->etype;
  tagint *tag = element->tag;
  int nlocal = element->nlocal;

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = etype[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = i;
    }
  }

}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_node_noscale_noimage(tagint *ids)
{

  int m, n;
  int i, j, k;
  int *mask = element->mask;
  double ****nodex = element->nodex;
  int **ctype = element->ctype;
  int *etype = element->etype;
  tagint *tag = element->tag;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;

  m = n = 0;
  for (i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) 
      for (j = 0; j < apc[etype[i]]; j++) 
        for (k = 0; k < npe[etype[i]]; k++) {
          buf[m++] = tag[i];
          buf[m++] = j+1;
          buf[m++] = k+1;
          buf[m++] = ctype[i][j];
          buf[m++] = nodex[i][j][k][0];
          buf[m++] = nodex[i][j][k][1];
          buf[m++] = nodex[i][j][k][2];
          if (ids) ids[n++] = i;
        }
}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_atom_noscale_reference(tagint *ids)
{
  double **xoriginal = fix->astore;
  int m, n;
  int i;
  int *mask = atom->mask;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double coord[3];

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = xoriginal[i][0];
      buf[m++] = xoriginal[i][1];
      buf[m++] = xoriginal[i][2];
      coord[0] = x[i][0];
      coord[1] = x[i][1];
      coord[2] = x[i][2];
      domain->unmap(coord, image[i]);
      buf[m++] = coord[0];
      buf[m++] = coord[1];
      buf[m++] = coord[2];
      if (ids) ids[n++] = i;
    }
  }
}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_elem_noscale_reference(tagint *ids)
{
  int m, n;
  int i, j, k, iapc, inpe;
  double xtmp, ytmp, ztmp;
  int *mask = element->mask;
  double ****nodexoriginal = fix->a4store;
  int **ctype = element->ctype;
  int *etype = element->etype;
  tagint *tag = element->tag;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = etype[i];
      xtmp = ytmp = ztmp = 0;
      inpe = npe[etype[i]]; 
      iapc = apc[etype[i]];
      for (j = 0; j < iapc; j++) 
        for (k = 0; k < inpe; k++) {
          xtmp += nodexoriginal[i][j][k][0]; 
          ytmp += nodexoriginal[i][j][k][1]; 
          ztmp += nodexoriginal[i][j][k][2]; 
        }
      buf[m++] = xtmp / inpe / iapc;
      buf[m++] = ytmp / inpe / iapc;
      buf[m++] = ztmp / inpe / iapc;
      if (ids) ids[n++] = i;
    }
  }
}

/* ------------------------------------------------------------------------------------------ */

void DumpCACData::pack_node_noscale_reference(tagint *ids)
{
  int m, n;
  int i, j, k;
  int *mask = element->mask;
  double ****nodex = element->nodex;
  double ****nodexoriginal = fix->a4store;
  int **ctype = element->ctype;
  int *etype = element->etype;
  tagint *tag = element->tag;
  imageint *image = element->image;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *apc = element->apc;
  double coord[3];

  m = n = 0;
  for (i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) 
      for (j = 0; j < apc[etype[i]]; j++) 
        for (k = 0; k < npe[etype[i]]; k++) {
          buf[m++] = tag[i];
          buf[m++] = j+1;
          buf[m++] = k+1;
          buf[m++] = ctype[i][j];
          buf[m++] = nodexoriginal[i][j][k][0];
          buf[m++] = nodexoriginal[i][j][k][1];
          buf[m++] = nodexoriginal[i][j][k][2];
          coord[0] = nodex[i][j][k][0];
          coord[1] = nodex[i][j][k][1];
          coord[2] = nodex[i][j][k][2];
          domain->unmap(coord, image[i]);
          buf[m++] = coord[0];
          buf[m++] = coord[1];
          buf[m++] = coord[2];
          if (ids) ids[n++] = i;
        }
}


/* -------------------------------------------------------------------------------------------- */
int DumpCACData::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n, size_one, mybuf);
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_elem_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_elem_choice)(n, size_one, mybuf);
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_node_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_node_choice)(n, size_one, mybuf);
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_atom_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_atom, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
  return offset;
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_elem_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_elem, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
  return offset;
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_node_noimage(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_node, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        static_cast<int> (mybuf[m+2]), 
        static_cast<int> (mybuf[m+3]), 
        mybuf[m+4], mybuf[m+5], mybuf[m+6]);
    m += 7;
  }
  return offset;
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_atom_reference(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;

  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_atom, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4], 
        mybuf[m+5], mybuf[m+6], mybuf[m+7]);
    m += 8;
  }

  return offset;
}

/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_elem_reference(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_elem, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
  return offset;
}


/* -------------------------------------------------------------------------------------------- */

int DumpCACData::convert_node_reference(int n, int size_one, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    offset += sprintf(&sbuf[offset], format_node, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        static_cast<int> (mybuf[m+2]), 
        static_cast<int> (mybuf[m+3]), 
        mybuf[m+4], mybuf[m+5], mybuf[m+6], 
        mybuf[m+7], mybuf[m+8], mybuf[m+9]);
    m += 10;
  }
  return offset;
}



/* ---------------------------------------------------------------------------------------------------- */

void DumpCACData::write_atom_data(int n, double *mybuf)
{
  (this->*write_atom_choice) (n, mybuf);
}

/* ---------------------------------------------------------------------------------------------------- */

void DumpCACData::write_elem_data(int n, double *mybuf)
{
  (this->*write_elem_choice) (n, mybuf);
}

/* ---------------------------------------------------------------------------------------------------- */

void DumpCACData::write_node_data(int n, double *mybuf)
{
  (this->*write_node_choice) (n, mybuf);
}

/* --------------------------------------------------------------------------------------------------- */

void DumpCACData::write_string(int n, double *mybuf)
{
  fwrite(mybuf, sizeof(char), n, fp);
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_atom_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_atom, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_elem_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_elem, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_node_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_node, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        static_cast<int> (mybuf[m+2]), 
        static_cast<int> (mybuf[m+3]), 
        mybuf[m+4], mybuf[m+5], mybuf[m+6]);
    m += 7;
  }
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_atom_lines_reference(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_atom, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4], 
        mybuf[m+5], mybuf[m+6], mybuf[m+7]);
    m += 8;
  }
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_elem_lines_reference(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_elem, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
  }
}

/* ------------------------------------------------------------------------------------------------- */

void DumpCACData::write_node_lines_reference(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, format_node, 
        static_cast<tagint> (mybuf[m]), 
        static_cast<int> (mybuf[m+1]), 
        static_cast<int> (mybuf[m+2]), 
        static_cast<int> (mybuf[m+3]), 
        mybuf[m+4], mybuf[m+5], mybuf[m+6], 
        mybuf[m+7], mybuf[m+8], mybuf[m+9]);
    m += 10;
  }
}

/*  ----------------------------------------------------------------------  */

int DumpCACData::count_atoms()
{
  int m = 0;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == 0 || igroup == 2) {
    return nlocal;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) m++;
  }

  return m;
}


/*  ----------------------------------------------------------------------  */

int DumpCACData::count_elements()
{
  int m = 0;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *etype = element->etype;
  int *npe = element->npe;
  int *apc = element->apc;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) 
      m += npe[etype[i]] * apc[etype[i]];

  return m;
}

/*  ----------------------------------------------------------------------  */

int DumpCACData::count_nodes()
{
  int m = 0;
  int *mask = element->mask;
  int nlocal = element->nlocal;

  if (igroup == 0 || igroup == 1) {
    return nlocal;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) m++;
  }

  return m;

}

/*  ----------------------------------------------------------------------
    initialize one atom's storage values, called when atom is created
    -------------------------------------------------------------------------  */

void DumpCACData::set_atom_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/*  ----------------------------------------------------------------------
    initialize one element's storage values, called when element is created
    -------------------------------------------------------------------------  */

void DumpCACData::set_elem_arrays(int i)
{
  double ***inodexoriginal = fix->a4store[i];
  double ***inodex = element->nodex[i];
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  for (int j = 0; j < apc[etype[i]]; j++) 
    for (int k = 0; k < npe[etype[i]]; k++) {
      inodexoriginal[j][k][0] = inodex[j][k][0];
      inodexoriginal[j][k][1] = inodex[j][k][1];
      inodexoriginal[j][k][2] = inodex[j][k][2];
    }
}

/*  ----------------------------------------------------------------------
    initialize one atom's storage values passed on from another element, called when atom is created from an element
    -------------------------------------------------------------------------  */

void DumpCACData::pass_on_atom_arrays(int i, int i_old, int iindex)
{
  if (reference_flag) {
    double **xoriginal = fix->astore;
    double ****nodexoriginal = fix->a4store;
    int iapc = element->apc[element->etype[i_old]];
    int iucell = iindex / iapc;
    int ibasis = iindex % iapc;
    element->evec->interpolate(xoriginal[i], nodexoriginal, i_old, ibasis, iucell, 3);
  }
}

/*  ----------------------------------------------------------------------
    initialize one element's storage values passed on from another element, called when element is created from another element
    -------------------------------------------------------------------------  */

void DumpCACData::pass_on_elem_arrays(int i, int i_old, int *iindexlist)
{
  if (reference_flag) {
    double ****nodexoriginal = fix->a4store;
    int *etype = element->etype;
    int inpe = element->npe[etype[i]];
    int iapc = element->apc[etype[i]];
    for (int j = 0; j < inpe; j++) {
      int iucell = iindexlist[j] / iapc;
      int ibasis = iindexlist[j] % iapc;
      element->evec->interpolate(nodexoriginal[i][ibasis][j], nodexoriginal, i_old, ibasis, iucell, 3);
    }
  }
}

