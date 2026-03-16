#include <string.h>
#include "dump_cac_data.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
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
#define NSECTIONS 7
#define CHUNK 1024
#define MAXLINE 256

enum{INTERPOLATE, INTEGRATION, NODE, CENTER};
enum{OFF, DATAFILE, CURRENT};

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
  reference_triclinic = 0;
  ref_fp = NULL;
  precision_flag = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "reference") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "current") == 0) reference_flag = CURRENT;
      else if (strcmp(arg[iarg+1], "data/file") == 0) reference_flag = DATAFILE;
      else if (strcmp(arg[iarg+1], "off") == 0) reference_flag = OFF;
      else error->all(FLERR, "Illegal dump atom command");
      if (reference_flag == DATAFILE) {
        if (iarg+3 > narg) error->all(FLERR, "Illegal dump atom command");
        if (comm->me == 0) {
          ref_fp = fopen(arg[iarg+2], "r");
          if (ref_fp == NULL) {
            char str[128];
            sprintf(str, "Cannot open file %s", arg[iarg+2]);
            error->one(FLERR, str);
          }
        }
        iarg += 3;
      } else 
        iarg += 2;
    } else if (strcmp(arg[iarg], "precision") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "high") == 0) precision_flag = 1;
      else if (strcmp(arg[iarg+1], "normal") == 0) precision_flag = 0;
      else error->all(FLERR, "Illegal dump atom command");
      iarg += 2;
    } else error->all(FLERR, "Illegal dump atom command");
  }

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

    // store current configuration as reference configuration

    if (reference_flag == CURRENT) {

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
        if (amask[i] & groupbit)
          domain->unmap(x[i], aimage[i], xoriginal[i]);
        else
          xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
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

      // store simulation box bounds

      boxxlo_ref = domain->boxlo[0];
      boxxhi_ref = domain->boxhi[0];
      boxylo_ref = domain->boxlo[1];
      boxyhi_ref = domain->boxhi[1];
      boxzlo_ref = domain->boxlo[2];
      boxzhi_ref = domain->boxhi[2];
      if (domain->triclinic) {
        reference_triclinic = 1;
        boxxy_ref = domain->xy;
        boxxz_ref = domain->xz;
        boxyz_ref = domain->yz;
      }
    } 

    // store configuration from external data file as reference configuration

    else {
      memory->create(atom_read_flag, atom->nlocal, "dump:atom_read_flag");
      memory->create(node_read_flag, element->nlocal, 
          element->maxapc, element->maxnpe, "dump:atom_read_flag");
      memset(&atom_read_flag[0], 0, atom->nlocal * sizeof(int));
      memset(&node_read_flag[0][0][0], 0, 
          element->nlocal * element->maxapc * element->maxnpe * sizeof(int));
      read_file();

      if (comm->me == 0)
        fclose(ref_fp); 
      
      // check if  all atoms/nodes have reference coords  
      
      int *etype = element->etype;
      int *npe = element->npe;
      int *apc = element->apc;
      for (int i = 0; i < atom->nlocal; i++)
        if (!atom_read_flag[i]) error->one(FLERR,"Missing atoms in external reference file");
      for (int i = 0; i < element->nlocal; i++)
        for (int ibasis = 0; ibasis < apc[etype[i]]; ibasis++)
          for (int inode = 0; inode < npe[etype[i]]; inode++)
            if (!node_read_flag[i][ibasis][inode]) {
              printf("i = %d %d %d etype = %d npe= %d\n",element->tag[i],inode,ibasis,etype[i],npe[etype[i]]);
              error->one(FLERR,"Missing nodes in external reference file");
            }

      memory->destroy(atom_read_flag);
      memory->destroy(node_read_flag);
    }
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
  if (reference_flag) {
    atom_size_one += 3;
    node_size_one += 3;
    if (image_flag) error->all(FLERR, "Cannot dump image and reference at the same time for dump cac/data");
  }

  max_size_one = MAX(atom_size_one, elem_size_one);
  max_size_one = MAX(max_size_one, node_size_one);

  int n;

  if (reference_flag) {
    if (precision_flag) {
      delete [] format_atom;
      char *str;
      str = (char *) TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_atom = new char[n];
      strcpy(format_atom, str);

      delete [] format_elem;
      str = (char *) TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_elem = new char[n];
      strcpy(format_elem, str);

      delete [] format_node;
      str = (char *) TAGINT_FORMAT " %d %d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n";
      n = strlen(str) + 2;
      format_node = new char[n];
      strcpy(format_node, str);
    } else {
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
    }
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


  fprintf(fp, "%-1.16e %-1.16e xlo xhi\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e ylo yhi\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e zlo zhi\n", domain->boxlo[2], domain->boxhi[2]);
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
    fprintf(fp, BIGINT_FORMAT " elements\n", nnode_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xloref xhiref\n", boxxlo_ref, boxxhi_ref);
  fprintf(fp, "%-2.16e %-1.16e yloref yhiref\n", boxylo_ref, boxyhi_ref);
  fprintf(fp, "%-1.16e %-1.16e zloref zhiref\n", boxzlo_ref, boxzhi_ref);
  if (reference_triclinic)
    fprintf(fp, "%-1.16e %-1.16e %-1.16e xyref xzref yzref\n", 
        boxxy_ref, boxxz_ref, boxyz_ref);
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
    fprintf(fp, BIGINT_FORMAT " elements\n", nnode_total);
  fprintf(fp, "%d atom types\n", atom->ntypes); 
  if (nelem_total) 
    fprintf(fp, "%d element types\n", element->netypes); 


  fprintf(fp, "%-1.16e %-1.16e xloref xhiref\n", boxxlo_ref, boxxhi_ref);
  fprintf(fp, "%-1.16e %-1.16e yloref yhiref\n", boxylo_ref, boxyhi_ref);
  fprintf(fp, "%-1.16e %-1.16e zloref zhiref\n", boxzlo_ref, boxzhi_ref);
  if (reference_triclinic)
    fprintf(fp, "%-1.16e %-1.16e %-1.16e xyref xzref yzref\n", 
        boxxy_ref, boxxz_ref, boxyz_ref);
  fprintf(fp, "%-1.16e %-1.16e xlocur xhicur\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e ylocur yhicur\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e zlocur zhicur\n", domain->boxlo[2], domain->boxhi[2]);
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
          buf[m++] = k+1;
          buf[m++] = j+1;
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
          buf[m++] = k+1;
          buf[m++] = j+1;
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

/*  ----------------------------------------------------------------------
    read free-format header of data file
    1st line and blank lines are skipped
    non-blank lines are checked for header keywords and leading value is read
    header ends with EOF or non-blank line containing no header keyword
    if EOF, line is set to blank line
    else line has first keyword line for rest of file
    some logic differs if adding atoms
    -------------------------------------------------------------------------  */

void DumpCACData::header()
{

  int n;
  char *ptr;

  // customize for new sections
  // When adding new sections, remember to change NSECTIONS
  const char *section_keywords[NSECTIONS] =
  {"Atoms", "Velocities", "Masses", "Element Types", 
    "Elements", "Nodes", "Node Velocities"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line, MAXLINE, ref_fp);
    if (eof == NULL) error->one(FLERR, "Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line, MAXLINE, ref_fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    MPI_Bcast(line, n, MPI_CHAR, 0, world); 

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    if (strspn(line, " \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable
    // customize for new header lines

    bigint natoms, nelements;
    if (strstr(line, "atoms")) {
      sscanf(line, BIGINT_FORMAT, &natoms);
      if (atom->natoms != natoms) 
        error->all(FLERR,"Number of atoms in the reference file does not match with current system");
    } else if (strstr(line, "elements")) {
      sscanf(line, BIGINT_FORMAT, &nelements);
      if (element->nelements != nelements) 
        error->all(FLERR,"Number of elements in the reference file does not match with current system");
    }

    // Atom/Element class type and grain settings are skipped

    else if (strstr(line, "atom types")) {
    } else if (strstr(line, "element types")) {
    } else if (strstr(line, "atom grains")) {
    }

    // local copy of box info
    // if file contain reference frame, store reference box and skip current box

    else if (strstr(line,"xlocur xhicur")) {
    } else if (strstr(line,"ylocur yhicur")) {
    } else if (strstr(line,"zlocur zhicur")) {
    } else if (strstr(line,"xycur xzcur yzcur")) {
    } else if (strstr(line,"xloref xhiref")) {
      external_file_reference = 1;
      sscanf(line,"%lg %lg",&boxxlo_ref,&boxxhi_ref);
    } else if (strstr(line,"yloref yhiref")) {
      sscanf(line,"%lg %lg",&boxylo_ref,&boxyhi_ref);
    } else if (strstr(line,"zloref zhiref")) {
      sscanf(line,"%lg %lg",&boxzlo_ref,&boxzhi_ref);
    } else if (strstr(line,"xyref xzref yzref")) {
      reference_triclinic = 1;
      sscanf(line,"%lg %lg %lg",&boxxy_ref,&boxxz_ref,&boxyz_ref);
    } else if (strstr(line, "xlo xhi")) {
      external_file_reference = 0;
      sscanf(line,"%lg %lg",&boxxlo_ref,&boxxhi_ref);
    } else if (strstr(line, "ylo yhi")) {
      sscanf(line,"%lg %lg",&boxylo_ref,&boxyhi_ref);
    } else if (strstr(line, "zlo zhi")) {
      sscanf(line,"%lg %lg",&boxzlo_ref,&boxzhi_ref);
    } else if (strstr(line, "xy xz yz")) {
      reference_triclinic = 1;
      sscanf(line,"%lg %lg %lg",&boxxy_ref,&boxxz_ref,&boxyz_ref);
    } else break;
  }

  // check that existing string is a valid section keyword

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword, section_keywords[n]) == 0) break;

  if (n == NSECTIONS) {
    char str[128];
    sprintf(str, "Unknown identifier in data file: %s", keyword);
    error->all(FLERR, str);
  } 

}

/*  ----------------------------------------------------------------------
    grab next keyword
    read lines until one is non-blank
    keyword is all text on line w/out leading & trailing white space
    optional style can be appended after comment char '#'
    read one additional line (assumed blank)
    if any read hits EOF, set keyword to empty
    if first = 1, line variable holds non-blank line that ended header
    -------------------------------------------------------------------------  */

void DumpCACData::parse_keyword(int first)
{
  int eof = 0;
  int done = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line, MAXLINE, ref_fp) == NULL) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line, " \t\n\r");
      if ((blank == strlen(line)) || (line[blank] == '#')) {
        if (fgets(line, MAXLINE, ref_fp) == NULL) eof = 1;
      } else done = 1;
    }
    if (fgets(buffer, MAXLINE, ref_fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof, 1, MPI_INT, 0, world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  MPI_Bcast(line, n, MPI_CHAR, 0, world);

  // store optional "style" following comment char '#' after keyword

  char *ptr;
  if ((ptr = strchr(line, '#'))) {
    *ptr++ = '\0';
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    int stop = strlen(ptr) - 1;
    while (ptr[stop] == ' ' || ptr[stop] == '\t'
        || ptr[stop] == '\n' || ptr[stop] == '\r') stop--;
    ptr[stop+1] = '\0';
    strcpy(style, ptr);
  } else style[0] = '\0';

  // copy non-whitespace portion of line into keyword

  int start = strspn(line, " \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
      || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword, &line[start]);
}

/*  ----------------------------------------------------------------------
    read external data file and store coords as reference coords in fix store
    -------------------------------------------------------------------------  */

void DumpCACData::read_file()
{

  int nchunk, eof;
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK * MAXLINE];
  bigint nread;

  // read header info

  header();

  // count total nodes to ensure element->nnodes is updated

  element->count_nodes(1);

  // customize for new sections
  // read rest of file in free format
  // skip Masses, Element Types, Velocities, Elements, and Node Velocities sections

  while (strlen(keyword)) {
    if (strcmp(keyword, "Masses") == 0) {
      eof = comm->read_lines_from_file(ref_fp, atom->ntypes, MAXLINE, buffer);
      if (eof) error->all(FLERR, "Unexpected end of data file");
    } else if (strcmp(keyword, "Element Types") == 0) {
      eof = comm->read_lines_from_file(ref_fp, element->netypes, MAXLINE, buffer);
      if (eof) error->all(FLERR, "Unexpected end of data file");
    } else if (strcmp(keyword, "Atoms") == 0) {
      atoms();
    } else if (strcmp(keyword, "Velocities") == 0) {
      nread = 0;
      while (nread < atom->natoms) {
        nchunk = MIN(atom->natoms - nread, CHUNK);
        eof = comm->read_lines_from_file(ref_fp, nchunk, MAXLINE, buffer);
        if (eof) error->all(FLERR, "Unexpected end of data file");
        nread += nchunk;
      }
    } else if (strcmp(keyword, "Elements") == 0) {
      nread = 0;
      while (nread < element->nelements) {
        nchunk = MIN(element->nelements - nread, CHUNK);
        eof = comm->read_lines_from_file(ref_fp, nchunk, MAXLINE, buffer);
        if (eof) error->all(FLERR, "Unexpected end of data file");
        nread += nchunk;
      }
    } else if (strcmp(keyword, "Nodes") == 0) {
      nodes();
    } else if (strcmp(keyword, "Node Velocities") == 0) {
      nread = 0;
      while (nread < element->nnodes) {
        nchunk = MIN(element->nnodes - nread, CHUNK);
        eof = comm->read_lines_from_file(ref_fp, nchunk, MAXLINE, buffer);
        if (eof) error->all(FLERR, "Unexpected end of data file");
        nread += nchunk;
      }
    } else {
      char str[128];
      sprintf(str, "Unknown identifier in data file: %s", keyword);
      error->all(FLERR, str);
    }
    parse_keyword(0);
  }
  delete [] line;
  delete [] keyword;
  delete [] buffer;
}

/*  ----------------------------------------------------------------------
    read all atoms
    -------------------------------------------------------------------------  */

void DumpCACData::atoms()
{
  int nchunk, eof;
  int i, m, xptr;
  double *coord;
  char *next;
  char **values = NULL;
  tagint itag;
  double **xoriginal = fix->astore;

  bigint nread = 0;
  bigint natoms = atom->natoms;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }
  while (nread < natoms) {
    nchunk = MIN(natoms - nread, CHUNK);
    eof = comm->read_lines_from_file(ref_fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");

    char *buf = buffer;
    next = strchr(buf, '\n');
    *next = '\0';
    int nwords = universe->count_words(buf);
    *next = '\n';

    if (values == NULL)
      values = new char*[nwords];

    // xptr = which word in line starts xyz coords
    // if reference file contain reference frame, 
    // it still start at the same position

    xptr = atom->avec->xcol_data - 1;

    // loop over lines of atom data
    // tokenize the line into values
    // extract xyz coords and image flags
    // remap atom into simulation box

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf, '\n');

      values[0] = strtok(buf, " \t\n\r\f");
      if (values[0] == NULL)
        error->all(FLERR, "Incorrect atom format in data file");
      for (m = 1; m < nwords; m++) {
        values[m] = strtok(NULL, " \t\n\r\f");
        if (values[m] == NULL)
          error->all(FLERR, "Incorrect atom format in data file");
      }

      itag = ATOTAGINT(values[0]);
      if ((m = atom->map(itag)) >= 0) {
        if (atom_read_flag[m]) error->one(FLERR,"Duplicate atoms in external reference file");
        atom_read_flag[m] = 1;
        xoriginal[m][0] = atof(values[xptr]);
        xoriginal[m][1] = atof(values[xptr + 1]);
        xoriginal[m][2] = atof(values[xptr + 2]);

      }
      buf = next + 1;

    }
    nread += nchunk;
  }

  delete [] values;
  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }
}

/*  ----------------------------------------------------------------------
    read all nodes
    to find elements, must build element map
    -------------------------------------------------------------------------  */

void DumpCACData::nodes()
{
  int nchunk, eof;
  int i, m, inode, ibasis;
  tagint itag;
  char *next;
  char **values = NULL;
  double ****nodexoriginal = fix->a4store;

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  bigint nread = 0;
  bigint nnodes = element->nnodes;
  while (nread < nnodes) {
    nchunk = MIN(nnodes - nread, CHUNK);
    eof = comm->read_lines_from_file(ref_fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");

    char *buf = buffer;
    next = strchr(buf, '\n');
    *next = '\0';
    int nwords = universe->count_words(buf);
    *next = '\n';

    if (values == NULL)
      values = new char*[nwords];

    // loop over lines of node 
    // tokenize the line into values
    // if I own element tag, unpack its values

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf, '\n');

      values[0] = strtok(buf, " \t\n\r\f");
      if (values[0] == NULL)
        error->all(FLERR, "Incorrect node format in data file");
      for (m = 1; m < nwords; m++) {
        values[m] = strtok(NULL, " \t\n\r\f");
        if (values[m] == NULL)
          error->all(FLERR, "Incorrect node format in data file");
      }

      itag = ATOTAGINT(values[0]);

      if ((m = element->map(itag)) >= 0) {
        inode = atoi(values[1]) - 1; 
        ibasis = atoi(values[2]) - 1; 
        if (node_read_flag[m][ibasis][inode]) 
          error->one(FLERR,"Duplicate nodes in external reference file");
        node_read_flag[m][ibasis][inode] = 1;
        nodexoriginal[m][ibasis][inode][0] = atof(values[4]);
        nodexoriginal[m][ibasis][inode][1] = atof(values[5]);
        nodexoriginal[m][ibasis][inode][2] = atof(values[6]);
      }

      buf = next + 1; 
    }

    nread += nchunk;
  }  

  delete [] values;
  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

}



