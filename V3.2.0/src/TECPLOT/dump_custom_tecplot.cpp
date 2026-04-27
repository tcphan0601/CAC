#include "TECIO.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "dump_custom_tecplot.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "modify.h"
#include "update.h"
#include "compute.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA   1048576

enum { PLT = 0, SZPLT = 1, ASCII = 2 };

/* ---------------------------------------------------------------------- */

DumpCustomTecplot::DumpCustomTecplot(CAC *cac, int narg, char **arg) :
  DumpCustom(cac, narg, arg)
{
  if (narg < 6) error->all(FLERR, "Illegal dump custom/tecplot command");

  npe_connect = (domain->dimension == 3) ? 8 : 4;
  tec_title = nullptr;
  tec_values = nullptr;
  valuelocation = nullptr;
  dummy = 0;
  dummy1 = 1;
  soltime = 0.0;
  debug_tec = 0;
  maxelem_buf = 0;
  nodecell_ids = nullptr;
  node_ids = nullptr;

  average_flag = 1;
  fileformat = PLT;
  // parse optional keywords, then columns
  // scan for optional keyword(s); everything else is a column keyword
  int cols_start;
  for (int i = 5; i < narg; i += 2) {
    if (strcmp(arg[i], "format") == 0) {
      if (i + 2 > narg) error->all(FLERR, "Illegal dump custom/tecplot command");
      if      (strcmp(arg[i+1], "plt")   == 0) fileformat = PLT;
      else if (strcmp(arg[i+1], "szplt") == 0) fileformat = SZPLT;
      else if (strcmp(arg[i+1], "ascii") == 0) fileformat = ASCII;
      else error->all(FLERR, "Illegal dump custom/tecplot command: unknown format");
    } else if (strcmp(arg[i], "average") == 0) {
      if (i + 2 > narg) error->all(FLERR, "Illegal dump custom/tecplot command");
      if      (strcmp(arg[i+1], "yes") == 0) average_flag = 1;
      else if (strcmp(arg[i+1], "no")  == 0) average_flag = 0;
      else error->all(FLERR, "Illegal dump custom/tecplot command: average yes/no expected");
    } else {
      // no other option keywords found, column keywords start from here
      cols_start = i;
      if (i + 1 > narg) error->all(FLERR, "No column keywords found for dump custom/tecplot command");
      break;
    }
  }

  parse_columns(cols_start, narg, arg);

  // detect format from file extension if not explicitly set
  char *suffix = filename + strlen(filename) - strlen(".plt");
  if (suffix > filename && strcmp(suffix, ".plt")   == 0) fileformat = PLT;
  suffix = filename + strlen(filename) - strlen(".szplt");
  if (suffix > filename && strcmp(suffix, ".szplt") == 0) fileformat = SZPLT;
  suffix = filename + strlen(filename) - strlen(".dat");
  if (suffix > filename && strcmp(suffix, ".dat")   == 0) fileformat = ASCII;

  buffer_allow = 1;
  buffer_flag  = 1;
  clearstep    = 1;
}

/* ---------------------------------------------------------------------- */

DumpCustomTecplot::~DumpCustomTecplot()
{
  delete [] tec_title;
  delete [] tec_values;
  delete [] valuelocation;
  memory->destroy(nodecell_ids);
  memory->destroy(node_ids);
  // ibuf freed by Dump base destructor

  if (multifile == 0 && filewriter && fileformat != ASCII)
    TECEND142();
}

/* ---------------------------------------------------------------------- */

void DumpCustomTecplot::init_style()
{
  resolve_computes();

  // size_one for ASCII path
  atom_size_one = node_size_one = ncols;
  elem_size_one = npe_connect * (average_flag ? 1 : element->maxapc);
  max_size_one  = MAX(atom_size_one, elem_size_one);

  // build TECIO variable string: "col0 col1 col2 ..."
  delete [] tec_values;
  int totlen = 0;
  for (int k = 0; k < ncols; k++) totlen += strlen(col_labels[k]) + 1;
  tec_values = new char[totlen + 2];
  tec_values[0] = '\0';
  for (int k = 0; k < ncols; k++) {
    if (k) strcat(tec_values, " ");
    strcat(tec_values, col_labels[k]);
  }

  // build ASCII columns string (same content, used in write_header)
  delete [] columns;
  int totlen2 = 0;
  for (int k = 0; k < ncols; k++) totlen2 += 2 + strlen(col_labels[k]) + 3;
  columns = new char[totlen2 + 2];
  columns[0] = '\0';
  for (int k = 0; k < ncols; k++) {
    if (k) strcat(columns, ", ");
    strcat(columns, "\"");
    strcat(columns, col_labels[k]);
    strcat(columns, "\"");
  }

  delete [] tec_title;
  tec_title = new char[32];
  strcpy(tec_title, "Output file from CAC");

  delete [] valuelocation;
  valuelocation = new int[ncols];
  for (int k = 0; k < ncols; k++) valuelocation[k] = 1;  // all nodal

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpCustomTecplot::openfile()
{
  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  char *filecurrent = filename;
  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar, '*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent, "%s" BIGINT_FORMAT "%s", filestar, update->ntimestep, ptr+1);
    else {
      char bif[8], pad[16];
      strcpy(bif, BIGINT_FORMAT);
      sprintf(pad, "%%s%%0%d%s%%s", padflag, &bif[1]);
      sprintf(filecurrent, pad, filestar, update->ntimestep, ptr+1);
    }
    *ptr = '*';
  }

  if (filewriter) {
    if (fileformat == ASCII) {
      fp = fopen(filecurrent, "w");
      if (!fp) {
        char str[128];
        sprintf(str, "Cannot open dump file %s", filecurrent);
        error->one(FLERR, str);
      }
    } else {
      int success = TECINI142(tec_title, tec_values, filecurrent,
                               (char *)".", &fileformat, &dummy, &debug_tec, &dummy1);
      if (success == -1) {
        char str[128];
        sprintf(str, "Cannot open dump file %s", filecurrent);
        error->one(FLERR, str);
      }
    }
  } else fp = nullptr;

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */
// Override write() to handle both ASCII and binary TECIO paths.

void DumpCustomTecplot::write()
{
  if (multifile) openfile();

  bigint bnatom_me  = count_atoms();
  bigint bnnode_me  = count_nodes();
  bigint bnelem_me  = count_elements();
  natom_me  = (int)bnatom_me;
  nnode_me  = (int)bnnode_me;
  nelem_me  = (int)bnelem_me;

  MPI_Allreduce(&bnatom_me, &natom_total, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  MPI_Allreduce(&bnnode_me, &nnode_total, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  MPI_Allreduce(&bnelem_me, &nelem_total, 1, MPI_CAC_BIGINT, MPI_SUM, world);

  bigint natom_hdr = natom_total, nnode_hdr = nnode_total, nelem_hdr = nelem_total;
  if (multiproc) {
    MPI_Allreduce(&bnatom_me, &natom_hdr, 1, MPI_CAC_BIGINT, MPI_SUM, clustercomm);
    MPI_Allreduce(&bnnode_me, &nnode_hdr, 1, MPI_CAC_BIGINT, MPI_SUM, clustercomm);
    MPI_Allreduce(&bnelem_me, &nelem_hdr, 1, MPI_CAC_BIGINT, MPI_SUM, clustercomm);
  }

  // invoke computes once for this timestep
  for (int ic = 0; ic < ncompute; ic++) {
    Compute *c = computes[ic];
    if (c->invoked_peratom != update->ntimestep) {
      c->compute_peratom();
      c->addstep(update->ntimestep + nevery);
    }
  }

  int nmax, tmp, nlines, nchars, nsmin, nsmax;
  MPI_Status status;
  MPI_Request request;

  if (fileformat == ASCII) {
    // ---- ASCII path: same pattern as base Dump::write() ----

    if (multifile) write_header(0);

    // atom zone
    if (natom_total) {
      if (filewriter) write_atom_header(natom_hdr);
      if (multiproc != nprocs) MPI_Allreduce(&natom_me, &nmax, 1, MPI_INT, MPI_MAX, world);
      else nmax = natom_me;
      if (nmax > maxbuf) {
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf, maxbuf * max_size_one, "dump:buf");
      }
      pack_atom(nullptr);
      int nsatom_me = convert_atom_string(natom_me, atom_size_one, buf);
      MPI_Allreduce(&nsatom_me, &nsmin, 1, MPI_INT, MPI_MIN, world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs) MPI_Allreduce(&nsatom_me, &nsmax, 1, MPI_INT, MPI_MAX, world);
      else nsmax = nsatom_me;
      if (nsmax > maxsbuf) { maxsbuf = nsmax; memory->grow(sbuf, maxsbuf, "dump:sbuf"); }
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(sbuf, maxsbuf, MPI_CHAR, me+iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_CHAR, &nchars);
          } else nchars = nsatom_me;
          write_atom_data(nchars, (double *)sbuf);
        }
        if (flush_flag && fp) fflush(fp);
      } else {
        MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(sbuf, nsatom_me, MPI_CHAR, fileproc, 0, world);
      }
    }

    // node zone
    if (element->nelements) {
      if (filewriter) write_node_header(nnode_hdr);
      if (multiproc != nprocs) MPI_Allreduce(&nnode_me, &nmax, 1, MPI_INT, MPI_MAX, world);
      else nmax = nnode_me;
      if (nmax > maxbuf) {
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf, maxbuf * max_size_one, "dump:buf");
      }
      pack_node(nullptr);
      int nsnode_me = convert_node_string(nnode_me, node_size_one, buf);
      MPI_Allreduce(&nsnode_me, &nsmin, 1, MPI_INT, MPI_MIN, world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs) MPI_Allreduce(&nsnode_me, &nsmax, 1, MPI_INT, MPI_MAX, world);
      else nsmax = nsnode_me;
      if (nsmax > maxsbuf) { maxsbuf = nsmax; memory->grow(sbuf, maxsbuf, "dump:sbuf"); }
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(sbuf, maxsbuf, MPI_CHAR, me+iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_CHAR, &nchars);
          } else nchars = nsnode_me;
          write_node_data(nchars, (double *)sbuf);
        }
        if (flush_flag && fp) fflush(fp);
      } else {
        MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(sbuf, nsnode_me, MPI_CHAR, fileproc, 0, world);
      }

      // element connectivity
      if (filewriter) write_elem_header(nelem_hdr);
      int nelem_vals = nelem_me * elem_size_one;
      if (multiproc != nprocs) MPI_Allreduce(&nelem_vals, &nmax, 1, MPI_INT, MPI_MAX, world);
      else nmax = nelem_vals;
      if (nmax > maxbuf) {
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf, maxbuf, "dump:buf");
      }
      pack_elem(nullptr);
      int nselem_me = convert_elem_string(nelem_me, elem_size_one, buf);
      MPI_Allreduce(&nselem_me, &nsmin, 1, MPI_INT, MPI_MIN, world);
      if (nsmin < 0) error->all(FLERR, "Too much buffered per-proc info for dump");
      if (multiproc != nprocs) MPI_Allreduce(&nselem_me, &nsmax, 1, MPI_INT, MPI_MAX, world);
      else nsmax = nselem_me;
      if (nsmax > maxsbuf) { maxsbuf = nsmax; memory->grow(sbuf, maxsbuf, "dump:sbuf"); }
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(sbuf, maxsbuf, MPI_CHAR, me+iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_CHAR, &nchars);
          } else nchars = nselem_me;
          write_elem_data(nchars, (double *)sbuf);
        }
        if (flush_flag && fp) fflush(fp);
      } else {
        MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(sbuf, nselem_me, MPI_CHAR, fileproc, 0, world);
      }
    }

    write_box_zone();

    if (multifile) {
      if (filewriter && fp) fclose(fp);
      fp = nullptr;
    }

  } else {
    // ---- Binary TECIO path ----

    int success;

    // --- atom zone ---
    if (natom_total) {
      if (filewriter) {
        zonetype = 0;   // ordered
        int imax = (int)natom_hdr, jmax = 1, kmax = 1;
        success = TECZNE142((char *)"Discrete Atoms", &zonetype,
            &imax, &jmax, &kmax,
            &dummy, &dummy, &dummy, &soltime,
            &dummy, &dummy, &dummy1,
            &dummy, &dummy, &dummy1, &dummy1, &dummy1,
            nullptr, nullptr, nullptr, &dummy);
        if (success == -1) error->one(FLERR, "Cannot write Discrete Atoms Zone");
      }
      grow_buf(natom_me);
      for (int k = 0; k < ncols; k++) {
        pack_atom_col(k, nullptr);
        comm_buf_tecio(natom_me);
      }
    }

    // --- CG element zone (nodes + connectivity) ---
    if (element->nelements) {
      if (filewriter) {
        zonetype = (domain->dimension == 3) ? 5 : 3;  // BRICK or QUAD
        int nnodes   = (int)nnode_hdr;
        int nelems_z = (int)nelem_hdr;
        success = TECZNE142((char *)"Coarse Elements", &zonetype,
            &nnodes, &nelems_z,
            &dummy, &dummy, &dummy, &dummy, &soltime,
            &dummy, &dummy, &dummy1,
            &dummy, &dummy,
            0, 0, 0,
            nullptr, valuelocation, nullptr, &dummy);
        if (success == -1) error->one(FLERR, "Cannot write Coarse Element zone");
      }
      grow_buf(nnode_me);
      for (int k = 0; k < ncols; k++) {
        pack_node_col(k, nullptr);
        comm_buf_tecio(nnode_me);
      }

      // element connectivity — uses int buffer (TECNODE142)
      int nelem_ints = nelem_me * npe_connect * (average_flag ? 1 : element->maxapc);
      if (multiproc != nprocs) MPI_Allreduce(&nelem_ints, &nmax, 1, MPI_INT, MPI_MAX, world);
      else nmax = nelem_ints;
      if (nmax > maxibuf) {
        maxibuf = nmax;
        memory->destroy(ibuf);
        memory->create(ibuf, maxibuf, "dump:ibuf");
      }
      pack_elem(nullptr);   // fills ibuf
      if (filewriter) {
        for (int iproc = 0; iproc < nclusterprocs; iproc++) {
          if (iproc) {
            MPI_Irecv(ibuf, maxibuf, MPI_INT, me+iproc, 0, world, &request);
            MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
            MPI_Wait(&request, &status);
            MPI_Get_count(&status, MPI_INT, &nlines);
          } else nlines = nelem_ints;
          success = TECNODE142(&nlines, ibuf);
          if (success == -1) error->one(FLERR, "Cannot write element connectivity");
        }
      } else {
        MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
        MPI_Rsend(ibuf, nelem_ints, MPI_INT, fileproc, 0, world);
      }
    }

    write_box_zone();

    if (multifile && filewriter) {
      TECEND142();
    }
  }
}

/* ---------------------------------------------------------------------- */
// Atom zone: real atoms only (no vatoms).

void DumpCustomTecplot::pack_atom(tagint *ids)
{
  int *mask   = atom->mask;
  int nlocal  = atom->nlocal;
  int m = 0, n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (int k = 0; k < ncols; k++)
      buf[m++] = get_atom_val(i, coldefs[k]);
    if (ids) ids[n++] = i;
  }
}

// Binary: column k for all real atoms.
void DumpCustomTecplot::pack_atom_col(int k, tagint *ids)
{
  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    buf[m++] = get_atom_val(i, coldefs[k]);
  }
}

int DumpCustomTecplot::count_atoms()
{
  int cnt = 0;
  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) cnt++;
  return cnt;
}

/* ---------------------------------------------------------------------- */
// FE node accessor: gets column value for node (element i, basis j, node inode).

double DumpCustomTecplot::get_node_val(int i, int ibasis, int inode, const ColDef &cd)
{
  int ietype = element->etype[i];
  switch (cd.type) {
    case COL_X:    return element->nodex[i][ibasis][inode][0];
    case COL_Y:    return element->nodex[i][ibasis][inode][1];
    case COL_Z:    return element->nodex[i][ibasis][inode][2];
    case COL_ID:   return element->tag[i];
    case COL_TYPE: return element->ctype[i][ibasis];
    case COL_VX:   return element->nodev[i][ibasis][inode][0];
    case COL_VY:   return element->nodev[i][ibasis][inode][1];
    case COL_VZ:   return element->nodev[i][ibasis][inode][2];
    case COL_FX:   return element->nodef[i][ibasis][inode][0];
    case COL_FY:   return element->nodef[i][ibasis][inode][1];
    case COL_FZ:   return element->nodef[i][ibasis][inode][2];
    case COL_COMPUTE_SCALAR:
      if (cd.compute->element_data_style == Compute::NODE)
        return cd.compute->vector_node[i][ibasis][inode];
      return cd.compute->vector_vatom[i][ibasis][element->n2u[ietype][inode]];
    case COL_COMPUTE_ARRAY:
      if (cd.compute->element_data_style == Compute::NODE)
        return cd.compute->array_node[i][ibasis][inode][cd.col];
      return cd.compute->array_vatom[i][ibasis][element->n2u[ietype][inode]][cd.col];
    case COL_FIX_SCALAR:
      if (cd.fix->element_data_style == Compute::NODE)
        return cd.fix->vector_node[i][ibasis][inode];
      return cd.compute->vector_vatom[i][ibasis][element->n2u[ietype][inode]];
    case COL_FIX_ARRAY:
      if (cd.fix->element_data_style == Compute::NODE)
        return cd.fix->array_node[i][ibasis][inode][cd.col];
      return cd.fix->array_vatom[i][ibasis][element->n2u[ietype][inode]][cd.col];

    // debug columns for FE nodes
    case COL_EID:        return element->tag[i];
    case COL_ETYPE:      return ietype;
    case COL_IBASIS:     return ibasis;
    case COL_INODE:      return inode;
    case COL_IGAUSS:     return element->n2g[ietype][inode];
    case COL_IUCELL:     return element->n2u[ietype][inode];
    case COL_OWNER_PROC: return me;
    default:
      error->all(FLERR, "Dump custom/tecplot: unrecognized column type in get_node_val");
      return 0.0;  // unreachable
  }
}

/* ---------------------------------------------------------------------- */
// Node zone: FE nodes of CG elements.

void DumpCustomTecplot::pack_node(tagint *ids)
{
  int *mask    = element->mask;
  int nelocal  = element->nlocal;
  int *apc     = element->apc;
  int *npe     = element->npe;
  int *etype   = element->etype;
  int m = 0, n = 0;

  if (average_flag) {
    for (int i = 0; i < nelocal; i++) {
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int inode = 0; inode < inpe; inode++) {
        for (int k = 0; k < ncols; k++) {
          double val = 0.0;
          for (int j = 0; j < iapc; j++)
            val += get_node_val(i, j, inode, coldefs[k]);
          buf[m++] = val / iapc;
        }
      }
      if (ids) ids[n++] = i;
    }
  } else {
    for (int i = 0; i < nelocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int j = 0; j < iapc; j++) {
        for (int inode = 0; inode < inpe; inode++) {
          for (int k = 0; k < ncols; k++)
            buf[m++] = get_node_val(i, j, inode, coldefs[k]);
        }
      }
      if (ids) ids[n++] = i;
    }
  }
}

// Binary: column k for FE nodes.
void DumpCustomTecplot::pack_node_col(int k, tagint *ids)
{
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *npe    = element->npe;
  int *etype  = element->etype;
  int m = 0;

  if (average_flag) {
    for (int i = 0; i < nelocal; i++) {
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int inode = 0; inode < inpe; inode++) {
        double val = 0.0;
        for (int j = 0; j < iapc; j++)
          val += get_node_val(i, j, inode, coldefs[k]);
        buf[m++] = val / iapc;
      }
    }
  } else {
    for (int i = 0; i < nelocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int j = 0; j < iapc; j++)
        for (int inode = 0; inode < inpe; inode++)
          buf[m++] = get_node_val(i, j, inode, coldefs[k]);
    }
  }
}

int DumpCustomTecplot::count_nodes()
{
  int cnt = 0;
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *npe    = element->npe;
  int *etype  = element->etype;
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    cnt += npe[etype[i]] * (average_flag ? 1 : apc[etype[i]]);
  }
  return cnt;
}

/* ---------------------------------------------------------------------- */
// Element connectivity — identical to DumpTecplot.

void DumpCustomTecplot::pack_elem(tagint *ids)
{
  int *mask  = element->mask;
  int nelocal = element->nlocal;
  int *apc   = element->apc;
  int *etype = element->etype;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;

  if (maxelem_buf < nelocal) {
    maxelem_buf = element->nmax;
    if (average_flag) {
      memory->destroy(nodecell_ids);
      memory->grow(nodecell_ids, maxelem_buf, npe_connect, "dump:nodecell_ids");
    } else {
      memory->destroy(node_ids);
      memory->grow(node_ids, maxelem_buf, element->maxapc, npe_connect, "dump:node_ids");
    }
  }

  if (average_flag)
    element->set_node_connectivities(igroup, nodecell_ids);
  else
    element->set_node_connectivities(igroup, node_ids);

  int *conn = new int[npe_connect];
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    int ietype = etype[i];
    if (average_flag) {
      element->node_connectivity(conn, element_shape_ids[ietype], nodecell_ids[i]);
      for (int j = 0; j < npe_connect; j++) buf[m++] = conn[j];
    } else {
      for (int j = 0; j < apc[ietype]; j++) {
        element->node_connectivity(conn, element_shape_ids[ietype], node_ids[i][j]);
        for (int k2 = 0; k2 < npe_connect; k2++) buf[m++] = conn[k2];
      }
    }
  }
  // binary path fills ibuf separately; for ASCII path buf holds the data
  // For binary, reuse buf as int buffer (only called from the int-buf branch)
  if (fileformat != ASCII) {
    // repack into ibuf
    for (int i = 0; i < m; i++) ibuf[i] = (int)buf[i];
  }
  delete [] conn;
}

int DumpCustomTecplot::count_elements()
{
  int cnt = 0;
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *etype  = element->etype;
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    cnt += (average_flag ? 1 : apc[etype[i]]);
  }
  return cnt;
}

/* ---------------------------------------------------------------------- */
// ASCII headers and write methods.

void DumpCustomTecplot::write_header(bigint /*ndump*/)
{
  if (filewriter)
    fprintf(fp, "title = \"CAC Simulation\"\nvariables = %s\n", columns);
}

void DumpCustomTecplot::write_atom_header(bigint /*ndump*/)
{
  if (filewriter)
    fprintf(fp, "zone t = \"Discrete Atoms\", datapacking = point\n");
}

void DumpCustomTecplot::write_node_header(bigint ndump)
{
  if (filewriter) {
    fprintf(fp, "zone t = \"Coarse Element\" n = " BIGINT_FORMAT " e = " BIGINT_FORMAT,
            ndump, nelem_total);
    if (domain->dimension == 3)
      fprintf(fp, " datapacking = point, zonetype = febrick\n");
    else
      fprintf(fp, " datapacking = point, zonetype = fequadrilateral\n");
  }
}

int DumpCustomTecplot::convert_atom_string(int n, int size_one, double *mybuf)
{
  int offset = 0, m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint)maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    for (int k = 0; k < size_one; k++) {
      if (k) offset += sprintf(&sbuf[offset], " ");
      offset += sprintf(&sbuf[offset], "%g", mybuf[m++]);
    }
    offset += sprintf(&sbuf[offset], "\n");
  }
  return offset;
}

int DumpCustomTecplot::convert_node_string(int n, int size_one, double *mybuf)
{
  return convert_atom_string(n, size_one, mybuf);
}

int DumpCustomTecplot::convert_elem_string(int n, int size_one, double *mybuf)
{
  int offset = 0, m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint)maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }
    for (int k = 0; k < size_one; k++) {
      if (k) offset += sprintf(&sbuf[offset], " ");
      offset += sprintf(&sbuf[offset], "%d", static_cast<int>(mybuf[m++]));
    }
    offset += sprintf(&sbuf[offset], "\n");
  }
  return offset;
}

void DumpCustomTecplot::write_atom_data(int n, double *mybuf)
{ fwrite(mybuf, sizeof(char), n, fp); }

void DumpCustomTecplot::write_node_data(int n, double *mybuf)
{ fwrite(mybuf, sizeof(char), n, fp); }

void DumpCustomTecplot::write_elem_data(int n, double *mybuf)
{ fwrite(mybuf, sizeof(char), n, fp); }

/* ---------------------------------------------------------------------- */
// TECIO helpers — copied from DumpTecplot.

void DumpCustomTecplot::grow_buf(int nme)
{
  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme, &nmax, 1, MPI_INT, MPI_MAX, world);
  else nmax = nme;
  if (nmax > maxbuf) {
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf, maxbuf, "dump:buf");
  }
}

void DumpCustomTecplot::comm_buf_tecio(int nme)
{
  int nmax, tmp;
  MPI_Status status;
  MPI_Request request;

  if (multiproc != nprocs) MPI_Allreduce(&nme, &nmax, 1, MPI_INT, MPI_MAX, world);
  else nmax = nme;

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      int nvals;
      if (iproc) {
        MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &nvals);
      } else nvals = nme;
      int success = TECDAT142(&nvals, buf, &dummy1);
      if (success == -1) error->one(FLERR, "Cannot write to TecPlot dump file");
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(buf, nme, MPI_DOUBLE, fileproc, 0, world);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomTecplot::write_box_zone()
{
  if (!filewriter) return;

  double xlo = domain->boxlo[0], ylo = domain->boxlo[1], zlo = domain->boxlo[2];
  double lx  = domain->prd[0],   ly  = domain->prd[1],   lz  = domain->prd[2];
  double xy  = domain->xy,       xz  = domain->xz,       yz  = domain->yz;

  double cx[8] = {
    xlo,            xlo+lx,
    xlo+lx+xy,      xlo+xy,
    xlo+xz,         xlo+lx+xz,
    xlo+lx+xy+xz,   xlo+xy+xz
  };
  double cy[8] = {
    ylo, ylo, ylo+ly, ylo+ly,
    ylo+yz, ylo+yz, ylo+ly+yz, ylo+ly+yz
  };
  double cz[8] = {
    zlo, zlo, zlo, zlo,
    zlo+lz, zlo+lz, zlo+lz, zlo+lz
  };

  if (fileformat == ASCII) {
    fprintf(fp, "zone t=\"Simulation Cell\", n=8, e=1, datapacking=point");
    fprintf(fp, (domain->dimension == 3) ? ", zonetype=febrick\n"
                                         : ", zonetype=fequadrilateral\n");
    int nextra = ncols - 3;  // zeros for all columns after x y z
    for (int i = 0; i < 8; i++) {
      fprintf(fp, "%g %g %g", cx[i], cy[i], cz[i]);
      for (int j = 0; j < nextra; j++) fprintf(fp, " 0");
      fprintf(fp, "\n");
    }
    fprintf(fp, "1 2 3 4 5 6 7 8\n");
  } else {
    int nvars = ncols;
    int *passive = new int[nvars];
    for (int i = 0; i < nvars; i++) passive[i] = (i < 3) ? 0 : 1;  // only x,y,z active
    int nnodes = 8, nelems = 1;
    int zt = (domain->dimension == 3) ? 5 : 3;
    int success = TECZNE142((char *)"Simulation Cell",
        &zt, &nnodes, &nelems,
        &dummy, &dummy, &dummy, &dummy, &soltime,
        &dummy, &dummy, &dummy1,
        &dummy, &dummy, 0, 0, 0,
        passive, nullptr, nullptr, &dummy);
    delete [] passive;
    if (success == -1) error->one(FLERR, "Cannot write Simulation Cell zone");
    int nvals = 8;
    TECDAT142(&nvals, cx, &dummy1);
    TECDAT142(&nvals, cy, &dummy1);
    TECDAT142(&nvals, cz, &dummy1);
    int connect[8] = {1,2,3,4,5,6,7,8};
    int nconnect = 8;
    TECNODE142(&nconnect, connect);
  }
}
