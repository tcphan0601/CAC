#ifdef TECIO_ENABLED        
#include "TECIO.h"
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "dump_custom.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "domain.h"
#include "modify.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "compute.h"
#include "fix.h"
#include "universe.h"

using namespace CAC_NS;

#define ONELINE 256
#define DELTA   1048576

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(CAC *cac, int narg, char **arg) :
  Dump(cac, narg, arg)
{
  if (narg < 6) error->all(FLERR, "Illegal dump custom/tecplot command");

  // tecplot format
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

  dump_format = LAMMPS;
  dump_style = FULLMAP;

  ncols = 0;
  coldefs = nullptr;
  ncompute = 0;
  id_compute = nullptr;
  computes = nullptr;
  nfix = 0;
  id_fix = nullptr;
  fixes = nullptr;
  wrap_flag = 0; 
  precision_flag = 0;

  int cols_start;

  for (int iarg = 5; iarg < narg; iarg += 2) {
    if (strcmp(arg[iarg], "format") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump custom/atom command: missing style value");
      if (strcmp(arg[iarg+1], "lammps/full/map") == 0) {
        dump_format = LAMMPS;
        dump_style = FULLMAP;
      } else if (strcmp(arg[iarg+1], "lammps/surface") == 0) {
        dump_format = LAMMPS;
        dump_style = SURFACE;
      } else if (strcmp(arg[iarg+1], "lammps/gauss/point") == 0) {
        dump_format = LAMMPS;
        dump_style = GAUSSPOINT;
      } else if (strcmp(arg[iarg+1], "lammps/node") == 0) {
        dump_format = LAMMPS;
        dump_style = NODE;
      } else if (strcmp(arg[iarg+1], "lammps/center") == 0) {
        dump_format = LAMMPS;
        dump_style = CENTER;
      } else if (strcmp(arg[iarg+1], "tecplot/ascii") == 0) {
        dump_format = TECPLOT;
        dump_style = ASCII;
#ifdef TECIO_ENABLED        
      } else if (strcmp(arg[iarg+1], "tecplot/plt") == 0) {
        dump_format = TECPLOT;
        dump_style = PLT;
      } else if (strcmp(arg[iarg+1], "tecplot/szplt") == 0) {
        dump_format = TECPLOT;
        dump_style = SZPLT;
#else
      } else if (strcmp(arg[iarg+1], "tecplot/plt") == 0 || strcmp(arg[iarg+1], "tecplot/szplt") == 0) {
        error->all(FLERR, "Install TECPLOT package for writing .plt or .szplt files by running 'make yes-tecplot'");
#endif
      } else error->all(FLERR, "Illegal dump custom/atom command: unknown style");
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump custom/tecplot command");
      if      (strcmp(arg[iarg+1], "yes") == 0) average_flag = 1;
      else if (strcmp(arg[iarg+1], "no")  == 0) average_flag = 0;
      else error->all(FLERR, "Illegal dump custom/tecplot command: average yes/no expected");
    } else if (strcmp(arg[iarg], "wrap") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "yes") == 0) wrap_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) wrap_flag = 0;
      else error->all(FLERR, "Illegal dump custom command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "precision") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dump atom command");
      if (strcmp(arg[iarg+1], "high") == 0) precision_flag = 1;
      else if (strcmp(arg[iarg+1], "normal") == 0) precision_flag = 0;
      else error->all(FLERR, "Illegal dump custom command");
      iarg += 2;
    } else {
      // no other option keywords found, column keywords start from here
      cols_start = iarg;
      if (iarg + 1 > narg) error->all(FLERR, "No column keywords found for dump custom/atom command");
      break;
    }
  }

  parse_columns(cols_start, narg, arg);

  // Tecplot style must have x y z columns 
  // Also move x y z columns to the front (columns 1 2 3 respectively)
  if (dump_format == TECPLOT) {
    int col_index;
    ColDef cdtmp;

    // check for column z and move it to the front
    col_index = -1;
    for (int k = 0; k < ncols; k++)
      if (coldefs[k].type == COL_Z) col_index = k;
    if (col_index < 0)
      error->all(FLERR, "Dump custom tecplot format must include z columns");
    cdtmp = coldefs[col_index];
    for (int k = col_index; k > 0; k--)
      coldefs[k] = coldefs[k-1];
    coldefs[0] = cdtmp;

    // check for column y and move it to the front
    col_index = -1;
    for (int k = 0; k < ncols; k++)
      if (coldefs[k].type == COL_Y) col_index = k;
    if (col_index < 0)
      error->all(FLERR, "Dump custom tecplot format must include y columns");
    cdtmp = coldefs[col_index];
    for (int k = col_index; k > 0; k--)
      coldefs[k] = coldefs[k-1];
    coldefs[0] = cdtmp;

    // check for column x and move it to the front
    col_index = -1;
    for (int k = 0; k < ncols; k++)
      if (coldefs[k].type == COL_X) col_index = k;
    if (col_index < 0)
      error->all(FLERR, "Dump custom tecplot format must include x columns");
    cdtmp = coldefs[col_index];
    for (int k = col_index; k > 0; k--)
      coldefs[k] = coldefs[k-1];
    coldefs[0] = cdtmp;

  }

  // check file extension for tecplot format
  if (dump_format == TECPLOT) {
    char *suffix = filename + strlen(filename) - strlen(".dat");
    if (suffix > filename && strcmp(suffix, ".dat") != 0 && dump_style == ASCII)
      error->all(FLERR, "Dump custom format tecplot/ascii requires .dat file extension");
#ifdef TECIO_ENABLED        
    suffix = filename + strlen(filename) - strlen(".plt");
    if (suffix > filename && strcmp(suffix, ".plt") != 0 && dump_style == PLT)
      error->all(FLERR, "Dump custom format tecplot/plt requires .plt file extension");
    suffix = filename + strlen(filename) - strlen(".szplt");
    if (suffix > filename && strcmp(suffix, ".szplt") != 0 && dump_style == SZPLT)
      error->all(FLERR, "Dump custom format tecplot/plt requires .szplt file extension");
#endif
  }

  buffer_allow = 1;
  buffer_flag  = 1;
  clearstep    = 1;

  nevery = universe->inumeric(FLERR, arg[3]);
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  for (int i = 0; i < ncols; i++)
    delete [] coldefs[i].col_label;
  delete [] coldefs;
  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  delete [] id_compute;
  delete [] computes;
  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  delete [] id_fix;
  delete [] fixes;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::init_style()
{
  // store all computes requested
  for (int icompute = 0; icompute < ncompute; icompute++) {
    int index_compute = modify->find_compute(id_compute[icompute]);
    if (index_compute < 0)
      error->all(FLERR, "Dump custom: compute ID does not exist");
    Compute *compute = modify->compute[index_compute];
    if (!compute->peratom_flag)
      error->all(FLERR, "Dump custom: compute does not produce per-atom output");
    computes[icompute] = compute;
  }

  // resolve packed col fields: separate ic and colidx from the packed int
  for (int k = 0; k < ncols; k++) {
    ColDef &cd = coldefs[k];
    if (cd.type == COL_COMPUTE_SCALAR) {
      int ic = -(cd.col + 2);
      cd.compute = computes[ic];
      cd.col = 0;
    } else if (cd.type == COL_COMPUTE_ARRAY) {
      int ic     = cd.col >> 16;
      int colidx = cd.col & 0xFFFF;
      cd.compute = computes[ic];
      cd.col     = colidx;
      if (cd.compute->size_peratom_cols == 0)
        error->all(FLERR, "Dump custom: compute produces scalar, not array; use c_ID not c_ID[N]");
      if (colidx >= cd.compute->size_peratom_cols)
        error->all(FLERR, "Dump custom: column index out of range for compute");
    }
  }

  for (int ifix = 0; ifix < nfix; ifix++) {
    int idx = modify->find_fix(id_fix[ifix]);
    if (idx < 0)
      error->all(FLERR, "Dump custom: fix ID does not exist");
    Fix *f = modify->fix[idx];
    if (!f->peratom_flag)
      error->all(FLERR, "Dump custom: fix does not produce per-atom output");
    fixes[ifix] = f;
  }

  for (int k = 0; k < ncols; k++) {
    ColDef &cd = coldefs[k];
    if (cd.type == COL_FIX_SCALAR) {
      int ifix = -(cd.col + 2);
      cd.fix = fixes[ifix];
      cd.col = 0;
    } else if (cd.type == COL_FIX_ARRAY) {
      int ifix   = cd.col >> 16;
      int colidx = cd.col & 0xFFFF;
      cd.fix = fixes[ifix];
      cd.col = colidx;
      if (cd.fix->size_peratom_cols == 0)
        error->all(FLERR, "Dump custom: fix produces scalar, not array; use f_ID not f_ID[N]");
      if (colidx >= cd.fix->size_peratom_cols)
        error->all(FLERR, "Dump custom: column index out of range for fix");
    }
  }

  if (dump_format == LAMMPS) {
    atom_size_one = ncols;
    elem_size_one = node_size_one = 0;
    max_size_one  = atom_size_one;

    domain->boundary_string(boundstr);
    if (domain->triclinic == 0)
      header_choice = &DumpCustom::header_item_lammps;
    else
      header_choice = &DumpCustom::header_item_lammps_triclinic;
    atom_header_choice = &DumpCustom::empty_header;
    node_header_choice = &DumpCustom::empty_header;
    elem_header_choice = &DumpCustom::empty_header;

    pack_atom_choice = &DumpCustom::pack_atom_lammps;
    pack_elem_choice = &DumpCustom::pack_empty;
    pack_node_choice = &DumpCustom::pack_empty;

    convert_atom_choice = &DumpCustom::convert_string;
    convert_elem_choice = &DumpCustom::convert_empty;
    convert_node_choice = &DumpCustom::convert_empty;

    if (buffer_flag == 1) 
      write_atom_choice = &DumpCustom::write_string;
    else 
      write_atom_choice = &DumpCustom::write_lines;
    write_elem_choice = &DumpCustom::write_empty;
    write_node_choice = &DumpCustom::write_empty;

    // build column string for header
    delete [] columns;
    int totlen = 0;
    for (int k = 0; k < ncols; k++) totlen += strlen(coldefs[k].col_label) + 1;
    columns = new char[totlen + 2];
    columns[0] = '\0';
    for (int k = 0; k < ncols; k++) {
      if (k) strcat(columns, " ");
      strcat(columns, coldefs[k].col_label);
    }
  } else if (dump_format == TECPLOT) {
    atom_size_one = node_size_one = ncols;
    elem_size_one = npe_connect * (average_flag ? 1 : element->maxapc);
    max_size_one  = MAX(atom_size_one, elem_size_one);

    if (dump_style == ASCII) {
      // build ASCII columns string (same content, used in write_header)
      delete [] columns;
      int totlen2 = 0;
      for (int k = 0; k < ncols; k++) totlen2 += 2 + strlen(coldefs[k].col_label) + 3;
      columns = new char[totlen2 + 2];
      columns[0] = '\0';
      for (int k = 0; k < ncols; k++) {
        if (k) strcat(columns, ", ");
        strcat(columns, "\"");
        strcat(columns, coldefs[k].col_label);
        strcat(columns, "\"");
      }
      header_choice = &DumpCustom::header_item_tecplot;
      atom_header_choice = &DumpCustom::atom_header_item_tecplot;
      elem_header_choice = &DumpCustom::empty_header;
      node_header_choice = &DumpCustom::node_header_item_tecplot;

      pack_atom_choice = &DumpCustom::pack_atom_tecplot;
      pack_elem_choice = &DumpCustom::pack_elem_tecplot;
      pack_node_choice = &DumpCustom::pack_node_tecplot;

      convert_atom_choice = &DumpCustom::convert_string;
      convert_node_choice = &DumpCustom::convert_string;
      convert_elem_choice = &DumpCustom::convert_string_node_connect;

      if (buffer_flag == 1) {
        write_atom_choice = &DumpCustom::write_string;
        write_node_choice = &DumpCustom::write_string;
        write_elem_choice = &DumpCustom::write_string;
      } else {
        write_atom_choice = &DumpCustom::write_lines;
        write_elem_choice = &DumpCustom::write_lines_node_connect;
        write_node_choice = &DumpCustom::write_lines;
      }
#ifdef TECIO_ENABLED        
    } else {

      header_choice = &DumpCustom::empty_header;
      atom_header_choice = &DumpCustom::empty_header;
      elem_header_choice = &DumpCustom::empty_header;
      node_header_choice = &DumpCustom::empty_header;

      pack_atom_choice = &DumpCustom::pack_atom_tecplot_binary;
      pack_node_choice = &DumpCustom::pack_node_tecplot_binary;
      pack_elem_choice = &DumpCustom::pack_elem_tecplot_binary;

      // writing handled by tecio library

      // build TECIO variable string: "col0 col1 col2 ..."
      delete [] tec_values;
      int totlen = 0;
      for (int k = 0; k < ncols; k++) totlen += strlen(coldefs[k].col_label) + 1;
      tec_values = new char[totlen + 2];
      tec_values[0] = '\0';
      for (int k = 0; k < ncols; k++) {
        if (k) strcat(tec_values, " ");
        strcat(tec_values, coldefs[k].col_label);
      }
      delete [] tec_title;
      tec_title = new char[32];
      strcpy(tec_title, "Output file from CAC");

      delete [] valuelocation;
      valuelocation = new int[ncols];
      for (int k = 0; k < ncols; k++) valuelocation[k] = 1;  // all nodal
#endif
    }
  }

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpCustom::openfile()
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
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent) + 2];
        strcpy(nameslist[numfiles], filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent) + 1];
        strcpy(nameslist[fileidx], filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  if (filewriter) {
    if (dump_style == ASCII || dump_format == LAMMPS) {
      fp = fopen(filecurrent, "w");
      if (fp == nullptr) {
        char str[128];
        sprintf(str, "Cannot open dump file %s", filecurrent);
        error->one(FLERR, str);
      }
#ifdef TECIO_ENABLED        
    } else {
      int success = TECINI142(tec_title, tec_values, filecurrent,
          (char *)".", &dump_style, &dummy, &debug_tec, &dummy1);
      if (success == -1) {
        char str[128];
        sprintf(str, "Cannot open dump file %s", filecurrent);
        error->one(FLERR, str);
      }
#endif
    }

  } else fp = nullptr;

  if (multifile) delete [] filecurrent;
}
/* ---------------------------------------------------------------------- */

// Parse column keywords from arg[istart..narg-1].
// Called by each subclass constructor after consuming its own keywords.

void DumpCustom::parse_columns(int istart, int narg, char **arg)
{
  // Phase 1: expand wildcards — c_ID[*], c_ID[N*M], f_ID[*], etc.
  // Each raw keyword is passed to expand_wildcard_col; wildcarded entries are
  // replaced with explicit c_ID[1]...c_ID[N] sequences.
  std::vector<std::string> keyword_list;
  for (int i = istart; i < narg; i++)
    universe->expand_args(arg[i], keyword_list);

  // Phase 2: allocate storage, then parse the expanded keyword list.
  int nkeywords = (int) keyword_list.size();
  if (nkeywords < 1) error->all(FLERR, "Dump custom: no columns specified");

  // pre-count compute and fix references so we can allocate exactly once
  int ncomp_alloc = 0, nfix_alloc = 0;
  for (const auto &keyword : keyword_list) {
    if (keyword.size() >= 2 && keyword[0] == 'c' && keyword[1] == '_') ncomp_alloc++;
    if (keyword.size() >= 2 && keyword[0] == 'f' && keyword[1] == '_') nfix_alloc++;
  }
  if (nkeywords < 1) error->all(FLERR, "Dump custom: no columns specified");

  coldefs    = new ColDef[nkeywords];
  ncols = 0;

  id_compute = new char*[ncomp_alloc + 1];
  computes   = new Compute*[ncomp_alloc + 1];
  ncompute   = 0;

  id_fix = new char*[nfix_alloc + 1];
  fixes  = new Fix*[nfix_alloc + 1];
  nfix   = 0;

  for (const auto &keyword_str : keyword_list) {
    const char *keyword = keyword_str.c_str();
    ColDef cd;
    cd.compute = nullptr;
    cd.fix     = nullptr;
    cd.col     = 0;

    if (strcmp(keyword, "id") == 0) {
      cd.type = COL_ID;
      if (sizeof(tagint) == sizeof(smallint)) cd.vtype = Dump::INT;
      else cd.vtype = Dump::BIGINT;
    } else if (strcmp(keyword, "type") == 0) {
      cd.type = COL_TYPE;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "x") == 0) {
      cd.type = COL_X;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "y") == 0) {
      cd.type = COL_Y;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "z") == 0) {
      cd.type = COL_Z;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "vx") == 0) {
      cd.type = COL_VX;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "vy") == 0) {
      cd.type = COL_VY;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "vz") == 0) {
      cd.type = COL_VZ;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "fx") == 0) {
      cd.type = COL_FX;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "fy") == 0) {
      cd.type = COL_FY;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "fz") == 0) {
      cd.type = COL_FZ;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "proc") == 0) {
      cd.type = COL_PROC;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "procp1") == 0) {
      cd.type = COL_PROCP1;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "mass") == 0) {
      cd.type = COL_MASS;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "ix") == 0) {
      cd.type = COL_IX;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "iy") == 0) {
      cd.type = COL_IY;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "iz") == 0) {
      cd.type = COL_IZ;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "xu") == 0) {
      if (domain->triclinic) cd.type = COL_XU_TRI;
      else cd.type = COL_XU;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "yu") == 0) {
      if (domain->triclinic) cd.type = COL_YU_TRI;
      else cd.type = COL_YU;
      cd.vtype = Dump::DOUBLE;
    } else if (strcmp(keyword, "zu") == 0) {
      if (domain->triclinic) cd.type = COL_ZU_TRI;
      else cd.type = COL_ZU;
      cd.vtype = Dump::DOUBLE;
    } // compute columns: c_ID or c_ID[N] 
    else if (strncmp(keyword, "c_", 2) == 0) {
      char idbuf[256];
      int col_idx = -1;
      const char *bracket = strchr(keyword + 2, '[');
      if (bracket) {
        int len = bracket - (keyword + 2);
        strncpy(idbuf, keyword + 2, len);
        idbuf[len] = '\0';
        col_idx = atoi(bracket + 1) - 1;
        if (col_idx < 0) error->all(FLERR, "Dump custom: column index must be >= 1");
      } else {
        strcpy(idbuf, keyword + 2);
      }

      // reuse existing entry if same ID was already seen
      int ic = -1;
      for (int k = 0; k < ncompute; k++)
        if (strcmp(id_compute[k], idbuf) == 0) { ic = k; break; }
      if (ic < 0) {
        int n = strlen(idbuf) + 1;
        id_compute[ncompute] = new char[n];
        strcpy(id_compute[ncompute], idbuf);
        computes[ncompute] = nullptr;
        ic = ncompute++;
      }
      if (col_idx < 0) {
        cd.type = COL_COMPUTE_SCALAR;
        cd.col  = -(ic + 2);
      } else {
        cd.type = COL_COMPUTE_ARRAY;
        cd.col  = col_idx | (ic << 16);
      }
      cd.vtype = Dump::DOUBLE;
      // fix columns: f_ID or f_ID[N]
    } else if (strncmp(keyword, "f_", 2) == 0) {
      char idbuf[256];
      int col_idx = -1;
      const char *bracket = strchr(keyword + 2, '[');
      if (bracket) {
        int len = bracket - (keyword + 2);
        strncpy(idbuf, keyword + 2, len);
        idbuf[len] = '\0';
        col_idx = atoi(bracket + 1) - 1;
        if (col_idx < 0) error->all(FLERR, "Dump custom: fix column index must be >= 1");
      } else {
        strcpy(idbuf, keyword + 2);
      }
      int ifix = -1;
      for (int k = 0; k < nfix; k++)
        if (strcmp(id_fix[k], idbuf) == 0) { ifix = k; break; }
      if (ifix < 0) {
        int n = strlen(idbuf) + 1;
        id_fix[nfix] = new char[n];
        strcpy(id_fix[nfix], idbuf);
        fixes[nfix] = nullptr;
        ifix = nfix++;
      }
      if (col_idx < 0) {
        cd.type = COL_FIX_SCALAR;
        cd.col  = -(ifix + 2);
      } else {
        cd.type = COL_FIX_ARRAY;
        cd.col  = col_idx | (ifix << 16);
      }
      cd.vtype = Dump::DOUBLE;
    } 
    // debug columns
    else if (strcmp(keyword, "eid") == 0) {
      cd.type = COL_EID;
      if (sizeof(tagint) == sizeof(smallint)) cd.vtype = Dump::INT;
      else cd.vtype = Dump::BIGINT;
    } else if (strcmp(keyword, "etype") == 0) {
      cd.type = COL_ETYPE;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "ibasis") == 0) {
      cd.type = COL_IBASIS;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "igauss") == 0) {
      cd.type = COL_IGAUSS;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "inode") == 0) {
      cd.type = COL_INODE;
      cd.vtype = Dump::INT;
    } else if (strcmp(keyword, "iucell") == 0) {
      cd.type = COL_IUCELL;
      cd.vtype = Dump::INT;
    } else {
      error->all(FLERR, "Dump custom: unknown column keyword");
    }

    int n = strlen(keyword) + 1;
    cd.col_label = new char[n];
    strcpy(cd.col_label, keyword);
    coldefs[ncols] = cd;
    ncols++;
  }

  if (precision_flag) 
    for (int icol = 0; icol < ncols; icol++) 
      if (coldefs[icol].vtype == Dump::DOUBLE)
        coldefs[icol].vtype == Dump::DOUBLE16;
}

/* ---------------------------------------------------------------------- */

// Value accessor for a real atom.  All column types are handled here so the
// caller never needs to switch on ColType for real atoms.

double DumpCustom::get_atom_val(int i, const ColDef &cd)
{
  double *x = atom->x[i];
  double *v = atom->v[i];
  double *f = atom->f[i];
  tagint tag = atom->tag[i];
  int itype = atom->type[i];
  int image = atom->image[i];
  switch (cd.type) {
    case COL_ID:   return tag;
    case COL_TYPE: return itype;
    case COL_X:    return x[0];
    case COL_Y:    return x[1];
    case COL_Z:    return x[2];
    case COL_VX:   return v[0];
    case COL_VY:   return v[1];
    case COL_VZ:   return v[2];
    case COL_FX:   return f[0];
    case COL_FY:   return f[1];
    case COL_FZ:   return f[2];
    case COL_IX: return (image & IMGMASK) - IMGMAX;             
    case COL_IY: return (image >> IMGBITS & IMGMASK) - IMGMAX;    
    case COL_IZ: return (image >> IMG2BITS) - IMGMAX;         
    case COL_XU: return x[0] + ((image & IMGMASK) - IMGMAX) * domain->xprd;
    case COL_YU: return x[1] + ((image >> IMGBITS & IMGMASK) - IMGMAX) * domain->yprd;
    case COL_ZU: return x[2] + ((image >> IMG2BITS) - IMGMAX) * domain->zprd;
    case COL_XU_TRI: 
                 {
                   double *h = domain->h;
                   int xbox = (image & IMGMASK) - IMGMAX;
                   int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
                 }
    case COL_YU_TRI: 
                 {
                   double *h = domain->h;
                   int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[1] + h[1]*ybox + h[3]*zbox;
                 }
    case COL_ZU_TRI: 
                 {
                   double *h = domain->h;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[2] + h[2]*zbox;
                 }
    case COL_MASS: return atom->mass[itype];
    case COL_PROC: return me;
    case COL_PROCP1: return (me + 1);
    case COL_COMPUTE_SCALAR: return cd.compute->vector_atom[i];
    case COL_COMPUTE_ARRAY:  return cd.compute->array_atom[i][cd.col];
    case COL_FIX_SCALAR:     return cd.fix->vector_atom[i];
    case COL_FIX_ARRAY:      return cd.fix->array_atom[i][cd.col];

                             // real atoms return -1 for all debug columns
    case COL_EID: case COL_ETYPE: case COL_IBASIS:
    case COL_IGAUSS: case COL_INODE: case COL_IUCELL:
                             return -1.0;
    default:
                             {
                               error->all(FLERR, "Dump custom: unrecognized column type in get_atom_val");
                               return 0.0;
                             }
  }
}

/* ---------------------------------------------------------------------- */

// Value accessor for a virtual atom (one ucell point of a CG element).
//
//   ibasis   = basis index within the element
//   iucell   = ucell (interpolation point) index
//   igcell   = gauss-cell index for this point (-1 if not a gauss point)
//   inode    = node index for this point      (-1 if not a node)
//   vatom_id = virtual atom ID to emit for COL_ID
//   coord    = interpolated position [3] to emit for COL_X/Y/Z
//
// vx/vy/vz and fx/fy/fz are interpolated from the element nodal arrays.


double DumpCustom::get_vatom_val(int i, int ibasis, int iucell, int igcell, int inode,
    tagint vatom_id, const double *x, const ColDef &cd)
{
  int image = element->image[i];
  switch (cd.type) {
    // position and identity — supplied by caller
    case COL_ID:   return vatom_id;
    case COL_TYPE: return element->ctype[i][ibasis];
    case COL_X:    return x[0];
    case COL_Y:    return x[1];
    case COL_Z:    return x[2];
                   // velocity/force interpolated from FE nodes
    case COL_VX: return element->evec->interpolate(element->nodev, i, ibasis, iucell, 0);
    case COL_VY: return element->evec->interpolate(element->nodev, i, ibasis, iucell, 1);
    case COL_VZ: return element->evec->interpolate(element->nodev, i, ibasis, iucell, 2);
    case COL_FX: return element->evec->interpolate(element->nodef, i, ibasis, iucell, 0);
    case COL_FY: return element->evec->interpolate(element->nodef, i, ibasis, iucell, 1);
    case COL_FZ: return element->evec->interpolate(element->nodef, i, ibasis, iucell, 2);
    case COL_IX: return (image & IMGMASK) - IMGMAX;             
    case COL_IY: return (image >> IMGBITS & IMGMASK) - IMGMAX;    
    case COL_IZ: return (image >> IMG2BITS) - IMGMAX;         
    case COL_XU: return x[0] + ((image & IMGMASK) - IMGMAX) * domain->xprd;
    case COL_YU: return x[1] + ((image >> IMGBITS & IMGMASK) - IMGMAX) * domain->yprd;
    case COL_ZU: return x[2] + ((image >> IMG2BITS) - IMGMAX) * domain->zprd;
    case COL_XU_TRI: 
                 {
                   double *h = domain->h;
                   int xbox = (image & IMGMASK) - IMGMAX;
                   int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
                 }
    case COL_YU_TRI: 
                 {
                   double *h = domain->h;
                   int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[1] + h[1]*ybox + h[3]*zbox;
                 }
    case COL_ZU_TRI: 
                 {
                   double *h = domain->h;
                   int zbox = (image >> IMG2BITS) - IMGMAX;
                   return x[2] + h[2]*zbox;
                 }

    case COL_MASS: return atom->mass[element->ctype[i][ibasis]];
    case COL_PROC: return me;
    case COL_PROCP1: return (me + 1);
                     // compute data
    case COL_COMPUTE_SCALAR:
                     {
                       if (cd.compute->element_data_style == Compute::NODE)
                         return element->evec->interpolate(cd.compute->vector_node, i, ibasis, iucell);
                       return cd.compute->vector_vatom[i][ibasis][iucell];
                     }
    case COL_COMPUTE_ARRAY:
                     {
                       if (cd.compute->element_data_style == Compute::NODE)
                         return element->evec->interpolate(cd.compute->array_node, i, ibasis, iucell, cd.col);
                       return cd.compute->array_vatom[i][ibasis][iucell][cd.col];
                     }
    case COL_FIX_SCALAR: 
                     {
                       if (cd.fix->element_data_style == Compute::NODE)
                         return element->evec->interpolate(cd.fix->vector_node, i, ibasis, iucell);
                       return cd.fix->vector_vatom[i][ibasis][iucell];
                     }
    case COL_FIX_ARRAY:
                     {
                       if (cd.fix->element_data_style == Compute::NODE)
                         return element->evec->interpolate(cd.fix->array_node, i, ibasis, iucell, cd.col);
                       return cd.fix->array_vatom[i][ibasis][iucell][cd.col];
                     }

                     // debug columns
    case COL_EID:        return element->tag[i];
    case COL_ETYPE:      return element->etype[i];
    case COL_IBASIS:     return ibasis;
    case COL_IGAUSS:     return igcell;  
    case COL_INODE:      return inode;  
    case COL_IUCELL:     return iucell;
    default:
                         {
                           error->all(FLERR, "Dump custom: unrecognized column type in get_vatom_val");
                           return 0.0;  // unreachable
                         }
  }
}

/* ---------------------------------------------------------------------- */
// FE node accessor: gets column value for node (element i, basis j, node inode).

double DumpCustom::get_node_val(int i, int ibasis, int inode, const ColDef &cd)
{
  int ietype = element->etype[i];
  double *x = element->nodex[i][ibasis][inode];
  int image = element->image[i];
  switch (cd.type) {
    case COL_ID:   return element->tag[i];
    case COL_TYPE: return element->ctype[i][ibasis];
    case COL_X:    return x[0];
    case COL_Y:    return x[1];
    case COL_Z:    return x[2];
    case COL_VX:   return element->nodev[i][ibasis][inode][0];
    case COL_VY:   return element->nodev[i][ibasis][inode][1];
    case COL_VZ:   return element->nodev[i][ibasis][inode][2];
    case COL_FX:   return element->nodef[i][ibasis][inode][0];
    case COL_FY:   return element->nodef[i][ibasis][inode][1];
    case COL_FZ:   return element->nodef[i][ibasis][inode][2];
    case COL_IX:   return (image & IMGMASK) - IMGMAX;             
    case COL_IY:   return (image >> IMGBITS & IMGMASK) - IMGMAX;    
    case COL_IZ:   return (image >> IMG2BITS) - IMGMAX;         
    case COL_XU:   return x[0] + ((image & IMGMASK) - IMGMAX) * domain->xprd;
    case COL_YU:   return x[1] + ((image >> IMGBITS & IMGMASK) - IMGMAX) * domain->yprd;
    case COL_ZU:   return x[2] + ((image >> IMG2BITS) - IMGMAX) * domain->zprd;
    case COL_XU_TRI: 
                   {
                     double *h = domain->h;
                     int xbox = (image & IMGMASK) - IMGMAX;
                     int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                     int zbox = (image >> IMG2BITS) - IMGMAX;
                     return x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
                   }
    case COL_YU_TRI: 
                   {
                     double *h = domain->h;
                     int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
                     int zbox = (image >> IMG2BITS) - IMGMAX;
                     return x[1] + h[1]*ybox + h[3]*zbox;
                   }
    case COL_ZU_TRI: 
                   {
                     double *h = domain->h;
                     int zbox = (image >> IMG2BITS) - IMGMAX;
                     return x[2] + h[2]*zbox;
                   }
    case COL_MASS:   return atom->mass[element->ctype[i][ibasis]];
    case COL_PROC:   return me;
    case COL_PROCP1: return (me + 1);

    case COL_COMPUTE_SCALAR:
                     {
                       if (cd.compute->element_data_style == Compute::NODE)
                         return cd.compute->vector_node[i][ibasis][inode];
                       return cd.compute->vector_vatom[i][ibasis][element->n2u[ietype][inode]];
                     }
    case COL_COMPUTE_ARRAY:
                     {
                       if (cd.compute->element_data_style == Compute::NODE)
                         return cd.compute->array_node[i][ibasis][inode][cd.col];
                       return cd.compute->array_vatom[i][ibasis][element->n2u[ietype][inode]][cd.col];
                     }
    case COL_FIX_SCALAR:
                     {
                       if (cd.fix->element_data_style == Compute::NODE)
                         return cd.fix->vector_node[i][ibasis][inode];
                       return cd.compute->vector_vatom[i][ibasis][element->n2u[ietype][inode]];
                     }
    case COL_FIX_ARRAY:
                     {
                       if (cd.fix->element_data_style == Compute::NODE)
                         return cd.fix->array_node[i][ibasis][inode][cd.col];
                       return cd.fix->array_vatom[i][ibasis][element->n2u[ietype][inode]][cd.col];
                     }

                     // debug columns for FE nodes
    case COL_EID:        return element->tag[i];
    case COL_ETYPE:      return ietype;
    case COL_IBASIS:     return ibasis;
    case COL_INODE:      return inode;
    case COL_IGAUSS:     return element->n2g[ietype][inode];
    case COL_IUCELL:     return element->n2u[ietype][inode];
    default:
                         {
                           error->all(FLERR, "Dump custom/tecplot: unrecognized column type in get_node_val");
                           return 0.0;  // unreachable
                         }
  }
}

/* ---------------------------------------------------------------------- */
// Override write() to handle Tecplot format
// Call Dump::write() to handle LAMMPS format

void DumpCustom::write()
{
  // invoke computes once for this timestep
  for (int icompute = 0; icompute < ncompute; icompute++) {
    Compute *compute = computes[icompute];
    if (compute->invoked_peratom != update->ntimestep) {
      compute->compute_peratom();
      compute->addstep(update->ntimestep + nevery);
    }
  }

  if (dump_format == LAMMPS) {
    Dump::write();
  } else {
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


    int nmax, tmp, nlines, nchars, nsmin, nsmax;
    MPI_Status status;
    MPI_Request request;

    if (dump_style == ASCII) {
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
            write_atom_data(nchars, atom_size_one, (double *)sbuf);
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
            write_node_data(nchars, node_size_one, (double *)sbuf);
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
            write_elem_data(nchars, elem_size_one, (double *)sbuf);
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
    }
#ifdef TECIO_ENABLED        
    else {
      // ---- Binary TECIO ----

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
          current_col = k;
          pack_atom(nullptr);
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
          current_col = k;
          pack_node(nullptr);
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

      if (multifile && filewriter) 
        TECEND142();
    }
#endif

  }
}

/*-------------------------------------------------------------------------------------------*/

void DumpCustom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (filewriter) (this->*header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpCustom::write_atom_header(bigint ndump)
{
  if (multiproc) (this->*atom_header_choice)(ndump);
  else if (filewriter) (this->*atom_header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpCustom::write_elem_header(bigint ndump)
{
  if (multiproc) (this->*elem_header_choice)(ndump);
  else if (filewriter) (this->*elem_header_choice)(ndump);
}

/*-------------------------------------------------------------------------------------------*/

void DumpCustom::write_node_header(bigint ndump)
{
  if (multiproc) (this->*node_header_choice)(ndump);
  else if (filewriter) (this->*node_header_choice)(ndump);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCustom::pack_atom(tagint *ids)
{
  (this->*pack_atom_choice)(ids);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCustom::pack_elem(tagint *ids)
{
  (this->*pack_elem_choice)(ids);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCustom::pack_node(tagint *ids)
{
  (this->*pack_node_choice)(ids);
}

/*--------------------------------------------------------------------------------------------*/

int DumpCustom::convert_atom_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_atom_choice)(n, size_one, mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpCustom::convert_node_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_node_choice)(n, size_one, mybuf);
}

/*--------------------------------------------------------------------------------------------*/

int DumpCustom::convert_elem_string(int n, int size_one, double *mybuf)
{
  return (this->*convert_elem_choice)(n, size_one, mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpCustom::write_node_data(int n, int size_one, double *mybuf)
{
  (this->*write_node_choice) (n, size_one, mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpCustom::write_atom_data(int n, int size_one, double *mybuf)
{
  (this->*write_atom_choice) (n, size_one, mybuf);
}

/*----------------------------------------------------------------------------------------------------*/

void DumpCustom::write_elem_data(int n, int size_one, double *mybuf)
{
  (this->*write_elem_choice) (n, size_one, mybuf);
}

/* ----------------------------------------------------------------------------------------- */

void DumpCustom::pack_atom_lammps(tagint *ids)
{
  int *mask   = atom->mask;
  double **x  = atom->x;
  int *type   = atom->type;
  tagint *tag = atom->tag;
  int nlocal  = atom->nlocal;

  // real atoms — all columns go through get_atom_val
  tagint maxtag = 0;
  int m = 0, n = 0;
  for (int i = 0; i < nlocal; i++) {
    maxtag = MAX(maxtag, tag[i]);
    if (!(mask[i] & groupbit)) continue;
    for (int k = 0; k < ncols; k++)
      buf[m++] = get_atom_val(i, coldefs[k]);
    if (ids) ids[n++] = i;
  }

  tagint maxtag_all = 0;
  if (atom->natoms)
    MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_CAC_TAGINT, MPI_MAX, world);

  // virtual atoms from CG elements
  int *emask       = element->mask;
  int *etype       = element->etype;
  int *apc         = element->apc;
  int *nucell      = element->nucell;
  int **ctype      = element->ctype;
  tagint *etag     = element->tag;
  double **ex      = element->x;
  double ****nodex = element->nodex;
  int maxucell     = element->maxucell;
  int maxapc       = element->maxapc;
  ElementVec *evec = element->evec;
  int nelocal      = element->nlocal;

  int **is_outer = (dump_style == SURFACE)      ? element->is_outer : nullptr;
  int *ngcell    = (dump_style == GAUSSPOINT) ? element->ngcell   : nullptr;
  int **g2u      = (dump_style == GAUSSPOINT) ? element->g2u      : nullptr;
  int maxgcell   = (dump_style == GAUSSPOINT) ? element->maxgcell : 0;
  int *npe       = element->npe;
  int **n2u      = element->n2u;
  int maxnpe     = element->maxnpe;
  int **n2g      = element->n2g;   // node  → gauss-cell index (-1 if not a gauss point)
  int **g2n      = element->g2n;   // gauss → node index       (-1 if not at a node)
  int **u2g      = element->u2g;   // ucell → gauss-cell index (-1 if not a gauss point)

  double coord[3];
  tagint itag = maxtag_all + 1;

  for (int i = 0; i < nelocal; i++) {
    if (!(emask[i] & groupbit)) continue;
    int ietype = etype[i];
    int iapc   = apc[ietype];

    if (dump_style == CENTER) {
      // one virtual atom per element, placed at the element centroid
      tagint vatom_id = itag + etag[i] - 1;
      int image = element->image[i];
      for (int k = 0; k < ncols; k++) {
        ColType ct = coldefs[k].type;
        if      (ct == COL_ID)          buf[m++] = vatom_id;
        else if (ct == COL_TYPE)        buf[m++] = 0;
        else if (ct == COL_X)           buf[m++] = ex[i][0];
        else if (ct == COL_Y)           buf[m++] = ex[i][1];
        else if (ct == COL_Z)           buf[m++] = ex[i][2];
        else if (ct == COL_IX)          buf[m++] = (image & IMGMASK) - IMGMAX;
        else if (ct == COL_IY)          buf[m++] = (image >> IMGBITS & IMGMASK) - IMGMAX;
        else if (ct == COL_IZ)          buf[m++] = (image >> IMG2BITS) - IMGMAX;
        else if (ct == COL_XU)          buf[m++] = ex[i][0] + ((image & IMGMASK) - IMGMAX) * domain->xprd;
        else if (ct == COL_YU)          buf[m++] = ex[i][1] + ((image >> IMGBITS & IMGMASK) - IMGMAX) * domain->yprd;
        else if (ct == COL_ZU)          buf[m++] = ex[i][2] + ((image >> IMG2BITS) - IMGMAX) * domain->zprd;
        else if (ct == COL_XU_TRI) {
          double *h = domain->h;
          int xbox = (image & IMGMASK) - IMGMAX;
          int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
          int zbox = (image >> IMG2BITS) - IMGMAX;
          buf[m++] = ex[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        } else if (ct == COL_YU_TRI) {
          double *h = domain->h;
          int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
          int zbox = (image >> IMG2BITS) - IMGMAX;
          buf[m++] = ex[i][1] + h[1]*ybox + h[3]*zbox;
        } else if (ct == COL_ZU_TRI) {
          double *h = domain->h;
          int zbox = (image >> IMG2BITS) - IMGMAX;
          buf[m++] =  ex[i][2] + h[2]*zbox;
        }
        else if (ct == COL_PROC)        buf[m++] = me;
        else if (ct == COL_PROCP1)      buf[m++] = me + 1;
        else if (ct == COL_EID)         buf[m++] = etag[i];
        else if (ct == COL_ETYPE)       buf[m++] = ietype;
        else if (ct == COL_IBASIS)      buf[m++] = -1;
        else if (ct == COL_IGAUSS)      buf[m++] = -1;
        else if (ct == COL_INODE)       buf[m++] = -1;
        else if (ct == COL_IUCELL)      buf[m++] = -1;
        else buf[m++] = 0.0;
      }
      if (ids) ids[n++] = i;

    } else if (dump_style == GAUSSPOINT) {
      int ingcell = ngcell[ietype];
      for (int ibasis = 0; ibasis < iapc; ibasis++) {
        for (int igcell = 0; igcell < ingcell; igcell++) {
          int iucell = g2u[ietype][igcell];
          int inode  = g2n[ietype][igcell]; 
          tagint vatom_id = itag + ((etag[i]-1)*maxgcell + igcell)*maxapc + ibasis;
          evec->interpolate(coord, nodex, i, ibasis, iucell, 3);
          if (wrap_flag) domain->remap(coord);
          for (int k = 0; k < ncols; k++)
            buf[m++] = get_vatom_val(i, ibasis, iucell, igcell, inode,
                vatom_id, coord, coldefs[k]);
          if (ids) ids[n++] = i;
        }
      }

    } else if (dump_style == NODE) {
      int inpe = npe[ietype];
      for (int ibasis = 0; ibasis < iapc; ibasis++) {
        for (int inode = 0; inode < inpe; inode++) {
          int iucell  = n2u[ietype][inode];
          int igcell  = n2g[ietype][inode];  // gauss at this node (-1 if none)
          tagint vatom_id = itag + ((etag[i]-1)*maxnpe + inode)*maxapc + ibasis;
          coord[0] = nodex[i][ibasis][inode][0];
          coord[1] = nodex[i][ibasis][inode][1];
          coord[2] = nodex[i][ibasis][inode][2];
          if (wrap_flag) domain->remap(coord);
          for (int k = 0; k < ncols; k++)
            buf[m++] = get_vatom_val(i, ibasis, iucell, igcell, inode,
                vatom_id, coord, coldefs[k]);
          if (ids) ids[n++] = i;
        }
      }

    } else {
      // FULLMAP or SURFACE
      int inpe = npe[ietype];
      for (int ibasis = 0; ibasis < iapc; ibasis++) {
        for (int iucell = 0; iucell < nucell[ietype]; iucell++) {
          if (dump_style == SURFACE && !is_outer[ietype][iucell]) continue;
          int inode = -1;
          for (int k = 0; k < inpe; k++)
            if (n2u[ietype][k] == iucell) { 
              inode = k; 
              break; 
            }
          int igcell = u2g[ietype][iucell];  // gauss at this ucell (-1 if none)
          tagint vatom_id = itag + ((etag[i]-1)*maxucell + iucell)*maxapc + ibasis;
          evec->interpolate(coord, nodex, i, ibasis, iucell, 3);
          if (wrap_flag) domain->remap(coord);
          for (int k = 0; k < ncols; k++)
            buf[m++] = get_vatom_val(i, ibasis, iucell, igcell, inode,
                vatom_id, coord, coldefs[k]);
          if (ids) ids[n++] = i;
        }
      }
    }
  }
}

void DumpCustom::pack_atom_tecplot(tagint *ids)
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

#ifdef TECIO_ENABLED        
void DumpCustom::pack_atom_tecplot_binary(tagint *ids)
{
  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    buf[m++] = get_atom_val(i, coldefs[current_col]);
  }
}

#endif

/*------------------------------------------------------------------------------------------*/

void DumpCustom::pack_elem_tecplot(tagint *ids)
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

  int *node_connect = new int[npe_connect];
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    int ietype = etype[i];
    if (average_flag) {
      element->node_connectivity(node_connect, element_shape_ids[ietype], nodecell_ids[i]);
      for (int j = 0; j < npe_connect; j++) 
        buf[m++] = node_connect[j];
    } else {
      for (int j = 0; j < apc[ietype]; j++) {
        element->node_connectivity(node_connect, element_shape_ids[ietype], node_ids[i][j]);
        for (int k = 0; k < npe_connect; k++) 
          buf[m++] = node_connect[k];
      }
    }
  }
  delete [] node_connect;
}

/*------------------------------------------------------------------------------------------*/

#ifdef TECIO_ENABLED        
void DumpCustom::pack_elem_tecplot_binary(tagint *ids)
{
  int *mask = element->mask;
  int nlocal = element->nlocal;

  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;
  int itype;

  if (maxelem_buf < nlocal) {
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

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      itype = etype[i];
      if (average_flag) {
        m += element->node_connectivity(&ibuf[m], element_shape_ids[itype], nodecell_ids[i]);
      } else {
        for (int j = 0; j < apc[itype]; j++) {
          m += element->node_connectivity(&ibuf[m], element_shape_ids[itype], node_ids[i][j]);
        }
      }
    }
  }
} 

#endif
/*------------------------------------------------------------------------------------------*/

void DumpCustom::pack_node_tecplot(tagint *ids)
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
          for (int ibasis = 0; ibasis < iapc; ibasis++)
            val += get_node_val(i, ibasis, inode, coldefs[k]);
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
      for (int ibasis = 0; ibasis < iapc; ibasis++) 
        for (int inode = 0; inode < inpe; inode++) 
          for (int k = 0; k < ncols; k++)
            buf[m++] = get_node_val(i, ibasis, inode, coldefs[k]);
      if (ids) ids[n++] = i;
    }
  }
}

/* ------------------------------------------------------------------------------------------ */

#ifdef TECIO_ENABLED        
void DumpCustom::pack_node_tecplot_binary(tagint *ids)
{
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *npe    = element->npe;
  int *etype  = element->etype;
  int m = 0, n = 0;

  if (average_flag) {
    for (int i = 0; i < nelocal; i++) {
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int inode = 0; inode < inpe; inode++) {
        double val = 0.0;
        for (int ibasis = 0; ibasis < iapc; ibasis++)
          val += get_node_val(i, ibasis, inode, coldefs[current_col]);
        buf[m++] = val / iapc;
      }
      if (ids) ids[n++] = i;
    }
  } else {
    for (int i = 0; i < nelocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      int iapc = apc[etype[i]];
      int inpe = npe[etype[i]];
      for (int ibasis = 0; ibasis < iapc; ibasis++)
        for (int inode = 0; inode < inpe; inode++)
          buf[m++] = get_node_val(i, ibasis, inode, coldefs[current_col]);
      if (ids) ids[n++] = i;
    }
  }
}
#endif

/*--------------------------------------------------------------------------------------------*/

int DumpCustom::convert_string(int n, int size_one, double *mybuf)
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
      if (coldefs[k].vtype == Dump::INT)
        offset += sprintf(&sbuf[offset], "%d", static_cast<int>(mybuf[m++]));
      else if (coldefs[k].vtype == Dump::BIGINT)
        offset += sprintf(&sbuf[offset], BIGINT_FORMAT, static_cast<bigint>(mybuf[m++]));
      else if (coldefs[k].vtype == Dump::DOUBLE)
        offset += sprintf(&sbuf[offset], "%g", mybuf[m++]);
      else if (coldefs[k].vtype == Dump::DOUBLE16)
        offset += sprintf(&sbuf[offset], "%-1.16e", mybuf[m++]);
    }
    offset += sprintf(&sbuf[offset], "\n");
  }
  return offset;
}

/*--------------------------------------------------------------------------------------------*/

int DumpCustom::convert_string_node_connect(int n, int size_one, double *mybuf)
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
      offset += sprintf(&sbuf[offset], TAGINT_FORMAT, static_cast<tagint>(mybuf[m++]));
    }
    offset += sprintf(&sbuf[offset], "\n");
  }
  return offset;
}

/*---------------------------------------------------------------------------------------------------*/

void DumpCustom::write_string(int n, int size_one, double *mybuf)
{
  fwrite(mybuf, sizeof(char), n, fp);
}
/*-------------------------------------------------------------------------------------------------*/

void DumpCustom::write_lines(int n, int size_one, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < size_one; k++) {
      if (k) fprintf(fp, " ");
      if (coldefs[k].vtype == Dump::INT)
        fprintf(fp, "%d", static_cast<int>(mybuf[m++]));
      else if (coldefs[k].vtype == Dump::BIGINT)
        fprintf(fp, BIGINT_FORMAT, static_cast<bigint>(mybuf[m++]));
      else if (coldefs[k].vtype == Dump::DOUBLE)
        fprintf(fp, "%g", (mybuf[m++]));
      else if (coldefs[k].vtype == Dump::DOUBLE16)
        fprintf(fp, "%-1.16e", (mybuf[m++]));
    }

    fprintf(fp, "\n");
  }
}

/*----------------------------------------------------------------------------------------------------*/

void DumpCustom::write_lines_node_connect(int n, int size_one, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < size_one; k++) {
      if (k) fprintf(fp, " ");
      fprintf(fp, TAGINT_FORMAT, static_cast<tagint>(mybuf[m++]));
    }
    fprintf(fp, "\n");
  }
}


/* ------------------------------------------------------------------------------------------ */

void DumpCustom::header_item_lammps(bigint ndump)
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

void DumpCustom::header_item_lammps_triclinic(bigint ndump)
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

/*------------------------------------------------------------------------------------------*/

void DumpCustom::header_item_tecplot(bigint ndump)
{
  fprintf(fp, "title = \"CAC Simulation\", ");
  fprintf(fp, "variables = %s\n", columns);
}

/*------------------------------------------------------------------------------------------*/

void DumpCustom::node_header_item_tecplot(bigint ndump)
{
  fprintf(fp, "zone t = \"Coarse Element\"");
  fprintf(fp, " n = " BIGINT_FORMAT " e = " BIGINT_FORMAT, ndump, nelem_total);
  if (domain->dimension == 3) 
    fprintf(fp, " datapacking = point, zonetype = febrick\n");
  else 
    fprintf(fp, " datapacking = point, zonetype = fequadrilateral\n");
}

/*------------------------------------------------------------------------------------------*/

void DumpCustom::atom_header_item_tecplot(bigint ndump)
{
  fprintf(fp, "zone t = \"Discrete Atoms\", datapacking = point\n");
}

/*------------------------------------------------------------------------------------------*/
void DumpCustom::grow_buf(int num_me)
{
  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&num_me, &nmax, 1, MPI_INT, MPI_MAX, world);
  else nmax = num_me;
  if (nmax > maxbuf) {
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf, maxbuf, "dump:buf");
  }
}

/*------------------------------------------------------------------------------------------*/
#ifdef TECIO_ENABLED        
void DumpCustom::comm_buf_tecio(int num_me)
{
  int nmax, tmp;
  MPI_Status status;
  MPI_Request request;

  if (multiproc != nprocs) MPI_Allreduce(&num_me, &nmax, 1, MPI_INT, MPI_MAX, world);
  else nmax = num_me;

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      int nvals;
      if (iproc) {
        MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &nvals);
      } else nvals = num_me;
      int success = TECDAT142(&nvals, buf, &dummy1);
      if (success == -1) error->one(FLERR, "Cannot write to TecPlot dump file");
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, fileproc, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(buf, num_me, MPI_DOUBLE, fileproc, 0, world);
  }
}
#endif

/* ---------------------------------------------------------------------- */

void DumpCustom::write_box_zone()
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

  if (dump_style == ASCII) {
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
  }
#ifdef TECIO_ENABLED        
  else {
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
#endif
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count_atoms()
{
  int count = 0;
  int *amask  = atom->mask;
  int nalocal = atom->nlocal;
  for (int i = 0; i < nalocal; i++)
    if (amask[i] & groupbit) count++;

  if (dump_format == LAMMPS) {
    int *emask  = element->mask;
    int nelocal = element->nlocal;
    int *etype  = element->etype;
    int *apc_e  = element->apc;
    int *nucell = element->nucell;

    if (dump_style == FULLMAP) {
      for (int i = 0; i < nelocal; i++)
        if (emask[i] & groupbit)
          count += nucell[etype[i]] * apc_e[etype[i]];
    } else if (dump_style == SURFACE) {
      int **is_outer = element->is_outer;
      for (int i = 0; i < nelocal; i++) {
        if (!(emask[i] & groupbit)) continue;
        int ietype = etype[i];
        int iapc   = apc_e[ietype];
        for (int iucell = 0; iucell < nucell[ietype]; iucell++)
          if (is_outer[ietype][iucell]) count += iapc;
      }
    } else if (dump_style == GAUSSPOINT) {
      int *ngcell = element->ngcell;
      for (int i = 0; i < nelocal; i++)
        if (emask[i] & groupbit)
          count += ngcell[etype[i]] * apc_e[etype[i]];
    } else if (dump_style == NODE) {
      int *npe = element->npe;
      for (int i = 0; i < nelocal; i++)
        if (emask[i] & groupbit)
          count += npe[etype[i]] * apc_e[etype[i]];
    } else {  // CENTER
      for (int i = 0; i < nelocal; i++)
        if (emask[i] & groupbit) count++;
    }
  }
  return count;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count_nodes()
{
  int count = 0;
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *npe    = element->npe;
  int *etype  = element->etype;
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    count += npe[etype[i]] * (average_flag ? 1 : apc[etype[i]]);
  }
  return count;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count_elements()
{
  int count = 0;
  int *mask   = element->mask;
  int nelocal = element->nlocal;
  int *apc    = element->apc;
  int *etype  = element->etype;
  for (int i = 0; i < nelocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    count += (average_flag ? 1 : apc[etype[i]]);
  }
  return count;
}

