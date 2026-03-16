#include "cactype.h" 
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "read_data.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "element_vec.h"
#include "universe.h"
#include "group.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "pair.h"
#include "universe.h"
#include "error.h"
#include "memory.h"

using namespace CAC_NS;

#define MAXLINE 256
#define LB_FACTOR 1.2
#define CHUNK 1024
#define DELTA 4            // must be 2 or larger
#define MAXBODY 20         // max # of lines in one body, also in Atom class
#define NSECTIONS 8       // change when add to header::section_keywords
#define EPSILON 1e-6

enum{NONE,APPEND,VALUE,MERGE};

/* ---------------------------------------------------------------------- */

ReadData::ReadData(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

}

/* ---------------------------------------------------------------------- */

ReadData::~ReadData()
{
  delete [] line;
  delete [] keyword;
  delete [] style;
  delete [] buffer;
  memory->sfree(arg);

  for (int i = 0; i < nfix; i++) {
    delete [] fix_header[i];
    delete [] fix_section[i];
  }
  memory->destroy(fix_index);
  memory->sfree(fix_header);
  memory->sfree(fix_section);
}

/* ---------------------------------------------------------------------- */

void ReadData::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal read_data command");

  addflag = NONE;
  aid_offset = eid_offset = 0;
  offsetflag = shiftflag = 0;
  toffset = eoffset = 0;
  extra_atom_types = extra_element_types = 0;
  shift[0] = shift[1] = shift[2] = 0.0;

  nelements = natoms = 0;
  groupbit = 0;

  nfix = 0;
  fix_index = NULL;
  fix_header = NULL;
  fix_section = NULL;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"add") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      if (strcmp(arg[iarg+1],"append") == 0) addflag = APPEND;
      else if (strcmp(arg[iarg+1],"merge") == 0) addflag = MERGE;
      else {
        error->all(FLERR,"Illegal read_data command");

        //if (atom->molecule_flag && (iarg+3 > narg))
        //  error->all(FLERR,"Illegal read_data command");
        //addflag = VALUE;
        //bigint offset = universe->bnumeric(FLERR,arg[iarg+1]);
        //if (offset > MAXTAGINT)
        //  error->all(FLERR,"Read data add atomID offset is too big");
        //id_offset = offset;

        //if (atom->molecule_flag) {
        //  offset = universe->bnumeric(FLERR,arg[iarg+2]);
        //  if (offset > MAXTAGINT)
        //    error->all(FLERR,"Read data add molID offset is too big");
        //  mol_offset = offset;
        //  iarg++;
        //}
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal read_data command");
      offsetflag = 1;
      toffset = universe->inumeric(FLERR,arg[iarg+1]);
      eoffset = universe->inumeric(FLERR,arg[iarg+2]);
      if (toffset < 0 || eoffset < 0)
        error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal read_data command");
      shiftflag = 1;
      shift[0] = universe->numeric(FLERR,arg[iarg+1]);
      shift[1] = universe->numeric(FLERR,arg[iarg+2]);
      shift[2] = universe->numeric(FLERR,arg[iarg+3]);
      if (domain->dimension == 2 && shift[2] != 0.0)
        error->all(FLERR,"Non-zero read_data shift z value for 2d simulation");
      iarg += 4;
    //} else if (strcmp(arg[iarg],"nocoeff") == 0) {
    //  coeffflag = 0;
    //  iarg ++;
    } else if (strcmp(arg[iarg],"extra/atom/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      extra_atom_types = universe->inumeric(FLERR,arg[iarg+1]);
      if (extra_atom_types < 0) error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra/element/types") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      extra_element_types = universe->inumeric(FLERR,arg[iarg+1]);
      if (extra_element_types < 0) error->all(FLERR,"Illegal read_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal read_data command");
      int igroup = group->find_or_create(arg[iarg+1]);
      groupbit = group->bitmask[igroup];
      iarg += 2;
    //} else if (strcmp(arg[iarg],"fix") == 0) {
    //  if (iarg+4 > narg)
    //    error->all(FLERR,"Illegal read_data command");
    //  memory->grow(fix_index,nfix+1,"read_data:fix_index");
    //  fix_header = (char **)
    //    memory->srealloc(fix_header,(nfix+1)*sizeof(char *),
    //                     "read_data:fix_header");
    //  fix_section = (char **)
    //    memory->srealloc(fix_section,(nfix+1)*sizeof(char *),
    //                     "read_data:fix_section");
    //  fix_index[nfix] = modify->find_fix(arg[iarg+1]);
    //  if (fix_index[nfix] < 0)
    //    error->all(FLERR,"Fix ID for read_data does not exist");
    //  if (strcmp(arg[iarg+2],"NULL") == 0) fix_header[nfix] = NULL;
    //  else {
    //    int n = strlen(arg[iarg+2]) + 1;
    //    fix_header[nfix] = new char[n];
    //    strcpy(fix_header[nfix],arg[iarg+2]);
    //  }
    //  int n = strlen(arg[iarg+3]) + 1;
    //  fix_section[nfix] = new char[n];
    //  strcpy(fix_section[nfix],arg[iarg+3]);
    //  nfix++;
    //  iarg += 4;

    } else error->all(FLERR,"Illegal read_data command");
  }

  // error checks

  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (domain->box_exist && !addflag)
    error->all(FLERR,"Cannot read_data without add keyword "
               "after simulation box is defined");
  if (!domain->box_exist && addflag)
    error->all(FLERR,"Cannot use read_data add before "
               "simulation box is defined");
  if (offsetflag && addflag == NONE)
    error->all(FLERR,"Cannot use read_data offset without add flag");
  if (shiftflag && addflag == NONE)
    error->all(FLERR,"Cannot use read_data shift without add flag");
  if (addflag != NONE && (extra_atom_types || extra_element_types))
    error->all(FLERR,"Cannot use read_data extra with add flag");

  // first time system initialization

  if (addflag == NONE) {
    domain->box_exist = 1;
    update->ntimestep = 0;
  }

  // compute atomID and elementID offset for addflag = APPEND

  if (addflag == APPEND) {
    tagint *atag = atom->tag;
    tagint *etag = element->tag;
    int nalocal = atom->nlocal;
    int nelocal = element->nlocal;
    tagint maxaid = 0, maxeid = 0;
    for (int i = 0; i < nalocal; i++) maxaid = MAX(maxaid,atag[i]);
    for (int i = 0; i < nelocal; i++) maxeid = MAX(maxeid,etag[i]);
    MPI_Allreduce(&maxaid,&aid_offset,1,MPI_CAC_TAGINT,MPI_MAX,world);
    MPI_Allreduce(&maxeid,&eid_offset,1,MPI_CAC_TAGINT,MPI_MAX,world);
  }


  // -----------------------------------------------------------------
 
  // flags for this data file
  
  int atomflag,elementflag,nodeflag;
  int element_type_flag,element_cluster_flag;

  elementflag = atomflag = nodeflag = 0;
  element_type_flag = element_cluster_flag = 0;

  // values in this data file

  natoms = nelements = nnodes = 0;
  netypes = natypes = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  triclinic = 0;

  nalocal_previous = atom->nlocal;
  nelocal_previous = element->nlocal;

  if (me == 0) {
    if (screen) fprintf(screen,"Reading data file '%s'...\n",arg[0]);
    if (logfile) fprintf(logfile,"Reading data file '%s'...\n",arg[0]);
    open(arg[0]);
  } else fp = NULL;

  // read header info

  header();

  // problem setup using info from header
  // only done once, if first data file

  if (addflag == NONE) {
    int n;
    if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
    else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);
    atom->allocate_type_arrays();
    atom->avec->grow(n);

    if (comm->nprocs == 1) n = static_cast<int> (element->nelements);
    else n = static_cast<int> (LB_FACTOR * element->nelements / comm->nprocs);
    element->evec->grow(n);

    domain->boxlo[0] = boxlo[0]; domain->boxhi[0] = boxhi[0];
    domain->boxlo[1] = boxlo[1]; domain->boxhi[1] = boxhi[1];
  	domain->boxlo[2] = boxlo[2]; domain->boxhi[2] = boxhi[2];

    if (triclinic) {
      domain->triclinic = 1;
      domain->xy = xy; domain->xz = xz; domain->yz = yz;
    }

    domain->print_box("  ");
    domain->set_initial_box();
    domain->set_global_box();
    comm->set_proc_grid();
    domain->set_local_box(); 
  
  // change simulation box to be union of existing box and new box + shift
  // only done if not first data file

  } else {
    domain->boxlo[0] = MIN(domain->boxlo[0],boxlo[0]+shift[0]);
    domain->boxhi[0] = MAX(domain->boxhi[0],boxhi[0]+shift[0]);
    domain->boxlo[1] = MIN(domain->boxlo[1],boxlo[1]+shift[1]);
    domain->boxhi[1] = MAX(domain->boxhi[1],boxhi[1]+shift[1]);
    domain->boxlo[2] = MIN(domain->boxlo[2],boxlo[2]+shift[2]);
    domain->boxhi[2] = MAX(domain->boxhi[2],boxhi[2]+shift[2]);

    // NOTE: not sure what to do about tilt value in subsequent data files
    //if (triclinic) {
    //  domain->xy = xy; domain->xz = xz; domain->yz = yz;
    // }

    domain->print_box("  ");
    domain->set_initial_box();
    domain->set_global_box();
    comm->set_proc_grid();
    domain->set_local_box();
  } 

  // customize for new sections
  // read rest of file in free format

  while (strlen(keyword)) {
    if (strcmp(keyword,"Masses") == 0) {
      mass();
    } else if (strcmp(keyword,"Element Types") == 0) {
      element_types();
      element_type_flag = 1;
      element_cluster_flag = element->element_cluster_flag;

    } else if (strcmp(keyword,"Atoms") == 0) {
      atomflag = 1;
      if (me == 0 && !style_match(style,atom->atom_style))
        error->warning(FLERR,"Atom style in data file differs "
            "from currently defined atom style");
      atoms();
    } else if (strcmp(keyword,"Velocities") == 0) {
      if (atomflag == 0)
        error->all(FLERR,"Atoms section should be read before Velocity section");
      atom_velocities();
    } else if (strcmp(keyword,"Elements") == 0) {
      if (element_type_flag == 0) 
        error->all(FLERR,"Element types must be defined before Elements section");
      elementflag = 1;
      elements();
    } else if (strcmp(keyword,"Element Clusters") == 0) {
      if (element_cluster_flag == 0) {
        error->warning(FLERR,"No element with multiple atoms per unit cell, Element Clusters section is ignored");
      }
      else if (elementflag == 0)
        error->all(FLERR,"Elements section should be read before Element Clusters section");
      element_clusters();
    } else if (strcmp(keyword,"Nodes") == 0) {
      if (elementflag == 0)
        error->all(FLERR,"Elements section should be read before Nodes section");
      nodeflag = 1;
      nodes();
    } else if (strcmp(keyword,"Node Velocities") == 0) {
      if (elementflag == 0)
        error->all(FLERR,"Elements section should be read before Nodes Velocity section");
      node_velocities();
    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(FLERR,str);
    }
    parse_keyword(0);
  }

  // error if natoms/nelements > 0 yet no atoms/elements were read

  if (natoms > 0 && atomflag == 0)
    error->all(FLERR,"No atoms in data file");
  if (nelements > 0 && (elementflag == 0 || nodeflag == 0))
    error->all(FLERR,"No elements/nodes in data file");

  // close file

  if (me == 0) {
    fclose(fp);
    fp = NULL;
  }


  // init per-atom fix/compute/variable values for created atoms/elements

  atom->data_fix_compute_variable(nalocal_previous,atom->nlocal);
  element->data_fix_compute_variable(nelocal_previous,element->nlocal);

  // assign atoms/elements/nodes added by this data file to specified group

  if (groupbit) {
    int *amask = atom->mask;
    int *emask = element->mask;
    int **nodemask = element->nodemask;
    int *etype = element->etype;
    int nalocal = atom->nlocal;
    int nelocal = element->nlocal;
    for (int i = nalocal_previous; i < nalocal; i++)
      amask[i] |= groupbit;
    for (int i = nelocal_previous; i < nelocal; i++) {
      emask[i] |= groupbit;
      for (int j = 0; j < element->npe[etype[i]]; j++)
        nodemask[i][j] != groupbit;
    }
  }


}


/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
   ------------------------------------------------------------------------- */

void ReadData::open(char *file)
{
  compressed = 0;

  if (!compressed) fp = fopen(file,"r");

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
   if EOF, line is set to blank line
   else line has first keyword line for rest of file
   some logic differs if adding atoms
   ------------------------------------------------------------------------- */

void ReadData::header()
{

  int n;
  char *ptr;

  // customize for new sections
  // When adding new sections, remember to change NSECTIONS
  const char *section_keywords[NSECTIONS] =
  {"Atoms","Velocity","Masses","Element Types",
    //"Ellipsoids","Lines","Triangles","Bodies","Bonds","Angles","Dihedrals","Impropers",
    //"Pair Coeffs","PairIJ Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs",
    //"BondBond Coeffs","BondAngle Coeffs","MiddleBondTorsion Coeffs",
    //"EndBondTorsion Coeffs","AngleTorsion Coeffs",
    //"AngleAngleTorsion Coeffs","BondBond13 Coeffs","AngleAngle Coeffs",
    "Elements","Nodes","Node Velocity","Element Clusters"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world); 

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable
    // customize for new header lines

    if (strstr(line,"atoms")) {
      sscanf(line,BIGINT_FORMAT,&natoms);
      if (addflag == NONE) atom->natoms = natoms;
      else atom->natoms += natoms;

    } else if (strstr(line,"elements")) {
      sscanf(line, BIGINT_FORMAT,&nelements);
      if (addflag == NONE) 
        element->nelements = nelements;
      else 
        element->nelements += nelements;

      // Atom/Element class type settings are only set by first data file

    } else if (strstr(line,"atom types")) {
      sscanf(line,"%d",&natypes);
      if (addflag == NONE) atom->ntypes = natypes + extra_atom_types;

    } else if (strstr(line,"element types")) {
      sscanf(line,"%d",&netypes);
      if (addflag == NONE) element->evec->grow_etype_arrays(netypes + extra_element_types);
      // local copy of box info
      // so can treat differently for first vs subsequent data files

    } else if (strstr(line,"xlo xhi")) {
      sscanf(line,"%lg %lg",&boxlo[0],&boxhi[0]);
    } else if (strstr(line,"ylo yhi")) {
      sscanf(line,"%lg %lg",&boxlo[1],&boxhi[1]); 
    } else if (strstr(line,"zlo zhi")) {
      sscanf(line,"%lg %lg",&boxlo[2],&boxhi[2]);
    } else if (strstr(line,"xy xz yz")) {
      triclinic = 1;
      sscanf(line,"%lg %lg %lg",&xy,&xz,&yz);

    } else break;
  }

  // error check on total system size

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT ||
      element->nelements < 0 || element->nelements >= MAXBIGINT)
    error->all(FLERR,"System in data file is too big");

  // check that existing string is a valid section keyword

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;

  if (n == NSECTIONS) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->all(FLERR,str);
  } 

}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   optional style can be appended after comment char '#'
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
   ------------------------------------------------------------------------- */

void ReadData::parse_keyword(int first)
{
  int eof = 0;
  int done = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line," \t\n\r");
      if ((blank == strlen(line)) || (line[blank] == '#')) {
        if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      } else done = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  // store optional "style" following comment char '#' after keyword

  char *ptr;
  if ((ptr = strchr(line,'#'))) {
    *ptr++ = '\0';
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    int stop = strlen(ptr) - 1;
    while (ptr[stop] == ' ' || ptr[stop] == '\t'
        || ptr[stop] == '\n' || ptr[stop] == '\r') stop--;
    ptr[stop+1] = '\0';
    strcpy(style,ptr);
  } else style[0] = '\0';

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
      || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   compare two style strings if they both exist
   one = comment in data file section, two = currently-defined style
   ignore suffixes listed in suffixes array at top of file
   ------------------------------------------------------------------------- */

int ReadData::style_match(const char *one, const char *two)
{
  int i,delta,len,len1,len2;

  if ((one == NULL) || (two == NULL)) return 1;

  len1 = strlen(one);
  len2 = strlen(two);

  //  for (i = 0; suffixes[i] != NULL; i++) {
  //    len = strlen(suffixes[i]);
  //    if ((delta = len1 - len) > 0)
  //      if (strcmp(one+delta,suffixes[i]) == 0) len1 = delta;
  //    if ((delta = len2 - len) > 0)
  //      if (strcmp(two+delta,suffixes[i]) == 0) len2 = delta;
  //  }

  if ((len1 == 0) || (len1 == len2) || (strncmp(one,two,len1) == 0)) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   read all atoms
   ------------------------------------------------------------------------- */

void ReadData::atoms()
{

  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading atoms ...\n");
    if (logfile) fprintf(logfile,"  reading atoms ...\n");
  }

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_atoms(nchunk,buffer,aid_offset,toffset,shiftflag,shift);
    nread += nchunk;
  }

  // check that all atoms were assigned correctly

  bigint n = atom->nlocal;
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
  bigint nassign = sum - (atom->natoms - natoms);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",nassign);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",nassign);
  }

  if (sum != atom->natoms)
    error->all(FLERR,"Did not assign all atoms correctly");
}

/* ---------------------------------------------------------------------- */

void ReadData::mass()
{
  char *next;
  char *buf = new char[natypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,natypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < natypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    atom->set_mass(FLERR,buf,toffset);
    buf = next + 1;
  }
  delete [] original;
}

/*-----------------------------------------------------------------------*/

void ReadData::element_types()
{
  char *next;
  char *buf = new char[netypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,netypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_element_types(buf,eoffset);
    buf = next + 1; 
  }
  delete [] original;
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
   ------------------------------------------------------------------------- */

void ReadData::skip_lines(bigint n)
{
  if (me) return;
  char *eof;
  for (bigint i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
}

/* ----------------------------------------------------------------------
   read all elements
   ------------------------------------------------------------------------- */

void ReadData::elements()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading elements ...\n");
    if (logfile) fprintf(logfile,"  reading elements ...\n");
  }

  bigint nread = 0;

  while (nread < nelements) {
    nchunk = MIN((nelements-nread), CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    element->data_elements(nchunk,buffer,eid_offset,toffset,eoffset,shiftflag,shift);
    nread += nchunk;
  }

  // count total number of nodes

  // check that all elements were assigned correctly 
  // count total number of nodes
  // count and check total number of clusters
  // scale up clusters count by max_apc factorial to make it integer

  bigint n = element->nlocal;

  bigint nnode_local = 0;
  int product = 1;
  for (int i = 2; i <= element->max_apc; i++)
    product *= i;

  bigint nclusters_local = 0;
  for (int i = 0; i < n; i++) {
    int itype = element->etype[i];
    nnode_local += element->npe[itype];
    nclusters_local += product/element->apc[itype];
  }

  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&nnode_local,&nnodes,1,MPI_CAC_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&nclusters_local,&nclusters,1,MPI_DOUBLE,MPI_SUM,world);
  if (nclusters%product && element->element_cluster_flag)
    error->all(FLERR,"Element might be missing for element clusters");
  if (addflag == NONE) 
    element->nnodes = nnodes;
  else 
    element->nnodes += nnodes;

  element->nclusters = nclusters = nclusters/product;

  bigint nassign = sum - (element->nelements - nelements);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " elements\n",nassign);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " elements\n",nassign);
  }

  if (sum != element->nelements)
    error->all(FLERR,"Did not assign all elements correctly");

  n = element->nintpls;
  MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
  element->nintpls = sum;

  // check that element IDs are valid

  element->tag_check();

  // create global mapping of elements

  if (element->map_style) {
    element->map_init();
    element->map_set();
  }
}

/* ----------------------------------------------------------------------
   read all nodes
   to find elements, must build element map
   ------------------------------------------------------------------------- */

void ReadData::nodes()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading nodes ...\n");
    if (logfile) fprintf(logfile,"  reading nodes ...\n");
  }

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  bigint nread = 0;

  while (nread < nnodes) {
    nchunk = MIN(nnodes-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    element->data_nodes(nchunk,buffer,eid_offset,shiftflag,shift);
    nread += nchunk;
  }  

  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  element->check_node_coords(1);
  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " nodes\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " nodes\n",nnodes);
  }

}

/* ----------------------------------------------------------------------
   read all element clusters
   to find elements, must build element map
   ------------------------------------------------------------------------- */

void ReadData::element_clusters()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading element clusters ...\n");
    if (logfile) fprintf(logfile,"  reading element clusters ...\n");
  }

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  // set element_clusters array to 0 for read elements

  tagint **element_clusters = element->element_clusters;
  for (int i = nelocal_previous; i < element->nlocal; i++) 
    for (int j = 0; j < element->max_apc; j++) 
      element_clusters[i][j] = 0;

  bigint nread = 0;

  while (nread < nclusters) {
    nchunk = MIN(nclusters-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    element->data_element_clusters(nchunk,buffer,eid_offset);
    nread += nchunk;
  }  

  // sort element clusters array so that the first ID is the smallest ID (swap the smallest with the first item)

  int *apc = element->apc;
  int *etype = element->etype;
  int min;
  int min_index;
  for (int i = nelocal_previous; i < element->nlocal; i++) {
    min = element_clusters[i][0];
    min_index = 0;
    for (int j = 1; j < apc[etype[i]]; j++) {
      if (min > element_clusters[i][j]) {
        min = element_clusters[i][j];
        min_index = j;
      }
    }
    element_clusters[i][min_index] = element_clusters[i][0];
    element_clusters[i][0] = min;
  }

  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  element->check_node_coords(1);
  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " nodes\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " nodes\n",nnodes);
  }
}

/* ----------------------------------------------------------------------
   read all atom velocities
   to find atoms, must build atom map if not a molecular system
   ------------------------------------------------------------------------- */

void ReadData::atom_velocities()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading atom velocities ...\n");
    if (logfile) fprintf(logfile,"  reading atom velocities ...\n");
  }

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  bigint nread = 0;

  while (nread < natoms) {
    nchunk = MIN(natoms-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    atom->data_vels(nchunk,buffer,aid_offset);
    nread += nchunk;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atom velocities\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atom velocities\n",natoms);
  }
}

/* ----------------------------------------------------------------------
   read all node velocities
   to find elements, must build element map
   ------------------------------------------------------------------------- */

void ReadData::node_velocities()
{
  int nchunk,eof;

  if (me == 0) {
    if (screen) fprintf(screen,"  reading node velocities ...\n");
    if (logfile) fprintf(logfile,"  reading node velocities ...\n");
  }

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  bigint nread = 0;

  while (nread < nnodes) {
    nchunk = MIN(nnodes-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    element->data_vels(nchunk,buffer,eid_offset);
    nread += nchunk;
  }

  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " node velocities\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " node velocities\n",nnodes);
  }
}


