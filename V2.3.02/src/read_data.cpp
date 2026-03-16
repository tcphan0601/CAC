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
#include "force.h"
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
#define NSECTIONS 31       // change when add to header::section_keywords
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
  id_offset = 0;
  offsetflag = shiftflag = 0;
  toffset = 0;
  extra_atom_types = 0;
  shift[0] = shift[1] = shift[2] = 0.0;

  groupbit = 0;

  nfix = 0;
  fix_index = NULL;
  fix_header = NULL;
  fix_section = NULL;

  if (addflag == NONE) {
    domain->box_exist = 1;
    update->ntimestep = 0;
  }

  int atomflag,elementflag;

  elementflag = atomflag = 0;

  // values in this data file

  natoms = 0;
  nelements = 0;
  nnodes = 0;
  netypes = natypes = 0;

  triclinic = 0;

  npe = element->npe;

  if (me == 0) {
    if (screen) fprintf(screen,"Reading data file ...\n");
    if (logfile) fprintf(logfile,"Reading data file ...\n");
    open(arg[0]);
  } else fp = NULL;

  // read header info

  header();

  /*-------------------------
   * need fix
   * -----------------------*/
  if (addflag == NONE) {
    int n;
    if (comm->nprocs == 1) n = static_cast<int> (atom->natoms);
    else n = static_cast<int> (LB_FACTOR * atom->natoms / comm->nprocs);
    atom->allocate_type_arrays();

    // allocate the array to store the number of atoms along each direction of the element

    atom->avec->grow(n);


    if (comm->nprocs == 1) n = static_cast<int> (element->nelements);
    else n = static_cast<int> (LB_FACTOR * element->nelements / comm->nprocs);

    element->evec->grow_etype_arrays(netypes);
    element->evec->grow(n);
    /*-------------------------*/

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
  }

  // customize for new sections
  // read rest of file in free format

  while (strlen(keyword)) {
    if (strcmp(keyword,"Atoms") == 0) {
      atomflag = 1;
      if (me == 0 && !style_match(style,atom->atom_style))
        error->warning(FLERR,"Atom style in data file differs "
            "from currently defined atom style");
      atoms();
    } else if (strcmp(keyword,"Atom Velocities") == 0) {
      if (atomflag == 0)
        error->all(FLERR,"Atoms section should be read before Atoms Velocity section");
      atom_velocities();
    } else if (strcmp(keyword,"Elements") == 0) {
      elementflag = 1;
      elements();
    } else if (strcmp(keyword, "Nodes") == 0){
      if (elementflag == 0)
        error->all(FLERR,"Elements section should be read before Nodes section");
      nodes();
    } else if (strcmp(keyword,"Node Velocities") == 0) {
      if (elementflag == 0)
        error->all(FLERR,"Elements section should be read before Nodes Velocity section");
      node_velocities();
    } else if (strcmp(keyword,"Masses") == 0) {
      mass();
    } else if (strcmp(keyword,"Interpolate") == 0) {
      if (elementflag) error->all(FLERR,"Interpolate section should be define before Element section");
      interpolate();
    } else if (strcmp(keyword,"Integration") == 0) {
      if (elementflag) error->all(FLERR,"Integration section should be define before Element section");
      integration();
      element->intg_input_style = 1;
    } else if (strcmp(keyword,"Integration List") == 0) {
      if (elementflag) error->all(FLERR,"Integration section should be define before Element section");
      integration_user();
      element->intg_input_style = 2;
    } else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(FLERR,str);
    }
    parse_keyword(0);
  }

  if (natoms > 0 && atomflag == 0)
    error->all(FLERR,"No atoms in data file");
  if (nelements > 0 && elementflag == 0)
    error->all(FLERR,"No elements in data file");

  if (me == 0) fclose(fp);

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
  {"Atoms","Atom Velocities","Ellipsoids","Lines","Triangles","Bodies",
    "Bonds","Angles","Dihedrals","Impropers",
    "Masses","Pair Coeffs","PairIJ Coeffs","Bond Coeffs","Angle Coeffs",
    "Dihedral Coeffs","Improper Coeffs", "Interpolate", 
    "Intergration List", "Integration",
    "BondBond Coeffs","BondAngle Coeffs","MiddleBondTorsion Coeffs",
    "EndBondTorsion Coeffs","AngleTorsion Coeffs","Elements","Nodes", "Node Velocity",
    "AngleAngleTorsion Coeffs","BondBond13 Coeffs","AngleAngle Coeffs"};

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
      element->nelements = nelements;
      nnodes = nelements*npe;
      element->nnodes = nnodes;
    } else if (strstr(line,"xlo xhi")) {
      sscanf(line,"%lg %lg",&boxlo[0],&boxhi[0]);
    } else if (strstr(line,"ylo yhi")) {
      sscanf(line,"%lg %lg",&boxlo[1],&boxhi[1]); 
    } else if (strstr(line,"zlo zhi")) {
      sscanf(line,"%lg %lg",&boxlo[2],&boxhi[2]);
    } else if (strstr(line,"xy xz yz")) {
      triclinic = 1;
      sscanf(line,"%lg %lg %lg",&xy,&xz,&yz);
    } else if (strstr(line,"atom types")) {
      sscanf(line,"%d",&natypes);
      if (addflag == NONE) atom->ntypes = natypes + extra_atom_types;
      element->nctypes = natypes + extra_atom_types;
    } else if (strstr(line,"element types")) {
      sscanf(line,"%d",&netypes);
    }  else break;
  }

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
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
    atom->data_atoms(nchunk,buffer,id_offset,toffset,shiftflag,shift);
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
    atom->set_mass(buf,toffset);
    buf = next + 1;
  }
  delete [] original;
}

/*-----------------------------------------------------------------------*/

void ReadData::interpolate()
{
  char *next;
  char *buf = new char[netypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,netypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_interpolate(buf,toffset);
    buf = next + 1; 
  }
  delete [] original;
}
/*-----------------------------------------------------------------------*/

void ReadData::integration()
{
  char *next;
  char *buf = new char[netypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,netypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_integration(buf,toffset);
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
    element->data_elements(nchunk,buffer);
    nread += nchunk;
  }
  
  bigint n = element->nlocal;
  bigint sum;
  MPI_Allreduce(&n,&sum,1,MPI_CAC_BIGINT,MPI_SUM,world);
  bigint nassign = sum - (element->nelements - nelements);

  if (me == 0) {

    if (sum != element->nelements)
      error->all(FLERR,"Did not assign all elements correctly");

  }
}

/* ----------------------------------------------------------------------
   read all nodes
   to find elements, must build element map
   since the node is taken to be the 
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
    element->data_nodes(nchunk,buffer);
    nread += nchunk;
  }  

  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }
  element->check_node_coords();
  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " nodes\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " nodes\n",nnodes);
  }
}

/*-------------------------------------------------------------------------------*/

void ReadData::integration_user()
{

  int n,itype,num;
  char *next;
  char *buf = new char[netypes*MAXLINE];
  int *nintg = element->nintg;
  int *nintpl = element->nintpl;

  int eof = comm->read_lines_from_file(fp,netypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_nintg_user(buf,toffset);
    buf = next+1;
  }
  delete [] original;

  // Count total number of integration points
  int totalnintg = 0;

  for (int i = 1; i <= netypes; i++)
    totalnintg += nintg[i];

  buf = new char[totalnintg*MAXLINE];

  eof = comm->read_lines_from_file(fp,totalnintg,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  original = buf;

  for (int i = 0; i < totalnintg; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_intg_user(buf,toffset);
    buf = next+1;
  }

  delete [] original;

  // check if integration points has been assigned properly

  int **i2ia = element->i2ia;
  int **n2i = element->n2i;
  double **weight = element->weight;
  double total_weight;
  for (itype = 1; itype <= netypes; itype++){
    itype += toffset;
    total_weight = 0.0;

    // look for duplication, missing integration points from user input

    for (int iintg = 0; iintg < nintg[itype]; iintg++){
      if (i2ia[itype][iintg] < 0) 
        error->all(FLERR,"Intergration points for this element type has not all been set");
      for (int jintg = iintg + 1; jintg < nintg[itype]; jintg++)
        if (i2ia[itype][iintg] == i2ia[itype][jintg]) {
          char *erstr = new char[100];
          sprintf(erstr,"Duplicate integration points in data file: i2ia[%d] = %d, i2ia[%d] = %d",iintg,i2ia[itype][iintg],jintg,i2ia[itype][jintg]);
          error->all(FLERR,erstr);
          delete [] erstr;
        }
      total_weight += weight[itype][iintg];
    }

    // check weight 

    if (fabs(total_weight - nintpl[itype]) > EPSILON) {
      char *erstr = new char[100];
      sprintf(erstr,"Total weight = %f does not match number of atoms = %d in this element type",total_weight,nintpl[itype]);
      error->all(FLERR,erstr);
      delete [] erstr;
    }

    // check node to integration point mapping

    for (int l = 0; l < npe; l++)
      if (n2i[itype][l] < 0) 
        error->all(FLERR,"Node to integration point mapping is not complete");
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
    atom->data_vels(nchunk,buffer,id_offset);
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
    element->data_vels(nchunk,buffer);
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


