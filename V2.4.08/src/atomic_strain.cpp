#include "TECIO.h"
#include "cactype.h" 
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "atomic_strain.h"
#include "domain.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "element_vec.h"
#include "universe.h"
#include "comm.h"
#include "balance.h"
#include "irregular_comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "universe.h"
#include "error.h"
#include "memory.h"
#include "math_extra.h"

using namespace CAC_NS;
using namespace MathExtra;

#define MAXLINE 256
#define LB_FACTOR 1.2
#define CHUNK 1024
#define DELTA 4            // must be 2 or larger
#define MAXBODY 20         // max # of lines in one body, also in Atom class
#define NSECTIONS 5       // change when add to header::section_keywords
#define EPSILON 1e-6

enum{PLT = 0, SZPLT = 1, ASCII = 2, ATOM = 3};        // output file format ASCII (.dat), binary (.plt), or SZL (.szplt)
enum{COORDINATE,ID,ETYPE,CTYPE,NUMNEIGHA,NUMNEIGHIA,VALIDITY,VONMISES,VOLSTRAIN,DEFGRAD,ROTATION,NSD,STRAINTENSOR,STRETCHTENSOR};
enum{REFERENCE,CURRENT};
/* ---------------------------------------------------------------------- */

AtomicStrain::AtomicStrain(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  balance = NULL;
  debug = 0;
  atomValidity = NULL;
  atomVonMisesStrain = NULL;
  atomVolStrain = NULL;
  atomF = NULL;
  atomQ = NULL;
  atomStretchTensor = NULL;
  atomStrainTensor = NULL;
  atomNSD = NULL;

  nodeValidity = NULL;
  nodeVonMisesStrain = NULL;
  nodeVolStrain = NULL;
  nodeF = NULL;
  nodeQ = NULL;
  nodeStretchTensor = NULL;
  nodeStrainTensor = NULL;
  nodeNSD = NULL;

  atom_displacement = NULL;
  node_displacement = NULL;

  valuelocation = NULL;
  buf = NULL;
  ibuf = NULL;
  maxbuf = maxibuf = 0;
}

/* ---------------------------------------------------------------------- */

AtomicStrain::~AtomicStrain()
{
  delete [] line;
  delete [] keyword;
  delete [] style;
  delete [] buffer;
  memory->sfree(arg);

  memory->destroy(atomValidity);
  memory->destroy(atomVonMisesStrain);
  memory->destroy(atomVolStrain);
  memory->destroy(atomF);
  memory->destroy(atomQ);
  memory->destroy(atomStretchTensor);
  memory->destroy(atomStrainTensor);
  memory->destroy(atomNSD);

  memory->destroy(nodeValidity);
  memory->destroy(nodeVonMisesStrain);
  memory->destroy(nodeVolStrain);
  memory->destroy(nodeF);
  memory->destroy(nodeQ);
  memory->destroy(nodeStretchTensor);
  memory->destroy(nodeStrainTensor);
  memory->destroy(nodeNSD);

  memory->destroy(atom_displacement);
  memory->destroy(node_displacement);

  memory->destroy(buf);
  memory->destroy(ibuf);
  delete [] valuelocation;
  delete balance;
}

/* ---------------------------------------------------------------------- */

void AtomicStrain::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal atomic_strain command");

  nelements = natoms = 0;

  atom->atom_strain_flag = 1;

  output_format = PLT;
  output_configuration = CURRENT;
  cut = 6;
  element_info_flag = 0;
  deformation_gradient_flag = 0;
  rotation_flag = 0;
  nsd_flag = 0;
  strain_tensor_flag = 0;
  stretch_tensor_flag = 0;
  wrap_flag = 1;
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"out/format") == 0) {
      if (strcmp(arg[iarg+1],"plt") == 0) output_format = PLT;
      else if (strcmp(arg[iarg+1],"szplt") == 0) output_format = SZPLT;
      else if (strcmp(arg[iarg+1],"ascii") == 0) output_format = ASCII;
      else if (strcmp(arg[iarg+1],"atom") == 0) output_format = ATOM;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      cut = universe->numeric(FLERR,arg[iarg+1]);
      iarg += 2; 
    } else if (strcmp(arg[iarg],"out/config") == 0) {
      if (strcmp(arg[iarg+1],"current") == 0) output_configuration = CURRENT;
      else if (strcmp(arg[iarg+1],"reference") == 0) output_configuration = REFERENCE;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"def/gradient") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) deformation_gradient_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) deformation_gradient_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"element/info") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) element_info_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) element_info_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rotation") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) rotation_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) rotation_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nonaff/squared/disp") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) nsd_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) nsd_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"stretch/tensor") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) stretch_tensor_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) stretch_tensor_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"strain/tensor") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) strain_tensor_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) strain_tensor_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"debug") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) debug = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) debug = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"wrap") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) wrap_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) wrap_flag = 0;
      else error->all(FLERR,"Illegal atomic_strain command");
      iarg += 2;
    } else error->all(FLERR,"Illegal atomic_strain command");
  }

  // error checks

  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (comm->style != 1) 
    error->all(FLERR,"atomic_strain command requires comm_style tiled for workload balancing");

  // box domain must exist

  //if (!domain->box_exist)
  //  error->all(FLERR,"Box domain not yet exist");

  // remove all atoms/elements

  atom->nlocal = 0;
  atom->natoms = 0;
  element->nlocal = 0;
  element->nelements = 0;

  // read input data
  // reference coordinates are wrapped
  // current coordinates are unwrapped

  read_input(arg[0]);
  domain->box_exist = 1;

  // reset border comm sizes

  element->evec->set_data_size();
  atom->avec->size_border = 11;

  // temporary set masses to avoid error (will not be used)
  
  double *masses = new double[atom->ntypes];
  for (int i = 0; i < atom->ntypes; i++) 
    masses[i] = 1.0;
  atom->set_mass(masses);
  delete [] masses;

  // do a balance work load
  
  balance = new Balance(cac);
  char **newarg = new char*[4];
  newarg[0] = (char *) "1.0";
  newarg[1] = (char *) "rcb";
  newarg[2] = (char *) "eweight";
  newarg[3] = (char *) "100";
  balance->command(4,newarg);

  // request a full atomic neighbor list for use by this command 
  // also sort out neighbors based on distance

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 2;
  neighbor->requests[irequest]->nodelist = 1;
  neighbor->requests[irequest]->intglist = 0;
  neighbor->requests[irequest]->cut = 1;
  neighbor->requests[irequest]->cutoff = cut;

  neighbor->init();
  comm->init();

  // setup domain, communication and neighboring
  // exchange and acquire ghosts
  // build standard neighbor list to setup bins & stencils

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
    domain->x2lamda_current(atom->nlocal,atom->x_current);
    domain->nodex2lamda_current(element->nlocal,element->nodex_current);
  }
  comm->setup_exchange();
  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 

  comm->borders();

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost,atom->x);
    domain->lamda2x(element->nlocal+element->nghost,element->x);
    domain->lamda2nodex(element->nlocal+element->nghost,element->nodex);
    domain->lamda2x_current(atom->nlocal+atom->nghost,atom->x_current);
    domain->lamda2nodex_current(element->nlocal+element->nghost,element->nodex_current);

  }
  neighbor->build_bins();

  // build neighbor list this command needs based on earlier request

  list = neighbor->lists[irequest];
  neighbor->build_one(list,1);

  compute();

  // write output file

  write_output(arg[1]);

  // reset number of atoms/elements

  atom->atom_strain_flag = 0;
  element->evec->set_data_size();
  atom->avec->size_border = 7;
  domain->box_exist = 0; 

  atom->natoms = 0;
  atom->nlocal = 0;
  atom->nghost = 0;
  element->nelements = 0;
  element->nlocal = 0;
  element->nghost = 0;
  element->netypes = 0;
  atom->atom_strain_flag = 0;
}


/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
   ------------------------------------------------------------------------- */

void AtomicStrain::open_input(char *file)
{
  fp = fopen(file,"r");

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open input file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens output file
   test if gzipped
   ------------------------------------------------------------------------- */

void AtomicStrain::open_output(char *file)
{
  fp = fopen(file,"w");

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open output file %s",file);
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

void AtomicStrain::header()
{

  int n;
  char *ptr;

  // customize for new sections
  // When adding new sections, remember to change NSECTIONS
  const char *section_keywords[NSECTIONS] =
  {"Atoms","Element Types",
    //"Ellipsoids","Lines","Triangles","Bodies","Bonds","Angles","Dihedrals","Impropers",
    //"Pair Coeffs","PairIJ Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs",
    //"BondBond Coeffs","BondAngle Coeffs","MiddleBondTorsion Coeffs",
    //"EndBondTorsion Coeffs","AngleTorsion Coeffs",
    //"AngleAngleTorsion Coeffs","BondBond13 Coeffs","AngleAngle Coeffs",
    "Elements","Nodes","Element Clusters"};

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
      atom->natoms = natoms;

    } else if (strstr(line,"elements")) {
      sscanf(line,BIGINT_FORMAT,&nelements);
      element->nelements = nelements;
    }

    // Atom/Element class type and grain settings are only set by first data file

    else if (strstr(line,"atom types")) {
      sscanf(line,"%d",&natypes);
      atom->ntypes = natypes;
    } else if (strstr(line,"element types")) {
      sscanf(line,"%d",&netypes);
      element->evec->grow_etype_arrays(0);
      element->evec->grow_etype_arrays(netypes);
    }

    // local copy of box info
    // so can treat differently for first vs subsequent data files

    else if (strstr(line,"xloref xhiref")) {
      sscanf(line,"%lg %lg",&boxlo[0],&boxhi[0]);
    } else if (strstr(line,"yloref yhiref")) {
      sscanf(line,"%lg %lg",&boxlo[1],&boxhi[1]); 
    } else if (strstr(line,"zloref zhiref")) {
      sscanf(line,"%lg %lg",&boxlo[2],&boxhi[2]);
    } else if (strstr(line,"xyref xzref yzref")) {
      triclinic = 1;
      sscanf(line,"%lg %lg %lg",&xy,&xz,&yz);
    }

    // local copy of current box info

    else if (strstr(line,"xlocur xhicur")) {
      sscanf(line,"%lg %lg",&boxlo_current[0],&boxhi_current[0]);
    } else if (strstr(line,"ylocur yhicur")) {
      sscanf(line,"%lg %lg",&boxlo_current[1],&boxhi_current[1]); 
    } else if (strstr(line,"zlocur zhicur")) {
      sscanf(line,"%lg %lg",&boxlo_current[2],&boxhi_current[2]);
    } else if (strstr(line,"xycur xzcur yzcur")) {
      triclinic = 1;
      sscanf(line,"%lg %lg %lg",&xy_current,&xz_current,&yz_current);

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

void AtomicStrain::parse_keyword(int first)
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

int AtomicStrain::style_match(const char *one, const char *two)
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

void AtomicStrain::atoms()
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
    atom->data_atoms_strain(nchunk,buffer);
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

/*-----------------------------------------------------------------------*/

void AtomicStrain::element_types()
{
  char *next;
  char *buf = new char[netypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,netypes,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    element->evec->set_element_types(buf,0);
    buf = next + 1; 
  }
  delete [] original;
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
   ------------------------------------------------------------------------- */

void AtomicStrain::skip_lines(bigint n)
{
  if (me) return;
  char *eof;
  for (bigint i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
}

/* ----------------------------------------------------------------------
   read all elements
   ------------------------------------------------------------------------- */

void AtomicStrain::elements()
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
    element->data_elements(nchunk,buffer,0,0,0,0,NULL);
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
  element->nnodes = nnodes;
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

void AtomicStrain::nodes()
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
    element->data_nodes_strain(nchunk,buffer);
    nread += nchunk;
  }  

  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  element->evec->update_node_coord();
  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " nodes\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " nodes\n",nnodes);
  }

}

/* ----------------------------------------------------------------------
   read all element clusters
   to find elements, must build element map
   ------------------------------------------------------------------------- */

void AtomicStrain::element_clusters()
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
  for (int i = 0; i < element->nlocal; i++) 
    for (int j = 0; j < element->max_apc; j++) 
      element_clusters[i][j] = 0;

  bigint nread = 0;

  while (nread < nclusters) {
    nchunk = MIN(nclusters-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) error->all(FLERR,"Unexpected end of data file");
    element->data_element_clusters(nchunk,buffer,0);
    nread += nchunk;
  }  

  // sort element clusters array so that the first ID is the smallest ID (swap the smallest with the first item)

  int *apc = element->apc;
  int *etype = element->etype;
  int min;
  int min_index;
  for (int i = 0; i < element->nlocal; i++) {
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

  element->evec->update_node_coord();
  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " nodes\n",nnodes);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " nodes\n",nnodes);
  }
}

/* ----------------------------------------------------------------------
   read input data
   ------------------------------------------------------------------------- */

void AtomicStrain::read_input(char *filename)
{

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


  if (me == 0) {
    if (screen) fprintf(screen,"Reading data file for atomic_strain command '%s'...\n",filename);
    if (logfile) fprintf(logfile,"Reading data file for atomic_strain command '%s'...\n",filename);
    open_input(filename);
  } else fp = NULL;

  // read header info

  header();

  // problem setup using info from header

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

  domain->boxlo_current[0] = boxlo_current[0]; domain->boxhi_current[0] = boxhi_current[0];
  domain->boxlo_current[1] = boxlo_current[1]; domain->boxhi_current[1] = boxhi_current[1];
  domain->boxlo_current[2] = boxlo_current[2]; domain->boxhi_current[2] = boxhi_current[2];

  if (triclinic) {
    domain->triclinic = 1;
    domain->xy = xy; domain->xz = xz; domain->yz = yz;
    domain->xy_current = xy_current; 
    domain->xz_current = xz_current; 
    domain->yz_current = yz_current;
  }

  domain->print_box("  ");
  domain->set_initial_box();
  domain->set_global_box();
  domain->set_current_global_box();
  comm->set_proc_grid();
  domain->set_local_box(); 

  // customize for new sections
  // read rest of file in free format

  while (strlen(keyword)) {
    if (strcmp(keyword,"Element Types") == 0) {
      element_types();
      element_type_flag = 1;
      element_cluster_flag = element->element_cluster_flag;

    } else if (strcmp(keyword,"Atoms") == 0) {
      atomflag = 1;
      if (me == 0 && !style_match(style,atom->atom_style))
        error->warning(FLERR,"Atom style in data file differs "
            "from currently defined atom style");
      atoms();
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

  // close input file

  if (me == 0) {
    fclose(fp);
    fp = NULL;
  }
}

/* ----------------------------------------------------------------------
   perform atomic strain calculation
   ------------------------------------------------------------------------- */

void AtomicStrain::compute()
{
  int nalocal = atom->nlocal;
  int naghost = atom->nghost;
  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int max_npe = element->max_npe;

  // grow arrays if needed based on flags

  memory->grow(atomValidity,nalocal,"atomic_strain:atomValidity");
  memory->grow(nodeValidity,nelocal,max_npe,"atomic_strain:nodeValidity");
  memory->grow(atomVonMisesStrain,nalocal,"atomic_strain:atomVonMisesStrain");
  memory->grow(nodeVonMisesStrain,nelocal,max_npe,"atomic_strain:nodeVonMisesStrain");
  memory->grow(atomVolStrain,nalocal,"atomic_strain:atomVolStrain");
  memory->grow(nodeVolStrain,nelocal,max_npe,"atomic_strain:nodeVolStrain");
  if (deformation_gradient_flag) {
    memory->grow(atomF,nalocal,9,"atomic_strain:atomF");  
    memory->grow(nodeF,nelocal,max_npe,9,"atomic_strain:nodeF");  
  }
  if (rotation_flag) {
    memory->grow(atomQ,nalocal,4,"atomic_strain:atomQ");  
    memory->grow(nodeQ,nelocal,max_npe,4,"atomic_strain:nodeQ");  
  }
  if (stretch_tensor_flag) {
    memory->grow(atomStretchTensor,nalocal,6,"atomic_strain:atomStretchTensor");  
    memory->grow(nodeStretchTensor,nelocal,max_npe,6,"atomic_strain:nodeStretchTensor");  
  }
  if (strain_tensor_flag) {
    memory->grow(atomStrainTensor,nalocal,6,"atomic_strain:atomStrainTensor");  
    memory->grow(nodeStrainTensor,nelocal,max_npe,6,"atomic_strain:nodeStrainTensor");  
  }
  if (nsd_flag) {
    memory->grow(atomNSD,nalocal,"atomic_strain:atomNSD");  
    memory->grow(nodeNSD,nelocal,max_npe,"atomic_strain:atomNSD");  
  }

  memory->grow(atom_displacement,atom->nmax,3,"atomic_strain:atom_displacement");
  memory->grow(node_displacement,element->nmax,max_npe,3,"atomic_strain:atom_displacement");

  double **x_ref = atom->x;
  double **x_cur = atom->x_current;
  double ***nodex_ref = element->nodex;
  double ***nodex_cur = element->nodex_current;
  int *npe = element->npe;
  int *etype = element->etype;
  double ***shape_array = element->shape_array;
  imageint *aimage = atom->image;
  imageint *eimage = element->image;
  double coord[3];

  // compute displacements vector for all atoms/nodes

  for (int i = 0; i < nalocal+naghost; i++) 
    sub3(x_cur[i],x_ref[i],atom_displacement[i]);     
  

  for (int i = 0; i < nelocal+neghost; i++) 
    for (int j = 0; j < npe[etype[i]]; j++) 
      sub3(nodex_cur[i][j],nodex_ref[i][j],node_displacement[i][j]);

//      if (element->tag[i] == 12818 && j == 0) {
//  int xbox = (eimage[i] & IMGMASK) - IMGMAX;
//  int ybox = (eimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
//  int zbox = (eimage[i] >> IMG2BITS) - IMGMAX;
//
//        printf("me = %d i = %d nlocal = %d j = %d image = %d %d %d\n ref_remap = %g %g %g \n       ref = %g %g %g\n       cur = %g %g %g\n       del = %g %g %g\n",
//            me,i,nelocal,j,xbox,ybox,zbox
//            ,nodex_ref[i][j][0],nodex_ref[i][j][1],nodex_ref[i][j][2]
//            ,coord[0],coord[1],coord[2]
//            ,nodex_cur[i][j][0],nodex_cur[i][j][1],nodex_cur[i][j][2]
//            ,nodex_cur[i][j][0]-coord[0] 
//            ,nodex_cur[i][j][1]-coord[1] 
//            ,nodex_cur[i][j][2]-coord[2] 
//            );
//      }
    
  int ainum = list->ainum; 
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int *numneighn2a = list->numneighn2a;
  int *numneighn2ia = list->numneighn2ia;
  int **firstneighn2a = list->firstneighn2a;
  int **firstneighn2ia = list->firstneighn2ia;
  int **firstneighn2ia_index = list->firstneighn2ia_index;
  int *e2nlist = list->e2nlist;

  double *defGradient,*Q,*stretchTensor,*strainTensor,NSD,vonMisesStrain,volStrain;

  // atomic strain calculation for atoms

  for (int i = 0; i < nalocal; i++) {
    if (deformation_gradient_flag) defGradient = atomF[i];
    else defGradient = NULL;
    if (rotation_flag) Q = atomQ[i];
    else Q = NULL;
    if (stretch_tensor_flag) stretchTensor = atomStretchTensor[i];
    else stretchTensor = NULL;
    if (strain_tensor_flag) strainTensor = atomStrainTensor[i];
    else strainTensor = NULL;

    atomValidity[i] = compute_strain(numneigha2a[i],firstneigha2a[i],numneigha2ia[i],
        firstneigha2ia[i],firstneigha2ia_index[i],x_ref[i],atom_displacement[i],
        defGradient,Q,stretchTensor,strainTensor,NSD,vonMisesStrain,volStrain);

    if (nsd_flag) atomNSD[i] = NSD;
    atomVonMisesStrain[i] = vonMisesStrain;
    atomVolStrain[i] = volStrain;
  }
  for (int i = 0; i < nelocal; i++) 
    for (int j = 0; j < npe[etype[i]]; j++) {
      int jj = e2nlist[i] + j;
      if (deformation_gradient_flag) defGradient = nodeF[i][j];
      else defGradient = NULL;
      if (rotation_flag) Q = nodeQ[i][j];
      else Q = NULL;
      if (stretch_tensor_flag) stretchTensor = nodeStretchTensor[i][j];
      else stretchTensor = NULL;
      if (strain_tensor_flag) strainTensor = nodeStrainTensor[i][j];
      else strainTensor = NULL;
      
      if (element->tag[i] == 178)
        debug = j+1;
      nodeValidity[i][j] = compute_strain(numneighn2a[jj],firstneighn2a[jj],numneighn2ia[jj],
          firstneighn2ia[jj],firstneighn2ia_index[jj],nodex_ref[i][j],node_displacement[i][j],
          defGradient,Q,stretchTensor,strainTensor,NSD,vonMisesStrain,volStrain);
      if (element->tag[i] == 178)
        debug = 0;

      if (nodeValidity[i][j]) {
        if (nsd_flag) nodeNSD[i][j] = NSD;
        nodeVonMisesStrain[i][j] = vonMisesStrain;
        nodeVolStrain[i][j] = volStrain;
      }
    }
}

/* ----------------------------------------------------------------------
   perform atomic strain calculation for an atom/node
   ------------------------------------------------------------------------- */

int AtomicStrain::compute_strain(int ajnum, int *ajlist, int iajnum, int *iajlist,
    int *iajindexlist, double *coord_ref, double *idisp, double *defGradient, double *quad, 
    double *stretchTensor, double *strainTensor, double &NSD, double &vonMisesStrain, double &volStrain)
{
  double del_ref[3],del_cur[3],*jdisp;
  double V[3][3],W[3][3];
  int i,j,k,l,node,jj;
  int ietype,iintpl,jetype,jintpl;

  double **x_ref = atom->x;
  double **x_cur = atom->x_current;
  double ***nodex_ref = element->nodex;
  double ***nodex_cur = element->nodex_current;
  int *npe = element->npe;
  int *etype = element->etype;
  double ***shape_array = element->shape_array;

  V[0][0] = W[0][0] = V[0][1] = W[0][1] = V[0][2] = W[0][2] = 0.0; 
  V[1][0] = W[1][0] = V[1][1] = W[1][1] = V[1][2] = W[1][2] = 0.0; 
  V[2][0] = W[2][0] = V[2][1] = W[2][1] = V[2][2] = W[2][2] = 0.0; 

  double sumSqDist = 0;
  // loop through neighbors to tally V and W
  // calculate current seperation vector from displacement vectors 
  // instead of current positions (due to wrapping)
  // del_cur = del_ref + jdisp - idisp

  // loop through neighboring atoms

  for (jj = 0; jj < ajnum; jj++) {

    j = ajlist[jj];
    jdisp = atom_displacement[j]; 

    for (k = 0; k < 3; k++) {
      del_ref[k] = x_ref[j][k] - coord_ref[k];
      del_cur[k] = del_ref[k] + jdisp[k] - idisp[k];
    }
    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++) {
        V[k][l] += del_ref[l] * del_ref[k];
        W[k][l] += del_ref[l] * del_cur[k];
      }
    sumSqDist += del_ref[0]*del_ref[0] 
      + del_ref[1]*del_ref[1] 
      + del_ref[2]*del_ref[2];
  }

  // loop through neighboring interpolated atoms

  for (jj = 0; jj < iajnum; jj++) {
    j = iajlist[jj];
    jetype = etype[j];
    jintpl = iajindexlist[jj];

    // interpolate positions at site J

    del_ref[0] = -coord_ref[0]; del_ref[1] = -coord_ref[1]; del_ref[2] = -coord_ref[2];
    for (node = 0; node < npe[jetype]; node++) 
      for (k = 0; k < 3; k++) 
        del_ref[k] += shape_array[jetype][jintpl][node]*nodex_ref[j][node][k];
    sub3(del_ref,idisp,del_cur);
    for (node = 0; node < npe[jetype]; node++) 
      for (k = 0; k < 3; k++) 
        del_cur[k] += shape_array[jetype][jintpl][node]*node_displacement[j][node][k];
//    if (debug) // && (fabs(del_cur[2]) >10 || fabs(del_cur[1]) >10 || fabs(del_cur[0]) >10))
//      printf("me = %d inode = %d jtag = %d j = %d jintpl = %d del_ref = %g %g %g del_cur = %g %g %g jdisp = %g %g %g\n",me,debug-1,element->tag[j],j,jintpl
//          ,del_ref[0],del_ref[1],del_ref[2]
//          ,del_cur[0],del_cur[1],del_cur[2]
//          ,del_cur[0]-del_ref[0]+idisp[0]
//          ,del_cur[1]-del_ref[1]+idisp[1]
//          ,del_cur[2]-del_ref[2]+idisp[2]
//          );

    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++) {
        V[k][l] += del_ref[l] * del_ref[k];
        W[k][l] += del_ref[l] * del_cur[k];
      }
    sumSqDist += del_ref[0]*del_ref[0] 
      + del_ref[1]*del_ref[1] 
      + del_ref[2]*del_ref[2];
  }

  int numneigh = ajnum + iajnum;
  double detThreshold = sumSqDist * 1e-12;

  if (domain->dimension == 2) {
    V[2][2] = W[2][2] = 1;
    V[0][2] = V[1][2] = V[2][0] = V[2][1] = 0;
    W[0][2] = W[1][2] = W[2][0] = W[2][1] = 0;
  }

  // check for invalid atoms (too little neighbors or V/W is non-invertible)
  // set all values to zero and skip to next atom

  double detV = det3(V);
  double detW = det3(W);


  if (numneigh < 2 || (domain->dimension == 3 && numneigh < 3) || 
      fabs(detV) < detThreshold || fabs(detW) < detThreshold) {
    vonMisesStrain = 0;
    volStrain = 0;
    if (deformation_gradient_flag) {
      defGradient[0] = defGradient[1] = defGradient[2] = 0;
      defGradient[3] = defGradient[4] = defGradient[5] = 0;
      defGradient[6] = defGradient[7] = defGradient[8] = 0;
    }
    if (rotation_flag) 
      quad[0] = quad[1] = quad[2] = quad[3] = 0;
    if (nsd_flag) NSD = 0;
    if (strain_tensor_flag) {
      strainTensor[0] = strainTensor[1] = 0;
      strainTensor[2] = strainTensor[3] = 0;
      strainTensor[4] = strainTensor[5] = 0;
    }
    if (stretch_tensor_flag) {
      stretchTensor[0] = stretchTensor[1] = 0;
      stretchTensor[2] = stretchTensor[3] = 0;
      stretchTensor[4] = stretchTensor[5] = 0;
    }
    return 0;
  }

  // Calculate deformation gradient tensor F

  double invV[3][3],F[3][3];
  invert3(V,invV); 
  times3(W,invV,F);

  if (deformation_gradient_flag) 
    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++)
        defGradient[k*3+l] = F[l][k];

  // Polar decomposition F = RU

  if (rotation_flag || stretch_tensor_flag) {
    double R[9],U[9],Farr[9];

    // convert F to array (row by row)

    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++)
        Farr[k*3+l] = F[k][l];

    polar_decomposition_3x3(Farr,false,R,U);
    if (rotation_flag) {

      // If F contains a reflection, R will not be a pure rotation matrix and the
      // conversion to a quaternion below will fail.
      // Thus, in the rather unlikely case that F contains a reflection, we simply flip the
      // R matrix to make it a pure rotation.

      if (matrix_determinant_3x3(R) < 0) {
        for (k = 0; k < 9; k++)
          R[k] = -R[k];
      }
      rotation_matrix_to_quaternion(R,quad);
    }
    if (stretch_tensor_flag) {
      stretchTensor[0] = U[0];
      stretchTensor[1] = U[4];
      stretchTensor[2] = U[8];
      stretchTensor[3] = U[1];
      stretchTensor[4] = U[2];
      stretchTensor[5] = U[5];
    }
  }

  // Calculate strain tensor

  double E[3][3];
  transpose_times3(F,F,E);
  E[0][0] -= 1;
  E[1][1] -= 1;
  E[2][2] -= 1;
  scalar_times3(0.5,E);
  double Sxx = E[0][0];
  double Syy = E[1][1];
  double Szz = E[2][2];
  double Sxy = E[0][1];
  double Sxz = E[0][2];
  double Syz = E[1][2];
  if (strain_tensor_flag) {
    strainTensor[0] = Sxx;
    strainTensor[1] = Syy;
    strainTensor[2] = Szz;
    strainTensor[3] = Sxy;
    strainTensor[4] = Sxz;
    strainTensor[5] = Syz;
  }

  // Calculate nonaffine squared displacement

  if (nsd_flag) {
    NSD = 0;
    double nonaffDisp[3];

    // loop through neighboring atoms

    for (jj = 0; jj < ajnum; jj++) {
      j = ajlist[jj];
      jdisp = atom_displacement[j]; 
      for (k = 0; k < 3; k++) {
        del_ref[k] = x_ref[j][k] - coord_ref[k];
        del_cur[k] = del_ref[k] + jdisp[k] - idisp[k];
      }
      matvec(F,del_ref,nonaffDisp);
      sub3(nonaffDisp,del_cur,nonaffDisp);
      NSD += lensq3(nonaffDisp);
    }

    // loop through neighboring interpolated atoms

    for (jj = 0; jj < iajnum; jj++) {
      j = iajlist[jj];
      jetype = etype[j];
      jintpl = iajindexlist[jj];

      // interpolate positions at site J

      del_ref[0] = -coord_ref[0]; del_ref[1] = -coord_ref[1]; del_ref[2] = -coord_ref[2];
      for (node = 0; node < npe[jetype]; node++) 
        for (k = 0; k < 3; k++) 
          del_ref[k] += shape_array[jetype][jintpl][node]*nodex_ref[j][node][k];
      sub3(del_ref,idisp,del_cur);
      for (node = 0; node < npe[jetype]; node++) 
        for (k = 0; k < 3; k++) 
          del_cur[k] += shape_array[jetype][jintpl][node]*node_displacement[j][node][k];
      matvec(F,del_ref,nonaffDisp);
      sub3(nonaffDisp,del_cur,nonaffDisp);
      NSD += lensq3(nonaffDisp);
    }
  }

  // Calculate von Mises shear strain

  double xydiff = Sxx - Syy;
  if (domain->dimension == 3) {
    double xzdiff = Sxx - Szz;
    double yzdiff = Syy - Szz;
    vonMisesStrain = sqrt(Sxy*Sxy + Sxz*Sxz + Syz*Syz +
        (xydiff*xydiff + xzdiff*xzdiff + yzdiff*yzdiff)/6.0);
  } else vonMisesStrain = sqrt(Sxy*Sxy + xydiff*xydiff/2.0);
  if (!ISFINITE(vonMisesStrain))
    vonMisesStrain = 0;

//  if (debug)
//    printf("vonMises = %g\n",vonMisesStrain);
  // Calculate volumetric component

  if (domain->dimension == 3)
    volStrain = (Sxx + Syy + Szz)/3.0;
  else
    volStrain = (Sxx + Syy)/2.0;
  if (!ISFINITE(volStrain))
    volStrain = 0;
  return 1;
}

/* ----------------------------------------------------------------------
   write output file
   ------------------------------------------------------------------------- */

void AtomicStrain::write_output(char *filename)
{
  if (me == 0) {
    if (screen) fprintf(screen,"Writing output file for atomic_strain command '%s'...\n",filename);
    if (logfile) fprintf(logfile,"Writing output file for atomic_strain command '%s'...\n",filename);
  }

  if (output_format != ATOM) {

    if (wrap_flag) {
      for (int i = 0; i < atom->nlocal; i++) {
        if (output_configuration == REFERENCE) domain->remap(atom->x[i]);
        else domain->remap_current(atom->x_current[i]);
      }
      int *npe = element->npe;
      int *etype = element->etype;
      double coord[3];
      for (int i = 0; i < element->nlocal; i++) {
        if (output_configuration == REFERENCE) 
          domain->remap(element->x[i],element->nodex[i],npe[etype[i]]);
        else {
          coord[0] = coord[1] = coord[2] = 0;
          for (int j = 0; j < npe[etype[i]]; j++) {
            coord[0] += element->nodex_current[i][j][0];
            coord[1] += element->nodex_current[i][j][1];
            coord[2] += element->nodex_current[i][j][2];
          }
          coord[0] /= npe[etype[i]];
          coord[1] /= npe[etype[i]];
          coord[2] /= npe[etype[i]];
          domain->remap_current(coord,element->nodex_current[i],npe[etype[i]]);
        }
      }
    }
    size_one = 6 + 2*debug + 3*element_info_flag + 9*deformation_gradient_flag + 4*rotation_flag + nsd_flag + 6*(strain_tensor_flag + stretch_tensor_flag);
    write_header_tecplot(filename);
    if (atom->natoms)
      write_tecplot_atoms();
    if (element->nelements) {
      write_tecplot_nodes();
      write_tecplot_element_connectivity();
    }
  } else {
    size_one = 8 + 2*debug + 4*element_info_flag + 9*deformation_gradient_flag + 4*rotation_flag + nsd_flag + 6*(strain_tensor_flag + stretch_tensor_flag);
    int nalocal = atom->nlocal;
    int nelocal = element->nlocal;
    int *npe = element->npe;
    int *etype = element->etype;
    int *nintpl = element->nintpl;
    bigint ntotal_local = nalocal;
    bigint ntotal;
    for (int i = 0; i < nelocal; i++)
      ntotal_local += nintpl[etype[i]];
    MPI_Allreduce(&ntotal_local,&ntotal,1,MPI_CAC_BIGINT,MPI_SUM,world);
    //printf(" me = %d ntotal_local = " BIGINT_FORMAT " ntotal = " BIGINT_FORMAT"\n",me,ntotal_local,ntotal);

    write_header_atom(filename,ntotal);

    int tmp,nlines;
    grow_buf(ntotal_local*size_one);
    pack_atom_all();
    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = ntotal_local;
        write_lines_atom(nlines,buf);
      }
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(buf,ntotal_local*size_one,MPI_DOUBLE,0,0,world);
    }
  }

  // close output file

  if (me == 0) {
    if (output_format == ASCII || output_format == ATOM) {
      fclose(fp);
      fp = NULL;
    } else {
      int success = TECEND142();
      if (success == -1)
        error->all(FLERR,"Cannot close output file");
    }
  }
}

/* ----------------------------------------------------------------------
   write header for atom output (LAMMPS dump style)
   ------------------------------------------------------------------------- */

void AtomicStrain::write_header_atom(char *file, bigint ndump) 
{
  // open data file and write headers

  if (me == 0) {

    char *columns = new char[500];
    char boundstr[9];
    double boxxlo,boxxhi;      // local copies of domain values
    double boxylo,boxyhi;      // lo/hi are bounding box for triclinic
    double boxzlo,boxzhi;
    double boxxy,boxxz,boxyz;

    if (domain->triclinic == 0) {
      if (output_configuration == REFERENCE) {
        boxxlo = domain->boxlo[0];
        boxxhi = domain->boxhi[0];
        boxylo = domain->boxlo[1];
        boxyhi = domain->boxhi[1];
        boxzlo = domain->boxlo[2];
        boxzhi = domain->boxhi[2];
      } else {
        boxxlo = domain->boxlo_current[0];
        boxxhi = domain->boxhi_current[0];
        boxylo = domain->boxlo_current[1];
        boxyhi = domain->boxhi_current[1];
        boxzlo = domain->boxlo_current[2];
        boxzhi = domain->boxhi_current[2];
      }
    } else {
      if (output_configuration == REFERENCE) {
        boxxlo = domain->boxlo_bound[0];
        boxxhi = domain->boxhi_bound[0];
        boxylo = domain->boxlo_bound[1];
        boxyhi = domain->boxhi_bound[1];
        boxzlo = domain->boxlo_bound[2];
        boxzhi = domain->boxhi_bound[2];
        boxxy = domain->xy;
        boxxz = domain->xz;
        boxyz = domain->yz;
      } else {
        boxxlo = domain->boxlo_current_bound[0];
        boxxhi = domain->boxhi_current_bound[0];
        boxylo = domain->boxlo_current_bound[1];
        boxyhi = domain->boxhi_current_bound[1];
        boxzlo = domain->boxlo_current_bound[2];
        boxzhi = domain->boxhi_current_bound[2];
        boxxy = domain->xy_current;
        boxxz = domain->xz_current;
        boxyz = domain->yz_current;
      }
    }

    fp = fopen(file, "w");
    if (fp == NULL) {
      char str[128];
      sprintf(str, "Cannot open data file %s",file);
      error->one(FLERR,str);
    }
    columns = strcpy(columns,"ID type x y z");
    if (element_info_flag)
      strcat(columns, " etype eid inode iintpl");
    if (debug)
      strcat(columns, " numneigha numneighia");
    strcat(columns," validity vonMises volStrain");
    if (deformation_gradient_flag)
      strcat(columns, " F11 F21 F31 F12 F22 F32 F13 F23 F33");
    if (rotation_flag)
      strcat(columns, " Rw Rx Ry Rz");
    if (nsd_flag)
      strcat(columns, " nsd");
    if (strain_tensor_flag)
      strcat(columns, " StrainXX StrainYY StrainZZ StrainXY StrainXZ StrainYZ");
    if (stretch_tensor_flag)
      strcat(columns, " StretchXX StretchYY StretchZZ StretchXY StretchXZ StretchYZ");

    domain->boundary_string(boundstr);
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,BIGINT_FORMAT "\n",ndump);
    if (domain->triclinic) {
      fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
      fprintf(fp,"%-1.16e %-1.16e\n",boxxlo,boxxhi);
      fprintf(fp,"%-1.16e %-1.16e\n",boxylo,boxyhi);
      fprintf(fp,"%-1.16e %-1.16e\n",boxzlo,boxzhi);
    } else {
      fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxxlo,boxxhi,boxxy);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxylo,boxyhi,boxxz);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxzlo,boxzhi,boxyz);
    }
    fprintf(fp,"ITEM: ATOMS %s\n",columns);
    delete [] columns;
  }
}


/* ----------------------------------------------------------------------
   write header for tecplot output
   ------------------------------------------------------------------------- */

void AtomicStrain::write_header_tecplot(char *file) 
{
  // open data file and write headers

  if (me == 0) {

    char *columns = new char[500];

    if (output_format == ASCII) {
      fp = fopen(file, "w");
      if (fp == NULL) {
        char str[128];
        sprintf(str, "Cannot open data file %s",file);
        error->one(FLERR,str);
      }
      columns = strcpy(columns,"\"x\", \"y\", \"z\"");
      if (element_info_flag)
        strcat(columns, ", \"ID\", \"etype\", \"ctype\"");
      if (debug)
        strcat(columns, ", \"numneigha\", \"numneighia\"");
      strcat(columns, ", \"validity\", \"vonMises\", \"volStrain\"");
      if (deformation_gradient_flag)
        strcat(columns, ", \"F11\", \"F21\", \"F31\", \"F12\", \"F22\", \"F32\", \"F13\", \"F23\", \"F33\"");
      if (rotation_flag)
        strcat(columns, ", \"Rw\", \"Rx\", \"Ry\", \"Rz\"");
      if (nsd_flag)
        strcat(columns, ", \"nsd\"");
      if (strain_tensor_flag)
        strcat(columns, ", \"StrainXX\", \"StrainYY\", \"StrainZZ\", \"StrainXY\", \"StrainXZ\", \"StrainYZ\"");
      if (stretch_tensor_flag)
        strcat(columns, ", \"StretchXX\", \"StretchYY\", \"StretchZZ\", \"StretchXY\", \"StretchXZ\", \"StretchYZ\"");

      fprintf(fp, "title=\"Strain Calculation from CAC\"\n");
      fprintf(fp, "variables = %s\n",columns);
    } else {
      soltime = update->ntimestep;
      dummy = 0;
      dummy1 = 1;
      columns = strcpy(columns,"x y z");
      if (element_info_flag)
        strcat(columns, " ID etype ctype");
      if (debug)
        strcat(columns, " numneigha numneighia");
      strcat(columns," validity vonMises volStrain");
      if (deformation_gradient_flag)
        strcat(columns, " F11 F21 F31 F12 F22 F32 F13 F23 F33");
      if (rotation_flag)
        strcat(columns, " Rw Rx Ry Rz");
      if (nsd_flag)
        strcat(columns, " nsd");
      if (strain_tensor_flag)
        strcat(columns, " StrainXX StrainYY StrainZZ StrainXY StrainXZ StrainYZ");
      if (stretch_tensor_flag)
        strcat(columns, " StretchXX StretchYY StretchZZ StretchXY StretchXZ StretchYZ");

      valuelocation = new int[size_one];
      for (int i = 0; i < size_one; i++) 
        valuelocation[i] = 1;

      int success = TECINI142((char *) "Strain Calculation from CAC",
          columns,
          file,
          (char *) ".",
          &output_format,
          &dummy,
          &debug,
          &dummy1);
      if (success == -1) {
        char str[128];
        sprintf(str,"Cannot open data file %s",file);
        error->one(FLERR,str);
      }
    }
    delete [] columns;
  }
}

/* ----------------------------------------------------------------------
   write atom section for tecplot
   ------------------------------------------------------------------------- */

void AtomicStrain::write_tecplot_atoms() 
{
  int nlocal = atom->nlocal;
  if (output_format == ASCII) {
    if (me == 0) 
      fprintf(fp,"zone t = \"Discrete Atoms\", datapacking = point\n");
    int tmp,nlines;
    grow_buf(nlocal*size_one);
    pack_atom_tecplot();
    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = nlocal;
        write_lines_tecplot(nlines,buf);
      }
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(buf,nlocal*size_one,MPI_DOUBLE,0,0,world);
    }

  } else {
    if (me == 0) {
      zonetype = 0;       // Ordered zone
      int imax = atom->natoms;
      int jmax = 1;
      int kmax = 1;
      int success = TECZNE142((char*)"Discrete Atoms",&zonetype,&imax,&jmax,&kmax,
          &dummy,&dummy,&dummy,
          &soltime,
          &dummy,&dummy,
          &dummy1,
          &dummy,&dummy,
          &dummy1,&dummy1,&dummy1,
          NULL,NULL,NULL,
          &dummy);
      if (success == -1) error->one(FLERR,"Cannot write Discrete Atoms Zone");
    }

    grow_buf(nlocal);

    // communicate buf and dump out each value one by one

    valueflag = COORDINATE;
    for (ival = 0; ival < 3; ival++) {
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);
    }

    if (element_info_flag) {
      valueflag = ID;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);

      valueflag = ETYPE;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);

      valueflag = CTYPE;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);
    }

    if (debug) {
      valueflag = NUMNEIGHA;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);

      valueflag = NUMNEIGHIA;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);
    }

    valueflag = VALIDITY;
    pack_atom_tecplot();
    comm_buf_tecio(nlocal);

    valueflag = VONMISES;
    pack_atom_tecplot();
    comm_buf_tecio(nlocal);

    valueflag = VOLSTRAIN;
    pack_atom_tecplot();
    comm_buf_tecio(nlocal);

    if (deformation_gradient_flag) {
      valueflag = DEFGRAD;
      for (ival = 0; ival < 9; ival++) {
        pack_atom_tecplot();
        comm_buf_tecio(nlocal);
      }
    }

    if (rotation_flag) {
      valueflag = ROTATION;
      for (ival = 0; ival < 4; ival++) {
        pack_atom_tecplot();
        comm_buf_tecio(nlocal);
      }
    }

    if (nsd_flag) {
      valueflag = NSD;
      pack_atom_tecplot();
      comm_buf_tecio(nlocal);
    }

    if (strain_tensor_flag) {
      valueflag = STRAINTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_atom_tecplot();
        comm_buf_tecio(nlocal);
      }
    }

    if (stretch_tensor_flag) {
      valueflag = STRETCHTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_atom_tecplot();
        comm_buf_tecio(nlocal);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   write element section for tecplot
   ------------------------------------------------------------------------- */

void AtomicStrain::write_tecplot_nodes() 
{
  int nlocal = element->nlocal;
  int nnode_local = 0;
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < nlocal; i++)
    nnode_local += npe[etype[i]];
  int nnodes;
  MPI_Allreduce(&nnode_local,&nnodes,1,MPI_INT,MPI_SUM,world);
  int nelements = element->nelements;

  if (output_format == ASCII) {
    int tmp,nlines;
    if (me == 0) {
      fprintf(fp,"zone t = \"Coarse Element\"");
      fprintf(fp," n = %d e = %d",nnodes,nelements);
      if (domain->dimension == 3) 
        fprintf(fp," datapacking = point, zonetype = febrick\n");
      else 
        fprintf(fp," datapacking = point, zonetype = fequadrilateral\n");
    }
    grow_buf(nnode_local*size_one);
    pack_node_tecplot();
    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
        } else nlines = nnode_local;
        write_lines_tecplot(nlines,buf);
      }
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(buf,nnode_local*size_one,MPI_DOUBLE,0,0,world);
    }
  } else {
    if (me == 0) {
      if (domain->dimension == 3) zonetype = 5;    // Brick zone
      else zonetype = 3;                           // Quadrilateral zone
      int success = TECZNE142((char *) "Coarse Elements",&zonetype,&nnodes,&nelements,
          &dummy,&dummy,&dummy,&dummy,
          &soltime,
          &dummy,&dummy,
          &dummy1,
          &dummy,&dummy,
          0,0,0,
          NULL,valuelocation,NULL,
          &dummy);
      if (success == -1) error->one(FLERR,"Cannot write Coarse Element zone");
    }

    grow_buf(nnode_local);

    // communicate buf and dump out each value one by one

    valueflag = COORDINATE;
    for (ival = 0; ival < 3; ival++) {
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);
    }

    if (element_info_flag) {
      valueflag = ID;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);

      valueflag = ETYPE;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);

      valueflag = CTYPE;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);
    }

    if (debug) {
      valueflag = NUMNEIGHA;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);

      valueflag = NUMNEIGHIA;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);
    }
    valueflag = VALIDITY;
    pack_node_tecplot();
    comm_buf_tecio(nnode_local);

    valueflag = VONMISES;
    pack_node_tecplot();
    comm_buf_tecio(nnode_local);

    valueflag = VOLSTRAIN;
    pack_node_tecplot();
    comm_buf_tecio(nnode_local);

    if (deformation_gradient_flag) {
      valueflag = DEFGRAD;
      for (ival = 0; ival < 9; ival++) {
        pack_node_tecplot();
        comm_buf_tecio(nnode_local);
      }
    }

    if (rotation_flag) {
      valueflag = ROTATION;
      for (ival = 0; ival < 4; ival++) {
        pack_node_tecplot();
        comm_buf_tecio(nnode_local);
      }
    }

    if (nsd_flag) {
      valueflag = NSD;
      pack_node_tecplot();
      comm_buf_tecio(nnode_local);
    }

    if (strain_tensor_flag) {
      valueflag = STRAINTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_node_tecplot();
        comm_buf_tecio(nnode_local);
      }
    }

    if (stretch_tensor_flag) {
      valueflag = STRETCHTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_node_tecplot();
        comm_buf_tecio(nnode_local);
      }
    }
  }
}
/* ----------------------------------------------------------------------
   write element section for tecplot
   ------------------------------------------------------------------------- */

void AtomicStrain::write_tecplot_element_connectivity()
{
  // dump node connectivity
  // use int buf for node connectivity 
  // due to tecio function requirement
  // nlines for TECNODE142 fuctions = # of values to write

  int nlocal = element->nlocal;
  int num_me,nlines,tmp;
  int size_one_line;
  char *format;
  if (domain->dimension == 3) {
    size_one_line = 8;
  } else {
    size_one_line = 4;
  }
  num_me = nlocal*size_one_line;
  grow_ibuf(num_me);

  pack_element_connectivity();
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(ibuf,maxibuf,MPI_INT,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_INT,&nlines);
      } else nlines = num_me;
      if (output_format == ASCII) {
        nlines /= size_one_line;
        int m = 0;
        for (int i = 0; i < nlines; i++) {
          for (int j = 0; j < size_one_line; j++)
            fprintf(fp,"%d ",ibuf[m++]);
          fprintf(fp,"\n");
        }
      } else {
        int success = TECNODE142(&nlines,ibuf);
        if (success == -1) error->one(FLERR,"Cannot write to dump file");
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(ibuf,num_me,MPI_INT,0,0,world);
  }

}

/*------------------------------------------------------------------------------------------*/

void AtomicStrain::pack_element_connectivity()
{
  int nlocal = element->nlocal;
  tagint *tag = element->tag;
  tagint **nodetag = element->nodetag;

  int *npe = element->npe;
  int *etype = element->etype;
  int **element_clusters = element->element_clusters;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;
  int inpe,itype;

  element->reset_node_tags(0,0);

  for (int i = 0; i < nlocal; i++) {
    itype = etype[i];
    if (element_shape_ids[itype] == Element::QUADRILATERAL) {
      for (int j = 0; j < 4; j++) 
        ibuf[m++] = nodetag[i][j];
      if (domain->dimension == 3) {
        for (int j = 0; j < 4; j++) 
          ibuf[m++] = nodetag[i][j];
      }

    } else if (element_shape_ids[itype] == Element::TRIANGLE) {

      // node connectivity scheme: 1 2 3 3

      for (int j = 0; j < 3; j++) 
        ibuf[m++] = nodetag[i][j];
      ibuf[m++] = nodetag[i][2];

    } else if (element_shape_ids[itype] == Element::HEXAHEDRON) {

      for (int j = 0; j < 8; j++) 
        ibuf[m++] = nodetag[i][j];

    } else if (element_shape_ids[itype] == Element::PYRAMID) {

      // node connectivity scheme: 1 2 3 4 5 5 5 5

      for (int j = 0; j < 4; j++) 
        ibuf[m++] = nodetag[i][j];
      for (int j = 4; j < 8; j++)
        ibuf[m++] = nodetag[i][4];

    } else if (element_shape_ids[itype] == Element::TETRAHEDRON) {

      // node connectivity scheme: 1 2 3 3 4 4 4 4

      for (int j = 0; j < 3; j++) 
        ibuf[m++] = nodetag[i][j];
      ibuf[m++] = nodetag[i][2];
      for (int j = 4; j < 8; j++) 
        ibuf[m++] = nodetag[i][3];

    } else if (element_shape_ids[itype] == Element::WEDGE) {

      // node connectivity scheme: 1 2 3 3 4 5 6 6

      for (int j = 0; j < 3; j++) 
        ibuf[m++] = nodetag[i][j];
      for (int j = 3; j < 7; j++) 
        ibuf[m++] = nodetag[i][j-1];
      ibuf[m++] = nodetag[i][5];

    }
  }
} 


/*----------------------------------------------------------------------------------------------------*/

void AtomicStrain::grow_buf(int num_me)
{
  int nmax;
  MPI_Allreduce(&num_me,&nmax,1,MPI_INT,MPI_MAX,world);

  if (nmax > maxbuf) {
    if ((bigint) nmax > MAXSMALLINT) {
      char errstr[100];  
      sprintf(errstr,"Too much per-proc info for atomic_strain command nmax = %d",nmax);
      error->all(FLERR,errstr);
    }

    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf,"dump:buf");
  }
}

/*----------------------------------------------------------------------------------------------------*/

void AtomicStrain::grow_ibuf(int num_me)
{
  int nmax;
  MPI_Allreduce(&num_me,&nmax,1,MPI_INT,MPI_MAX,world);

  if (nmax > maxibuf) {
    if ((bigint) nmax > MAXSMALLINT) {
      char errstr[100];  
      sprintf(errstr,"Too much per-proc info for atomic_strain command nmax = %d",nmax);
      error->all(FLERR,errstr);
    }

    maxibuf = nmax;
    memory->destroy(ibuf);
    memory->create(ibuf,maxibuf,"dump:ibuf");
  }
}

/*----------------------------------------------------------------------------------------------------*/

void AtomicStrain::comm_buf_tecio(int num_me)
{
  int nlines,tmp,success;
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
      } else nlines = num_me;
      success = TECDAT142(&nlines,buf,&dummy1);
      if (success == -1) error->one(FLERR,"Cannot write to dump file");
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,num_me,MPI_DOUBLE,0,0,world);
  }
}

/*----------------------------------------------------------------------------------------------------*/

void AtomicStrain::pack_atom_tecplot()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x;
  if (output_configuration == REFERENCE) x = atom->x;
  else x = atom->x_current;

  int m = 0; 
  if (output_format == ASCII) {
    for (int i = 0; i < nlocal; i++) {
      if (wrap_flag) {
        if (output_configuration == REFERENCE) domain->remap(x[i]);
        else domain->remap_current(x[i]);
      }
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (element_info_flag) {
        buf[m++] = tag[i];
        buf[m++] = 0;
        buf[m++] = type[i];
      }
      if (debug) {
        buf[m++] = list->numneigha2a[i];
        buf[m++] = list->numneigha2ia[i];
      }
      buf[m++] = atomValidity[i];
      buf[m++] = atomVonMisesStrain[i];
      buf[m++] = atomVolStrain[i];
      if (deformation_gradient_flag) {
        buf[m++] = atomF[i][0];
        buf[m++] = atomF[i][1];
        buf[m++] = atomF[i][2];
        buf[m++] = atomF[i][3];
        buf[m++] = atomF[i][4];
        buf[m++] = atomF[i][5];
        buf[m++] = atomF[i][6];
        buf[m++] = atomF[i][7];
        buf[m++] = atomF[i][8];
      }
      if (rotation_flag) {
        buf[m++] = atomQ[i][0];
        buf[m++] = atomQ[i][1];
        buf[m++] = atomQ[i][2];
        buf[m++] = atomQ[i][3];
      }
      if (nsd_flag)
        buf[m++] = atomNSD[i];
      if (strain_tensor_flag) {
        buf[m++] = atomStrainTensor[i][0];
        buf[m++] = atomStrainTensor[i][1];
        buf[m++] = atomStrainTensor[i][2];
        buf[m++] = atomStrainTensor[i][3];
        buf[m++] = atomStrainTensor[i][4];
        buf[m++] = atomStrainTensor[i][5];
      }
      if (stretch_tensor_flag) {
        buf[m++] = atomStretchTensor[i][0];
        buf[m++] = atomStretchTensor[i][1];
        buf[m++] = atomStretchTensor[i][2];
        buf[m++] = atomStretchTensor[i][3];
        buf[m++] = atomStretchTensor[i][4];
        buf[m++] = atomStretchTensor[i][5];
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (valueflag == COORDINATE)
        buf[m++] = x[i][ival];
      else if (valueflag == ID) 
        buf[m++] = tag[i];
      else if (valueflag == ETYPE) 
        buf[m++] = 0;
      else if (valueflag == CTYPE) 
        buf[m++] = type[i];
      else if (valueflag == NUMNEIGHA) 
        buf[m++] = list->numneigha2a[i];
      else if (valueflag == NUMNEIGHIA) 
        buf[m++] = list->numneigha2ia[i];
      else if (valueflag == VALIDITY) 
        buf[m++] = atomValidity[i];
      else if (valueflag == VONMISES) 
        buf[m++] = atomVonMisesStrain[i];
      else if (valueflag == VOLSTRAIN)
        buf[m++] = atomVolStrain[i];
      else if (valueflag == DEFGRAD)
        buf[m++] = atomF[i][ival];
      else if (valueflag == ROTATION)
        buf[m++] = atomQ[i][ival];
      else if (valueflag == NSD) 
        buf[m++] = atomNSD[i];
      else if (valueflag == STRAINTENSOR)
        buf[m++] = atomStrainTensor[i][ival];
      else if (valueflag == STRETCHTENSOR)
        buf[m++] = atomStretchTensor[i][ival];
    }
  }
}

/*----------------------------------------------------------------------------------------------------*/

void AtomicStrain::pack_node_tecplot()
{
  int nlocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tagint *tag = element->tag;
  int *npe = element->npe;
  double ***nodex;
  double coord[3];
  int inpe;
  if (output_configuration == REFERENCE) nodex = element->nodex;
  else nodex = element->nodex_current;

  int m = 0; 
  int n = 0;
  if (output_format == ASCII) {
    for (int i = 0; i < nlocal; i++) {
      inpe = npe[etype[i]]; 
      for (int j = 0; j < inpe; j++) {
        buf[m++] = nodex[i][j][0];
        buf[m++] = nodex[i][j][1];
        buf[m++] = nodex[i][j][2];
        if (element_info_flag) {
          buf[m++] = tag[i];
          buf[m++] = etype[i];
          buf[m++] = ctype[i];
        }
        if (debug) {
          buf[m++] = list->numneighn2a[n];
          buf[m++] = list->numneighn2ia[n];
          n++;
        }
        buf[m++] = nodeValidity[i][j];
        buf[m++] = nodeVonMisesStrain[i][j];
        buf[m++] = nodeVolStrain[i][j];
        if (deformation_gradient_flag) {
          buf[m++] = nodeF[i][j][0];
          buf[m++] = nodeF[i][j][1];
          buf[m++] = nodeF[i][j][2];
          buf[m++] = nodeF[i][j][3];
          buf[m++] = nodeF[i][j][4];
          buf[m++] = nodeF[i][j][5];
          buf[m++] = nodeF[i][j][6];
          buf[m++] = nodeF[i][j][7];
          buf[m++] = nodeF[i][j][8];
        }
        if (rotation_flag) {
          buf[m++] = nodeQ[i][j][0];
          buf[m++] = nodeQ[i][j][1];
          buf[m++] = nodeQ[i][j][2];
          buf[m++] = nodeQ[i][j][3];
        }
        if (nsd_flag)
          buf[m++] = nodeNSD[i][j];
        if (strain_tensor_flag) {
          buf[m++] = nodeStrainTensor[i][j][0];
          buf[m++] = nodeStrainTensor[i][j][1];
          buf[m++] = nodeStrainTensor[i][j][2];
          buf[m++] = nodeStrainTensor[i][j][3];
          buf[m++] = nodeStrainTensor[i][j][4];
          buf[m++] = nodeStrainTensor[i][j][5];
        }
        if (stretch_tensor_flag) {
          buf[m++] = nodeStretchTensor[i][j][0];
          buf[m++] = nodeStretchTensor[i][j][1];
          buf[m++] = nodeStretchTensor[i][j][2];
          buf[m++] = nodeStretchTensor[i][j][3];
          buf[m++] = nodeStretchTensor[i][j][4];
          buf[m++] = nodeStretchTensor[i][j][5];
        }
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      inpe = npe[etype[i]]; 
      for (int j = 0; j < npe[etype[i]]; j++) {
        if (valueflag == COORDINATE)
          buf[m++] = nodex[i][j][ival];
        else if (valueflag == ID) 
          buf[m++] = tag[i];
        else if (valueflag == ETYPE) 
          buf[m++] = etype[i];
        else if (valueflag == CTYPE) 
          buf[m++] = ctype[i];
        else if (valueflag == NUMNEIGHA)
          buf[m++] = list->numneighn2a[n++];
        else if (valueflag == NUMNEIGHIA)
          buf[m++] = list->numneighn2ia[n++];
        else if (valueflag == VALIDITY) 
          buf[m++] = nodeValidity[i][j];
        else if (valueflag == VONMISES) 
          buf[m++] = nodeVonMisesStrain[i][j];
        else if (valueflag == VOLSTRAIN)
          buf[m++] = nodeVolStrain[i][j];
        else if (valueflag == DEFGRAD)
          buf[m++] = nodeF[i][j][ival];
        else if (valueflag == ROTATION)
          buf[m++] = nodeQ[i][j][ival];
        else if (valueflag == NSD) 
          buf[m++] = nodeNSD[i][j];
        else if (valueflag == STRAINTENSOR)
          buf[m++] = nodeStrainTensor[i][j][ival];
        else if (valueflag == STRETCHTENSOR)
          buf[m++] = nodeStretchTensor[i][j][ival];
      }
    }
  }
}

/*-------------------------------------------------------------------------------------------------*/

void AtomicStrain::write_lines_tecplot(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,"%g %g %g",mybuf[m],mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (element_info_flag) {
      fprintf(fp, TAGINT_FORMAT " %d %d",
          static_cast<tagint> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          static_cast<int> (mybuf[m+2]));
      m += 3;
    }
    if (debug) {
      fprintf(fp," %d %d",static_cast<int> (mybuf[m]),static_cast<int> (mybuf[m+1]));
      m += 2;
    }
    fprintf(fp," %d %g %g",static_cast<int> (mybuf[m]),mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (deformation_gradient_flag) {
      fprintf(fp," %g %g %g %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5],
          mybuf[m+6],mybuf[m+7],mybuf[m+8]);
      m += 9;
    }
    if (rotation_flag) {
      fprintf(fp," %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }

    if (nsd_flag) {
      fprintf(fp," %g",mybuf[m]);
      m += 1;
    }
    if (strain_tensor_flag) {
      fprintf(fp," %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (stretch_tensor_flag) {
      fprintf(fp," %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    fprintf(fp,"\n");
  }
}

/*------------------------------------------------------------------------------------------*/

void AtomicStrain::pack_atom_all()
{
  int m,n;
  int i,j,k,iintpl,inode,inode_local,ietype;
  double **x;
  if (output_configuration == REFERENCE) x = atom->x;
  else x = atom->x_current;

  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  tagint maxtag,maxtag_all;
  maxtag = 0;

  m = n = 0;
  for (i = 0; i < nlocal; i++) {
    if (wrap_flag) {
      if (output_configuration == REFERENCE) domain->remap(x[i]);
      else domain->remap_current(x[i]);
    }

    maxtag = MAX(maxtag,tag[i]);
    buf[m++] = tag[i];
    buf[m++] = type[i];
    buf[m++] = x[i][0];
    buf[m++] = x[i][1];
    buf[m++] = x[i][2];
    if (element_info_flag) {
      buf[m++] = 0;
      buf[m++] = 0;
      buf[m++] = -1;
      buf[m++] = -1;
    }
    if (debug) {
      buf[m++] = list->numneigha2a[i];
      buf[m++] = list->numneigha2ia[i];
    }
    buf[m++] = atomValidity[i];
    buf[m++] = atomVonMisesStrain[i];
    buf[m++] = atomVolStrain[i];
    if (deformation_gradient_flag) {
      buf[m++] = atomF[i][0];
      buf[m++] = atomF[i][1];
      buf[m++] = atomF[i][2];
      buf[m++] = atomF[i][3];
      buf[m++] = atomF[i][4];
      buf[m++] = atomF[i][5];
      buf[m++] = atomF[i][6];
      buf[m++] = atomF[i][7];
      buf[m++] = atomF[i][8];
    }
    if (rotation_flag) {
      buf[m++] = atomQ[i][0];
      buf[m++] = atomQ[i][1];
      buf[m++] = atomQ[i][2];
      buf[m++] = atomQ[i][3];
    }
    if (nsd_flag)
      buf[m++] = atomNSD[i];
    if (strain_tensor_flag) {
      buf[m++] = atomStrainTensor[i][0];
      buf[m++] = atomStrainTensor[i][1];
      buf[m++] = atomStrainTensor[i][2];
      buf[m++] = atomStrainTensor[i][3];
      buf[m++] = atomStrainTensor[i][4];
      buf[m++] = atomStrainTensor[i][5];
    }
    if (stretch_tensor_flag) {
      buf[m++] = atomStretchTensor[i][0];
      buf[m++] = atomStretchTensor[i][1];
      buf[m++] = atomStretchTensor[i][2];
      buf[m++] = atomStretchTensor[i][3];
      buf[m++] = atomStretchTensor[i][4];
      buf[m++] = atomStretchTensor[i][5];
    }
  }

  if (atom->natoms) {
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);
  } else maxtag_all = 0;

  // pack atoms from elements 

  nlocal = element->nlocal;
  double ***nodex;
  if (output_configuration == REFERENCE) nodex = element->nodex;
  else nodex = element->nodex_current;
  int *nintpl = element->nintpl;
  int *etype = element->etype;
  int *ctype = element->ctype;
  tag = element->tag;
  int *npe = element->npe;
  double ***shape_array = element->shape_array;
  double tmp;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int *e2nlist = list->e2nlist;

  // set IDs for interpolated atoms from element tag IDs to have consistent tag IDs across different simulations of the same system

  tagint itag = maxtag_all + 1;
  int maxintpl = element->maxintpl;
  for (i = 0; i < nlocal; i++) {
    ietype = etype[i];
    for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {

      buf[m++] = itag + (tag[i]-1)*maxintpl + iintpl;
      buf[m++] = ctype[i];
      for (j = 0; j < 3; j++) {
        tmp = 0.0;
        for (inode = 0; inode < npe[ietype]; inode++) 
          tmp += shape_array[ietype][iintpl][inode]*nodex[i][inode][j];
        buf[m++] = tmp;
      }
      if (wrap_flag) {
        if (output_configuration == REFERENCE) domain->remap(&buf[m-3]);
        else domain->remap_current(&buf[m-3]);
      }

      if (ia2i[ietype][iintpl] < 0) inode = -1;
      else inode = i2n[ietype][ia2i[ietype][iintpl]];
      if (element_info_flag) {
        buf[m++] = ietype;
        buf[m++] = tag[i];
        buf[m++] = inode;
        buf[m++] = iintpl;
      }
      if (debug) {
        if (inode >= 0) {
          inode_local = e2nlist[i]+inode;
          buf[m++] = list->numneighn2a[inode_local];
          buf[m++] = list->numneighn2ia[inode_local];
        } else {
          buf[m++] = -1;
          buf[m++] = -1;
        }
      }

      tmp = 0.0;
      for (inode = 0; inode < npe[ietype]; inode++) 
        tmp += shape_array[ietype][iintpl][inode]*nodeValidity[i][inode];
      buf[m++] = tmp;

      tmp = 0.0;
      for (inode = 0; inode < npe[ietype]; inode++) 
        tmp += shape_array[ietype][iintpl][inode]*nodeVonMisesStrain[i][inode];
      buf[m++] = tmp;

      tmp = 0.0;
      for (inode = 0; inode < npe[ietype]; inode++) 
        tmp += shape_array[ietype][iintpl][inode]*nodeVolStrain[i][inode];
      buf[m++] = tmp;

      if (deformation_gradient_flag) {
        for (j = 0; j < 9; j++) {
          tmp = 0.0;
          for (inode = 0; inode < npe[ietype]; inode++) 
            tmp += shape_array[ietype][iintpl][inode]*nodeF[i][inode][j];
          buf[m++] = tmp;
        }
      }
      if (rotation_flag) {
        for (j = 0; j < 4; j++) {
          tmp = 0.0;
          for (inode = 0; inode < npe[ietype]; inode++) 
            tmp += shape_array[ietype][iintpl][inode]*nodeQ[i][inode][j];
          buf[m++] = tmp;
        }
      }
      if (nsd_flag) {
        tmp = 0.0;
        for (inode = 0; inode < npe[ietype]; inode++) 
          tmp += shape_array[ietype][iintpl][inode]*nodeNSD[i][inode];
        buf[m++] = tmp;
      }
      if (strain_tensor_flag) {
        for (j = 0; j < 6; j++) {
          tmp = 0.0;
          for (inode = 0; inode < npe[ietype]; inode++) 
            tmp += shape_array[ietype][iintpl][inode]*nodeStrainTensor[i][inode][j];
          buf[m++] = tmp;
        }
      }
      if (stretch_tensor_flag) {
        for (j = 0; j < 6; j++) {
          tmp = 0.0;
          for (inode = 0; inode < npe[ietype]; inode++) 
            tmp += shape_array[ietype][iintpl][inode]*nodeStretchTensor[i][inode][j];
          buf[m++] = tmp;
        }
      }
    }
  }  
}

/*-------------------------------------------------------------------------------------------------*/

void AtomicStrain::write_lines_atom(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, TAGINT_FORMAT " %d %g %g %g",
        static_cast<tagint> (mybuf[m]),static_cast<int> (mybuf[m+1]),
        mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += 5;
    if (element_info_flag) {
      fprintf(fp," %d %d %d %d",
          static_cast<int> (mybuf[m]),
          static_cast<int> (mybuf[m+1]),
          static_cast<int> (mybuf[m+2]),
          static_cast<int> (mybuf[m+3]));
      m += 4;
    }
    if (debug) {
      fprintf(fp," %d %d",static_cast<int> (mybuf[m]),static_cast<int> (mybuf[m+1]));
      m += 2;
    }
    fprintf(fp," %d %g %g",static_cast<int> (mybuf[m]),mybuf[m+1],mybuf[m+2]);
    m += 3;
    if (deformation_gradient_flag) {
      fprintf(fp," %g %g %g %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5],
          mybuf[m+6],mybuf[m+7],mybuf[m+8]);
      m += 9;
    }
    if (rotation_flag) {
      fprintf(fp," %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3]);
      m += 4;
    }

    if (nsd_flag) {
      fprintf(fp," %g",mybuf[m]);
      m += 1;
    }
    if (strain_tensor_flag) {
      fprintf(fp," %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    if (stretch_tensor_flag) {
      fprintf(fp," %g %g %g %g %g %g",
          mybuf[m],mybuf[m+1],mybuf[m+2],
          mybuf[m+3],mybuf[m+4],mybuf[m+5]);
      m += 6;
    }
    fprintf(fp,"\n");
  }
}


