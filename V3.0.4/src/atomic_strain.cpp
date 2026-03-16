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
#include "group.h"
#include "region.h"
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
#define NSECTIONS 6       // change when add to header::section_keywords
#define EPSILON 1e-6

enum{PLT = 0, SZPLT = 1, ASCII = 2, ATOM = 3, NODE = 4, FULLATOM = 5};        // output file format ASCII (.dat), binary (.plt), or SZL (.szplt)
enum{COORDINATE, ID, ETYPE, CTYPE, VALIDITY, VONMISES, VOLSTRAIN, DEFGRAD, 
  ROTATIONQUAT, ROTATIONVECT, DISPLACEMENT, NSD, STRAINTENSOR, STRETCHTENSOR, NUMNEIGH};
enum{REFERENCE, CURRENT};
/*  ----------------------------------------------------------------------  */

AtomicStrain::AtomicStrain(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK * MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  balance = NULL;
  atomValidity = NULL;
  atomVonMisesStrain = NULL;
  atomVolStrain = NULL;
  atomF = NULL;
  atomQ = NULL;
  atomW = NULL;
  atomD = NULL;
  atomStretchTensor = NULL;
  atomStrainTensor = NULL;
  atomNSD = NULL;

  nodeValidity = NULL;
  nodeVonMisesStrain = NULL;
  nodeVolStrain = NULL;
  nodeF = NULL;
  nodeQ = NULL;
  nodeW = NULL;
  nodeD = NULL;
  nodeStretchTensor = NULL;
  nodeStrainTensor = NULL;
  nodeNSD = NULL;


  valuelocation = NULL;
  buf = NULL;
  ibuf = NULL;
  maxbuf = maxibuf = 0;
}

/*  ----------------------------------------------------------------------  */

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
  memory->destroy(atomW);
  memory->destroy(atomD);
  memory->destroy(atomStretchTensor);
  memory->destroy(atomStrainTensor);
  memory->destroy(atomNSD);

  memory->destroy(nodeValidity);
  memory->destroy(nodeVonMisesStrain);
  memory->destroy(nodeVolStrain);
  memory->destroy(nodeF);
  memory->destroy(nodeQ);
  memory->destroy(nodeW);
  memory->destroy(nodeD);
  memory->destroy(nodeStretchTensor);
  memory->destroy(nodeStrainTensor);
  memory->destroy(nodeNSD);


  memory->destroy(buf);
  memory->destroy(ibuf);
  delete [] valuelocation;
  delete balance;
}

/*  ----------------------------------------------------------------------  */

void AtomicStrain::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal atomic_strain command");

  nelements = natoms = 0;

  atom->atom_strain_flag = 1;

  noutputs = universe->inumeric(FLERR,arg[1]);
  if (noutputs < 1) error->all(FLERR, "Illegal atomic_strain command");
  output_format_list = new int[noutputs];
  output_configuration_list = new int[noutputs];
  output_filename_list = new char*[noutputs];

  if (narg < 2 + 3 * noutputs) error->all(FLERR, "Illegal atomic_strain command");
  int iarg = 2;
  for (int i = 0; i < noutputs; i++) {
    output_filename_list[i] = arg[iarg];  
    if (strcmp(arg[iarg+1], "plt") == 0) output_format_list[i] = PLT;
    else if (strcmp(arg[iarg+1], "szplt") == 0) output_format_list[i] = SZPLT;
    else if (strcmp(arg[iarg+1], "ascii") == 0) output_format_list[i] = ASCII;
    else if (strcmp(arg[iarg+1], "atom") == 0) output_format_list[i] = ATOM;
    else if (strcmp(arg[iarg+1], "node") == 0) output_format_list[i] = NODE;
    else if (strcmp(arg[iarg+1], "full/atom") == 0) output_format_list[i] = FULLATOM;
    else error->all(FLERR, "Illegal atomic_strain command");
    if (strcmp(arg[iarg+2], "current") == 0) output_configuration_list[i] = CURRENT;
    else if (strcmp(arg[iarg+2], "reference") == 0) output_configuration_list[i] = REFERENCE;
    else error->all(FLERR, "Illegal atomic_strain command");
    iarg += 3;
  }
  cut = 6;
  element_info_flag = 0;
  deformation_gradient_flag = 0;
  rotation_quat_flag = 0;
  rotation_vect_flag = 0;
  displacement_flag = 0;
  nsd_flag = 0;
  strain_tensor_flag = 0;
  stretch_tensor_flag = 0;
  wrap_flag = 1;
  debug = 0;
  average_flag = 1;
  region_id = NULL;
  no_elem_flag = no_atom_flag = 0;
  ntimestep = 0;

  double eweight = 0;
  int ref_file = 0;
  char *ref_file_name;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      cut = universe->numeric(FLERR, arg[iarg+1]);
      iarg += 2; 
    } else if (strcmp(arg[iarg], "def/grad") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) deformation_gradient_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) deformation_gradient_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "elem/info") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) element_info_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) element_info_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rotation/quat") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) rotation_quat_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) rotation_quat_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rotation/vect") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) rotation_vect_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) rotation_vect_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "displacement") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) displacement_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) displacement_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "nsd") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) nsd_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) nsd_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "stretch/tensor") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) stretch_tensor_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) stretch_tensor_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "strain/tensor") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) strain_tensor_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) strain_tensor_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      if (strcmp(arg[iarg+1], "yes") == 0) average_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) average_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ref/file") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      ref_file = 1;
      ref_file_name = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "eweight") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal atomic_strain command");
      eweight = universe->numeric(FLERR, arg[iarg+1]);
      if (eweight <= 0) error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "wrap") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) wrap_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) wrap_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "no/atom") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) no_atom_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) no_atom_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "no/element") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) no_elem_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) no_elem_flag = 0;
      else error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "region") == 0) {
      region_id = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "step") == 0) {
      ntimestep = universe->inumeric(FLERR,arg[iarg+1]); 
      if (ntimestep < 0) error->all(FLERR, "Illegal atomic_strain command");
      iarg += 2;
    } else error->all(FLERR, "Illegal atomic_strain command");

  }

  // error checks

  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR, "Cannot run 2d simulation with nonperiodic Z dimension");
  if (comm->style != 1) 
    error->all(FLERR, "atomic_strain command requires comm_style tiled for workload balancing");

  // box domain must exist

  //if (!domain->box_exist)
  //  error->all(FLERR, "Box domain not yet exist");

  // remove all atoms/elements

  atom->nlocal = 0;
  atom->natoms = 0;
  element->nlocal = 0;
  element->nelements = 0;

  // read input data
  // reference coordinates are wrapped
  // current coordinates are unwrapped

  read_input(arg[0]);
  if (ref_file) read_ref_file(ref_file_name);
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
  char *eweight_str = new char[20];
  newarg[0] = (char *) "1.0";
  newarg[1] = (char *) "rcb";
  newarg[2] = (char *) "eweight";
  if (eweight == 0) eweight = 100;
  sprintf(eweight_str, "%g", eweight); 
  newarg[3] = eweight_str;
  balance->command(4, newarg);

  delete [] newarg;
  delete [] eweight_str;

  // request a full atomic neighbor list for use by this command 
  // also sort out neighbors based on distance

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 2;
  neighbor->requests[irequest]->nodelist = 1;
  neighbor->requests[irequest]->gausslist = 0;
  neighbor->requests[irequest]->cut = 1;
  neighbor->requests[irequest]->cutoff = cut;

  neighbor->init();

  // manually set cutneighmax for border
  // use a 1.5 factor just to be safe

  neighbor->cutneighmax = cut * 1.5;
  comm->init();

  // setup domain, communication and neighboring
  // exchange and acquire ghosts
  // build standard neighbor list to setup bins & stencils

  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal, atom->x);
    domain->x2lamda(element->nlocal, element->x);
    domain->nodex2lamda(element->nlocal, element->nodex);
    domain->x2lamda_current(atom->nlocal, atom->x_current);
    domain->nodex2lamda_current(element->nlocal, element->nodex_current);
  }
  comm->setup_exchange();
  comm->exchange();
  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 

  comm->borders();

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost, atom->x);
    domain->lamda2x(element->nlocal+element->nghost, element->x);
    domain->lamda2nodex(element->nlocal+element->nghost, element->nodex);
    domain->lamda2x_current(atom->nlocal+atom->nghost, atom->x_current);
    domain->lamda2nodex_current(element->nlocal+element->nghost, element->nodex_current);

  }
  neighbor->build_bins();

  // build neighbor list this command needs based on earlier request

  list = neighbor->lists[irequest];
  neighbor->build_one(list, 1);

  compute();


  // write output file

  for (int i = 0; i < noutputs; i++)
    write_output(output_filename_list[i], output_format_list[i], output_configuration_list[i]);


  // reset number of atoms/elements and other settings

  atom->atom_strain_flag = 0;
  element->evec->set_data_size();
  atom->avec->size_border = 7;
  domain->box_exist = 0; 
  neighbor->cutneighmax = 0.0;

  atom->natoms = 0;
  atom->nlocal = 0;
  atom->nghost = 0;
  element->nelements = 0;
  element->nlocal = 0;
  element->nghost = 0;
  element->netypes = 0;
  delete [] output_format_list;
  delete [] output_configuration_list;
  delete [] output_filename_list;
}


/*  ----------------------------------------------------------------------
    proc 0 opens data file
    test if gzipped
    -------------------------------------------------------------------------  */

void AtomicStrain::open_input(char *file)
{
  fp = fopen(file, "r");

  if (fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open input file %s", file);
    error->one(FLERR, str);
  }
}

/*  ----------------------------------------------------------------------
    proc 0 opens output file
    test if gzipped
    -------------------------------------------------------------------------  */

void AtomicStrain::open_output(char *file)
{
  fp = fopen(file, "w");

  if (fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open output file %s", file);
    error->one(FLERR, str);
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

void AtomicStrain::header_ref()
{

  int n;
  char *ptr;

  // customize for new sections
  // When adding new sections, remember to change NSECTIONS
  const char *section_keywords[NSECTIONS] =
  {"Atoms", "Element Types", 
    //"Ellipsoids", "Lines", "Triangles", "Bodies", "Bonds", "Angles", "Dihedrals", "Impropers", 
    //"Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", 
    //"BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs", 
    //"EndBondTorsion Coeffs", "AngleTorsion Coeffs", 
    //"AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs", 
    "Elements", "Nodes", "Velocities", "Node Velocities"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line, MAXLINE, fp);
    if (eof == NULL) error->one(FLERR, "Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line, MAXLINE, fp) == NULL) n = 0;
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

    if (strstr(line, "atoms")) {
      sscanf(line, BIGINT_FORMAT, &natoms);
      if (atom->natoms != natoms) error->all(FLERR, "Number of atoms in ref/file does not match input file");

    } else if (strstr(line, "elements")) {
      sscanf(line, BIGINT_FORMAT, &nelements);
      if (element->nelements != nelements) error->all(FLERR, "Number of elements in ref/file does not match input file");
    }

    // Atom/Element class type and grain settings are only set by first data file

    else if (strstr(line, "atom types")) {
    } else if (strstr(line, "element types")) {
    }

    // local copy of box info
    // so can treat differently for first vs subsequent data files

    else if (strstr(line, "xlo xhi")) {
      sscanf(line, "%lg %lg", &boxlo[0], &boxhi[0]);
    } else if (strstr(line, "ylo yhi")) {
      sscanf(line, "%lg %lg", &boxlo[1], &boxhi[1]); 
    } else if (strstr(line, "zlo zhi")) {
      sscanf(line, "%lg %lg", &boxlo[2], &boxhi[2]);
    } else if (strstr(line, "xy xz yz")) {
      triclinic = 1;
      sscanf(line, "%lg %lg %lg", &xy, &xz, &yz);

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
    read free-format header of data file
    1st line and blank lines are skipped
    non-blank lines are checked for header keywords and leading value is read
    header ends with EOF or non-blank line containing no header keyword
    if EOF, line is set to blank line
    else line has first keyword line for rest of file
    some logic differs if adding atoms
    -------------------------------------------------------------------------  */

void AtomicStrain::header()
{

  int n;
  char *ptr;

  // customize for new sections
  // When adding new sections, remember to change NSECTIONS
  const char *section_keywords[NSECTIONS] =
  {"Atoms", "Element Types", 
    //"Ellipsoids", "Lines", "Triangles", "Bodies", "Bonds", "Angles", "Dihedrals", "Impropers", 
    //"Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", 
    //"BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs", 
    //"EndBondTorsion Coeffs", "AngleTorsion Coeffs", 
    //"AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs", 
    "Elements", "Nodes"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line, MAXLINE, fp);
    if (eof == NULL) error->one(FLERR, "Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line, MAXLINE, fp) == NULL) n = 0;
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

    if (strstr(line, "atoms")) {
      sscanf(line, BIGINT_FORMAT, &natoms);
      atom->natoms = natoms;

    } else if (strstr(line, "elements")) {
      sscanf(line, BIGINT_FORMAT, &nelements);
      element->nelements = nelements;
    }

    // Atom/Element class type and grain settings are only set by first data file

    else if (strstr(line, "atom types")) {
      sscanf(line, "%d", &natypes);
      atom->ntypes = natypes;
    } else if (strstr(line, "element types")) {
      sscanf(line, "%d", &netypes);
      element->evec->grow_etype_arrays(0);
      element->evec->grow_etype_arrays(netypes);
    }

    // local copy of box info
    // so can treat differently for first vs subsequent data files

    else if (strstr(line, "xloref xhiref")) {
      sscanf(line, "%lg %lg", &boxlo[0], &boxhi[0]);
    } else if (strstr(line, "yloref yhiref")) {
      sscanf(line, "%lg %lg", &boxlo[1], &boxhi[1]); 
    } else if (strstr(line, "zloref zhiref")) {
      sscanf(line, "%lg %lg", &boxlo[2], &boxhi[2]);
    } else if (strstr(line, "xyref xzref yzref")) {
      triclinic = 1;
      sscanf(line, "%lg %lg %lg", &xy, &xz, &yz);
    }

    // local copy of current box info

    else if (strstr(line, "xlocur xhicur")) {
      sscanf(line, "%lg %lg", &boxlo_current[0], &boxhi_current[0]);
    } else if (strstr(line, "ylocur yhicur")) {
      sscanf(line, "%lg %lg", &boxlo_current[1], &boxhi_current[1]); 
    } else if (strstr(line, "zlocur zhicur")) {
      sscanf(line, "%lg %lg", &boxlo_current[2], &boxhi_current[2]);
    } else if (strstr(line, "xycur xzcur yzcur")) {
      triclinic = 1;
      sscanf(line, "%lg %lg %lg", &xy_current, &xz_current, &yz_current);

    } else break;
  }

  // error check on total system size

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT ||
      element->nelements < 0 || element->nelements >= MAXBIGINT)
    error->all(FLERR, "System in data file is too big");

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

void AtomicStrain::parse_keyword(int first)
{
  int eof = 0;
  int done = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line, MAXLINE, fp) == NULL) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line, " \t\n\r");
      if ((blank == strlen(line)) || (line[blank] == '#')) {
        if (fgets(line, MAXLINE, fp) == NULL) eof = 1;
      } else done = 1;
    }
    if (fgets(buffer, MAXLINE, fp) == NULL) eof = 1;
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
    compare two style strings if they both exist
    one = comment in data file section, two = currently-defined style
    ignore suffixes listed in suffixes array at top of file
    -------------------------------------------------------------------------  */

int AtomicStrain::style_match(const char *one, const char *two)
{
  int i, delta, len, len1, len2;

  if ((one == NULL) || (two == NULL)) return 1;

  len1 = strlen(one);
  len2 = strlen(two);

  //  for (i = 0; suffixes[i] != NULL; i++) {
  //    len = strlen(suffixes[i]);
  //    if ((delta = len1 - len) > 0)
  //      if (strcmp(one+delta, suffixes[i]) == 0) len1 = delta;
  //    if ((delta = len2 - len) > 0)
  //      if (strcmp(two+delta, suffixes[i]) == 0) len2 = delta;
  //  }

  if ((len1 == 0) || (len1 == len2) || (strncmp(one, two, len1) == 0)) return 1;
  return 0;
}

/*  ----------------------------------------------------------------------
    read all atoms
    -------------------------------------------------------------------------  */

void AtomicStrain::atoms()
{

  int nchunk, eof;

  if (me == 0) {
    if (screen) fprintf(screen, "  reading atoms ...\n");
    if (logfile) fprintf(logfile, "  reading atoms ...\n");
  }

  bigint nread = 0;
  while (nread < natoms) {
    nchunk = MIN(natoms-nread, CHUNK);
    eof = comm->read_lines_from_file(fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    atom->data_atoms_strain(nchunk, buffer);
    nread += nchunk;
  }

  // check that all atoms were assigned correctly

  bigint n = atom->nlocal;
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  bigint nassign = sum - (atom->natoms - natoms);

  if (me == 0) {
    if (screen) fprintf(screen, "  " BIGINT_FORMAT " atoms\n", nassign);
    if (logfile) fprintf(logfile, "  " BIGINT_FORMAT " atoms\n", nassign);
  }

  if (sum != atom->natoms)
    error->all(FLERR, "Did not assign all atoms correctly");
}

/*  ----------------------------------------------------------------------
    read all atoms
    -------------------------------------------------------------------------  */

void AtomicStrain::atoms_ref()
{

  int nchunk, eof;

  if (me == 0) {
    if (screen) fprintf(screen, "  reading atoms ...\n");
    if (logfile) fprintf(logfile, "  reading atoms ...\n");
  }

  bigint nread = 0;
  while (nread < natoms) {
    nchunk = MIN(natoms-nread, CHUNK);
    eof = comm->read_lines_from_file(fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    atom->data_atoms(nchunk, buffer, 0, 0, 0, NULL);
    nread += nchunk;
  }

  // check that all atoms were assigned correctly

  bigint n = atom->nlocal;
  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  bigint nassign = sum - (atom->natoms - natoms);

  if (me == 0) {
    if (screen) fprintf(screen, "  " BIGINT_FORMAT " atoms\n", nassign);
    if (logfile) fprintf(logfile, "  " BIGINT_FORMAT " atoms\n", nassign);
  }

  if (sum != atom->natoms)
    error->all(FLERR, "Did not assign all atoms correctly");
}

/* ----------------------------------------------------------------------- */

void AtomicStrain::element_types()
{
  char *next;
  char *buf = new char[netypes * MAXLINE];

  int eof = comm->read_lines_from_file(fp, netypes, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < netypes; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    element->evec->set_element_types(buf, 0);
    buf = next + 1; 
  }
  delete [] original;
}

/*  ----------------------------------------------------------------------
    proc 0 reads N lines from file
    could be skipping Natoms lines, so use bigints
    -------------------------------------------------------------------------  */

void AtomicStrain::skip_lines(bigint n)
{
  if (me) return;
  char *eof;
  for (bigint i = 0; i < n; i++) eof = fgets(line, MAXLINE, fp);
  if (eof == NULL) error->one(FLERR, "Unexpected end of data file");
}

/*  ----------------------------------------------------------------------
    read all elements
    -------------------------------------------------------------------------  */

void AtomicStrain::elements()
{
  int nchunk, eof;

  if (me == 0) {
    if (screen) fprintf(screen, "  reading elements ...\n");
    if (logfile) fprintf(logfile, "  reading elements ...\n");
  }

  bigint nread = 0;

  while (nread < nelements) {
    nchunk = MIN((nelements-nread), CHUNK);
    eof = comm->read_lines_from_file(fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    element->data_elements(nchunk, buffer, 0, 0, 0, NULL);
    nread += nchunk;
  }

  // count total number of nodes

  // check that all elements were assigned correctly 
  // count total number of nodes

  bigint n = element->nlocal;

  bigint nnode_local = 0;

  for (int i = 0; i < n; i++) {
    int itype = element->etype[i];
    nnode_local += element->npe[itype] * element->apc[itype];
  }

  bigint sum;
  MPI_Allreduce(&n, &sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  MPI_Allreduce(&nnode_local, &nnodes, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  element->nnodes = nnodes;

  bigint nassign = sum - (element->nelements - nelements);

  if (me == 0) {
    if (screen) fprintf(screen, "  " BIGINT_FORMAT " elements\n", nassign);
    if (logfile) fprintf(logfile, "  " BIGINT_FORMAT " elements\n", nassign);
  }

  if (sum != element->nelements)
    error->all(FLERR, "Did not assign all elements correctly");

  n = element->nucells;
  MPI_Allreduce(&n, &sum, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  element->nucells = sum;

  // check that element IDs are valid

  element->tag_check();

  // create global mapping of elements

  if (element->map_style) {
    element->map_init();
    element->map_set();
  }
}

/*  ----------------------------------------------------------------------
    read all nodes
    to find elements, must build element map
    -------------------------------------------------------------------------  */

void AtomicStrain::nodes()
{
  int nchunk, eof;

  if (me == 0) {
    if (screen) fprintf(screen, "  reading nodes ...\n");
    if (logfile) fprintf(logfile, "  reading nodes ...\n");
  }

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  bigint nread = 0;
  while (nread < nnodes) {
    nchunk = MIN(nnodes - nread, CHUNK);
    eof = comm->read_lines_from_file(fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    element->data_nodes_strain(nchunk, buffer);
    nread += nchunk;
  }  
  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  element->evec->update_node_coord();
  if (me == 0) {
    if (screen) fprintf(screen, "  " BIGINT_FORMAT " nodes\n", nnodes);
    if (logfile) fprintf(logfile, "  " BIGINT_FORMAT " nodes\n", nnodes);
  }

}

/*  ----------------------------------------------------------------------
    read all nodes
    to find elements, must build element map
    -------------------------------------------------------------------------  */

void AtomicStrain::nodes_ref()
{
  int nchunk, eof;

  if (me == 0) {
    if (screen) fprintf(screen, "  reading nodes ...\n");
    if (logfile) fprintf(logfile, "  reading nodes ...\n");
  }

  int mapflag = 0;
  if (element->map_style == 0) {
    mapflag = 1;
    element->map_init();
    element->map_set();
  }

  bigint nread = 0;
  while (nread < nnodes) {
    nchunk = MIN(nnodes - nread, CHUNK);
    eof = comm->read_lines_from_file(fp, nchunk, MAXLINE, buffer);
    if (eof) error->all(FLERR, "Unexpected end of data file");
    element->data_nodes(nchunk, buffer, 0, 0, 0, NULL);
    nread += nchunk;
  }  
  if (mapflag) {
    element->map_delete();
    element->map_style = 0;
  }

  element->evec->update_node_coord();
  if (me == 0) {
    if (screen) fprintf(screen, "  " BIGINT_FORMAT " nodes\n", nnodes);
    if (logfile) fprintf(logfile, "  " BIGINT_FORMAT " nodes\n", nnodes);
  }

}

/*  ----------------------------------------------------------------------
    read input data
    -------------------------------------------------------------------------  */

void AtomicStrain::read_input(char *filename)
{

  // flags for this data file

  int atomflag, elementflag, nodeflag;
  int element_type_flag;

  elementflag = atomflag = nodeflag = 0;
  element_type_flag = 0;

  // values in this data file

  natoms = nelements = nnodes = 0;
  netypes = natypes = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  triclinic = 0;


  if (me == 0) {
    if (screen) fprintf(screen, "Reading data file for atomic_strain command '%s'...\n", filename);
    if (logfile) fprintf(logfile, "Reading data file for atomic_strain command '%s'...\n", filename);
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
    if (strcmp(keyword, "Element Types") == 0) {
      element_types();
      element_type_flag = 1;

    } else if (strcmp(keyword, "Atoms") == 0) {
      atomflag = 1;
      if (me == 0 && !style_match(style, atom->atom_style))
        error->warning(FLERR, "Atom style in data file differs "
            "from currently defined atom style");
      atoms();
    } else if (strcmp(keyword, "Elements") == 0) {
      if (element_type_flag == 0) 
        error->all(FLERR, "Element types must be defined before Elements section");
      elementflag = 1;
      elements();
    } else if (strcmp(keyword, "Nodes") == 0) {
      if (elementflag == 0)
        error->all(FLERR, "Elements section should be read before Nodes section");
      nodeflag = 1;
      nodes();
    } else {
      char str[128];
      sprintf(str, "Unknown identifier in data file: %s", keyword);
      error->all(FLERR, str);
    }
    parse_keyword(0);
  }

  // error if natoms/nelements > 0 yet no atoms/elements were read

  if (natoms > 0 && atomflag == 0)
    error->all(FLERR, "No atoms in data file");
  if (nelements > 0 && (elementflag == 0 || nodeflag == 0))
    error->all(FLERR, "No elements/nodes in data file");

  // close input file

  if (me == 0) {
    fclose(fp);
    fp = NULL;
  }
}

/*  ----------------------------------------------------------------------
    read input data
    -------------------------------------------------------------------------  */

void AtomicStrain::read_ref_file(char *filename)
{

  // flags for this data file

  int atomflag, elementflag, nodeflag;
  int element_type_flag;

  elementflag = atomflag = nodeflag = 0;
  element_type_flag = 0;

  // values in this data file

  natoms = nelements = nnodes = 0;
  netypes = natypes = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  triclinic = 0;


  if (me == 0) {
    if (screen) fprintf(screen, "Reading reference data file for atomic_strain command '%s'...\n", filename);
    if (logfile) fprintf(logfile, "Reading reference data file for atomic_strain command '%s'...\n", filename);
    open_input(filename);
  } else fp = NULL;

  // read header info

  header_ref();

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
  domain->set_current_global_box();
  comm->set_proc_grid();
  domain->set_local_box(); 

  // customize for new sections
  // read rest of file in free format

  while (strlen(keyword)) {
    if (strcmp(keyword, "Element Types") == 0) {
      skip_lines(element->netypes);
      element_type_flag = 1;

    } else if (strcmp(keyword, "Atoms") == 0) {
      atomflag = 1;
      if (me == 0 && !style_match(style, atom->atom_style))
        error->warning(FLERR, "Atom style in data file differs "
            "from currently defined atom style");
      atom->nlocal = 0;
      atoms_ref();
    } else if (strcmp(keyword, "Elements") == 0) {
      if (element_type_flag == 0) 
        error->all(FLERR, "Element types must be defined before Elements section");
      elementflag = 1;
      element->nlocal = 0;
      elements();
    } else if (strcmp(keyword, "Nodes") == 0) {
      if (elementflag == 0)
        error->all(FLERR, "Elements section should be read before Nodes section");
      nodeflag = 1;
      nodes_ref();
    } else if (strcmp(keyword, "Velocities") == 0) {
      skip_lines(atom->natoms);
    } else if (strcmp(keyword, "Node Velocities") == 0) {
      element->count_nodes(1);
      skip_lines(element->nnodes);
    } else {
      char str[128];
      sprintf(str, "Unknown identifier in data file: %s", keyword);
      error->all(FLERR, str);
    }
    parse_keyword(0);
  }

  // error if natoms/nelements > 0 yet no atoms/elements were read

  if (natoms > 0 && atomflag == 0)
    error->all(FLERR, "No atoms in data file");
  if (nelements > 0 && (elementflag == 0 || nodeflag == 0))
    error->all(FLERR, "No elements/nodes in data file");

  // close input file

  if (me == 0) {
    fclose(fp);
    fp = NULL;
  }
}


/*  ----------------------------------------------------------------------
    perform atomic strain calculation
    -------------------------------------------------------------------------  */

void AtomicStrain::compute()
{
  int nalocal = atom->nlocal;
  int naghost = atom->nghost;
  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;

  // grow arrays if needed based on flags

  memory->grow(atomValidity, nalocal, "atomic_strain:atomValidity");
  memory->grow(nodeValidity, nelocal, maxapc, maxnpe, "atomic_strain:nodeValidity");
  memory->grow(atomVonMisesStrain, nalocal, "atomic_strain:atomVonMisesStrain");
  memory->grow(nodeVonMisesStrain, nelocal, maxapc, maxnpe, "atomic_strain:nodeVonMisesStrain");
  memory->grow(atomVolStrain, nalocal, "atomic_strain:atomVolStrain");
  memory->grow(nodeVolStrain, nelocal, maxapc, maxnpe, "atomic_strain:nodeVolStrain");
  if (deformation_gradient_flag) {
    memory->grow(atomF, nalocal, 9, "atomic_strain:atomF");  
    memory->grow(nodeF, nelocal, maxapc, maxnpe, 9, "atomic_strain:nodeF");  
  }
  if (rotation_quat_flag) {
    memory->grow(atomQ, nalocal, 4, "atomic_strain:atomQ");  
    memory->grow(nodeQ, nelocal, maxapc, maxnpe, 4, "atomic_strain:nodeQ");  
  }
  if (rotation_vect_flag) {
    memory->grow(atomW, nalocal, 3, "atomic_strain:atomW");  
    memory->grow(nodeW, nelocal, maxapc, maxnpe, 3, "atomic_strain:nodeW");  
  }
  if (stretch_tensor_flag) {
    memory->grow(atomStretchTensor, nalocal, 6, "atomic_strain:atomStretchTensor");  
    memory->grow(nodeStretchTensor, nelocal, maxapc, maxnpe, 6, "atomic_strain:nodeStretchTensor");  
  }
  if (strain_tensor_flag) {
    memory->grow(atomStrainTensor, nalocal, 6, "atomic_strain:atomStrainTensor");  
    memory->grow(nodeStrainTensor, nelocal, maxapc, maxnpe, 6, "atomic_strain:nodeStrainTensor");  
  }
  if (nsd_flag) {
    memory->grow(atomNSD, nalocal, "atomic_strain:atomNSD");  
    memory->grow(nodeNSD, nelocal, maxapc, maxnpe, "atomic_strain:nodeNSD");  
  }

  memory->grow(atomD, atom->nmax, 3, "atomic_strain:atomD");
  memory->grow(nodeD, element->nmax, maxapc, maxnpe, 3, "atomic_strain:nodeD");

  double **x_ref = atom->x;
  double **x_cur = atom->x_current;
  double ****nodex_ref = element->nodex;
  double ****nodex_cur = element->nodex_current;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  imageint *aimage = atom->image;
  imageint *eimage = element->image;
  double coord[3];

  // compute displacements vector for all atoms/nodes

  for (int i = 0; i < nalocal + naghost; i++) 
    sub3(x_cur[i], x_ref[i], atomD[i]);     


  for (int i = 0; i < nelocal + neghost; i++) 
    for (int j = 0; j < apc[etype[i]]; j++) 
      for (int k = 0; k < npe[etype[i]]; k++) 
        sub3(nodex_cur[i][j][k], nodex_ref[i][j][k], nodeD[i][j][k]);

  int inum = list->inum; 
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  double *defGradient, *Q, *W, *stretchTensor, *strainTensor;
  double *NSD, *vonMisesStrain, *volStrain;
  int *validity;
  double *ref_coord, *disp;
  int i, iindex, ibasis, inode, iapc;

  for (int ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];

    // i is atom

    if (iindex < 0) {
      if (deformation_gradient_flag) defGradient = atomF[i];
      else defGradient = NULL;
      if (rotation_quat_flag) Q = atomQ[i];
      else Q = NULL;
      if (rotation_vect_flag) W = atomW[i];
      else W = NULL;
      if (stretch_tensor_flag) stretchTensor = atomStretchTensor[i];
      else stretchTensor = NULL;
      if (strain_tensor_flag) strainTensor = atomStrainTensor[i];
      else strainTensor = NULL;
      if (nsd_flag) NSD = &atomNSD[i];
      vonMisesStrain = &atomVonMisesStrain[i];
      volStrain = &atomVolStrain[i];
      validity = &atomValidity[i];
      ref_coord = x_ref[i];
      disp = atomD[i];
    } 

    // i is node

    else {
      iapc = apc[etype[i]];
      inode = iindex / iapc;
      ibasis = iindex % iapc;
      if (deformation_gradient_flag) defGradient = nodeF[i][ibasis][inode];
      else defGradient = NULL;
      if (rotation_quat_flag) Q = nodeQ[i][ibasis][inode];
      else Q = NULL;
      if (rotation_vect_flag) W = nodeW[i][ibasis][inode];
      else W = NULL;
      if (stretch_tensor_flag) stretchTensor = nodeStretchTensor[i][ibasis][inode];
      else stretchTensor = NULL;
      if (strain_tensor_flag) strainTensor = nodeStrainTensor[i][ibasis][inode];
      else strainTensor = NULL;
      if (nsd_flag) NSD = &nodeNSD[i][ibasis][inode];
      vonMisesStrain = &nodeVonMisesStrain[i][ibasis][inode];
      volStrain = &nodeVolStrain[i][ibasis][inode];
      validity = &nodeValidity[i][ibasis][inode];
      ref_coord = nodex_ref[i][ibasis][inode];
      disp = nodeD[i][ibasis][inode];
    }

    *validity = compute_strain(numneigh[ii], firstneigh[ii], firstneighindex[ii], ref_coord, disp, 
        defGradient, Q, W, stretchTensor, strainTensor, NSD, vonMisesStrain, volStrain);
  }
}

/*  ----------------------------------------------------------------------
    perform atomic strain calculation for an atom/node
    -------------------------------------------------------------------------  */

int AtomicStrain::compute_strain(int numneigh, int *jlist, int *jindexlist, 
    double *coord_ref, double *idisp, double *defGradient, double *rot_quat, 
    double *rot_vect, double *stretchTensor, double *strainTensor, double *NSD, 
    double *vonMisesStrain, double *volStrain)
{
  double del_ref[3], del_cur[3], jdisp[3], jx_ref[3];
  double V[3][3], W[3][3];
  int j, k, l, node, jj, col, row;
  int jindex, jetype, jucell, jbasis, japc;

  double **x_ref = atom->x;
  double **x_cur = atom->x_current;
  double ****nodex_ref = element->nodex;
  double ****nodex_cur = element->nodex_current;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  ElementVec *evec = element->evec;

  V[0][0] = W[0][0] = V[0][1] = W[0][1] = V[0][2] = W[0][2] = 0.0; 
  V[1][0] = W[1][0] = V[1][1] = W[1][1] = V[1][2] = W[1][2] = 0.0; 
  V[2][0] = W[2][0] = V[2][1] = W[2][1] = V[2][2] = W[2][2] = 0.0; 

  double sumSqDist = 0;

  // loop through neighbors to tally V and W
  // calculate current seperation vector from displacement vectors 
  // instead of current positions (due to wrapping)
  // del_cur = del_ref + jdisp - idisp


  for (jj = 0; jj < numneigh; jj++) {

    j = jlist[jj];
    jindex = jindexlist[jj];

    // j is atom

    if (jindex < 0) {
      jdisp[0] = atomD[j][0]; 
      jdisp[1] = atomD[j][1]; 
      jdisp[2] = atomD[j][2]; 
      jx_ref[0] = x_ref[j][0];
      jx_ref[1] = x_ref[j][1];
      jx_ref[2] = x_ref[j][2];
    }

    // j is virtual atom

    else {
      japc = apc[etype[j]];
      jucell = jindex / japc;
      jbasis = jindex % japc;
      evec->interpolate(jdisp, nodeD, j, jbasis, jucell, 3);
      evec->interpolate(jx_ref, nodex_ref, j, jbasis, jucell, 3);
    }

    for (k = 0; k < 3; k++) {
      del_ref[k] = jx_ref[k] - coord_ref[k];
      del_cur[k] = del_ref[k] + jdisp[k] - idisp[k];
    }
    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++) {
        V[k][l] += del_ref[l] * del_ref[k];
        W[k][l] += del_ref[l] * del_cur[k];
      }
    sumSqDist += del_ref[0] * del_ref[0] 
      + del_ref[1] * del_ref[1] 
      + del_ref[2] * del_ref[2];
  }

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
    *vonMisesStrain = 0;
    *volStrain = 0;
    if (deformation_gradient_flag) {
      defGradient[0] = defGradient[1] = defGradient[2] = 0;
      defGradient[3] = defGradient[4] = defGradient[5] = 0;
      defGradient[6] = defGradient[7] = defGradient[8] = 0;
    }
    if (rotation_quat_flag) 
      rot_quat[0] = rot_quat[1] = rot_quat[2] = rot_quat[3] = 0;
    if (rotation_vect_flag) 
      rot_vect[0] = rot_vect[1] = rot_vect[2] = 0;
    if (nsd_flag) *NSD = 0;
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

  double invV[3][3], F[3][3];
  invert3(V, invV); 
  times3(W, invV, F);

  if (deformation_gradient_flag) 
    for (row = 0; row < 3; row++)
      for (col = 0; col < 3; col++)
        defGradient[row * 3 + col] = F[row][col];

  if (rotation_vect_flag) {
    rot_vect[0] = 0.5 * (F[2][1] - F[1][2]);
    rot_vect[1] = 0.5 * (F[0][2] - F[2][0]);
    rot_vect[2] = 0.5 * (F[1][0] - F[0][1]);
  }

  // Polar decomposition F = RU

  if (rotation_quat_flag || stretch_tensor_flag) {
    double R[9], U[9], Farr[9];

    // convert F to array (row by row)

    for (row = 0; row < 3; row++)
      for (col = 0; col < 3; col++)
        Farr[row * 3 + col] = F[row][col];

    polar_decomposition_3x3(Farr, false, R, U);
    if (rotation_quat_flag) {

      // If F contains a reflection, R will not be a pure rotation_quat matrix and the
      // conversion to a quaternion below will fail.
      // Thus, in the rather unlikely case that F contains a reflection, we simply flip the
      // R matrix to make it a pure rotation_quat.

      if (matrix_determinant_3x3(R) < 0) {
        for (k = 0; k < 9; k++)
          R[k] = -R[k];
      }
      rotation_matrix_to_quaternion(R, rot_quat);
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
  transpose_times3(F, F, E);
  E[0][0] -= 1;
  E[1][1] -= 1;
  E[2][2] -= 1;
  scalar_times3(0.5, E);
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
    *NSD = 0;
    double nonaffDisp[3];

    for (jj = 0; jj < numneigh; jj++) {

      j = jlist[jj];
      jindex = jindexlist[jj];

      // j is atom

      if (jindex < 0) {
        jdisp[0] = atomD[j][0]; 
        jdisp[1] = atomD[j][1]; 
        jdisp[2] = atomD[j][2]; 
        jx_ref[0] = x_ref[j][0];
        jx_ref[1] = x_ref[j][1];
        jx_ref[2] = x_ref[j][2];
      }

      // j is virtual atom

      else {
        japc = apc[etype[j]];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        evec->interpolate(jdisp, nodeD, j, japc, jucell, 3);
        evec->interpolate(jx_ref, nodex_ref, j, japc, jucell, 3);
      }

      for (k = 0; k < 3; k++) {
        del_ref[k] = jx_ref[k] - coord_ref[k];
        del_cur[k] = del_ref[k] + jdisp[k] - idisp[k];
      }
      matvec(F, del_ref, nonaffDisp);
      sub3(nonaffDisp, del_cur, nonaffDisp);
      *NSD += lensq3(nonaffDisp);
    }
  }

  // Calculate von Mises shear strain

  double xydiff = Sxx - Syy;
  if (domain->dimension == 3) {
    double xzdiff = Sxx - Szz;
    double yzdiff = Syy - Szz;
    *vonMisesStrain = sqrt(Sxy * Sxy + Sxz * Sxz + Syz * Syz +
        (xydiff * xydiff + xzdiff * xzdiff + yzdiff * yzdiff) / 6.0);
  } else *vonMisesStrain = sqrt(Sxy * Sxy + xydiff * xydiff / 2.0);

  if (!ISFINITE(*vonMisesStrain))
    *vonMisesStrain = 0;

  // Calculate volumetric component

  if (domain->dimension == 3)
    *volStrain = (Sxx + Syy + Szz)/3.0;
  else
    *volStrain = (Sxx + Syy)/2.0;
  if (!ISFINITE(*volStrain))
    *volStrain = 0;
  return 1;
}

/*  ----------------------------------------------------------------------
    write output file
    -------------------------------------------------------------------------  */

void AtomicStrain::write_output(char *filename, int output_format, int output_configuration)
{
  if (me == 0) {
    if (screen) fprintf(screen, "Writing output file for atomic_strain command '%s'...\n", filename);
    if (logfile) fprintf(logfile, "Writing output file for atomic_strain command '%s'...\n", filename);
  }

  if (region_id) {
    char *group_id = (char *) "ATOMIC_STRAIN_OUTPUT";
    char **newarg = new char*[5];
    newarg[0] = group_id;
    newarg[1] = (char *) "region";
    newarg[2] = region_id;
    newarg[3] = (char *) "style";
    newarg[4] = (char *) "oneatom";

    // swap current and reference x pointers to assign group
    // swap back after assigning group

    if (output_configuration == CURRENT) {
      double **tmp2 = atom->x;
      atom->x = atom->x_current;
      atom->x_current = tmp2;
      double ****tmp4 = element->nodex;
      element->nodex = element->nodex_current;
      element->nodex_current = tmp4;
    }
    group->assign(5, newarg);
    if (output_configuration == CURRENT) {
      double **tmp2 = atom->x;
      atom->x = atom->x_current;
      atom->x_current = tmp2;
      double ****tmp4 = element->nodex;
      element->nodex = element->nodex_current;
      element->nodex_current = tmp4;
    }

    delete [] newarg;

    igroup = group->find(group_id);
    groupbit = group->bitmask[igroup];
  } else {
    igroup = group->find((char *) "all");
    groupbit = group->bitmask[igroup];
  }

  if (group->count_atom(igroup) == 0) no_atom_flag == 1;
  if (group->count_elem(igroup) == 0) no_elem_flag == 1;

  if (output_format < ATOM) {
    if (wrap_flag) {
      for (int i = 0; i < atom->nlocal; i++) {
        if (output_configuration == REFERENCE) domain->remap(atom->x[i]);
        else domain->remap_current(atom->x_current[i]);
      }
      int *npe = element->npe;
      int *apc = element->apc;
      int *etype = element->etype;
      double coord[3];
      for (int i = 0; i < element->nlocal; i++) {
        if (output_configuration == REFERENCE) 
          domain->remap(element->x[i], element->nodex[i], apc[etype[i]], npe[etype[i]]);
        else {
          coord[0] = coord[1] = coord[2] = 0;
          for (int j = 0; j < apc[etype[i]]; j++) 
            for (int k = 0; k < npe[etype[i]]; k++) {
              coord[0] += element->nodex_current[i][j][k][0];
              coord[1] += element->nodex_current[i][j][k][1];
              coord[2] += element->nodex_current[i][j][k][2];
            }
          coord[0] /= npe[etype[i]] * apc[etype[i]];
          coord[1] /= npe[etype[i]] * apc[etype[i]];
          coord[2] /= npe[etype[i]] * apc[etype[i]];
          domain->remap_current(coord, element->nodex_current[i], apc[etype[i]], npe[etype[i]]);
        }
      }
    }
    size_one = 6 + 4 * (element_info_flag + rotation_quat_flag) + 9 * deformation_gradient_flag + nsd_flag 
      + 6 * (strain_tensor_flag + stretch_tensor_flag) + 3 * (rotation_vect_flag + displacement_flag);
    write_header_tecplot(filename, output_format);

    if (no_atom_flag == 0)
      write_tecplot_atoms(output_format, output_configuration);

    if (no_elem_flag == 0) {
      write_tecplot_nodes(output_format, output_configuration);
      write_tecplot_element_connectivity(output_format);
    }
  } else {
    size_one = 8 + 6 * element_info_flag + 9 * deformation_gradient_flag + 4 * rotation_quat_flag + nsd_flag 
      + 6 * (strain_tensor_flag + stretch_tensor_flag) + 3 * (rotation_vect_flag + displacement_flag);
    int nalocal = atom->nlocal;
    int nelocal = element->nlocal;
    int *apc = element->apc;
    int *npe = element->npe;
    int *etype = element->etype;
    int *nucell = element->nucell;
    int *amask = atom->mask;
    int *emask = element->mask;
    bigint ntotal_local = 0;
    bigint ntotal;

    if (no_atom_flag == 0)
      for (int i = 0; i < nalocal; i++)
        if (amask[i] & groupbit) ntotal_local++;
    if (no_elem_flag == 0)
      for (int i = 0; i < nelocal; i++)
        if (emask[i] & groupbit)
          if (output_format == FULLATOM) {
            ntotal_local += nucell[etype[i]] * apc[etype[i]];
          } else if (output_format == NODE) {
            ntotal_local += npe[etype[i]] * apc[etype[i]];
          }

    MPI_Allreduce(&ntotal_local, &ntotal, 1, MPI_CAC_BIGINT, MPI_SUM, world);

    write_header_atom(filename, output_configuration, ntotal);

    int tmp, nlines;
    grow_buf(ntotal_local * size_one);
    pack_atom_all(output_format, output_configuration);
    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
          MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
          MPI_Wait(&request, &status);
          MPI_Get_count(&status, MPI_DOUBLE, &nlines);
          nlines /= size_one;
        } else nlines = ntotal_local;
        write_lines_atom(nlines, buf);
      }
    } else {
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(buf, ntotal_local * size_one, MPI_DOUBLE, 0, 0, world);
    }
  }

  // close output file

  if (me == 0) {
    if (output_format == ASCII || output_format >= ATOM) {
      fclose(fp);
      fp = NULL;
    } else {
      int success = TECEND142();
      if (success == -1)
        error->all(FLERR, "Cannot close output file");
    }
  }
}

/*  ----------------------------------------------------------------------
    write header for atom output (LAMMPS dump style)
    -------------------------------------------------------------------------  */

void AtomicStrain::write_header_atom(char *file, int output_configuration, bigint ndump) 
{
  // open data file and write headers

  if (me == 0) {

    char *columns = new char[500];
    char boundstr[9];
    double boxxlo, boxxhi;      // local copies of domain values
    double boxylo, boxyhi;      // lo/hi are bounding box for triclinic
    double boxzlo, boxzhi;
    double boxxy, boxxz, boxyz;

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
      sprintf(str, "Cannot open data file %s", file);
      error->one(FLERR, str);
    }
    columns = strcpy(columns, "ID type x y z");
    if (element_info_flag)
      strcat(columns, " etype eid inode iucell ibasis numneigh");
    strcat(columns, " validity vonMises volStrain");
    if (deformation_gradient_flag)
      strcat(columns, " F11 F12 F13 F21 F22 F23 F31 F32 F33");
    if (rotation_quat_flag)
      strcat(columns, " Rw Rx Ry Rz");
    if (rotation_vect_flag)
      strcat(columns, " w1 w2 w3");
    if (displacement_flag)
      strcat(columns, " dx dy dz");
    if (nsd_flag)
      strcat(columns, " nsd");
    if (strain_tensor_flag)
      strcat(columns, " StrainXX StrainYY StrainZZ StrainXY StrainXZ StrainYZ");
    if (stretch_tensor_flag)
      strcat(columns, " StretchXX StretchYY StretchZZ StretchXY StretchXZ StretchYZ");

    domain->boundary_string(boundstr);
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, BIGINT_FORMAT "\n", ndump);
    if (domain->triclinic) {
      fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
      fprintf(fp, "%-1.16e %-1.16e\n", boxxlo, boxxhi);
      fprintf(fp, "%-1.16e %-1.16e\n", boxylo, boxyhi);
      fprintf(fp, "%-1.16e %-1.16e\n", boxzlo, boxzhi);
    } else {
      fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxxlo, boxxhi, boxxy);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxylo, boxyhi, boxxz);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", boxzlo, boxzhi, boxyz);
    }
    fprintf(fp, "ITEM: ATOMS %s\n", columns);
    delete [] columns;
  }
}


/*  ----------------------------------------------------------------------
    write header for tecplot output
    -------------------------------------------------------------------------  */

void AtomicStrain::write_header_tecplot(char *file, int output_format) 
{
  // open data file and write headers

  if (me == 0) {

    char *columns = new char[500];

    if (output_format == ASCII) {
      fp = fopen(file, "w");
      if (fp == NULL) {
        char str[128];
        sprintf(str, "Cannot open data file %s", file);
        error->one(FLERR, str);
      }
      columns = strcpy(columns, "\"x\", \"y\", \"z\"");
      if (element_info_flag)
        strcat(columns, ", \"ID\", \"etype\", \"ctype\", \"numneigh\"");
      strcat(columns, ", \"validity\", \"vonMises\", \"volStrain\"");
      if (deformation_gradient_flag)
        strcat(columns, ", \"F11\", \"F12\", \"F12\", \"F21\", \"F22\", \"F23\", \"F31\", \"F32\", \"F33\"");
      if (rotation_quat_flag)
        strcat(columns, ", \"Rw\", \"Rx\", \"Ry\", \"Rz\"");
      if (rotation_vect_flag)
        strcat(columns, ", \"w1\", \"w2\", \"w3\"");
      if (displacement_flag)
        strcat(columns, ", \"dx\", \"dy\", \"dz\"");
      if (nsd_flag)
        strcat(columns, ", \"nsd\"");
      if (strain_tensor_flag)
        strcat(columns, ", \"StrainXX\", \"StrainYY\", \"StrainZZ\", \"StrainXY\", \"StrainXZ\", \"StrainYZ\"");
      if (stretch_tensor_flag)
        strcat(columns, ", \"StretchXX\", \"StretchYY\", \"StretchZZ\", \"StretchXY\", \"StretchXZ\", \"StretchYZ\"");

      fprintf(fp, "title=\"Strain Calculation from CAC with cutoff = %g\"\n", cut);
      fprintf(fp, "variables = %s\n", columns);
    } else {
      soltime = ntimestep;
      dummy = 0;
      dummy1 = 1;
      columns = strcpy(columns, "x y z");
      if (element_info_flag)
        strcat(columns, " ID etype ctype");
      strcat(columns, " validity vonMises volStrain");
      if (deformation_gradient_flag)
        strcat(columns, " F11 F12 F13 F21 F22 F23 F31 F32 F33");
      if (rotation_quat_flag)
        strcat(columns, " Rw Rx Ry Rz");
      if (rotation_vect_flag)
        strcat(columns, " w1 w2 w3");
      if (displacement_flag)
        strcat(columns, " dx dy dz");
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
        sprintf(str, "Cannot open data file %s", file);
        error->one(FLERR, str);
      }
    }
    delete [] columns;
  }
}

/*  ----------------------------------------------------------------------
    write atom section for tecplot
    -------------------------------------------------------------------------  */

void AtomicStrain::write_tecplot_atoms(int output_format, int output_configuration) 
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;
  int ntotal; 
  MPI_Allreduce(&n, &ntotal, 1, MPI_INT, MPI_SUM, world);

  if (output_format == ASCII) {
    if (me == 0) 
      fprintf(fp, "zone t = \"Timestep: %d, Discrete Atoms\", datapacking = point\n", ntimestep);
    int tmp, nlines;
    grow_buf(n * size_one);
    pack_atom_tecplot(output_format, output_configuration);
    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
          MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
          MPI_Wait(&request, &status);
          MPI_Get_count(&status, MPI_DOUBLE, &nlines);
          nlines /= size_one;
        } else nlines = n;
        write_lines_tecplot(nlines, buf);
      }
    } else {
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(buf, n * size_one, MPI_DOUBLE, 0, 0, world);
    }

  } else {
    if (me == 0) {
      zonetype = 0;       // Ordered zone
      int imax = ntotal;
      int jmax = 1;
      int kmax = 1;
      char title[100];
      sprintf(title, "Timestep: %d, Discrete Atoms", ntimestep);
      int success = TECZNE142(title, &zonetype, &imax, &jmax, &kmax, 
          &dummy, &dummy, &dummy, 
          &soltime, 
          &dummy, &dummy, 
          &dummy1, 
          &dummy, &dummy, 
          &dummy1, &dummy1, &dummy1, 
          NULL, NULL, NULL, 
          &dummy);
      if (success == -1) error->one(FLERR, "Cannot write Discrete Atoms Zone");
    }

    grow_buf(n);

    // communicate buf and dump out each value one by one

    valueflag = COORDINATE;
    for (ival = 0; ival < 3; ival++) {
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);
    }

    if (element_info_flag) {
      valueflag = ID;
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);

      valueflag = ETYPE;
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);

      valueflag = CTYPE;
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);

      valueflag = NUMNEIGH;
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);
    }

    valueflag = VALIDITY;
    pack_atom_tecplot(output_format, output_configuration);
    comm_buf_tecio(n);

    valueflag = VONMISES;
    pack_atom_tecplot(output_format, output_configuration);
    comm_buf_tecio(n);

    valueflag = VOLSTRAIN;
    pack_atom_tecplot(output_format, output_configuration);
    comm_buf_tecio(n);

    if (deformation_gradient_flag) {
      valueflag = DEFGRAD;
      for (ival = 0; ival < 9; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }

    if (rotation_quat_flag) {
      valueflag = ROTATIONQUAT;
      for (ival = 0; ival < 4; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }

    if (rotation_vect_flag) {
      valueflag = ROTATIONVECT;
      for (ival = 0; ival < 3; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }

    if (displacement_flag) {
      valueflag = DISPLACEMENT;
      for (ival = 0; ival < 3; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }

    if (nsd_flag) {
      valueflag = NSD;
      pack_atom_tecplot(output_format, output_configuration);
      comm_buf_tecio(n);
    }

    if (strain_tensor_flag) {
      valueflag = STRAINTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }

    if (stretch_tensor_flag) {
      valueflag = STRETCHTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_atom_tecplot(output_format, output_configuration);
        comm_buf_tecio(n);
      }
    }
  }
}


/*  ----------------------------------------------------------------------
    write element section for tecplot
    -------------------------------------------------------------------------  */

void AtomicStrain::write_tecplot_nodes(int output_format, int output_configuration) 
{
  int nlocal = element->nlocal;
  int nnode_local = 0;
  int nelem_local = 0;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  int *mask = element->mask;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) 
      if (average_flag) {
        nnode_local += npe[etype[i]];
        nelem_local++;
      } else {
        nnode_local += npe[etype[i]] * apc[etype[i]];
        nelem_local += apc[etype[i]];
      }
  }

  int nnodes;
  MPI_Allreduce(&nnode_local, &nnodes, 1, MPI_INT, MPI_SUM, world);
  int nelements;
  MPI_Allreduce(&nelem_local, &nelements, 1, MPI_INT, MPI_SUM, world);

  if (output_format == ASCII) {
    int tmp, nlines;
    if (me == 0) {
      fprintf(fp, "zone t = \"Timestep: %d, Coarse Element\"", ntimestep);
      fprintf(fp, " n = %d e = %d", nnodes, nelements);
      if (domain->dimension == 3) 
        fprintf(fp, " datapacking = point, zonetype = febrick\n");
      else 
        fprintf(fp, " datapacking = point, zonetype = fequadrilateral\n");
    }
    grow_buf(nnode_local * size_one);
    pack_node_tecplot(output_format, output_configuration);

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
          MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
          MPI_Wait(&request, &status);
          MPI_Get_count(&status, MPI_DOUBLE, &nlines);
          nlines /= size_one;
        } else nlines = nnode_local;
        write_lines_tecplot(nlines, buf);
      }
    } else {
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(buf, nnode_local * size_one, MPI_DOUBLE, 0, 0, world);
    }
  } else {
    if (me == 0) {
      if (domain->dimension == 3) zonetype = 5;    // Brick zone
      else zonetype = 3;                           // Quadrilateral zone
      char title[100];
      sprintf(title, "Timestep: %d, Coarse Element", ntimestep);
      int success = TECZNE142(title, &zonetype, &nnodes, &nelements, 
          &dummy, &dummy, &dummy, &dummy, 
          &soltime, 
          &dummy, &dummy, 
          &dummy1, 
          &dummy, &dummy, 
          0, 0, 0, 
          NULL, valuelocation, NULL, 
          &dummy);
      if (success == -1) error->one(FLERR, "Cannot write Coarse Element zone");
    }

    grow_buf(nnode_local);

    // communicate buf and dump out each value one by one

    valueflag = COORDINATE;
    for (ival = 0; ival < 3; ival++) {
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);
    }

    if (element_info_flag) {
      valueflag = ID;
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);

      valueflag = ETYPE;
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);

      valueflag = CTYPE;
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);

      valueflag = NUMNEIGH;
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);
    }

    valueflag = VALIDITY;
    pack_node_tecplot(output_format, output_configuration);
    comm_buf_tecio(nnode_local);

    valueflag = VONMISES;
    pack_node_tecplot(output_format, output_configuration);
    comm_buf_tecio(nnode_local);

    valueflag = VOLSTRAIN;
    pack_node_tecplot(output_format, output_configuration);
    comm_buf_tecio(nnode_local);

    if (deformation_gradient_flag) {
      valueflag = DEFGRAD;
      for (ival = 0; ival < 9; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }

    if (rotation_quat_flag) {
      valueflag = ROTATIONQUAT;
      for (ival = 0; ival < 4; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }

    if (rotation_vect_flag) {
      valueflag = ROTATIONVECT;
      for (ival = 0; ival < 3; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }
    if (displacement_flag) {
      valueflag = DISPLACEMENT;
      for (ival = 0; ival < 3; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }

    if (nsd_flag) {
      valueflag = NSD;
      pack_node_tecplot(output_format, output_configuration);
      comm_buf_tecio(nnode_local);
    }

    if (strain_tensor_flag) {
      valueflag = STRAINTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }

    if (stretch_tensor_flag) {
      valueflag = STRETCHTENSOR;
      for (ival = 0; ival < 6; ival++) {
        pack_node_tecplot(output_format, output_configuration);
        comm_buf_tecio(nnode_local);
      }
    }
  }
}
/*  ----------------------------------------------------------------------
    write element section for tecplot
    -------------------------------------------------------------------------  */

void AtomicStrain::write_tecplot_element_connectivity(int output_format)
{
  // dump node connectivity
  // use int buf for node connectivity 
  // due to tecio function requirement
  // nlines for TECNODE142 fuctions = # of values to write

  int nlocal = element->nlocal;
  int n = 0;
  int *mask = element->mask;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;

  int num_me, nlines, tmp;
  int size_one_line;
  char *format;
  if (domain->dimension == 3) {
    size_one_line = 8;
  } else {
    size_one_line = 4;
  }
  num_me = n * size_one_line;
  if (!average_flag) num_me *= element->maxapc;
  grow_ibuf(num_me);

  pack_element_connectivity();
  if (me == 0) {
    int ntotal = 0;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(ibuf, maxibuf, MPI_INT, me+iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_INT, &nlines);
      } else nlines = num_me;

      if (output_format == ASCII) {
        nlines /= size_one_line;
        int m = 0;
        for (int i = 0; i < nlines; i++) {
          for (int j = 0; j < size_one_line; j++)
            fprintf(fp, "%d ", ibuf[m++]);
          fprintf(fp, "\n");
        }
      } else {
        ntotal += nlines;
        int success = TECNODE142(&nlines, ibuf);
        if (success == -1) error->one(FLERR, "Cannot write to dump file");
      }
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(ibuf, num_me, MPI_INT, 0, 0, world);
  }

}

/* ------------------------------------------------------------------------------------------ */

void AtomicStrain::pack_element_connectivity()
{
  int nlocal = element->nlocal;
  int *mask = element->mask;

  int *apc = element->apc;
  int *etype = element->etype;
  int *element_shape_ids = element->element_shape_ids;
  int m = 0;

  int npe_connect;
  if (domain->dimension == 3) npe_connect = 8;
  else npe_connect = 4;

  if (average_flag) {
    int **nodecell_ids;
    memory->create(nodecell_ids, nlocal, npe_connect, "atomic_strain:nodecell_ids");
    element->set_node_connectivities(igroup, nodecell_ids);
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
        int itype = etype[i];
        m += element->node_connectivity(&ibuf[m], element_shape_ids[itype], nodecell_ids[i]);
      }
    memory->destroy(nodecell_ids);
  } else {
    int ***node_ids;
    memory->create(node_ids, nlocal, element->maxapc, npe_connect, "atomic_strain:node_ids");
    element->set_node_connectivities(igroup, node_ids);
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
        int itype = etype[i];
        for (int j = 0; j < apc[itype]; j++) {
          m += element->node_connectivity(&ibuf[m], element_shape_ids[itype], node_ids[i][j]);
        }
      }
    memory->destroy(node_ids);
  }

} 


/* ---------------------------------------------------------------------------------------------------- */

void AtomicStrain::grow_buf(int num_me)
{
  int nmax;
  MPI_Allreduce(&num_me, &nmax, 1, MPI_INT, MPI_MAX, world);

  if (nmax > maxbuf) {
    if ((bigint) nmax > MAXSMALLINT) {
      char errstr[100];  
      sprintf(errstr, "Too much per-proc info for atomic_strain command nmax = %d", nmax);
      error->all(FLERR, errstr);
    }

    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf, maxbuf, "dump:buf");
  }
}

/* ---------------------------------------------------------------------------------------------------- */

void AtomicStrain::grow_ibuf(int num_me)
{
  int nmax;
  MPI_Allreduce(&num_me, &nmax, 1, MPI_INT, MPI_MAX, world);

  if (nmax > maxibuf) {
    if ((bigint) nmax > MAXSMALLINT) {
      char errstr[100];  
      sprintf(errstr, "Too much per-proc info for atomic_strain command nmax = %d", nmax);
      error->all(FLERR, errstr);
    }

    maxibuf = nmax;
    memory->destroy(ibuf);
    memory->create(ibuf, maxibuf, "dump:ibuf");
  }
}

/* ---------------------------------------------------------------------------------------------------- */

void AtomicStrain::comm_buf_tecio(int num_me)
{
  int nlines, tmp, success;
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf, maxbuf, MPI_DOUBLE, me+iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, me+iproc, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &nlines);
      } else nlines = num_me;
      success = TECDAT142(&nlines, buf, &dummy1);
      if (success == -1) error->one(FLERR, "Cannot write to dump file");
    }
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(buf, num_me, MPI_DOUBLE, 0, 0, world);
  }
}

/* ---------------------------------------------------------------------------------------------------- */

void AtomicStrain::pack_atom_tecplot(int output_format, int output_configuration)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x;
  if (output_configuration == REFERENCE) x = atom->x;
  else x = atom->x_current;

  int m = 0; 
  if (output_format == ASCII) {
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
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
          buf[m++] = list->numneigh[i];
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
        if (rotation_quat_flag) {
          buf[m++] = atomQ[i][0];
          buf[m++] = atomQ[i][1];
          buf[m++] = atomQ[i][2];
          buf[m++] = atomQ[i][3];
        }
        if (rotation_vect_flag) {
          buf[m++] = atomW[i][0];
          buf[m++] = atomW[i][1];
          buf[m++] = atomW[i][2];
        }
        if (displacement_flag) {
          buf[m++] = atomD[i][0];
          buf[m++] = atomD[i][1];
          buf[m++] = atomD[i][2];
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
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
        if (valueflag == COORDINATE)
          buf[m++] = x[i][ival];
        else if (valueflag == ID) 
          buf[m++] = tag[i];
        else if (valueflag == ETYPE) 
          buf[m++] = 0;
        else if (valueflag == CTYPE) 
          buf[m++] = type[i];
        else if (valueflag == NUMNEIGH) 
          buf[m++] = list->numneigh[i];
        else if (valueflag == VALIDITY) 
          buf[m++] = atomValidity[i];
        else if (valueflag == VONMISES) 
          buf[m++] = atomVonMisesStrain[i];
        else if (valueflag == VOLSTRAIN)
          buf[m++] = atomVolStrain[i];
        else if (valueflag == DEFGRAD)
          buf[m++] = atomF[i][ival];
        else if (valueflag == ROTATIONQUAT)
          buf[m++] = atomQ[i][ival];
        else if (valueflag == ROTATIONVECT)
          buf[m++] = atomW[i][ival];
        else if (valueflag == DISPLACEMENT)
          buf[m++] = atomD[i][ival];
        else if (valueflag == NSD) 
          buf[m++] = atomNSD[i];
        else if (valueflag == STRAINTENSOR)
          buf[m++] = atomStrainTensor[i][ival];
        else if (valueflag == STRETCHTENSOR)
          buf[m++] = atomStretchTensor[i][ival];
      }
  }
}

/* ---------------------------------------------------------------------------------------------------- */

void AtomicStrain::pack_node_tecplot(int output_format, int output_configuration)
{
  int nlocal = element->nlocal;
  int *mask = element->mask;
  int *etype = element->etype;
  int **ctype = element->ctype;
  tagint *tag = element->tag;
  int *npe = element->npe;
  int *apc = element->apc;
  double ****nodex;
  double coord[3];
  int inpe, iapc, ietype;
  if (output_configuration == REFERENCE) nodex = element->nodex;
  else nodex = element->nodex_current;

  int m = 0; 
  if (output_format == ASCII) {
    if (average_flag) {
      for (int i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          ietype = etype[i];
          inpe = npe[ietype]; 
          iapc = apc[ietype]; 
          for (int k = 0; k < inpe; k++) {
            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodex[i][j][k][0];
            buf[m++] /= iapc;
            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodex[i][j][k][1];
            buf[m++] /= iapc;
            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodex[i][j][k][2];
            buf[m++] /= iapc;
            if (element_info_flag) {
              buf[m++] = tag[i];
              buf[m++] = ietype;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += ctype[i][j];
              buf[m++] /= iapc;
              buf[m++] = 0.0;
            }
            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodeValidity[i][j][k];
            buf[m++] /= iapc;

            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodeVonMisesStrain[i][j][k];
            buf[m++] /= iapc;

            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) 
              buf[m] += nodeVolStrain[i][j][k];
            buf[m++] /= iapc;
            if (deformation_gradient_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][2];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][3];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][4];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][5];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][6];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][7];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeF[i][j][k][8];
              buf[m++] /= iapc;
            }
            if (rotation_quat_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeQ[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeQ[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeQ[i][j][k][2];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeQ[i][j][k][3];
              buf[m++] /= iapc;
            }
            if (rotation_vect_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeW[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeW[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeW[i][j][k][2];
              buf[m++] /= iapc;
            }
            if (displacement_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeD[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeD[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeD[i][j][k][2];
              buf[m++] /= iapc;
            }
            if (nsd_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeNSD[i][j][k];
              buf[m++] /= iapc;
            }
            if (strain_tensor_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][2];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][3];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][4];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStrainTensor[i][j][k][5];
              buf[m++] /= iapc;
            }
            if (stretch_tensor_flag) {
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][0];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][1];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][2];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][3];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][4];
              buf[m++] /= iapc;
              buf[m] = 0.0;
              for (int j = 0; j < iapc; j++) 
                buf[m] += nodeStretchTensor[i][j][k][5];
              buf[m++] /= iapc;
            }
          }
        }
    } else {
      for (int i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          ietype = etype[i];
          inpe = npe[ietype]; 
          iapc = apc[ietype]; 
          for (int j = 0; j < iapc; j++) 
            for (int k = 0; k < inpe; k++) {
              buf[m++] = nodex[i][j][k][0];
              buf[m++] = nodex[i][j][k][1];
              buf[m++] = nodex[i][j][k][2];
              if (element_info_flag) {
                buf[m++] = tag[i];
                buf[m++] = ietype;
                buf[m++] = ctype[i][j];
                buf[m++] = 0.0;
              }
              buf[m++] = nodeValidity[i][j][k];
              buf[m++] = nodeVonMisesStrain[i][j][k];
              buf[m++] = nodeVolStrain[i][j][k];
              if (deformation_gradient_flag) {
                buf[m++] = nodeF[i][j][k][0];
                buf[m++] = nodeF[i][j][k][1];
                buf[m++] = nodeF[i][j][k][2];
                buf[m++] = nodeF[i][j][k][3];
                buf[m++] = nodeF[i][j][k][4];
                buf[m++] = nodeF[i][j][k][5];
                buf[m++] = nodeF[i][j][k][6];
                buf[m++] = nodeF[i][j][k][7];
                buf[m++] = nodeF[i][j][k][8];
              }
              if (rotation_quat_flag) {
                buf[m++] = nodeQ[i][j][k][0];
                buf[m++] = nodeQ[i][j][k][1];
                buf[m++] = nodeQ[i][j][k][2];
                buf[m++] = nodeQ[i][j][k][3];
              }
              if (rotation_vect_flag) {
                buf[m++] = nodeW[i][j][k][0];
                buf[m++] = nodeW[i][j][k][1];
                buf[m++] = nodeW[i][j][k][2];
              }
              if (displacement_flag) {
                buf[m++] = nodeD[i][j][k][0];
                buf[m++] = nodeD[i][j][k][1];
                buf[m++] = nodeD[i][j][k][2];
              }

              if (nsd_flag)
                buf[m++] = nodeNSD[i][j][k];

              if (strain_tensor_flag) {
                buf[m++] = nodeStrainTensor[i][j][k][0];
                buf[m++] = nodeStrainTensor[i][j][k][1];
                buf[m++] = nodeStrainTensor[i][j][k][2];
                buf[m++] = nodeStrainTensor[i][j][k][3];
                buf[m++] = nodeStrainTensor[i][j][k][4];
                buf[m++] = nodeStrainTensor[i][j][k][5];
              }
              if (stretch_tensor_flag) {
                buf[m++] = nodeStretchTensor[i][j][k][0];
                buf[m++] = nodeStretchTensor[i][j][k][1];
                buf[m++] = nodeStretchTensor[i][j][k][2];
                buf[m++] = nodeStretchTensor[i][j][k][3];
                buf[m++] = nodeStretchTensor[i][j][k][4];
                buf[m++] = nodeStretchTensor[i][j][k][5];
              }
            }
        }
    }
  } else {
    if (average_flag) {
      for (int i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          inpe = npe[etype[i]]; 
          iapc = apc[etype[i]]; 
          for (int k = 0; k < inpe; k++) {
            buf[m] = 0.0;
            for (int j = 0; j < iapc; j++) {
              if (valueflag == COORDINATE)
                buf[m] += nodex[i][j][k][ival];
              else if (valueflag == ID) 
                buf[m] += tag[i];
              else if (valueflag == ETYPE) 
                buf[m] += etype[i];
              else if (valueflag == CTYPE) 
                buf[m] += ctype[i][j];
              else if (valueflag == VALIDITY) 
                buf[m] += nodeValidity[i][j][k];
              else if (valueflag == VONMISES) 
                buf[m] += nodeVonMisesStrain[i][j][k];
              else if (valueflag == VOLSTRAIN)
                buf[m] += nodeVolStrain[i][j][k];
              else if (valueflag == DEFGRAD)
                buf[m] += nodeF[i][j][k][ival];
              else if (valueflag == ROTATIONQUAT)
                buf[m] += nodeQ[i][j][k][ival];
              else if (valueflag == ROTATIONVECT)
                buf[m] += nodeW[i][j][k][ival];
              else if (valueflag == DISPLACEMENT)
                buf[m] += nodeD[i][j][k][ival];
              else if (valueflag == NSD) 
                buf[m] += nodeNSD[i][j][k];
              else if (valueflag == STRAINTENSOR)
                buf[m] += nodeStrainTensor[i][j][k][ival];
              else if (valueflag == STRETCHTENSOR)
                buf[m] += nodeStretchTensor[i][j][k][ival];
            }
            buf[m++] /= iapc; 
          }
        }

    } else {
      for (int i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          inpe = npe[etype[i]]; 
          iapc = apc[etype[i]]; 
          for (int j = 0; j < iapc; j++) 
            for (int k = 0; k < inpe; k++) {
              if (valueflag == COORDINATE)
                buf[m++] = nodex[i][j][k][ival];
              else if (valueflag == ID) 
                buf[m++] = tag[i];
              else if (valueflag == ETYPE) 
                buf[m++] = etype[i];
              else if (valueflag == CTYPE) 
                buf[m++] = ctype[i][j];
              else if (valueflag == VALIDITY) 
                buf[m++] = nodeValidity[i][j][k];
              else if (valueflag == VONMISES) 
                buf[m++] = nodeVonMisesStrain[i][j][k];
              else if (valueflag == VOLSTRAIN)
                buf[m++] = nodeVolStrain[i][j][k];
              else if (valueflag == DEFGRAD)
                buf[m++] = nodeF[i][j][k][ival];
              else if (valueflag == ROTATIONQUAT)
                buf[m++] = nodeQ[i][j][k][ival];
              else if (valueflag == ROTATIONVECT)
                buf[m++] = nodeW[i][j][k][ival];
              else if (valueflag == DISPLACEMENT)
                buf[m++] = nodeD[i][j][k][ival];
              else if (valueflag == NSD) 
                buf[m++] = nodeNSD[i][j][k];
              else if (valueflag == STRAINTENSOR)
                buf[m++] = nodeStrainTensor[i][j][k][ival];
              else if (valueflag == STRETCHTENSOR)
                buf[m++] = nodeStretchTensor[i][j][k][ival];
            }
        }
    }
  }
}

/* ------------------------------------------------------------------------------------------------- */

void AtomicStrain::write_lines_tecplot(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, "%g %g %g", mybuf[m], mybuf[m+1], mybuf[m+2]);
    m += 3;
    if (element_info_flag) {
      fprintf(fp, TAGINT_FORMAT " %d %d %d", 
          static_cast<tagint> (mybuf[m]), 
          static_cast<int> (mybuf[m+1]), 
          static_cast<int> (mybuf[m+2]),
          static_cast<int> (mybuf[m+3]));
      m += 4;
    }
    fprintf(fp, " %d %g %g", static_cast<int> (mybuf[m]), mybuf[m+1], mybuf[m+2]);
    m += 3;
    if (deformation_gradient_flag) {
      fprintf(fp, " %g %g %g %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5], 
          mybuf[m+6], mybuf[m+7], mybuf[m+8]);
      m += 9;
    }
    if (rotation_quat_flag) {
      fprintf(fp, " %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], mybuf[m+3]);
      m += 4;
    }
    if (rotation_vect_flag) {
      fprintf(fp, " %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2]);
      m += 3;
    }
    if (displacement_flag) {
      fprintf(fp, " %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2]);
      m += 3;
    }
    if (nsd_flag) {
      fprintf(fp, " %g", mybuf[m]);
      m += 1;
    }
    if (strain_tensor_flag) {
      fprintf(fp, " %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5]);
      m += 6;
    }
    if (stretch_tensor_flag) {
      fprintf(fp, " %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5]);
      m += 6;
    }
    fprintf(fp, "\n");
  }
}

/* ------------------------------------------------------------------------------------------ */

void AtomicStrain::pack_atom_all(int output_format, int output_configuration)
{
  int m;
  int i, j, k, iucell, inode, ibasis, iapc, inode_local, ietype;
  double **x;
  if (output_configuration == REFERENCE) x = atom->x;
  else x = atom->x_current;

  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint maxtag, maxtag_all;

  m = 0;
  if (no_atom_flag == 0)
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (wrap_flag) {
          if (output_configuration == REFERENCE) domain->remap(x[i]);
          else domain->remap_current(x[i]);
        }

        maxtag = MAX(maxtag, tag[i]);
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
          buf[m++] = -1;
          buf[m++] = list->numneigh[i];
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
        if (rotation_quat_flag) {
          buf[m++] = atomQ[i][0];
          buf[m++] = atomQ[i][1];
          buf[m++] = atomQ[i][2];
          buf[m++] = atomQ[i][3];
        }

        if (rotation_vect_flag) {
          buf[m++] = atomW[i][0];
          buf[m++] = atomW[i][1];
          buf[m++] = atomW[i][2];
        }

        if (displacement_flag) {
          buf[m++] = atomD[i][0];
          buf[m++] = atomD[i][1];
          buf[m++] = atomD[i][2];
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

  if (output_format == ATOM) return;

  if (no_atom_flag == 0) {
    maxtag = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) 
        maxtag = MAX(maxtag, tag[i]);

    MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_CAC_TAGINT, MPI_MAX, world);
  } else maxtag_all = 0;

  // pack atoms from elements 

  if (no_elem_flag == 0) { 
    nlocal = element->nlocal;
    mask = element->mask;
    double ****nodex;
    if (output_configuration == REFERENCE) nodex = element->nodex;
    else nodex = element->nodex_current;
    int *nucell = element->nucell;
    int *etype = element->etype;
    int **ctype = element->ctype;
    tag = element->tag;
    int *npe = element->npe;
    int *apc = element->apc;
    int maxapc = element->maxapc;
    int **u2g = element->u2g;
    int **n2u = element->n2u;
    int **g2n = element->g2n;
    ElementVec *evec = element->evec;

    // set IDs for virtual atoms from element tag IDs to have consistent tag IDs across different simulations of the same system

    tagint itag = maxtag_all + 1;
    if (output_format == FULLATOM) {
      int maxucell = element->maxucell;
      for (i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          ietype = etype[i];
          iapc = apc[ietype];
          for (ibasis = 0; ibasis < iapc; ibasis++) {
            for (iucell = 0; iucell < nucell[ietype]; iucell++) {

              buf[m++] = itag + 
                ((tag[i]-1) * maxucell + iucell) * maxapc + ibasis;
              buf[m++] = ctype[i][ibasis];
              for (j = 0; j < 3; j++) 
                buf[m++] = evec->interpolate(nodex, i, ibasis, iucell, j);
              if (wrap_flag) {
                if (output_configuration == REFERENCE) domain->remap(&buf[m-3]);
                else domain->remap_current(&buf[m-3]);
              }

              if (u2g[ietype][iucell] < 0) inode = -1;
              else inode = g2n[ietype][u2g[ietype][iucell]];
              if (element_info_flag) {
                buf[m++] = ietype;
                buf[m++] = tag[i];
                buf[m++] = inode;
                buf[m++] = iucell;
                buf[m++] = ibasis;
                buf[m++] = 0.0;
              }

              buf[m++] = evec->interpolate(nodeValidity, i, ibasis, iucell);
              buf[m++] = evec->interpolate(nodeVonMisesStrain, i, ibasis, iucell);
              buf[m++] = evec->interpolate(nodeVolStrain, i, ibasis, iucell);

              if (deformation_gradient_flag) 
                for (j = 0; j < 9; j++) 
                  buf[m++] = evec->interpolate(nodeF, i, ibasis, iucell, j);

              if (rotation_quat_flag) 
                for (j = 0; j < 4; j++) 
                  buf[m++] = evec->interpolate(nodeQ, i, ibasis, iucell, j);

              if (rotation_vect_flag) 
                for (j = 0; j < 3; j++) 
                  buf[m++] = evec->interpolate(nodeW, i, ibasis, iucell, j);

              if (displacement_flag) 
                for (j = 0; j < 3; j++) 
                  buf[m++] = evec->interpolate(nodeD, i, ibasis, iucell, j);

              if (nsd_flag) 
                buf[m++] = evec->interpolate(nodeNSD, i, ibasis, iucell);

              if (strain_tensor_flag) 
                for (j = 0; j < 6; j++) 
                  buf[m++] = evec->interpolate(nodeStrainTensor, i, ibasis, iucell, j);

              if (stretch_tensor_flag) 
                for (j = 0; j < 6; j++) 
                  buf[m++] = evec->interpolate(nodeStretchTensor, i, ibasis, iucell, j);
            }
          }
        }  
    } else if (output_format == NODE) {
      int maxnpe = element->maxnpe;
      for (i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          ietype = etype[i];
          iapc = apc[ietype];
          for (ibasis = 0; ibasis < iapc; ibasis++) {
            for (inode = 0; inode < npe[ietype]; inode++) {

              buf[m++] = itag + ((tag[i]-1) * maxnpe + inode) * maxapc + ibasis;
              buf[m++] = ctype[i][ibasis];
              buf[m++] = nodex[i][ibasis][inode][0];
              buf[m++] = nodex[i][ibasis][inode][1];
              buf[m++] = nodex[i][ibasis][inode][2];
              if (wrap_flag) {
                if (output_configuration == REFERENCE) domain->remap(&buf[m-3]);
                else domain->remap_current(&buf[m-3]);
              }

              if (element_info_flag) {
                buf[m++] = ietype;
                buf[m++] = tag[i];
                buf[m++] = inode;
                buf[m++] = n2u[ietype][inode];
                buf[m++] = ibasis;
                buf[m++] = 0.0;
              }

              buf[m++] = nodeValidity[i][ibasis][iucell];
              buf[m++] = nodeVonMisesStrain[i][ibasis][iucell];
              buf[m++] = nodeVolStrain[i][ibasis][iucell];

              if (deformation_gradient_flag) 
                for (j = 0; j < 9; j++) 
                  buf[m++] = nodeF[i][ibasis][iucell][j];

              if (rotation_quat_flag) 
                for (j = 0; j < 4; j++) 
                  buf[m++] = nodeQ[i][ibasis][iucell][j];

              if (rotation_vect_flag) 
                for (j = 0; j < 3; j++) 
                  buf[m++] = nodeW[i][ibasis][iucell][j];

              if (displacement_flag) 
                for (j = 0; j < 3; j++) 
                  buf[m++] = nodeD[i][ibasis][iucell][j];

              if (nsd_flag) 
                buf[m++] = nodeNSD[i][ibasis][iucell];

              if (strain_tensor_flag) 
                for (j = 0; j < 6; j++) 
                  buf[m++] = nodeStrainTensor[i][ibasis][iucell][j];

              if (stretch_tensor_flag) 
                for (j = 0; j < 6; j++) 
                  buf[m++] = nodeStretchTensor[i][ibasis][iucell][j];
            }
          }
        }  
    }
  }
}

/* ------------------------------------------------------------------------------------------------- */

void AtomicStrain::write_lines_atom(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp, TAGINT_FORMAT " %d %g %g %g", 
        static_cast<tagint> (mybuf[m]), static_cast<int> (mybuf[m+1]), 
        mybuf[m+2], mybuf[m+3], mybuf[m+4]);
    m += 5;
    if (element_info_flag) {
      fprintf(fp, " %d %d %d %d %d %d", 
          static_cast<int> (mybuf[m]), 
          static_cast<int> (mybuf[m+1]), 
          static_cast<int> (mybuf[m+2]), 
          static_cast<int> (mybuf[m+3]), 
          static_cast<int> (mybuf[m+4]), 
          static_cast<int> (mybuf[m+5]));
      m += 6;
    }
    fprintf(fp, " %d %g %g", static_cast<int> (mybuf[m]), mybuf[m+1], mybuf[m+2]);
    m += 3;
    if (deformation_gradient_flag) {
      fprintf(fp, " %g %g %g %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5], 
          mybuf[m+6], mybuf[m+7], mybuf[m+8]);
      m += 9;
    }
    if (rotation_quat_flag) {
      fprintf(fp, " %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], mybuf[m+3]);
      m += 4;
    }
    if (rotation_vect_flag) {
      fprintf(fp, " %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2]);
      m += 3;
    }
    if (displacement_flag) {
      fprintf(fp, " %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2]);
      m += 3;
    }
    if (nsd_flag) {
      fprintf(fp, " %g", mybuf[m]);
      m += 1;
    }
    if (strain_tensor_flag) {
      fprintf(fp, " %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5]);
      m += 6;
    }
    if (stretch_tensor_flag) {
      fprintf(fp, " %g %g %g %g %g %g", 
          mybuf[m], mybuf[m+1], mybuf[m+2], 
          mybuf[m+3], mybuf[m+4], mybuf[m+5]);
      m += 6;
    }
    fprintf(fp, "\n");
  }
}


