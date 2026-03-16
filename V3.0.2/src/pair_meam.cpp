#include "pair_meam.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "meam.h"
#include "atom.h"
#include "element.h"
#include "force.h"
#include "universe.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define MAXLINE 1024

static const int nkeywords = 22;
static const char *keywords[] = {
  "Ec", "alpha", "rho0", "delta", "lattce", 
  "attrac", "repuls", "nn2", "Cmin", "Cmax", "rc", "delr", 
  "augt1", "gsmooth_factor", "re", "ialloy", 
  "mixture_ref_t", "erose_form", "zbl", 
  "emb_lin_neg", "bkgd_dyn", "theta"};

/*  ----------------------------------------------------------------------  */

PairMEAM::PairMEAM(CAC *cac) : Pair(cac)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  threebody_flag = 1;
  allocated = 0;

  nelements = 0;
  elements = NULL;
  mass = NULL;
  meam_inst = new MEAM(memory, atom, element, error);
  meam_inst->cac = cac;
  scale = NULL;

  // set comm size needed by this Pair

  comm_atom_forward = 38;
  comm_elem_forward = comm_atom_forward * element->maxnpe * element->maxapc;
  comm_atom_reverse = 30;
  comm_elem_reverse = comm_atom_reverse * element->maxnpe * element->maxapc;
}

/*  ----------------------------------------------------------------------
    free all arrays
    check if allocated, since class can be destructed when incomplete
    -------------------------------------------------------------------------  */

PairMEAM::~PairMEAM()
{
  delete meam_inst;

  for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  delete [] mass;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
    delete [] map;
  }
}

/*  ----------------------------------------------------------------------  */

void PairMEAM::compute(int eflag, int vflag)
{
  int i, ii, n, errorflag;
  int ietype, inode, iindex, inpe, iapc;
  ev_init(eflag, vflag);

  // pointers and parameters from neighbor list

  int inum = listhalf->inum;
  int *ilist = listhalf->ilist;
  int *iindexlist = listhalf->iindexlist;
  int *numneigh_half = listhalf->numneigh;
  int **firstneigh_half = listhalf->firstneigh;
  int **firstneighindex_half = listhalf->firstneighindex;

  int *numneigh_full = listfull->numneigh;
  int **firstneigh_full = listfull->firstneigh;
  int **firstneighindex_full = listfull->firstneighindex;

  int *nvanumneigh = listfull->nvanumneigh;
  int **nvafirstneigh = listfull->nvafirstneigh;
  int **nvafirstneighindex = listfull->nvafirstneighindex;
  int ***va2nvalist = listfull->va2nvalist;

  // strip neighbor lists of any special bond flags before using with MEAM
  // necessary before doing neigh_f2c and neigh_c2f conversions each step

  //if (neighbor->ago == 0) {
  //  neigh_strip(inum_half, ilist_half, numneigh_half, firstneigh_half);
  //  neigh_strip(inum_half, ilist_half, numneigh_full, firstneigh_full);
  //}

  // check size of scrfcn based on half neighbor list

  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int naall = nalocal + atom->nghost;
  int neall = nelocal + element->nghost;
  int *etype = element->etype;
  int *npe = element->npe;
  int *apc = element->apc;

  n = 0;

  for (ii = 0; ii < inum; ii++) 
    n += numneigh_half[ii];

  meam_inst->meam_dens_setup(atom->nmax, element->nmax, naall, neall, n);

  // 3 stages of MEAM calculation
  // - 1st stage:
  //    + Compute screening function and derivatives (scrfcn, dscrfcn, and fcpair)
  //    + Compute intermediate density terms to be communicated (only for atoms' electron densities). Variables to be computed:
  //    t_ave, tsq_ave, rho0, arho2b, arho1, arho2, arho3, arho3b 
  //    + NOTE: All these quantities for elements are calculated at the nodes only, 
  //            for virtual atoms, simply interpolate the value from nodes

  // - 2nd stage:
  //    Calculate the remaining density terms
  // loop over my atoms followed by communication
  //
  // - Forward comm
  // - 3rd stage: force calculation

  int offset = 0;
  errorflag = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];
    meam_inst->meam_dens_init(i, iindex, map, offset,
        numneigh_half[ii], firstneigh_half[ii], firstneighindex_half[ii], 
        numneigh_full[ii], firstneigh_full[ii], firstneighindex_full[ii]);
    offset += numneigh_half[ii];
  }

  comm->reverse_comm_pair(this);

  meam_inst->meam_dens_final(eflag_either, eflag_global, eflag_atom, 
      &eng_vdwl, eatom, enode, map, scale, errorflag);

  //meam_inst->write_dens(); 
  if (errorflag) {
    char str[128];
    sprintf(str, "MEAM library error %d", errorflag);
    error->one(FLERR, str);
  }

  comm->forward_comm_pair(this);

  offset = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];
    meam_inst->meam_force(i, iindex, eflag_either, eflag_global, eflag_atom, 
        vflag_either, vflag_global, vflag_atom, virial, &eng_vdwl, eatom, enode, map, scale, 
        numneigh_half[ii], firstneigh_half[ii], firstneighindex_half[ii], 
        numneigh_full[ii], firstneigh_full[ii], firstneighindex_full[ii], 
        nvanumneigh, nvafirstneigh, nvafirstneighindex, va2nvalist,
        offset, vatom, vnode);
    offset += numneigh_half[ii];
  }

  if (vflag_fdotr) 
    virial_fdotr_compute();

}

/*  ----------------------------------------------------------------------  */

void PairMEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(scale, n+1, n+1, "pair:scale");

  map = new int[n+1];
}

/*  ----------------------------------------------------------------------
   global settings
   -------------------------------------------------------------------------  */

void PairMEAM::settings(int narg, char ** /* arg */)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/*  ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   -------------------------------------------------------------------------  */

void PairMEAM::coeff(int narg, char **arg)
{
  int m, n;

  if (!allocated) allocate();

  if (narg < 6) error->all(FLERR, "Incorrect args for pair coefficients");

  // insure I, J args are **

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // read MEAM element names between 2 filenames
  // nelements = # of MEAM elements
  // elements = list of unique element names

  if (nelements) {
    for (int i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
    delete [] mass;
  }
  nelements = narg - 4 - atom->ntypes;
  if (nelements < 1) {
    //printf("narg = %d nelements = %d ntypes = %d\n", narg, nelements, atom->ntypes);
    error->all(FLERR, "Incorrect args for pair coefficients");
  }
  if (nelements > maxelt)
    error->all(FLERR, "Too many elements extracted from MEAM library. "
        "Increase 'maxelt' in meam.h and recompile.");
  elements = new char*[nelements];
  mass = new double[nelements];

  for (int i = 0; i < nelements; i++) {
    n = strlen(arg[i+3]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i], arg[i+3]);
  }

  // read MEAM library and parameter files
  // pass all parameters to MEAM package
  // tell MEAM package that setup is done

  read_files(arg[2], arg[2+nelements+1]);
  meam_inst->meam_setup_done(&cutmax);

  // read args that map atom types to MEAM elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (int i = 4 + nelements; i < narg; i++) {
    m = i - (4+nelements) + 1;
    int j;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i], elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i], "NULL") == 0) map[m] = -1;
    else error->all(FLERR, "Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I, J = **

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i, j for type pairs where both are mapped to elements
  // set mass for i, i in atom class

  int count = 0;
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR, i, mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*  ----------------------------------------------------------------------
   init specific to this pair style
   -------------------------------------------------------------------------  */

void PairMEAM::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style MEAM requires newton pair on");

  // need full double list and half neighbor list
  int irequest_full = neighbor->request(this, instance_me);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
  neighbor->requests[irequest_full]->threebody = 1;
  int irequest_half = neighbor->request(this, instance_me);
  neighbor->requests[irequest_half]->id = 2;
  neighbor->requests[irequest_half]->halffull = 1;
  neighbor->requests[irequest_half]->halffulllist = irequest_full;

}

/*  ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
   -------------------------------------------------------------------------  */

void PairMEAM::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = list = ptr;
  else if (id == 2) listhalf = ptr;

}


/*  ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
   -------------------------------------------------------------------------  */

double PairMEAM::init_one(int i, int j)
{

  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  return cutmax;
}

/*  ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
   -------------------------------------------------------------------------  */

double PairMEAM::init_one_inner(int i, int j)
{
  return cutmax;
}

/*  ----------------------------------------------------------------------  */

void PairMEAM::read_files(char *globalfile, char *userfile)
{
  // open global meamf file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(globalfile);
    if (fp == NULL) {
      char str[128];
      snprintf(str, 128, "Cannot open MEAM potential file %s", globalfile);
      error->one(FLERR, str);
    }
  }

  // allocate parameter arrays

  int params_per_line = 19;

  lattice_t *lat = new lattice_t[nelements];
  int *ielement = new int[nelements];
  int *ibar = new int[nelements];
  double *z = new double[nelements];
  double *atwt = new double[nelements];
  double *alpha = new double[nelements];
  double *b0 = new double[nelements];
  double *b1 = new double[nelements];
  double *b2 = new double[nelements];
  double *b3 = new double[nelements];
  double *alat = new double[nelements];
  double *esub = new double[nelements];
  double *asub = new double[nelements];
  double *t0 = new double[nelements];
  double *t1 = new double[nelements];
  double *t2 = new double[nelements];
  double *t3 = new double[nelements];
  double *rozero = new double[nelements];

  bool *found = new bool[nelements];
  for (int i = 0; i < nelements; i++) found[i] = false;

  // read each set of params from global MEAM file
  // one set of params can span multiple lines
  // store params if element name is in element list
  // if element name appears multiple times, only store 1st entry

  int i, n, nwords;
  char **words = new char*[params_per_line+1];
  char line[MAXLINE], *ptr;
  int eof = 0;

  int nset = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    nwords = universe->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n], MAXLINE-n, fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof, 1, MPI_INT, 0, world);
      if (eof) break;
      MPI_Bcast(&n, 1, MPI_INT, 0, world);
      MPI_Bcast(line, n, MPI_CHAR, 0, world);
      if ((ptr = strchr(line, '#'))) *ptr = '\0';
      nwords = universe->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR, "Incorrect format in MEAM potential file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    nwords = 0;
    words[nwords++] = strtok(line, "' \t\n\r\f");
    while ((words[nwords++] = strtok(NULL, "' \t\n\r\f"))) continue;

    // skip if element name isn't in element list

    for (i = 0; i < nelements; i++)
      if (strcmp(words[0], elements[i]) == 0) break;
    if (i >= nelements) continue;

    // skip if element already appeared

    if (found[i] == true) continue;
    found[i] = true;

    // map lat string to an integer

    if (!MEAM::str_to_lat(words[1], true, lat[i]))
      error->all(FLERR, "Unrecognized lattice type in MEAM file 1");

    // store parameters

    z[i] = atof(words[2]);
    ielement[i] = atoi(words[3]);
    atwt[i] = atof(words[4]);
    alpha[i] = atof(words[5]);
    b0[i] = atof(words[6]);
    b1[i] = atof(words[7]);
    b2[i] = atof(words[8]);
    b3[i] = atof(words[9]);
    alat[i] = atof(words[10]);
    esub[i] = atof(words[11]);
    asub[i] = atof(words[12]);
    t0[i] = atof(words[13]);
    t1[i] = atof(words[14]);
    t2[i] = atof(words[15]);
    t3[i] = atof(words[16]);
    rozero[i] = atof(words[17]);
    ibar[i] = atoi(words[18]);

    if (!isone(t0[i]))
      error->all(FLERR, "Unsupported parameter in MEAM potential file: t0!=1");

    // z given is ignored: if this is mismatched, we definitely won't do what the user said -> fatal error
    if (z[i] != MEAM::get_Zij(lat[i]))
      error->all(FLERR, "Mismatched parameter in MEAM potential file: z!=lat");

    nset++;
  }

  // error if didn't find all elements in file

  if (nset != nelements)
    error->all(FLERR, "Did not find all elements in MEAM library file");

  // pass element parameters to MEAM package

  meam_inst->meam_setup_global(nelements, lat, ielement, atwt, alpha, b0, b1, b2, b3, 
      alat, esub, asub, t0, t1, t2, t3, rozero, ibar);

  // set element masses

  for (i = 0; i < nelements; i++) mass[i] = atwt[i];

  // clean-up memory

  delete [] words;

  delete [] lat;
  delete [] ielement;
  delete [] ibar;
  delete [] z;
  delete [] atwt;
  delete [] alpha;
  delete [] b0;
  delete [] b1;
  delete [] b2;
  delete [] b3;
  delete [] alat;
  delete [] esub;
  delete [] asub;
  delete [] t0;
  delete [] t1;
  delete [] t2;
  delete [] t3;
  delete [] rozero;
  delete [] found;

  // done if user param file is NULL

  if (strcmp(userfile, "NULL") == 0) return;

  // open user param file on proc 0

  if (comm->me == 0) {
    fp = force->open_potential(userfile);
    if (fp == NULL) {
      char str[128];
      snprintf(str, 128, "Cannot open MEAM potential file %s", userfile);
      error->one(FLERR, str);
    }
  }

  // read settings
  // pass them one at a time to MEAM package
  // match strings to list of corresponding ints

  int which;
  double value;
  lattice_t latt;
  int nindex, index[3];
  int maxparams = 6;
  char **params = new char*[maxparams];
  int nparams;

  eof = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    nparams = universe->count_words(line);
    if (nparams == 0) continue;

    // words = ptrs to all words in line

    nparams = 0;
    params[nparams++] = strtok(line, "=(), '\t\n\r\f");
    while (nparams < maxparams &&
        (params[nparams++] = strtok(NULL, "=(), '\t\n\r\f")))
      continue;
    nparams--;

    for (which = 0; which < nkeywords; which++)
      if (strcmp(params[0], keywords[which]) == 0) break;
    if (which == nkeywords) {
      char str[128];
      snprintf(str, 128, "Keyword %s in MEAM parameter file not recognized", 
          params[0]);
      error->all(FLERR, str);
    }
    nindex = nparams - 2;
    for (i = 0; i < nindex; i++) index[i] = atoi(params[i+1]) - 1;

    // map lattce_meam value to an integer

    if (which == 4) {
      if (!MEAM::str_to_lat(params[nparams-1], false, latt))
        error->all(FLERR, "Unrecognized lattice type in MEAM file 2");
      value = latt;
    }
    else value = atof(params[nparams-1]);

    // pass single setting to MEAM package

    int errorflag = 0;
    meam_inst->meam_setup_param(which, value, nindex, index, &errorflag);
    if (errorflag) {
      char str[128];
      sprintf(str, "MEAM library error %d", errorflag);
      error->all(FLERR, str);
    }
  }

  delete [] params;
}

/*  ----------------------------------------------------------------------  */

int PairMEAM::pack_atom_forward_comm(int n, int *list, double *buf, 
    int /* pbc_flag */, int * /* pbc */)
{
  int i, j, k, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = meam_inst->atomrho0[j];
    buf[m++] = meam_inst->atomrho1[j];
    buf[m++] = meam_inst->atomrho2[j];
    buf[m++] = meam_inst->atomrho3[j];
    buf[m++] = meam_inst->atomfrhop[j];
    buf[m++] = meam_inst->atomgamma[j];
    buf[m++] = meam_inst->atomdgamma1[j];
    buf[m++] = meam_inst->atomdgamma2[j];
    buf[m++] = meam_inst->atomdgamma3[j];
    buf[m++] = meam_inst->atomarho2b[j];
    for (k = 0; k < 3; k++) {
      buf[m++] = meam_inst->atomarho1[j][k];
      buf[m++] = meam_inst->atomarho3b[j][k];
      buf[m++] = meam_inst->atomt_ave[j][k];
      buf[m++] = meam_inst->atomtsq_ave[j][k];
    }
    for (k = 0; k < 6; k++) buf[m++] = meam_inst->atomarho2[j][k];
    for (k = 0; k < 10; k++) buf[m++] = meam_inst->atomarho3[j][k];
  }

  return m;
}

/*  ----------------------------------------------------------------------  */

void PairMEAM::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    meam_inst->atomrho0[i] = buf[m++];
    meam_inst->atomrho1[i] = buf[m++];
    meam_inst->atomrho2[i] = buf[m++];
    meam_inst->atomrho3[i] = buf[m++];
    meam_inst->atomfrhop[i] = buf[m++];
    meam_inst->atomgamma[i] = buf[m++];
    meam_inst->atomdgamma1[i] = buf[m++];
    meam_inst->atomdgamma2[i] = buf[m++];
    meam_inst->atomdgamma3[i] = buf[m++];
    meam_inst->atomarho2b[i] = buf[m++];
    for (k = 0; k < 3; k++) {
      meam_inst->atomarho1[i][k] = buf[m++];
      meam_inst->atomarho3b[i][k] = buf[m++];
      meam_inst->atomt_ave[i][k] = buf[m++];
      meam_inst->atomtsq_ave[i][k] = buf[m++];
    }
    for (k = 0; k < 6; k++) meam_inst->atomarho2[i][k] = buf[m++];
    for (k = 0; k < 10; k++) meam_inst->atomarho3[i][k] = buf[m++];
  }
}


/*  ----------------------------------------------------------------------  */

int PairMEAM::pack_elem_forward_comm(int n, int *list, double *buf, 
    int /* pbc_flag */, int * /* pbc */)
{
  int ii, i, j, k, l, m;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    i = list[ii];
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) {
        buf[m++] = meam_inst->noderho0[i][j][k];
        buf[m++] = meam_inst->noderho1[i][j][k];
        buf[m++] = meam_inst->noderho2[i][j][k];
        buf[m++] = meam_inst->noderho3[i][j][k];
        buf[m++] = meam_inst->nodefrhop[i][j][k];
        buf[m++] = meam_inst->nodegamma[i][j][k];
        buf[m++] = meam_inst->nodedgamma1[i][j][k];
        buf[m++] = meam_inst->nodedgamma2[i][j][k];
        buf[m++] = meam_inst->nodedgamma3[i][j][k];
        buf[m++] = meam_inst->nodearho2b[i][j][k];
        for (l = 0; l < 3; l++) {
          buf[m++] = meam_inst->nodearho1[i][j][k][l];
          buf[m++] = meam_inst->nodearho3b[i][j][k][l];
          buf[m++] = meam_inst->nodet_ave[i][j][k][l];
          buf[m++] = meam_inst->nodetsq_ave[i][j][k][l];
        }
        for (l = 0; l < 6; l++) buf[m++] = meam_inst->nodearho2[i][j][k][l];
        for (l = 0; l < 10; l++) buf[m++] = meam_inst->nodearho3[i][j][k][l];
      }
  }

  return m;

}

/*  ----------------------------------------------------------------------  */

void PairMEAM::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i, j, k, l, m, last;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) {
        meam_inst->noderho0[i][j][k] = buf[m++];
        meam_inst->noderho1[i][j][k] = buf[m++];
        meam_inst->noderho2[i][j][k] = buf[m++];
        meam_inst->noderho3[i][j][k] = buf[m++];
        meam_inst->nodefrhop[i][j][k] = buf[m++];
        meam_inst->nodegamma[i][j][k] = buf[m++];
        meam_inst->nodedgamma1[i][j][k] = buf[m++];
        meam_inst->nodedgamma2[i][j][k] = buf[m++];
        meam_inst->nodedgamma3[i][j][k] = buf[m++];
        meam_inst->nodearho2b[i][j][k] = buf[m++];
        for (l = 0; l < 3; l++) {
          meam_inst->nodearho1[i][j][k][l] = buf[m++];
          meam_inst->nodearho3b[i][j][k][l] = buf[m++];
          meam_inst->nodet_ave[i][j][k][l] = buf[m++];
          meam_inst->nodetsq_ave[i][j][k][l] = buf[m++];
        }
        for (l = 0; l < 6; l++) meam_inst->nodearho2[i][j][k][l] = buf[m++];
        for (l = 0; l < 10; l++) meam_inst->nodearho3[i][j][k][l] = buf[m++];
      }
  }
}

/*  ----------------------------------------------------------------------  */

int PairMEAM::pack_atom_reverse_comm(int n, int first, double *buf)
{
  int i, k, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = meam_inst->atomrho0[i];
    buf[m++] = meam_inst->atomarho2b[i];
    for (k = 0; k < 3; k++) {
      buf[m++] = meam_inst->atomarho1[i][k];
      buf[m++] = meam_inst->atomarho3b[i][k];
      buf[m++] = meam_inst->atomt_ave[i][k];
      buf[m++] = meam_inst->atomtsq_ave[i][k];
    }
    for (k = 0; k < 6; k++) buf[m++] = meam_inst->atomarho2[i][k];
    for (k = 0; k < 10; k++) buf[m++] = meam_inst->atomarho3[i][k];
  }

  return m;
}

/*  ----------------------------------------------------------------------  */

void PairMEAM::unpack_atom_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    meam_inst->atomrho0[j] += buf[m++];
    meam_inst->atomarho2b[j] += buf[m++];
    for (k = 0; k < 3; k++) {
      meam_inst->atomarho1[j][k] += buf[m++];
      meam_inst->atomarho3b[j][k] += buf[m++];
      meam_inst->atomt_ave[j][k] += buf[m++];
      meam_inst->atomtsq_ave[j][k] += buf[m++];
    }
    for (k = 0; k < 6; k++) meam_inst->atomarho2[j][k] += buf[m++];
    for (k = 0; k < 10; k++) meam_inst->atomarho3[j][k] += buf[m++];
  }
}

/*  ----------------------------------------------------------------------  */

int PairMEAM::pack_elem_reverse_comm(int n, int first, double *buf)
{
  int i, j, k, l, m, last;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) {
        buf[m++] = meam_inst->noderho0[i][j][k];
        buf[m++] = meam_inst->nodearho2b[i][j][k];
        for (l = 0; l < 3; l++) {
          buf[m++] = meam_inst->nodearho1[i][j][k][l];
          buf[m++] = meam_inst->nodearho3b[i][j][k][l];
          buf[m++] = meam_inst->nodet_ave[i][j][k][l];
          buf[m++] = meam_inst->nodetsq_ave[i][j][k][l];
        }
        for (l = 0; l < 6; l++) buf[m++] = meam_inst->nodearho2[i][j][k][l];
        for (l = 0; l < 10; l++) buf[m++] = meam_inst->nodearho3[i][j][k][l];
      }
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

void PairMEAM::unpack_elem_reverse_comm(int n, int *list, double *buf)
{
  int ii, i, j, k, l, m;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    i = list[ii];
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) {
        meam_inst->noderho0[i][j][k] += buf[m++];
        meam_inst->nodearho2b[i][j][k] += buf[m++];
        for (l = 0; l < 3; l++) {
          meam_inst->nodearho1[i][j][k][l] += buf[m++];
          meam_inst->nodearho3b[i][j][k][l] += buf[m++];
          meam_inst->nodet_ave[i][j][k][l] += buf[m++];
          meam_inst->nodetsq_ave[i][j][k][l] += buf[m++];
        }
        for (l = 0; l < 6; l++) meam_inst->nodearho2[i][j][k][l] += buf[m++];
        for (l = 0; l < 10; l++) meam_inst->nodearho3[i][j][k][l] += buf[m++];
      }
  }
}

/*  ----------------------------------------------------------------------
    memory usage of local atom-based arrays
    -------------------------------------------------------------------------  */

double PairMEAM::memory_usage()
{
  double bytes = 37 * (meam_inst->namax + meam_inst->nemax * 
      element->maxapc * element->maxnpe) * sizeof(double);
  bytes += 3 * meam_inst->maxneigh * sizeof(double);
  return bytes;
}

/*  ----------------------------------------------------------------------
    strip special bond flags from neighbor list entries
    are not used with MEAM
    need to do here so Fortran lib doesn't see them
    done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
    -------------------------------------------------------------------------  */

void PairMEAM::neigh_strip(int inum, int *ilist, 
    int *numneigh, int **firstneigh)
{
  int i, j, ii, jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j] &= NEIGHMASK;
  }
}

/*  ----------------------------------------------------------------------  */

void *PairMEAM::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  return NULL;
}
