#include <mpi.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "universe.h"
#include "version.h"
#include "error.h"
#include "memory.h"

using namespace CAC_NS;

#define MAXLINE 256

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
   define frequently used universial functions
   ------------------------------------------------------------------------- */

Universe::Universe(CAC *cac, MPI_Comm communicator) : Pointers(cac)
{
  version = (const char *) CAC_VERSION;
  uworld = uorig = communicator;
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  uscreen = stdout;
  ulogfile = nullptr;

  existflag = 0;
  nworlds = 0;
  procs_per_world = nullptr;
  root_proc = nullptr;

  memory->create(uni2orig,nprocs,"universe:uni2orig");
  for (int i = 0; i < nprocs; i++) uni2orig[i] = i;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  if (uworld != uorig) MPI_Comm_free(&uworld);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(uni2orig);
//  delete [] version;
}

/* ----------------------------------------------------------------------
   reorder universe processors
   create uni2orig as inverse mapping
   re-create uworld communicator with new ordering via Comm_split()
   style = "nth", arg = N
   move every Nth proc to end of rankings
   style = "custom", arg = filename
   file has nprocs lines with I J
   I = universe proc ID in original communicator uorig
   J = universe proc ID in reordered communicator uworld
   ------------------------------------------------------------------------- */

void Universe::reorder(char *style, char *arg)
{
  char line[MAXLINE];

  if (uworld != uorig) MPI_Comm_free(&uworld);

  if (strcmp(style,"nth") == 0) {
    int n = inumeric(FLERR,arg);
    if (n <= 0)
      error->universe_all(FLERR,"Invalid -reorder N value");
    if (nprocs % n)
      error->universe_all(FLERR,"Nprocs not a multiple of N for -reorder");
    for (int i = 0; i < nprocs; i++) {
      if (i < (n-1)*nprocs/n) uni2orig[i] = i/(n-1) * n + (i % (n-1));
      else uni2orig[i] = (i - (n-1)*nprocs/n) * n + n-1;
    }

  } else if (strcmp(style,"custom") == 0) {

    if (me == 0) {
      FILE *fp = fopen(arg,"r");
      if (fp == nullptr) error->universe_one(FLERR,"Cannot open -reorder file");

      // skip header = blank and comment lines

      char *ptr;
      if (!fgets(line,MAXLINE,fp))
        error->one(FLERR,"Unexpected end of -reorder file");
      while (1) {
        if ((ptr = strchr(line,'#'))) *ptr = '\0';
        if (strspn(line," \t\n\r") != strlen(line)) break;
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
      }

      // read nprocs lines
      // uni2orig = inverse mapping

      int me_orig,me_new;
      sscanf(line,"%d %d",&me_orig,&me_new);
      if (me_orig < 0 || me_orig >= nprocs ||
          me_new < 0 || me_new >= nprocs)
        error->one(FLERR,"Invalid entry in -reorder file");
      uni2orig[me_new] = me_orig;

      for (int i = 1; i < nprocs; i++) {
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
        sscanf(line,"%d %d",&me_orig,&me_new);
        if (me_orig < 0 || me_orig >= nprocs ||
            me_new < 0 || me_new >= nprocs)
          error->one(FLERR,"Invalid entry in -reorder file");
        uni2orig[me_new] = me_orig;
      }
      fclose(fp);
    }

    // bcast uni2org from proc 0 to all other universe procs

    MPI_Bcast(uni2orig,nprocs,MPI_INT,0,uorig);

  } else error->universe_all(FLERR,"Invalid command-line argument");

  // create new uworld communicator

  int ome,key;
  MPI_Comm_rank(uorig,&ome);
  for (int i = 0; i < nprocs; i++)
    if (uni2orig[i] == ome) key = i;

  MPI_Comm_split(uorig,0,key,&uworld);
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == nullptr -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
   ------------------------------------------------------------------------- */

void Universe::add_world(char *str)
{
  int n,nper;
  char *ptr;

  if (str == nullptr) {
    n = 1;
    nper = nprocs;
  } else if ((ptr = strchr(str,'x')) != nullptr) {
    *ptr = '\0';
    n = atoi(str);
    nper = atoi(ptr+1);
  } else {
    n = 1;
    nper = atoi(str);
  }

  memory->grow(procs_per_world,nworlds+n,"universe:procs_per_world");
  memory->grow(root_proc,(nworlds+n),"universe:root_proc");

  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0) root_proc[nworlds] = 0;
    else
      root_proc[nworlds] = root_proc[nworlds-1] + procs_per_world[nworlds-1];
    if (me >= root_proc[nworlds]) iworld = nworlds;
    nworlds++;
  }
}

/* ----------------------------------------------------------------------
   check if total procs in all worlds = procs in universe
   ------------------------------------------------------------------------- */

int Universe::consistent()
{
  int n = 0;
  for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
  if (n == nprocs) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
   ------------------------------------------------------------------------- */
double Universe::numeric(const char *file, int line, char *str)
{
  if (!str)
    error->all(file,line,"Expected floating point parameter "
        "in input script or data file");
  int n = strlen(str);
  if (n == 0)
    error->all(file,line,"Expected floating point parameter "
        "in input script or data file");
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    error->all(file,line,"Expected floating point parameter "
        "in input script or data file");
  }
  return atof(str);
}
/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
   ------------------------------------------------------------------------- */
int Universe::inumeric(const char *file, int line, char *str)
{
  if (!str) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  int n = strlen(str);
  if (n == 0) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  }
  return atoi(str);
}
/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
   ------------------------------------------------------------------------- */
bigint Universe::bnumeric(const char *file, int line, char *str)
{
  if (!str) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  int n = strlen(str);
  if (n == 0) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  }
  return ATOBIGINT(str);
}
/* ----------------------------------------------------------------------
   read a tag integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
   ------------------------------------------------------------------------- */
tagint Universe::tnumeric(const char *file, int line, char *str)
{
  if (!str) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  int n = strlen(str);
  if (n == 0) 
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    error->all(file,line,
        "Expected integer parameter in input script or data file");
  }
  return ATOTAGINT(str);
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
   ------------------------------------------------------------------------- */

int Universe::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy,n,"universe:copy");
  strcpy(copy,line);

  char *ptr;
  if ((ptr = strchr(copy,'#'))) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == nullptr) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(nullptr," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}

/* ----------------------------------------------------------------------
   return next prime larger than n
------------------------------------------------------------------------- */

int Universe::next_prime(int n)
{
  int factor;

  int nprime = n+1;
  if (nprime % 2 == 0) nprime++;
  int root = static_cast<int> (sqrt(1.0*n)) + 2;

  while (nprime < MAXSMALLINT) {
    for (factor = 3; factor < root; factor++)
      if (nprime % factor == 0) break;
    if (factor == root) return nprime;
    nprime += 2;
  }

  return MAXSMALLINT;
}

/* ----------------------------------------------------------------------
    compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = nmin to nmax,
     (3) i* = i to nmax, (4) *j = nmin to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Universe::bounds(const char *file, int line, char *str, int nmax, int &nlo, int &nhi, int nmin)
{
  char *ptr = strchr(str,'*');

  if (ptr == nullptr) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = nmin;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < nmin || nhi > nmax || nlo > nhi) 
    error->all(file,line,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = nmin to nmax,
     (3) i* = i to nmax, (4) *j = nmin to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Universe::boundsbig(const char *file, int line, char *str,
                      bigint nmax, bigint &nlo, bigint &nhi, bigint nmin)
{
  char *ptr = strchr(str,'*');

  if (ptr == nullptr) {
    nlo = nhi = ATOBIGINT(str);
  } else if (strlen(str) == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = nmin;
    nhi = ATOBIGINT(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = ATOBIGINT(str);
    nhi = nmax;
  } else {
    nlo = ATOBIGINT(str);
    nhi = ATOBIGINT(ptr+1);
  }

  if (nlo < nmin || nhi > nmax || nlo > nhi)
    error->all(file,line,"Numeric index is out of bounds");
}

/* specialization for the case of just a single string argument */

void Universe::logmesg(CAC *cac, const std::string &mesg)
{
  if (cac->screen) fputs(mesg.c_str(), cac->screen);
  if (cac->logfile) fputs(mesg.c_str(), cac->logfile);
}

