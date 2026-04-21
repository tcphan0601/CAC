#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "force.h"
#include "style_pair.h"
#include "comm.h"
#include "pair.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

#define MAXLINE 1024

/* ------------------------------------------------------------- */

Force::Force(CAC *cac) : Pointers(cac)
{
  newton = newton_pair = newton_bond = 1;

  special_lj[0] = special_coul[0] = 1.0;
  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;
  special_extra = 0;

  dielectric = 1.0;

  pair = nullptr;

  char *str = (char *) "none";
  int n = strlen(str) + 1;
  pair_style = new char[n];
  strcpy(pair_style, str);

  pair_map = new PairCreatorMap();

#define PAIR_CLASS
#define PairStyle(key, Class) \
  (*pair_map)[#key] = &pair_creator<Class>;
#include "style_pair.h"
#undef PairStyle
#undef PAIR_CLASS

}

/*  ----------------------------------------------------------------------  */

Force::~Force()
{
  delete [] pair_style;

  if (pair) delete pair;
  pair = nullptr;

  delete pair_map;
}


/*  ----------------------------------------------------------------------
 one instance per pair style in style_pair.h
 -------------------------------------------------------------------------  */

template <typename T>
Pair *Force::pair_creator(CAC *cac)
{
  return new T(cac);
}

/*  ----------------------------------------------------------------------  */

void Force::setup()
{
  if (pair) pair->setup();
}


/*  ----------------------------------------------------------------------
    create a pair style, called from input script or restart file
-------------------------------------------------------------------------  */

void Force::create_pair(const char *style, int trysuffix)
{
  delete [] pair_style;
  if (pair) delete pair;

  int sflag;
  pair = new_pair(style, trysuffix, sflag);
  store_style(pair_style, style, sflag);
}

/*  ----------------------------------------------------------------------
   generate a pair class
   if trysuffix = 1, try first with suffix1/2 appended
   return sflag = 0 for no suffix added, 1 or 2 for suffix1/2 added
-------------------------------------------------------------------------  */

Pair *Force::new_pair(const char *style, int trysuffix, int &sflag)
{
  if (trysuffix && cac->suffix_enable) {
    if (cac->suffix) {
      sflag = 1;
      char estyle[256];
      sprintf(estyle, "%s/%s", style, cac->suffix);
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator pair_creator = (*pair_map)[estyle];
        return pair_creator(cac);
      }
    }
    if (cac->suffix2) {
      sflag = 2;
      char estyle[256];
      sprintf(estyle, "%s/%s", style, cac->suffix2);
      if (pair_map->find(estyle) != pair_map->end()) {
        PairCreator pair_creator = (*pair_map)[estyle];
        return pair_creator(cac);
      }
    }
  }

  sflag = 0;
  if (strcmp(style, "none") == 0) return nullptr;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator pair_creator = (*pair_map)[style];
    return pair_creator(cac);
  }

  error->all(FLERR, "Unknown pair style");
  return nullptr;

}

/*  ----------------------------------------------------------------------
   store style name in str allocated here
   if sflag = 0, no suffix
   if sflag = 1/2, append suffix or suffix2 to style
-------------------------------------------------------------------------  */

void Force::store_style(char *&str, const char *style, int sflag)
{
  if (sflag) {
    char estyle[256];
    if (sflag == 1) sprintf(estyle, "%s/%s", style, cac->suffix);
    else sprintf(estyle, "%s/%s", style, cac->suffix2);
    int n = strlen(estyle) + 1;
    str = new char[n];
    strcpy(str, estyle);
  } else {
    int n = strlen(style) + 1;
    str = new char[n];
    strcpy(str, style);
  }
}

/*  ----------------------------------------------------------------------  */

void Force::init()
{
  qqrd2e = qqr2e/dielectric;
  if (pair) pair->init();             // so g_ewald is defined
}

/*  ----------------------------------------------------------------------  */

bigint Force::memory_usage()
{
  bigint bytes = 0;
  if (pair) bytes += static_cast<bigint> (pair->memory_usage());
  return bytes;
}

/* --------------------------------------------------------------------------
  open a potential file as specified by name
  if fails, search in dir specified by env variable CAC_POTENTIALS
----------------------------------------------------------------------------- */

FILE *Force::open_potential(const char *name)
{
  FILE *fp;

  if (name == nullptr) return nullptr;

  // attempt to open file directly
  // if successful, return ptr

  fp = fopen(name, "r");
  if (fp) {
    if (comm->me == 0) potential_date(fp, name);
    rewind(fp);
    return fp;
  }

  // try the enviroment variable directory

  const char *path = getenv("CAC_POTENTIALS");
  if (path == nullptr) return nullptr;

  const char *pot = potential_name(name);
  if (pot == nullptr) return nullptr;

  size_t len1 = strlen(path);
  size_t len2 = strlen(pot);
  char *newpath = new char[len1+len2+2];

  strcpy(newpath, path);
  newpath[len1] = '/';
  newpath[len1+1] = 0;
  strcat(newpath, pot);

  fp = fopen(newpath, "r");
  if (fp) {
    if (comm->me == 0) potential_date(fp, name);
    rewind(fp);
  }

  delete [] newpath;
  return fp;
}

/* --------------------------------------------------------------------------------
  strip off leading part of path, return just the filename
---------------------------------------------------------------------------------- */

const char *Force::potential_name(const char *path)
{
  const char *pot;

  if (path == nullptr) return nullptr;

  for (pot = path; *path != '\0'; ++path) {
  if (*path == '/') pot = path + 1;
  }
  return pot;
}

/* ----------------------------------------------------------------------------------
  read first line of potential file
  if has DATE field, print following word
----------------------------------------------------------------------------------- */

void Force::potential_date(FILE *fp, const char *name)
{
  char line[MAXLINE];
  char *ptr = fgets(line, MAXLINE, fp);
  if (ptr == nullptr) return;

  char *word;
  word = strtok(line, " \t\n\r\f");
  while(word) {
    if (strcmp(word, "DATE:") == 0) {
      word = strtok(nullptr, " \t\n\r\f");
      if (word == nullptr) return;
      if (screen)
	fprintf(screen, "Reading potential file %s with DATE: %s\n", name, word);
      if (logfile)
	fprintf(logfile, "Reading potential file %s with DATE: %s\n", name, word);
      return;
    }
    word = strtok(nullptr, " \t\n\r\f");
  }
}

/*  ----------------------------------------------------------------------
   return style name of Pair class that matches Pair ptr
   called by Neighbor::print_neigh_info()
   return nullptr if no match
-------------------------------------------------------------------------  */

char *Force::pair_match_ptr(Pair *ptr)
{
  if (ptr == pair) return pair_style;

  if (strstr(pair_style, "hybrid")) {
    PairHybrid *hybrid = (PairHybrid *) pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (ptr == hybrid->styles[i]) return hybrid->keywords[i];
  }

  return nullptr;
}
