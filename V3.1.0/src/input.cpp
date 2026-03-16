#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/stat.h>
#include "input.h"
#include "style_command.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "comm.h"
#include "comm_brick.h"
#include "comm_tiled.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "force.h"
#include "pair.h"
#include "min.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "neighbor.h"
#include "universe.h"
#include "error.h"
#include "variable.h"
#include "memory.h"
#include "fix.h"
#include "thermo.h"

using namespace CAC_NS;

#define DELTALINE 256
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(CAC *cac, int argc, char **argv) : Pointers(cac)
{
  MPI_Comm_rank(world, &me);

  maxline = maxcopy = maxwork = 0;
  line = copy = work = nullptr;
  narg = maxarg = 0;
  arg = nullptr;


  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = nullptr;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *), "input:infiles");
    infiles[0] = infile;
  } else infiles = nullptr;

  variable = new Variable(cac);

  // fill map with commands listed in style_command.h

  command_map = new CommandCreatorMap();

#define COMMAND_CLASS
#define CommandStyle(key, Class) \
  (*command_map)[#key] = &command_creator<Class>;
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[iarg], "-var") == 0 || strcmp(argv[iarg], "-v") == 0) {
      //int jarg = iarg+3;
      //while (jarg < argc && argv[jarg][0] != '-') jarg++;
      //variable->set(argv[iarg+1], jarg-iarg-2, &argv[iarg+2]);
      //iarg = jarg;
    } else if (strcmp(argv[iarg], "-echo") == 0 ||
        strcmp(argv[iarg], "-e") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 2;
    } else iarg++;
  }

}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  memory->sfree(line);
  memory->sfree(copy);
  memory->sfree(work);
  if (labelstr) delete [] labelstr;
  memory->sfree(arg);
  memory->sfree(infiles);
  delete variable;

  delete command_map;
}

/* ----------------------------------------------------------------------
   one instance per command in style_command.h
   ------------------------------------------------------------------------- */

  template <typename T>
void Input::command_creator(CAC *cac, int narg, char **arg)
{
  T cmd(cac);
  cmd.command(narg, arg);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
   ------------------------------------------------------------------------- */

void Input::file()
{
  int m, n;
  while (1) {

    // read a line from input script
    // n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line

    if (me == 0) {
      m = 0;
      while (1) {
        if (maxline-m < 2) reallocate(line, maxline, 0);

        // end of file reached, so break
        //  n == 0 if nothing read, else n = line with str terminator

        if (fgets(&line[m], maxline-m, infile) == nullptr) {
          if (m) n = strlen(line) + 1;
          else n = 0;
          break;
        }

        // continue if last char read was not a newline
        // could happen if line is very long

        m = strlen(line);
        if (line[m-1] != '\n') continue;

        // continue reading if final printable char is & char
        // or if odd number of triple quotes
        // else break with n = line with str terminator

        m--;
        while (m >= 0 && isspace(line[m])) m--;
        if (m < 0 || line[m] != '&') {
          if (numtriple(line) % 2) {
            m += 2;
            continue;
          }
          line[m+1] = '\0';
          n = m+2;
          break;
        }
      } 
    }
    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n == 0) {
      if (label_active) error->all(FLERR, "Label wasn't found in input script");
      if (me == 0) {
        if (infile != stdin) {
          fclose(infile);
          infile = nullptr;
        }
        nfile--;
      }
      MPI_Bcast(&nfile, 1, MPI_INT, 0, world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    if (n > maxline) reallocate(line, maxline, n);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen, "%s\n", line);
      if (echo_log && logfile) fprintf(logfile, "%s\n", line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();

    if (command == nullptr) continue;

    // test input script

    if (!strcmp(command, "flag")) error->all(FLERR, "\n*****FLAG*****\n");//

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command, "label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char *str = new char[maxline+32];
      sprintf(str, "Unknown command: %s", line);
      error->all(FLERR, str);
    }
  }
}

/* ----------------------------------------------------------------------
   invoke one command in single
   first copy to line, then parse, then execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(const char *single)
{
  int n = strlen(single) + 1;
  if (n > maxline) reallocate(line, maxline, n);
  strcpy(line, single);

  // echo the command unless scanning for label

  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen, "%s\n", line);
    if (echo_log && logfile) fprintf(logfile, "%s\n", line);
  }

  // parse the line
  // if no command, just return nullptr

  parse();
  if (command == nullptr) return nullptr;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command, "label") != 0) return nullptr;

  // execute the command and return its name

  if (execute_command()) {
    char *str = new char[maxline+32];
    sprintf(str, "Unknown command: %s", line);
    error->all(FLERR, str);
  }

  return command;
}


/* ----------------------------------------------------------------------
   parse copy of command line by inserting string terminators
   strip comment = all chars from # on
   replace all $ via variable substitution except within quotes
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double/triple quotes as one arg via nextword()
   ------------------------------------------------------------------------- */

void Input::parse()
{
  // duplicate line into copy string to break into words

  int n = strlen(line) + 1;
  if (n > maxcopy) reallocate(copy, maxcopy, n);
  strcpy(copy, line);

  // strip any # comment by replacing it with 0
  // do not strip from a # inside single/double/triple quotes
  // quoteflag = 1, 2, 3 when encounter first single/double, triple quote
  // quoteflag = 0 when encounter matching single/double, triple quote

  int quoteflag = 0;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && !quoteflag) {
      *ptr = '\0';
      break;
    }
    if (quoteflag == 0) {
      if (strstr(ptr, "\"\"\"") == ptr) {
        quoteflag = 3;
        ptr += 2;
      }
      else if (*ptr == '"') quoteflag = 2;
      else if (*ptr == '\'') quoteflag = 1;
    } else {
      if (quoteflag == 3 && strstr(ptr, "\"\"\"") == ptr) {
        quoteflag = 0;
        ptr += 2;
      }
      else if (quoteflag == 2 && *ptr == '"') quoteflag = 0;
      else if (quoteflag == 1 && *ptr == '\'') quoteflag = 0;
    }
    ptr++;
  }

  // perform $ variable substitution (print changes)
  // except if searching for a label since earlier variable may not be defined

  if (!label_active) substitute(copy, work, maxcopy, maxwork, 1);

  // command = 1st arg in copy string

  char *next;
  command = nextword(copy, &next);
  if (command == nullptr) return;

  // point arg[] at each subsequent arg in copy string
  // nextword() inserts string terminators into copy string to delimit args
  // nextword() treats text between single/double/triple quotes as one arg

  narg = 0;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg, maxarg*sizeof(char *), "input:arg");
    }
    arg[narg] = nextword(ptr, &next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double/triple quotes as one arg
   matching quote must be followed by whitespace char if not end of string
   strip quotes from returned word
   return ptr to start of word or nullptr if no word in string
   also return next = ptr after word
   ------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start, *stop;

  // start = first non-whitespace char

  start = &str[strspn(str, " \t\n\v\f\r")];
  if (*start == '\0') return nullptr;

  // if start is single/double/triple quote:
  // start = first char beyond quote
  // stop = first char of matching quote
  // next = first char beyond matching quote
  // next must be nullptr or whitespace
  // if start is not single/double/triple quote:
  // stop = first whitespace char after start
  // next = char after stop, or stop itself if stop is nullptr

  if (strstr(start, "\"\"\"") == start) {
    stop = strstr(&start[3], "\"\"\"");
    if (!stop) error->all(FLERR, "Unbalanced quotes in input line");
    start += 3;
    *next = stop+3;
    if (**next && !isspace(**next))
      error->all(FLERR, "Input line quote not followed by whitespace");
  } else if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1], *start);
    if (!stop) error->all(FLERR, "Unbalanced quotes in input line");
    start++;
    *next = stop+1;
    if (**next && !isspace(**next))
      error->all(FLERR, "Input line quote not followed by whitespace");
  } else {
    stop = &start[strcspn(start, " \t\n\v\f\r")];
    if (*stop == '\0') *next = stop;
    else *next = stop+1;
  }

  // set stop to nullptr to terminate word

  *stop = '\0';
  return start;
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str using work str2 and return it
   reallocate str/str2 to hold expanded version if necessary & reset max/max2
   print updated string if flag is set and not searching for label
   label_active will be 0 if called from external class
   ------------------------------------------------------------------------- */

void Input::substitute(char *&str, char *&str2, int &max, int &max2, int flag)
{
  // use str2 as scratch space to expand str, then copy back to str
  // reallocate str and str2 as necessary
  // do not replace $ inside single/double/triple quotes
  // var = pts at variable name, ended by nullptr
  //   if $ is followed by '{', trailing '}' becomes nullptr
  //   else $x becomes x followed by nullptr
  // beyond = points to text following variable

  int i, n, paren_count;
  char immediate[256];
  char *var, *value, *beyond;
  int quoteflag = 0;
  char *ptr = str;

  n = strlen(str) + 1;
  if (n > max2) reallocate(str2, max2, n);
  *str2 = '\0';
  char *ptr2 = str2;

  while (*ptr) {

    // variable substitution

    if (*ptr == '$' && !quoteflag) {

      // value = ptr to expanded variable
      // variable name between curly braces, e.g. ${a}

      if (*(ptr+1) == '{') {
        var = ptr+2;
        i = 0;

        while (var[i] != '\0' && var[i] != '}') i++;

        if (var[i] == '\0') error->one(FLERR, "Invalid variable name");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
        value = variable->retrieve(var);

        // immediate variable between parenthesis, e.g. $(1/2)

      } else if (*(ptr+1) == '(') {
        var = ptr+2;
        paren_count = 0;
        i = 0;

        while (var[i] != '\0' && !(var[i] == ')' && paren_count == 0)) {
          switch (var[i]) {
            case '(': paren_count++; break;
            case ')': paren_count--; break;
            default: ;
          }
          i++;
        }

        if (var[i] == '\0') error->one(FLERR, "Invalid immediate variable");
        var[i] = '\0';
        beyond = ptr + strlen(var) + 3;
        sprintf(immediate, "%.20g", variable->compute_equal(var));
        value = immediate;

        // single character variable name, e.g. $a

      } else {
        var = ptr;
        var[0] = var[1];
        var[1] = '\0';
        beyond = ptr + 2;
        value = variable->retrieve(var);
      }

      if (value == nullptr) {
        char errstr[200];
        sprintf(errstr, "Substitution for illegal variable name: %s", ptr+2);
        error->one(FLERR, errstr);
      }

      // check if storage in str2 needs to be expanded
      // re-initialize ptr and ptr2 to the point beyond the variable.

      n = strlen(str2) + strlen(value) + strlen(beyond) + 1;
      if (n > max2) reallocate(str2, max2, n);
      strcat(str2, value);
      ptr2 = str2 + strlen(str2);
      ptr = beyond;

      // output substitution progress if requested

      if (flag && me == 0 && label_active == 0) {
        if (echo_screen && screen) fprintf(screen, "%s%s\n", str2, beyond);
        if (echo_log && logfile) fprintf(logfile, "%s%s\n", str2, beyond);
      }

      continue;
    }

    // quoteflag = 1, 2, 3 when encounter first single/double, triple quote
    // quoteflag = 0 when encounter matching single/double, triple quote
    // copy 2 extra triple quote chars into str2

    if (quoteflag == 0) {
      if (strstr(ptr, "\"\"\"") == ptr) {
        quoteflag = 3;
        *ptr2++ = *ptr++;
        *ptr2++ = *ptr++;
      }
      else if (*ptr == '"') quoteflag = 2;
      else if (*ptr == '\'') quoteflag = 1;
    } else {
      if (quoteflag == 3 && strstr(ptr, "\"\"\"") == ptr) {
        quoteflag = 0;
        *ptr2++ = *ptr++;
        *ptr2++ = *ptr++;
      }
      else if (quoteflag == 2 && *ptr == '"') quoteflag = 0;
      else if (quoteflag == 1 && *ptr == '\'') quoteflag = 0;
    }

    // copy current character into str2

    *ptr2++ = *ptr++;
    *ptr2 = '\0';
  }

  // set length of input str to length of work str2
  // copy work string back to input str

  if (max2 > max) reallocate(str, max, max2);
  strcpy(str, str2);
}

/* ----------------------------------------------------------------------
   expand arg to earg, for arguments with syntax c_ID[*] or f_ID[*]
   fields to consider in input arg range from iarg to narg
   return new expanded # of values, and copy them w/out "*" into earg
   if any expansion occurs, earg is new allocation, must be freed by caller
   if no expansion occurs, earg just points to arg, caller need not free
   ------------------------------------------------------------------------- */

int Input::expand_args(int narg, char **arg, int mode, char **&earg)
{
  int n, iarg, index, nlo, nhi, nmax, expandflag, icompute, ifix;
  char *ptr1, *ptr2, *str;

  ptr1 = nullptr;
  for (iarg = 0; iarg < narg; iarg++) {
    ptr1 = strchr(arg[iarg], '*');
    if (ptr1) break;
  }

  if (!ptr1) {
    earg = arg;
    return narg;
  }

  // maxarg should always end up equal to newarg, so caller can free earg

  int maxarg = narg-iarg;
  earg = (char **) memory->smalloc(maxarg*sizeof(char *), "input:earg");

  int newarg = 0;
  for (iarg = 0; iarg < narg; iarg++) {
    expandflag = 0;

    if (strncmp(arg[iarg], "c_", 2) == 0 ||
        strncmp(arg[iarg], "f_", 2) == 0) {

      ptr1 = strchr(&arg[iarg][2], '[');
      if (ptr1) {
        ptr2 = strchr(ptr1, ']');
        if (ptr2) {
          *ptr2 = '\0';
          if (strchr(ptr1, '*')) {
            if (arg[iarg][0] == 'c') {
              *ptr1 = '\0';
              icompute = modify->find_compute(&arg[iarg][2]);
              *ptr1 = '[';

              // check for global vector/array, peratom array, local array

              if (icompute >= 0) {
                Compute *compute = modify->compute[icompute];

                if (mode == 0 && compute->vector_flag) {
                  nmax = compute->size_vector;
                  expandflag = 1;
                } else if (mode == 1 && compute->array_flag) {
                  nmax = compute->size_array_cols;
                  expandflag = 1;
                } else if (compute->peratom_flag && 
                    compute->size_peratom_cols) {
                  nmax = compute->size_peratom_cols;
                  expandflag = 1;
                } else if (compute->local_flag == 2) {
                  nmax = compute->size_vector_local;
                  expandflag = 1;
                } else if (compute->local_flag == 3) {
                  nmax = compute->size_array_local_cols;
                  expandflag = 1;
                }
              }	      
            } else if (arg[iarg][0] == 'f') {
              *ptr1 = '\0';
              ifix = modify->find_fix(&arg[iarg][2]);
              *ptr1 = '[';

              // check for global vector/array, peratom array, local array

              if (ifix >= 0) {
                if (mode == 0 && modify->fix[ifix]->vector_flag) {
                  nmax = modify->fix[ifix]->size_vector;
                  expandflag = 1;
                } else if (mode == 1 && modify->fix[ifix]->array_flag) {
                  nmax = modify->fix[ifix]->size_array_cols;
                  expandflag = 1;
                } else if (modify->fix[ifix]->peratom_flag && 
                    modify->fix[ifix]->size_peratom_cols) {
                  nmax = modify->fix[ifix]->size_peratom_cols;
                  expandflag = 1;
                } else if (modify->fix[ifix]->local_flag && 
                    modify->fix[ifix]->size_local_cols) {
                  nmax = modify->fix[ifix]->size_local_cols;
                  expandflag = 1;
                }
              }
            }
          }
          *ptr2 = ']';
        }
      }
    }

    if (expandflag) {
      *ptr2 = '\0';
      universe->bounds(FLERR, ptr1+1, nmax, nlo, nhi);
      *ptr2 = ']';
      if (newarg+nhi-nlo+1 > maxarg) {
        maxarg += nhi-nlo+1;
        earg = (char **) 
          memory->srealloc(earg, maxarg*sizeof(char *), "input:earg");
      }
      for (index = nlo; index <= nhi; index++) {
        n = strlen(arg[iarg]) + 16;   // 16 = space for large inserted integer
        str = earg[newarg] = new char[n];
        strncpy(str, arg[iarg], ptr1+1-arg[iarg]);
        sprintf(&str[ptr1+1-arg[iarg]], "%d", index);
        strcat(str, ptr2);
        newarg++;
      }

    } else {
      if (newarg == maxarg) {
        maxarg++;
        earg = (char **) 
          memory->srealloc(earg, maxarg*sizeof(char *), "input:earg");
      }
      n = strlen(arg[iarg]) + 1;
      earg[newarg] = new char[n];
      strcpy(earg[newarg], arg[iarg]);
      newarg++;
    }
  }

  return newarg;
}

/* ----------------------------------------------------------------------
   return number of triple quotes in line
   ------------------------------------------------------------------------- */

int Input::numtriple(char *line)
{
  int count = 0;
  char *ptr = line;
  while ((ptr = strstr(ptr, "\"\"\""))) {
    ptr += 3;
    count++;
  }
  return count;
}


/* ----------------------------------------------------------------------
   rellocate a string
   if n > 0: set max >= n in increments of DELTALINE
   if n = 0: just increment max by DELTALINE
   ------------------------------------------------------------------------- */

void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;

  str = (char *) memory->srealloc(str, max*sizeof(char), "input:str");
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
   ------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;
  if (!strcmp(command, "clear")) clear();
  else if (!strcmp(command, "echo")) echo();
  else if (!strcmp(command, "if")) ifthenelse();
  else if (!strcmp(command, "jump")) jump();
  else if (!strcmp(command, "label")) label();
  else if (!strcmp(command, "log")) log();
  else if (!strcmp(command, "next")) next_command();
  else if (!strcmp(command, "print")) print();
  else if (!strcmp(command, "quit")) quit();
  else if (!strcmp(command, "variable")) variable_command();

  else if (!strcmp(command, "add_atoms")) add_atoms();
  else if (!strcmp(command, "add_etype")) add_etype();
  else if (!strcmp(command, "atom_modify")) atom_modify();
  else if (!strcmp(command, "atom_style")) atom_style();
  else if (!strcmp(command, "boundary")) boundary();
  else if (!strcmp(command, "charge")) charge();
  else if (!strcmp(command, "comm_modify")) comm_modify();
  else if (!strcmp(command, "comm_style")) comm_style();
  else if (!strcmp(command, "compress_etype")) compress_etype();
  else if (!strcmp(command, "compute")) compute();
  else if (!strcmp(command, "dielectric")) dielectric();
  else if (!strcmp(command, "dimension")) dimension();
  else if (!strcmp(command, "dump")) dump();
  else if (!strcmp(command, "dump_modify")) dump_modify();
  else if (!strcmp(command, "element_modify")) element_modify();
  else if (!strcmp(command, "element_style")) element_style(); 
  else if (!strcmp(command, "fix")) fix();
  else if (!strcmp(command, "group")) group_command();
  else if (!strcmp(command, "lattice")) lattice();
  else if (!strcmp(command, "mass")) mass();
  else if (!strcmp(command, "min_modify")) min_modify();
  else if (!strcmp(command, "min_style")) min_style();
  else if (!strcmp(command, "neighbor")) neighbor_command();
  else if (!strcmp(command, "neigh_modify")) neigh_modify();
  else if (!strcmp(command, "newton")) newton();
  else if (!strcmp(command, "region")) region();
  else if (!strcmp(command, "reset_timestep")) reset_timestep();
  else if (!strcmp(command, "pair_style")) pair_style();
  else if (!strcmp(command, "pair_coeff")) pair_coeff();
  else if (!strcmp(command, "processors")) processors();
  else if (!strcmp(command, "thermo")) thermo();
  else if (!strcmp(command, "thermo_style")) thermo_style();
  else if (!strcmp(command, "thermo_modify")) thermo_modify();
  else if (!strcmp(command, "timestep")) timestep();
  else if (!strcmp(command, "uncompute")) uncompute();
  else if (!strcmp(command, "undump")) undump();
  else if (!strcmp(command, "unfix")) unfix();
  else if (!strcmp(command, "units")) units();

  // debug commands
  
  else if (!strcmp(command, "debug")) debug();
  else if (!strcmp(command, "debug_comm")) debug_comm();
  //else if (!strcmp(command, "debug_dump_force_integration")) debug_dump_force_integration();
  //else if (!strcmp(command, "debug_element_intg_off")) debug_element_intg_off();
  //else if (!strcmp(command, "debug_test_adaptive_mesh")) debug_test_adaptive_mesh();

  else flag = 0;

  // return if command was listed above

  if (flag) return 0;

  // invoke commands added via style_command.h

  if (command_map->find(command) != command_map->end()) {
    CommandCreator command_creator = (*command_map)[command];
    command_creator(cac, narg, arg);
    return 0;
  }

  // unrecognized command

  return -1;
}
/* ---------------------------------------------------------------------- */

void Input::region()
{
  domain->add_region(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::reset_timestep()
{
  update->reset_timestep(narg, arg);
}

/*------------------------------------------------------------------------*/

void Input::group_command()
{
  group->assign(narg, arg);
}	

/*------------------------------------------------------------------------*/

void Input::element_style()
{
  if (narg < 1) error->all(FLERR, "Illegal element_style command");
  if (domain->box_exist)
    error->all(FLERR, "Element_style command after simulation box is defined");
  element->create_evec(arg[0], narg-1, &arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::atom_modify()
{
  atom->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::atom_style()
{
  if (narg < 1) error->all(FLERR, "Illegal atom_style command");
  if (domain->box_exist)
    error->all(FLERR, "Atom_style command after simulation box is defined");
  atom->create_avec(arg[0], narg-1, &arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all(FLERR, "Illegal clear command");
  cac->destroy();
  cac->create();
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
  if (narg != 1) error->all(FLERR, "Illegal units command");
  if (domain->box_exist)
    error->all(FLERR, "Units command after simulation box is defined");
  update->set_units(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (narg != 1) error->all(FLERR, "Illegal dimension command");
  if (domain->box_exist)
    error->all(FLERR, "Dimension command after simulation box is defined");
  domain->dimension = universe->inumeric(FLERR, arg[0]);
  if (domain->dimension != 2 && domain->dimension != 3)
    error->all(FLERR, "Illegal dimension command");

  // must reset default extra_dof of all computes
  // since some were created before dimension command is encountered

  for (int i = 0; i < modify->ncompute; i++)
    modify->compute[i]->reset_extra_dof();
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  if (domain->box_exist)
    error->all(FLERR, "Boundary command after simulation box is defined");
  domain->set_boundary(narg, arg, 0);
}

/* ---------------------------------------------------------------------- */

void Input::newton()
{
  int newton_pair = 1, newton_bond = 1;

  if (narg == 1) {
    if (strcmp(arg[0], "off") == 0) newton_pair = newton_bond = 0;
    else if (strcmp(arg[0], "on") == 0) newton_pair = newton_bond = 1;
    else error->all(FLERR, "Illegal newton command");
  } else if (narg == 2) {
    if (strcmp(arg[0], "off") == 0) newton_pair = 0;
    else if (strcmp(arg[0], "on") == 0) newton_pair= 1;
    else error->all(FLERR, "Illegal newton command");
    if (strcmp(arg[1], "off") == 0) newton_bond = 0;
    else if (strcmp(arg[1], "on") == 0) newton_bond = 1;
    else error->all(FLERR, "Illegal newton command");
  } else error->all(FLERR, "Illegal newton command");

  force->newton_pair = newton_pair;

  if (domain->box_exist && (newton_bond != force->newton_bond))
    error->all(FLERR, "Newton bond change after simulation box is defined");
  force->newton_bond = newton_bond;

  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;

}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
  neighbor->set(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
  neighbor->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::element_modify()
{
  element->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::add_etype()
{
  element->add_etype(narg, arg);
}


/* ----------------------------------------------------------------------
   if old pair style exists and new style is same, just change settings
   else create new pair class
   ------------------------------------------------------------------------- */

void Input::pair_style()
{
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");
  if (force->pair) {
    int match = 0;
    if (strcmp(arg[0], force->pair_style) == 0) match = 1;
    if (!match && cac->suffix_enable) {
      char estyle[256];
      if (cac->suffix) {
        sprintf(estyle, "%s/%s", arg[0], cac->suffix);
        if (strcmp(estyle, force->pair_style) == 0) match = 1;
      }
      if (cac->suffix2) {
        sprintf(estyle, "%s/%s", arg[0], cac->suffix2);
        if (strcmp(estyle, force->pair_style) == 0) match = 1;
      }
    }
    if (match) {
      force->pair->settings(narg-1, &arg[1]);
      return;
    }
  }

  force->create_pair(arg[0], 1);
  if (force->pair) force->pair->settings(narg-1, &arg[1]);  
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Pair_coeff command before simulation box is defined");
  if (force->pair == nullptr)
    error->all(FLERR, "Pair_coeff command before pair_style is defined");
  force->pair->coeff(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
  modify->add_fix(narg, arg, 1);

}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
  output->set_thermo(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_modify()
{
  output->thermo->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
  output->create_thermo(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  output->add_dump(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  output->modify_dump(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all(FLERR, "Illegal timestep command");
  update->dt = universe->numeric(FLERR, arg[0]);
}

/*-----------------------------------------------------------------------*/

void Input::processors()
{
  if (domain->box_exist)
    error->all(FLERR, "Processors command after simulation box is defined");
  comm->set_processors(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::comm_modify()
{
  comm->modify_params(narg, arg);
}


/* ---------------------------------------------------------------------- */

void Input::comm_style()
{

  if (narg < 1) error->all(FLERR, "Illegal comm_style command");
  if (strcmp(arg[0], "brick") == 0) {
    if (comm->style == 0) return;
    Comm *oldcomm = comm;
    comm = new CommBrick(cac, oldcomm);
    delete oldcomm;
  } else if (strcmp(arg[0], "tiled") == 0) {
    if (comm->style == 1) return;
    Comm *oldcomm = comm;
    comm = new CommTiled(cac, oldcomm);
    delete oldcomm;
   } else error->all(FLERR, "Illegal comm_style command");
}


/* ---------------------------------------------------------------------- */

void Input::compute()
{
  modify->add_compute(narg, arg);
}


/* ---------------------------------------------------------------------- */

void Input::uncompute()
{
  if (narg != 1) error->all(FLERR, "Illegal uncompute command");
  modify->delete_compute(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
  if (narg != 1) error->all(FLERR, "Illegal undump command");
  output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
  //if (narg != 1) error->all(FLERR, "Illegal unfix command");

  // can delete multiple fix in one line

  for (int i = 0; i < narg; i++) 
    modify->delete_fix(arg[i]);
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg, arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg < 1) error->all(FLERR, "Illegal print command");

  // copy 1st arg back into line (copy is being used)
  // check maxline since arg[0] could have been exanded by variables
  // substitute for $ variables (no printing) and print arg

  int n = strlen(arg[0]) + 1;
  if (n > maxline) reallocate(line, maxline, n);
  strcpy(line, arg[0]);
  substitute(line, work, maxline, maxwork, 0);

  // parse optional args

  FILE *fp = nullptr;
  int screenflag = 1;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "file") == 0 || strcmp(arg[iarg], "append") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal print command");
      if (me == 0) {
        if (fp != nullptr) fclose(fp);
        if (strcmp(arg[iarg], "file") == 0) fp = fopen(arg[iarg+1], "w");
        else fp = fopen(arg[iarg+1], "a");
        if (fp == nullptr) {
          char str[128];
          sprintf(str, "Cannot open print file %s", arg[iarg+1]);
          error->one(FLERR, str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "screen") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal print command");
      if (strcmp(arg[iarg+1], "yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) screenflag = 0;
      else error->all(FLERR, "Illegal print command");
      iarg += 2;
    } else error->all(FLERR, "Illegal print command");
  }

  if (me == 0) {
    if (screenflag && screen) fprintf(screen, "%s\n", line);
    if (screenflag && logfile) fprintf(logfile, "%s\n", line);
    if (fp) {
      fprintf(fp, "%s\n", line);
      fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all(FLERR, "Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (strcmp(arg[0], "SELF") == 0) rewind(infile);
    else {
      if (infile && infile != stdin) fclose(infile);
      infile = fopen(arg[0], "r");
      if (infile == nullptr) {
        char str[128];
        sprintf(str, "Cannot open input script %s", arg[0]);
        error->one(FLERR, str);
      }
      infiles[nfile-1] = infile;
    }
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr, arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all(FLERR, "Illegal label command");
  if (label_active && strcmp(labelstr, arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::add_atoms()
{
  element->evec->add_atoms(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::lattice()
{
  domain->set_lattice(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::mass()
{
  if (narg != 2) error->all(FLERR, "Illegal mass command");
  if (domain->box_exist == 0)
    error->all(FLERR, "Mass command before simulation box is defined");
  atom->set_mass(FLERR, narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::charge()
{
  if (narg != 2) error->all(FLERR, "Illegal charge command");
  if (domain->box_exist == 0)
    error->all(FLERR, "Charge command before simulation box is defined");
  atom->set_charge(FLERR, narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_modify()
{
  update->minimize->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_style()
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Min_style command before simulation box is defined");
  update->create_minimize(narg, arg);
}


/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all(FLERR, "Illegal echo command");

  if (strcmp(arg[0], "none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0], "screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0], "log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0], "both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all(FLERR, "Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::ifthenelse()
{
  if (narg < 3) error->all(FLERR, "Illegal if command");

  // substitute for variables in Boolean expression for "if"
  // in case expression was enclosed in quotes
  // must substitute on copy of arg else will step on subsequent args

  int n = strlen(arg[0]) + 1;
  if (n > maxline) reallocate(line, maxline, n);
  strcpy(line, arg[0]);
  substitute(line, work, maxline, maxwork, 0);

  // evaluate Boolean expression for "if"

  double btest = variable->evaluate_boolean(line);

  // bound "then" commands

  if (strcmp(arg[1], "then") != 0) error->all(FLERR, "Illegal if command");

  int first = 2;
  int iarg = first;
  while (iarg < narg &&
         (strcmp(arg[iarg], "elif") != 0 && strcmp(arg[iarg], "else") != 0))
    iarg++;
  int last = iarg-1;

  // execute "then" commands
  // make copies of all arg string commands
  // required because re-parsing a command via one() will wipe out args

  if (btest != 0.0) {
    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR, "Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR, "Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands], arg[i]);
      ncommands++;
    }

    for (int i = 0; i < ncommands; i++) {
      one(commands[i]);
      delete [] commands[i];
    }
    delete [] commands;

    return;
  }

  // done if no "elif" or "else"

  if (iarg == narg) return;

  // check "elif" or "else" until find commands to execute
  // substitute for variables and evaluate Boolean expression for "elif"
  // must substitute on copy of arg else will step on subsequent args
  // bound and execute "elif" or "else" commands

  while (iarg != narg) {
    if (iarg+2 > narg) error->all(FLERR, "Illegal if command");
    if (strcmp(arg[iarg], "elif") == 0) {
      n = strlen(arg[iarg+1]) + 1;
      if (n > maxline) reallocate(line, maxline, n);
      strcpy(line, arg[iarg+1]);
      substitute(line, work, maxline, maxwork, 0);
      btest = variable->evaluate_boolean(line);
      first = iarg+2;
    } else {
      btest = 1.0;
      first = iarg+1;
    }

    iarg = first;
    while (iarg < narg &&
           (strcmp(arg[iarg], "elif") != 0 && strcmp(arg[iarg], "else") != 0))
      iarg++;
    last = iarg-1;

    if (btest == 0.0) continue;

    int ncommands = last-first + 1;
    if (ncommands <= 0) error->all(FLERR, "Illegal if command");

    char **commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      int n = strlen(arg[i]) + 1;
      if (n == 1) error->all(FLERR, "Illegal if command");
      commands[ncommands] = new char[n];
      strcpy(commands[ncommands], arg[i]);
      ncommands++;
    }

    // execute the list of commands

    for (int i = 0; i < ncommands; i++) {
      one(commands[i]);
      delete [] commands[i];
    }
    delete [] commands;

    return;
  }
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg > 2) error->all(FLERR, "Illegal log command");

  int appendflag = 0;
  if (narg == 2) {
    if (strcmp(arg[1], "append") == 0) appendflag = 1;
    else error->all(FLERR, "Illegal log command");
  }

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0], "none") == 0) logfile = nullptr;
    else {
      if (appendflag) logfile = fopen(arg[0], "a");
      else logfile = fopen(arg[0], "w");

      if (logfile == nullptr) {
        char str[128];
        snprintf(str, 128, "Cannot open logfile %s", arg[0]);
        error->one(FLERR, str);
      }
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::quit()
{

  if (narg == 0) error->done(0); // 1 would be fully backwards compatible
  if (narg == 1) error->done(universe->inumeric(FLERR, arg[0]));
  error->all(FLERR, "Illegal quit command");
}

/* ---------------------------------------------------------------------- */

void Input::compress_etype()
{
  if (narg) error->all(FLERR, "Invalid compress_etype command");
  element->compress_etype();
}

/* ---------------------------------------------------------------------- */

void Input::dielectric()
{
  if (narg != 1) error->all(FLERR, "Invalid dielectric command");
  force->dielectric = universe->numeric(FLERR, arg[0]);
}
/* ---------------------------------------------------------------------- */

void Input::debug_comm()
{
  if (narg != 0) error->all(FLERR, "Invalid debug_comm command");
  comm->debug = 1;
}

/* ---------------------------------------------------------------------- */

void Input::debug()
{
//  int a, b;
//  a = element->count_element_clusters();
//  MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_SUM, world);
//  printf("me = %d ncluster = %d total = %d %d\n", comm->me, element->count_element_clusters(), b, element->nelements/8);
}

/* ---------------------------------------------------------------------- */

  /*
void Input::debug_dump_force_integration()
{
  if (narg != 3) error->all(FLERR, "Invalid debug_dump_force_integration command");
  force->pair->debug_mode = 1;
  force->pair->debug_nevery = universe->inumeric(FLERR, arg[1]);
  int igroup = group->find(arg[0]);
  force->pair->debug_groupbit = group->bitmask[igroup];
  strcpy(arg[2], force->pair->debug_file);
}

void Input::debug_element_intg_off()
{
  element->debug_element_intg_off(narg, arg);
}

void Input::debug_test_adaptive_mesh()
{
  if (element->nmax == 0) element->evec->grow(1);

  int nsplit = 2;
  int **split_list;
  memory->create(split_list, nsplit, 3, "input:split_list");
  for (int i = 0; i < nsplit; i++) {
    split_list[i][0] = 0;
  }
  split_list[0][1] = 7;
  split_list[0][2] = 0;

  split_list[1][1] = 3;
  split_list[1][2] = 2;

  element->nlocal = 1;
  atom->nlocal = 0;
  element->nelements = comm->nprocs;
  element->x[0][0] = 0.5;
  element->x[0][1] = 0.5;
  element->x[0][2] = 0.5;
  element->nodex[0][0][0] = 0;
  element->nodex[0][0][1] = 0;
  element->nodex[0][0][2] = 0;
  element->nodex[0][1][0] = 1;
  element->nodex[0][1][1] = 0;
  element->nodex[0][1][2] = 0;
  element->nodex[0][2][0] = 1;
  element->nodex[0][2][1] = 1;
  element->nodex[0][2][2] = 0;
  element->nodex[0][3][0] = 0;
  element->nodex[0][3][1] = 1;
  element->nodex[0][3][2] = 0;
  element->nodex[0][4][0] = 0;
  element->nodex[0][4][1] = 0;
  element->nodex[0][4][2] = 1;
  element->nodex[0][5][0] = 1;
  element->nodex[0][5][1] = 0;
  element->nodex[0][5][2] = 1;
  element->nodex[0][6][0] = 1;
  element->nodex[0][6][1] = 1;
  element->nodex[0][6][2] = 1;
  element->nodex[0][7][0] = 0;
  element->nodex[0][7][1] = 1;
  element->nodex[0][7][2] = 1;
  element->etype[0] = 1;
  element->ctype[0] = 1;
  element->tag[0] = 1;
  element->evec->split_element(split_list, nsplit);
  int sum = atom->nlocal;
  for (int i = 0; i < element->nlocal; i++)
    sum += element->nintpl[element->etype[i]];

  printf("natom = %d sum = %d\n", atom->nlocal, sum);
//  if (sum != element->nintpl[1]) {
//    error->all(FLERR, "TEST"); 
//  }
  atom->tag_extend();
  element->tag_extend();

  int n = atom->nlocal;
  for (int i = 0; i < element->nlocal;i++){
    n += element->npe[element->etype[i]]+1;
  }
  FILE *fp = fopen("dump_test.atom", "w");
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, "%d \n", n);
  fprintf(fp, "ITEM: BOX BOUNDS p p p\n");
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
  fprintf(fp, "ITEM: ATOMS id type x y z eid inode etype\n");
  int count = atom->nlocal+1;
  for (int i = 0; i < count-1; i++)
    fprintf(fp, "%d %d %g %g %g %d %d %d\n", atom->tag[i], 1
        , atom->x[i][0]
        , atom->x[i][1]
        , atom->x[i][2]
        , -1, -1, -1 
        );
  for (int i = 0; i < element->nlocal; i++) {
    fprintf(fp, "%d %d %g %g %g %d %d %d\n", count++, 2
        , element->x[i][0]
        , element->x[i][1]
        , element->x[i][2], 
        element->tag[i], -1, element->etype[i]);
    for (int j = 0; j < element->npe[element->etype[i]]; j++) 
      fprintf(fp, "%d %d %g %g %g %d %d %d\n", count++, 3
          , element->nodex[i][j][0]
          , element->nodex[i][j][1]
          , element->nodex[i][j][2], 
          element->tag[i], j, element->etype[i]);
  }
  fclose(fp);
  
}
*/

