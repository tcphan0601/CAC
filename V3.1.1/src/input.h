#ifndef CAC_INPUT_H
#define CAC_INPUT_H

#include <stdio.h>
#include "pointers.h"
#include <map>
#include <string>

namespace CAC_NS {

class Input : protected Pointers {
  friend class Error;
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command

  class Variable *variable;    // defined variables

  Input(class CAC *, int, char **);
  ~Input();
  void file();
  char *one(const char *);       // process a single command

  void substitute(char *&, char *&, int &, int &, int);
                                 // substitute for variables in a string
  int expand_args(int, char **, int, char **&);  // expand args due to wildcard

  typedef void ( * CommandCreator)(CAC *, int, char **);
   typedef std::map<std::string, CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;

 protected:
  template <typename T> static void command_creator(CAC *, int, char **);

 private:

  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line, *copy, *work;      // input line & copy and work string
  int maxline, maxcopy, maxwork; // max lengths of char strings
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile, maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise

  FILE **infiles;              // list of open input files


  void parse();                          // parse an input text line
  char *nextword(char *, char **);       // find next word in string with quotes
  int numtriple(char *);                 // count number of triple quotes
  void reallocate(char *&, int &, int);  // reallocate a char string
  int execute_command();                 // execute a single command

  // input script commands
  
  void clear(); 
  void echo();
  void ifthenelse();
  void jump();
  void label();
  void log();
  void next_command();
  void print();
  void quit();
  void variable_command();

  // CAC commands
  
  void add_atoms();
  void add_etype();
  void atom_modify();
  void atom_style();
  void boundary();
  void charge();
  void comm_modify();
  void comm_style();
  void compress_etype();
  void compute();
  void dielectric();
  void dimension();
  void dump();
  void dump_modify();
  void element_modify();
  void element_style();
  void fix();
  void group_command();
  void lattice();
  void mass();
  void min_modify();
  void min_style();
  void neighbor_command();
  void neigh_modify();
  void newton(); 
  void region();
  void reset_timestep();
  void pair_style();
  void pair_coeff();
  void processors();
  void thermo();
  void thermo_style();
  void thermo_modify();
  void timestep();
  void uncompute();
  void undump();
  void unfix();
  void units();

  // debug commands
  
  void debug();
  void debug_comm();
  //void debug_dump_force_integration();
  //void debug_element_gcell_off();
  //void debug_test_adaptive_mesh();
};

}

#endif
