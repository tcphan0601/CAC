#ifndef CAC_VARIABLE_H
#define CAC_VARIABLE_H

#include <stdlib.h>
#include "pointers.h"

namespace CAC_NS {

class Variable : protected Pointers {
 friend class Info;
 public:
  Variable(class CAC *);
  ~Variable();
  void set(int, char **);
  void set(char *, int, char **);
  int set_string(char *, char *);
  int next(int, char **);

  int find(char *);
  void set_atom_arrays(int);
  void set_elem_arrays(int);
  void pass_on_atom_arrays(int, int, int) {}
  void pass_on_elem_arrays(int, int, int *) {}

  //void python_command(int, char **);

  int equalstyle(int);
  int atomstyle(int);
  int vectorstyle(int);
  //char *pythonstyle(char *, char *);
  int internalstyle(int);

  char *retrieve(char *);
  double compute_equal(int);
  double compute_equal(char *);
  void compute_atom(int, int, double *, int, int);
  int compute_vector(int, double **);
  void internal_set(int, double);

  tagint int_between_brackets(char *&, int);
  double evaluate_boolean(char *);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables following lists can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  //class VarReader **reader;   // variable that reads from file
  char ***data;            // str value of each variable's values
  double *dvalue;          // single numeric value for internal variables

  struct VecVar {
    int n,nmax;
    bigint currentstep;
    double *values;
  };
  VecVar *vecs;

  int *eval_in_progress;       // flag if evaluation of variable is in progress
  int treetype;                // ATOM or VECTOR flag for formula evaluation

  class RanMars *randomequal;   // random number generator for equal-style vars
  class RanMars *randomatom;    // random number generator for atom-style vars

  int precedence[18];      // precedence level of math operators
                           // set length to include up to XOR in enum

  struct Tree {            // parse tree for atom-style or vector-style vars
    double value;          // single scalar
    double *array;         // per-atom or per-type list of doubles
    int *iarray;           // per-atom list of ints
    bigint *barray;        // per-atom list of bigints
    int type;              // operation, see enum{} in variable.cpp
    int nvector;           // length of array for vector-style variable
    int nstride;           // stride between atoms if array is a 2d array
    int selfalloc;         // 1 if array is allocated here, else 0
    int ivalue1,ivalue2;   // extra values for needed for gmask,rmask,grmask
    int nextra;            // # of additional args beyond first 2
    Tree *first,*second;   // ptrs further down tree for first 2 args
    Tree **extra;          // ptrs further down tree for nextra args
  };

  //int compute_python(int);
  void remove(int);
  void grow();
  void copy(int, char **, char **);
  double evaluate(char *, Tree **, int);
  double collapse_tree(Tree *);
  double eval_tree(Tree *, int);
  int size_tree_vector(Tree *);
  int compare_tree_vector(int, int);
  void free_tree(Tree *);
  int find_matching_paren(char *, int, char *&, int);
  int math_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  int group_function(char *, char *, Tree **, Tree **, int &, double *, int &, int);
  int region_function(char *, int);
  int special_function(char *, char *, Tree **, Tree **,
                       int &, double *, int &, int);
  void peratom2global(int, char *, double *, int, tagint, int,
                      Tree **, Tree **, int &, double *, int &);
  int is_atom_vector(char *);
  void atom_vector(char *, Tree **, Tree **, int &);
  int is_constant(char *);
  double constant(char *);
  int parse_args(char *, char **);
  char *find_next_comma(char *);
  void print_var_error(const char *, int, const char *, int);
  void print_tree(Tree *, int);
};
/*
class VarReader : protected Pointers {
 public:
  class FixStore *fixstore;
  char *id_fix;

  VarReader(class CAC *, char *, char *, int);
  ~VarReader();
  int read_scalar(char *);
  int read_peratom();

 private:
  int me,style;
  FILE *fp;
  char *buffer;
};
*/
}

#endif
