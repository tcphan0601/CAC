/*-------------------------------------------------------------------------
 CAC Class:Similar to lammps.h in lammps
Hao Chen 
-------------------------------------------------------------------------- */

#ifndef CAC_CAC_H
#define CAC_CAC_H

#include <stdio.h>
#include <mpi.h>

namespace CAC_NS {

class CAC {
 public:
  // ptrs to fundamental CAC classes 
  class Memory *memory;     // memory allocation functions
  class Error *error;       // error handling
  class Universe *universe; // universe of processors
  class Input *input;       // input script processing
  // ptrs to top-level CAC-specific classes
  class Element *element;   // element-based quantities
  class Atom *atom;         // atom-based quantities
  class Update *update;     // integrators/minimizers
  class Neighbor *neighbor; // neighbor lists
  class Comm *comm;         // inter-processor communication
  class Domain *domain;     // simulation box
  class Force *force;       // inter-particle forces
  class Modify *modify;     // fixes and computes
  class Group *group;       // groups of atoms
  class Output *output;     // thermo/dump/restart
  class Timer *timer;       // CPU timing info

  MPI_Comm world;           // MPI communicator
  FILE *infile;             // infile
  FILE *screen;             // screen output
  FILE *logfile;            // logfile

  double initclock;         // wall clock at instantiation

  char *suffix,*suffix2;    // suffixes to add to input script style names
  int suffix_enable;        // 1 if suffixes are enabled, 0 if disabled
  int nclass;               // number of atomic class (atom and/or element), should be either 1 or 2

  CAC(int, char **, MPI_Comm);
  ~CAC();
  void create();
  void destroy();
  void init();

 private:
  CAC() {};                 //prohibit using the default constructor
  CAC(const CAC &) {};      //prohibit using the copy constructor
};

}
#endif
