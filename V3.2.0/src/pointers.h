#ifndef CAC_POINTERS_H
#define CAC_POINTERS_H

#include "cactype.h"
#include <mpi.h>        // IWYU pragma: export
#include <cstddef>      // IWYU pragme: export
#include <cstdio>       // IWYU pragma: export
#include <string>       // IWYU pragma: export
#include <vector>       // IWYU pragma: export

#include "cac.h"

namespace CAC_NS {

// universal defines inside namespace

#define FLERR __FILE__, __LINE__

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

class Pointers {
 public:
  Pointers(CAC *ptr) :
    cac(ptr), 
    memory(ptr->memory), 
    error(ptr->error), 
    universe(ptr->universe), 
    input(ptr->input), 
    atom(ptr->atom), 
    element(ptr->element), 
    update(ptr->update), 
    neighbor(ptr->neighbor), 
    comm(ptr->comm), 
    domain(ptr->domain), 
    force(ptr->force), 
    modify(ptr->modify), 
    group(ptr->group), 
    output(ptr->output), 
    timer(ptr->timer), 
    world(ptr->world), 
    infile(ptr->infile), 
    screen(ptr->screen), 
    logfile(ptr->logfile) {}
  virtual ~Pointers() {}

 protected:
  CAC *cac;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;
 
  Atom *&atom;
  Element *&element;
  Update *&update; 
  Neighbor *&neighbor;
  Comm *&comm;
  Domain *&domain;
  Force *&force;
  Modify *&modify;
  Group *&group;
  Output *&output;
  Timer *&timer;

  MPI_Comm &world;
  FILE *&infile;
  FILE *&screen;
  FILE *&logfile;
};

}

#endif
