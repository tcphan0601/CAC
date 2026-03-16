#ifndef CAC_COMM_BRICK_H
#define CAC_COMM_BRICK_H

#include "comm.h"

namespace CAC_NS {

class CommBrick : public Comm {
 public:
  CommBrick(class CAC *);
  CommBrick(class CAC *, class Comm *);
  virtual ~CommBrick();

  virtual void init();
  virtual void setup_borders();                // setup 3d borders comm pattern
  virtual void forward_comm(int dummy = 0);    // forward comm of atom coords
  virtual void reverse_comm();                 // reverse comm of forces
  virtual void exchange();                     // move atoms to new procs
  virtual void borders();                      // setup list of atoms to comm

  virtual void forward_comm_pair(class Pair *);    // forward comm from a Pair
  virtual void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  virtual void forward_comm_fix(class Fix *, int size = 0);
                                                   // forward comm from a Fix
  virtual void reverse_comm_fix(class Fix *, int size = 0);
                                                   // reverse comm from a Fix
  virtual void reverse_comm_fix_variable(class Fix *);
                                     // variable size reverse comm from a Fix
  virtual void forward_comm_compute(class Compute *);  // forward from a Compute
  virtual void reverse_comm_compute(class Compute *);  // reverse from a Compute
  virtual void forward_comm_dump(class Dump *);    // forward comm from a Dump
  virtual void reverse_comm_dump(class Dump *);    // reverse comm from a Dump

  void forward_comm_atom_array(int, double **);         // forward comm of per-atom array
  void forward_comm_elem_array(int, double **);         // forward comm of per-element array
  int exchange_variable(int, double *, double *&);  // exchange on neigh stencil
  virtual bigint memory_usage();

 protected:
  int nswap;                        // # of swaps to perform = sum of maxneed
  int recvneed[3][2];               // # of procs away I recv atoms from
  int sendneed[3][2];               // # of procs away I send atoms to
  int maxneed[3];                   // max procs away any proc needs, per dim
  int maxswap;                      // max # of swaps memory is allocated for
  int *sendproc,*recvproc;          // proc to send/recv to/from at each swap

  int *atom_sendnum, *atom_recvnum;          // # of atoms to send/recv in each swap
  int *atom_size_forward_recv;          // # of values to recv in each forward comm
  int *atom_size_reverse_send;          // # to send in each reverse comm
  int *atom_size_reverse_recv;          // # to recv in each reverse comm
  int *atom_firstrecv;                  // where to put 1st recv atom in each swap
  int **atom_sendlist;                  // list of atoms to send in each swap
  int *atom_maxsendlist;                // max size of send list for each swap

  int *elem_sendnum, *elem_recvnum;        // # of elements to send/recv in each swap 
  int *elem_size_forward_recv;         // # of values to recv in each forward comm for elements
  int *elem_size_reverse_send;         // # to send in each reverse comm for elements
  int *elem_size_reverse_recv;         // # to recv in each reverse comm for elements
  int *elem_firstrecv;                 // where to put 1st recv elements in each swap
  int **elem_sendlist;                 // list of elements to send in each swap 
  int *elem_maxsendlist;               // max size of send list of elments for each swap 


  double *slablo,*slabhi;           // bounds of slab to send at each swap
  double **multilo,**multihi;       // bounds of slabs for multi-type swap
  double **cutghostmulti;           // cutghost on a per-type basis
  int *pbc_flag;                    // general flag for sending atoms thru PBC
  int **pbc;                        // dimension flags for PBC adjustments

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int bufextra;                     // extra space beyond maxsend in send buffer
  int smax,rmax;             // max size in atoms of single borders send/recv

  // NOTE: init_buffers is called from a constructor and must not be made virtual
  void init_buffers();

  int updown(int, int, int, double, int, double *);
                                            // compare cutoff to procs
  virtual void grow_send(int, int);         // reallocate send buffer
  virtual void grow_recv(int);              // free/allocate recv buffer
  virtual void grow_atom_list(int, int);    // reallocate one atom_sendlist
  virtual void grow_elem_list(int, int);    // reallocate one elem_sendlist
  virtual void grow_swap(int);              // grow swap and multi arrays
  virtual void allocate_swap(int);          // allocate swap arrays
  virtual void allocate_multi(int);         // allocate multi arrays
  virtual void free_swap();                 // free swap arrays
  virtual void free_multi();                // free multi arrays
};

}

#endif

