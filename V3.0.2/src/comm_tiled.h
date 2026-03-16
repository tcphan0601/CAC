#ifndef CAC_COMM_TILED_H
#define CAC_COMM_TILED_H

#include "comm.h"

namespace CAC_NS {

class CommTiled : public Comm {
 public:
  CommTiled(class CAC *);
  CommTiled(class CAC *, class Comm *);
  virtual ~CommTiled();

  void init();
  void setup_exchange();                       // setup 3d exchange comm pattern
  void setup_borders();                        // setup 3d borders comm pattern
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

  virtual void forward_comm_atom_array(int, double **);         // forward comm of per-atom array
  virtual void forward_comm_elem_array(int, double **);         // forward comm of per-element array

  void coord2proc_setup();
  int coord2proc(double *, int &, int &, int &);

  virtual bigint memory_usage();

 private:

  // forward/reverse comm info, proc lists include self
  
  int noverlap;                // # of overlapping procs
  int maxoverlap;              // current max length of overlap
  int *overlap;                // list of overlapping procs
  int nprocmax;

  int overlapother;            // 1 if overlaping with other procs
  int overlapself;             // # of overlaping with self in overlap list
  int *atom_sendnum, *atom_recvnum; // # of atoms to send/recv per proc
  int *atom_size_forward_recv;      // # of atom values to recv in each forward proc
  int *atom_firstrecv;              // where to put 1st recv atom per proc
  int *atom_size_reverse_send;      // # of atom values to send in each reverse proc
  int *atom_size_reverse_recv;      // # of atom values to recv in each reverse proc
  int *atom_forward_recv_offset;    // forward comm offsets in buf_recv per proc
  int *atom_reverse_recv_offset;    // reverse comm offsets in buf_recv per proc
  int **atom_sendlist;              // list of atoms to send per proc
  int *atom_maxsendlist;            // max size of send list per proc

  int *elem_sendnum, *elem_recvnum; // # of elements to send/recv per proc
  int *elem_size_forward_recv;      // # of element values to recv in each forward proc
  int *elem_firstrecv;              // where to put 1st recv element per proc
  int *elem_size_reverse_send;      // # of element values to send in each reverse proc
  int *elem_size_reverse_recv;      // # of element values to recv in each reverse proc
  int *elem_forward_recv_offset;    // forward comm offsets in buf_recv per proc
  int *elem_reverse_recv_offset;    // reverse comm offsets in buf_recv per proc
  int **elem_sendlist;              // list of elements to send per proc
  int *elem_maxsendlist;            // max size of send list per proc
 
  int *pbc_flag;               // general flag for sending atoms/elements thru PBC
  int **pbc;                   // dimension flags for PBC adjustments

  double **minsendbox;         // minimum bounding box of atoms/elements to send per proc
  int *minsendbox_flag;        // 1 if minsendbox for this send proc is finite
  double **maxsendbox;         // maximum bounding box of atoms/elements to send per proc due to element extension

  // exchange foreign box info
  
  int *overlap_sendnum, *overlap_recvnum;
  int **overlap_sendlist;
  int *overlap_maxsendlist;
  int *overlap_firstrecv;

  // exchange comm info, proc lists do not include self

  int *nexchproc;              // # of procs to send/recv to/from in each dim
  int *nexchprocmax;           // current max # of exch procs for each dim
  int **exchproc;              // procs to exchange with per dim
  int **exchnum;               // # of values received per dim/proc

  double *buf_send;            // send buffer for all comm
  double *buf_recv;            // recv buffer for all comm
  int maxsend, maxrecv;         // current size of send/recv buffer
  int bufextra;                // extra space beyond maxsend in send buffer
  int smaxone, rmaxone;         // max size in atoms of single borders send/recv
  int smaxall, rmaxall;         // max size in atoms of any borders send/recv
                               //   for comm to all procs in one swap

  int maxreqstat;              // max size of Request and Status vectors
  MPI_Request *requests;

  struct RCBinfo {
    double mysplit[3][2];      // fractional RCB bounding box for one proc
    double cutfrac;    	       // fractional position of cut this proc owns
    int dim;	                 // dimension = 0/1/2 of cut
  };

  RCBinfo *rcbinfo;            // list of RCB info for all procs

  // data for setup_borders()
  
  // total 7x7x7 = 343 boxes so need 11 ints of 32 bits each (each bit is a flag)
  int **overlap_repeat;        // 1st flags (for boxes 0-->342) for each proc to determine if overlap array has repeats in O(P)
  int **overlap_pbc;           // list of image offsets in each dim for each overlap
  double **overlap2box;        // list of overlap domain boxes setup by overlap calculation for brick 
  double expand_subbox[6];     // subbox of me expanded by my element bounding boxes

  // foreign boxes to further expand ghost sendbox
  
  int *foreign_overlap;        // which overlap this foreign box is relevant to
  double **foreign_boxes;      // foreign box that I need to send ghost
  int nforeign;                // number of foreign boxes I received
  int maxforeign;              // max size of foreign arrays 
  double maxfbox[3];           // max size of foreign boxes I received

  // foreign box bin info
  
  int nbinx, nbiny, nbinz;
  int mbins;
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;
  double binsizex, binsizey, binsizez;
  double bininvx, bininvy, bininvz;
  int *fboxbins;
  int *fboxbinhead;
  int maxbin;

  // foreign box stencil info

  int maxstencil;
  int nstencil;
  int *stencil;

  double *prd;                 // local ptrs to Domain attributes
  double *boxlo, *boxhi;
  double *sublo, *subhi;
  int dimension;

  // NOTE: init_buffers is called from a constructor and must not be made virtual
  void init_buffers();
  void deallocate_buffers();

  // box drop and other functions

  typedef void (CommTiled:: * BoxDropPtr)(int, double *, double *, int *, int &, int);
  BoxDropPtr box_drop;
  void box_drop_brick(int, double *, double *, int *, int &, int);
  void box_drop_tiled(int, double *, double *, int *, int &, int);
  void box_drop_tiled_recurse(double *, double *, int, int, int *, int &, int);

  typedef void (CommTiled:: * BoxOtherPtr)(int, int, int, double *, double *, int);
  BoxOtherPtr box_other;
  void box_other_brick(int, int, int, double *, double *, int);
  void box_other_tiled(int, int, int, double *, double *, int);

  typedef int (CommTiled:: * BoxTouchPtr)(int, int, int);
  BoxTouchPtr box_touch;
  int box_touch_brick(int, int, int);
  int box_touch_tiled(int, int, int);

  typedef int (CommTiled:: * PointDropPtr)(int, double *);
  PointDropPtr point_drop;
  int point_drop_brick(int, double *);
  int point_drop_tiled(int, double *);
  int point_drop_tiled_recurse(double *, int, int);
  int closer_subbox_edge(int, double *);

  // functions for foreign box bins and stencil

  void setup_bins();
  void bin_foreign_boxes();
  int coord2bin(double, double, double);
  double bin_distance(int, int);
  void create_stencil();

  void grow_send(int, int);           // reallocate send buffer
  void grow_recv(int);                // free/allocate recv buffer
  void grow_atom_list(int, int);      // reallocate one atom_sendlist
  void grow_elem_list(int, int);      // reallocate one elem_sendlist
  void grow_overlap_list(int, int);   // reallocate one overlap_sendlist
  void grow_overlap(int, int);        // allocate arrays for send/recv
  void grow_foreign_box();            // reallocate foreign box arrays
  void grow_overlap_array();
  int pack_element_bound_box(int, int, int *, double *);
  int unpack_element_bound_box(int, int, double *, int);
 
  // debug info, flag, and function (can be removed later)
  
  int *foreign_tag;          // ID of element this foreign box belong to
  int debugflag;
  void debug_write_foreign_box(int);
  void debug_write_domain_box();
};

}

#endif
