#ifdef FIX_CLASS

FixStyle(move, FixMove)

#else

#ifndef CAC_FIX_MOVE_H
#define CAC_FIX_MOVE_H

#include "fix.h"

namespace CAC_NS {

class FixMove : public Fix {
 public:
  FixMove(class CAC *, int, char **);
  ~FixMove();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  //void initial_integrate_respa(int, int, int);
  //void final_integrate_respa(int, int);

  double memory_usage();
  //void write_restart(FILE *);
  //void restart(char *);
  void grow_atom_arrays(int);
  void grow_elem_arrays(int);
  void copy_atom_arrays(int, int, int);
  void copy_elem_arrays(int, int, int);
  void set_atom_arrays(int);
  void set_elem_arrays(int);
  int pack_atom_exchange(int, double *);
  int unpack_atom_exchange(int, double *);
  int pack_elem_exchange(int, double *);
  int unpack_elem_exchange(int, double *);

  //int pack_restart(int, double *);
  //void unpack_restart(int, int);
  //int maxsize_restart();
  //int size_restart(int);

  void reset_dt();

 private:
  char *xvarstr, *yvarstr, *zvarstr, *vxvarstr, *vyvarstr, *vzvarstr;
  int mstyle;
  int vxflag, vyflag, vzflag, axflag, ayflag, azflag;
  double vx, vy, vz, ax, ay, az;
  double period, omega_rotate;
  double point[3], axis[3], runit[3];
  double dt, dtv, dtf;
  int xvar, yvar, zvar, vxvar, vyvar, vzvar;
  int xvarstyle, yvarstyle, zvarstyle, vxvarstyle, vyvarstyle, vzvarstyle;
  //int extra_flag, omega_flag, angmom_flag;
  //int radius_flag, ellipsoid_flag, line_flag, tri_flag, body_flag;
  //int theta_flag, quat_flag;
  //int nlevels_respa, nrestart;
  int time_origin;

  int element_group_style;
  double **xoriginal;         // original coords of atoms
  double ****nodexoriginal;    // original coords of elements' nodes
  //double *toriginal;          // original theta of atoms
  //double **qoriginal;         // original quat of atoms
  int displaceflag, velocityflag;
  int maxatom;
  int maxelem;
  double **displace, **velocity;
  double ****nodedisplace, ****nodevelocity;

  //class AtomVecEllipsoid *avec_ellipsoid;
  //class AtomVecLine *avec_line;
  //class AtomVecTri *avec_tri;
  //class AtomVecBody *avec_body;
};

}

#endif
#endif
