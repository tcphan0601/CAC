#ifndef CAC_LATTICE_H
#define CAC_LATTICE_H

#include "pointers.h"

namespace CAC_NS {

class Lattice : protected Pointers {
 public:
  int style;                           // NONE,SC,FCC,etc
  double xlattice,ylattice,zlattice;   // lattice scale factors in 3 dims
  double a1[3],a2[3],a3[3];            // edge vectors of unit cell
  double scale;
  int nbasis;                          // # of basis atoms in unit cell
  double **basis;                      // fractional coords of each basis atom
                                       // within unit cell (0 <= coord < 1)

  Lattice(class CAC *, int, char **);
  ~Lattice();
  void lattice2box(double &, double &, double &);
  void box2lattice(double &, double &, double &);
  void bbox(int, double, double, double,
            double &, double &, double &, double &, double &, double &);

 private:

  double origin[3];                    // lattice origin
  int orientx[3];                      // lattice orientation vecs
  int orienty[3];                      // orientx = what lattice dir lies
  int orientz[3];                      //           along x dim in box

  double primitive[3][3];              // lattice <-> box transform matrices
  double priminv[3][3];
  double rotaterow[3][3];
  double rotatecol[3][3];

  int orthogonal();
  int right_handed();
  int collinear();
  void setup_transform();
  void add_basis(double, double, double);
  double dot(double *, double *);
  void cross(double *, double *, double *);
};

}

#endif

