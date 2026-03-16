#ifdef COMMAND_CLASS

CommandStyle(displace,Displace)

#else

#ifndef CAC_DISPLACE_H
#define CAC_DISPLACE_H

#include "pointers.h"

namespace CAC_NS {

class Displace : protected Pointers {
 public:
  Displace(class CAC *);
  ~Displace();
  void command(int, char **);

 private:
  double xscale[3];
  int igroup,groupbit;
  int scaleflag;
  int group_flag;
  double *mvec;

  void move(int, char *);
  void ramp(char **);
  void rotate(char **);
  void mirror(char **);
  void dislocation(char **);
  void multi_dislocations(char **);

  void move(int, double);
  void rotate_point(double *, double *, double *, double, double);
  void options(int, char **);

  void edge_dislocation_field(double &, double &, double, double, double, double);
  void screw_dislocation_field(double &, double, double, double);
  void ramp_coord(double, double, double, double, double *, int, int);
};

}

#endif
#endif

