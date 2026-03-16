#ifdef COMMAND_CLASS

CommandStyle(displace, Displace)

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
  int igroup, groupbit;
  int scaleflag;
  int group_flag;
  int sequence_flag;
  int tag_match_flag;
  int index_offset;
  int reverse_flag;
  int nskip;
  int style;
  double *mvec;

  void move(int, char *);
  void ramp(char **);
  void polar_ramp(char **);
  void rotate(char **);
  void mirror(char **);
  void multi_dislocations(char **);
  void disp_file(char **);

  void dislocations(char **);
  struct dislocation_info {
    int d_dim;
    int slip_dim;
    int dim1;
    int dim2;
    int direction;
    double center1;
    double center2;
    double b;           // Burgers vector
    double theta;       // rotation angle (in radians)
    double disl_angle;  // angle (in degrees)
    double v;           // Poisson's ratio
  };

  void move(int, double);
  void rotate_point(double *, double *, double *, double, double);
  void options(int, char **);

  void edge_dislocation_field(double &, double &, double, double, double, double, double);
  void screw_dislocation_field(double &, double, double, double, double);
  void ramp_coord(double, double, double, double, double *, int, int);
  void polar_ramp_coord(double, double, double, double, double *, int, int, double *, double *);
};

}

#endif
#endif

