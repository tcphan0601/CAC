#ifdef COMMAND_CLASS

CommandStyle(change_box, ChangeBox)

#else

#ifndef CAC_CHANGE_BOX_H
#define CAC_CHANGE_BOX_H

#include "pointers.h"

namespace CAC_NS {

class ChangeBox : protected Pointers {
 public:
  ChangeBox(class CAC *);
  void command(int, char **);

 private:
  int scaleflag;
  double scale[3];

  struct Operation {
    int style, flavor;
    int dim, boundindex;
    int vdim1, vdim2;
    double flo, fhi, ftilt;
    double dlo, dhi, dtilt;
    double scale;
  };

  Operation *ops;
  int nops;

  double boxlo[3], h_inv[6];

  void options(int, char **);
  void save_box_state();
  void volume_preserve(int, int, double);
};

}

#endif
#endif
