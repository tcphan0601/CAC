#ifdef COMMAND_CLASS

CommandStyle(fill_surface,FillSurface)

#else

#ifndef CAC_FILL_SURFACE_H
#define CAC_FILL_SURFACE_H

#include "pointers.h"

namespace CAC_NS {

class FillSurface : protected Pointers {
 public:
  FillSurface(class CAC *);
  void command(int, char **);
  int *ctypelist;
  double **xtemplist;

 private:
  int igroup,igroupfill,groupbit,groupfillbit,mygroup,iregion;
  int direction;
  int dim;
  double lo,hi;
  double width;

  int num_surface_element;
  int max_surface_element;
  int *surface_element_list;


};

}

#endif
#endif

