#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet, Verlet)

#else

#ifndef CAC_VERLET_H
#define CAC_VERLET_H

#include "integrate.h"

namespace CAC_NS {

class Verlet : public Integrate {

 public:
  Verlet(class CAC *, int, char **);
  virtual ~Verlet() {}
  virtual void init();
  virtual void setup(int, int flag = 1);
  virtual void run(int, int);
  void cleanup();

 protected:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int torqueflag, extraflag;

  virtual void force_clear();

};

}
#endif
#endif
