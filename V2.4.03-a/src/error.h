#ifndef CAC_ERROR_H
#define CAC_ERROR_H

#include "pointers.h"

namespace CAC_NS {

class Error : protected Pointers {
 public:
  Error(class CAC *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);
  void universe_warn(const char *, int, const char *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done();
};

}

#endif

/* ERROR/WARNING messages:

*/
