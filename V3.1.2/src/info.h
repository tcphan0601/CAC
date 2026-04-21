#ifdef COMMAND_CLASS

CommandStyle(info, Info)

#else

#ifndef CAC_INFO_H
#define CAC_INFO_H

#include "pointers.h"

namespace CAC_NS {

class Info : protected Pointers {
 public:
  Info(class CAC *cac) : Pointers(cac) {};
  void command(int, char **);

  bool is_active(const char *, const char *);
  bool is_defined(const char *, const char *);
  bool is_available(const char *, const char *);

  //bool has_gzip_support() const;
  //bool has_png_support() const;
  //bool has_jpeg_support() const;
  //bool has_ffmpeg_support() const;
  //bool has_exceptions() const;

  char **get_variable_names(int &num);

 private:
  void available_styles(FILE * out, int flags);

  void atom_styles(FILE * out);
  void element_styles(FILE * out);
  void integrate_styles(FILE * out);
  void minimize_styles(FILE * out);
  void pair_styles(FILE * out);
  void fix_styles(FILE * out);
  void compute_styles(FILE * out);
  void region_styles(FILE * out);
  void dump_styles(FILE * out);
  void command_styles(FILE * out);
};

}

#endif
#endif

