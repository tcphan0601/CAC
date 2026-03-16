#ifndef CAC_GROUP_H
#define CAC_GROUP_H

#include <stdio.h>
#include "pointers.h"
#include <map>

namespace CAC_NS {

class Group : protected Pointers {

 public:
  int ngroup;                  // # of defined groups
  char **names;                // name of each group
  int *bitmask;                // one-bit mask for each group

  int *inversemask;            // inverse mask for each group
  int *dynamic;                // 1 if dynamic, 0 if not
  
  Group(class CAC *);
  ~Group();

  void assign(int,char **); 
  int find(const char *);            // lookup name in list of groups
  int find_or_create(const char *);  // lookup name or create new group
  bigint count_atom(int);
  bigint count_atom(int, int);
  bigint count_elem(int);
  bigint count_elem(int, int);
  bigint count_node(int);
  bigint count_node(int, int);
  bigint count_intpl(int);
  bigint count_intpl(int, int);

 private:
  int me;
  std::map<tagint,int> *hash;
  
  int find_unused();
  void clear(int);

  static Group *cptr;
  int molbit;
  

};

}

#endif
