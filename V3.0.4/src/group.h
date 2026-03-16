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

  enum{ATOM, NODE, ELEMENT};

  void assign(int, char **); 
  int find(const char *);            // lookup name in list of groups
  int find_or_create(const char *);  // lookup name or create new group
  bigint count_atom(int);
  bigint count_atom(int, int);
  bigint count_elem(int);
  bigint count_elem(int, int);
  bigint count_element_node(int);
  bigint count_element_node(int, int);
  bigint count_node(int);
  bigint count_node(int, int);
  bigint count_vatom(int);
  bigint count_vatom(int, int);

  double mass(int);                              // total mass of atoms in group
  double mass(int, int);
  void xcm(int, double, double *);               // center-of-mass coords of group
  void xcm(int, double, double *, int);
  void vcm(int, double, double *);               // center-of-mass velocity of group
  void vcm(int, double, double *, int);
  //void fcm(int, double *);                     // total force on group
  //void fcm(int, double *, int);
  double ke(int);                                // kinetic energy of group
  double ke(int, int);
  void bounds(int, double *);
  void bounds(int, double *, int);
  void angmom(int, double *, double *);          // angular momentum of group
  void angmom(int, double *, double *, int);
  void inertia(int, double *, double [3][3]);    // inertia tensor
  void inertia(int, double *, double [3][3], int);
  void omega(double *, double [3][3], double *); // angular velocity
  void omega(double *, double [3][3], double *, int);

 private:
  int me;
  std::map<tagint, int> *hash;
  
  int find_unused();
  void clear(int);

  static Group *cptr;
  int molbit;
  

};

}

#endif
