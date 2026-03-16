#ifndef CAC_FORCE_H
#define CAC_FORCE_H

#include "pointers.h"
#include <stdio.h>
#include <map>
#include <string>

namespace CAC_NS {

class Force : protected Pointers {
 public:
  double boltz;                      // Boltzmann constant (eng/degree-K)
  double hplanck;                    // Planck's constant (energy-time)
  double mvv2e;                      // conversion of mv^2 to energy
  double ftm2v;                      // conversion of ft/m to velocity
  double mv2d;                       // conversion of mass/volume to density
  double nktv2p;                     // conversion of NkT/V to pressure
  double qqr2e;                      // conversion of q^2/r to energy
  double qe2f;                       // conversion of qE to force
  double vxmu2f;                     // conversion of vx dynamic-visc to force
  double xxt2kmu;                    // conversion of xx/t to kinematic-visc
  double dielectric;                 // dielectric constant
  double qqrd2e;                     // q^2/r to energy w/ dielectric constant
  double e_mass;                     // electron mass
  double hhmrr2e;                    // conversion of (hbar)^2/(mr^2) to energy
  double mvh2r;                      // conversion of mv/hbar to distance
                                     // hbar = h/(2*pi)

  double angstrom;                   // 1 angstrom in native units
  double femtosecond;                // 1 femtosecond in native units
  double qelectron;                  // 1 electron charge abs() in native units

  int newton,newton_pair,newton_bond;   // Newton's 3rd law settings

  class Pair *pair;
  char *pair_style;
  char *pair_match_ptr(Pair *);

  typedef Pair *(*PairCreator)(CAC *);
  typedef std::map<std::string,PairCreator> PairCreatorMap;
  PairCreatorMap *pair_map;
 
  double special_lj[4];      // 1-2, 1-3, 1-4 prefactors for LJ
  double special_coul[4];    // 1-2, 1-3, 1-4 prefactors for Coulombics
  int special_angle;         // 0 if defined angles are ignored
                             //  1 if only weight 1,3 atoms if in an angle
  int special_dihedral;      // 0 if defined dihedrals are ignored
                             // 1 if only weight 1,4 atoms if in a dihedral
  int special_extra;         // extra space for added bonds

  Force(class CAC *);
  ~Force();
  void init();
  void setup();
  void create_pair(const char *, int);
  class Pair *new_pair(const char *, int, int &);
  void store_style(char *&, const char *, int);
  void bounds(const char *, int, char*, int, int &, int &, int nmin=1);
  FILE *open_potential(const char *);
  void potential_date(FILE *, const char *);
  const char *potential_name(const char *);
  bigint memory_usage();
 
 private:
  template <typename T> static Pair *pair_creator(CAC *);

};

}

#endif
