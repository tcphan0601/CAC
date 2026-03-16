// SPDX-License-Identifier: LGPL-3.0-or-later
#ifdef PAIR_CLASS

PairStyle(deepmd, PairDeepMD)

#else

#ifndef CAC_PAIR_NNP_H
#define CAC_PAIR_NNP_H

//#ifdef DP_USE_CXX_API
//#ifdef CACPLUGIN
//#include "DeepPot.h"
//#else
#include "deepmd/DeepPot.h"
//#endif
namespace deepmd_compat = deepmd;
//#else
//#ifdef CACPLUGIN
//#include "deepmd.hpp"
//#else
//#include "deepmd/deepmd.hpp"
//#endif
//namespace deepmd_compat = deepmd::hpp;
//#endif

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include "comm_brick.h"
#include "comm_tiled.h"
#include "pair_base.h"
#include "universe.h"
#define FLOAT_PREC double

namespace CAC_NS {
class CommBrickDeepMD : public CommBrick {
  friend class PairDeepMD;
};
class CommTiledDeepMD : public CommTiled {
  friend class PairDeepMD;
};

class PairDeepMD : public PairDeepBaseModel {
 public:
  PairDeepMD(class CAC *);
  ~PairDeepMD() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void compute(int, int) override;
  int pack_atom_reverse_comm(int, int, double *) override;
  void unpack_atom_reverse_comm(int, int *, double *) override;

  virtual void print_summary(const std::string pre) const;
 protected:
  deepmd_compat::DeepPot deep_pot;
  deepmd_compat::DeepPotModelDevi deep_pot_model_devi;
  virtual void ev_setup(int, int, int alloc = 1);

 private:
  CommBrickDeepMD *commbrickdata_;
  CommTiledDeepMD *commtileddata_;

  void reconstruct_input_nlist(std::vector<double> &, std::vector<int> &,
      std::unordered_map<std::vector<int>, int, VectorHasher<int>> &,
      std::unordered_map<int, std::vector<int>> &
      );

  int insert_to_deepmd_arrays(int, int, double *, int, int,
      std::unordered_map<std::vector<int>, int, VectorHasher<int>> &,
      std::unordered_map<int, std::vector<int>> &,
      std::vector<double> &, std::vector<int> &
      );

  void dump_debug(int, std::vector<double>
      , std::unordered_map<int, std::vector<int>>, std::unordered_set<std::vector<int>, VectorHasher<int>>);
  void dump_debug_neighbor(int, std::vector<double>
      , std::unordered_map<int, std::vector<int>>, std::unordered_set<std::vector<int>, VectorHasher<int>>);
};

}  // namespace CAC_NS

#endif
#endif
