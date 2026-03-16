// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef CAC_PAIR_NNP_BASE_H
#define CAC_PAIR_NNP_BASE_H

#include "pair.h"
#include "my_page.h"
//#ifdef DP_USE_CXX_API
//#ifdef CACPLUGIN
//#include "DeepBaseModel.h"
//#else
#include "deepmd/DeepBaseModel.h"
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
#include <map>
#define FLOAT_PREC double

namespace CAC_NS {
class PairDeepBaseModel : public Pair {
 public:
  PairDeepBaseModel(class CAC *,
                    const char *,
                    deepmd_compat::DeepBaseModel *,
                    deepmd_compat::DeepBaseModelDevi *);
  virtual ~PairDeepBaseModel() override;
  //void *extract(const char *, int &) override;
  void init_style() override;
  //void write_restart(FILE *) override;
  //void read_restart(FILE *) override;
  double init_one(int i, int j) override;
  virtual void print_summary(const std::string pre) const;
  int get_node_rank();
  void cum_sum(std::map<int, int> &, std::map<int, int> &);

  std::string get_file_content(const std::string &model);
  std::vector<std::string> get_file_content(
      const std::vector<std::string> &models);
  std::vector<std::string> type_names;
  double ener_unit_cvt_factor, dist_unit_cvt_factor, force_unit_cvt_factor;

  //virtual void print_summary_deep_base(const std::string pre);
 protected:
  // these seems to caused segmentation fault in the constructor when passing deep_model and deep_model_devi from child class to deep_base and deep_base_model_devi
  // Not sure why this is the case in CAC and not in LAMMPS
  // Fix: changed them to pointers;

  deepmd_compat::DeepBaseModel *deep_base;
  deepmd_compat::DeepBaseModelDevi *deep_base_model_devi;
  virtual void allocate();
  double **scale;
  unsigned numb_models;
  double cutoff;
  int numb_types;
  int numb_types_spin;
  std::vector<std::vector<double> > all_force;
  std::ofstream fp;
  int out_freq;
  std::string out_file;
  int dim_fparam;
  int dim_aparam;
  int out_each;
  int out_rel;
  int out_rel_v;
  int stdf_comm_buff_size;
  bool single_model;
  bool multi_models_mod_devi;
  bool multi_models_no_mod_devi;
  bool is_restart;
  std::vector<double> virtual_len;
  std::vector<double> spin_norm;
  // for spin systems, search new index of atoms by their old index
  //std::map<int, int> new_idx_map;
  //std::map<int, int> old_idx_map;
  std::vector<double> fparam;
  std::vector<double> aparam;
  double eps;
  double eps_v;

  void make_fparam_from_compute(std::vector<double> &fparam);
  bool do_compute_fparam;
  char *compute_fparam_id;
  void make_aparam_from_compute(std::vector<double> &aparam);
  bool do_compute_aparam;
  char *compute_aparam_id;
  deepmd::DPBackend backend_type;

  //void make_ttm_fparam(std::vector<double> &fparam);

  //void make_ttm_aparam(std::vector<double> &dparam);
  //bool do_ttm;
  //std::string ttm_fix_id;
  int *counts, *displacements;
  tagint *tagsend, *tagrecv;
  double *stdfsend, *stdfrecv;
  std::vector<int> type_idx_map;

  // reconstructed neighbor list
  int nall;
  int inum;

  int nmax;
  int *ilist;
  int *numneigh;
  int **firstneigh;
  MyPage<int> *ipage;

  // reconstructed communication data
  int nswap, maxswap, *sendnum, *recvnum, *firstrecv;
  int **sendlist, *sendproc, *recvproc;

};

}  // namespace CAC_NS

void make_uniform_aparam(std::vector<double> &daparam,
                         const std::vector<double> &aparam,
                         const int &nlocal);
void ana_st(double &max,
            double &min,
            double &sum,
            const std::vector<double> &vec,
            const int &nloc);

#endif
