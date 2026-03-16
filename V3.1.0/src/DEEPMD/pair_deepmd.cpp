// SPDX-License-Identifier: LGPL-3.0-or-later
#include <string.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <sstream>

#include "atom.h"
#include "element.h"
#include "element_vec.h"
//#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "update.h"
#include "universe.h"
#include "my_page.h"
/*
#if CAC_VERSION_NUMBER >= 20210831
// in lammps #2902, fix_ttm members turns from private to protected
#define USE_TTM 1
#include "fix_ttm_dp.h"
#endif
*/
#include "deepmd_version.h"
#include "pair_deepmd.h"


#define STEP 11

using namespace CAC_NS;
using namespace std;

static const char cite_user_deepmd_package[] =
    "USER-DEEPMD package:\n\n"
    "@article{Wang_ComputPhysCommun_2018_v228_p178,\n"
    "  author = {Wang, Han and Zhang, Linfeng and Han, Jiequn and E, Weinan},\n"
    "  doi = {10.1016/j.cpc.2018.03.016},\n"
    "  url = {https://doi.org/10.1016/j.cpc.2018.03.016},\n"
    "  year = 2018,\n"
    "  month = {jul},\n"
    "  publisher = {Elsevier {BV}},\n"
    "  volume = 228,\n"
    "  journal = {Comput. Phys. Commun.},\n"
    "  title = {{DeePMD-kit: A deep learning package for many-body potential "
    "energy representation and molecular dynamics}},\n"
    "  pages = {178--184}\n"
    "}\n"
    "@article{Zeng_JChemPhys_2023_v159_p054801,\n"
    "  title  = {{DeePMD-kit v2: A software package for deep potential "
    "models}},\n"
    "  author =   {Jinzhe Zeng and Duo Zhang and Denghui Lu and Pinghui Mo and "
    "Zeyu Li\n"
    "         and Yixiao Chen and Mari{\\'a}n Rynik and Li'ang Huang and Ziyao "
    "Li and \n"
    "         Shaochen Shi and Yingze Wang and Haotian Ye and Ping Tuo and "
    "Jiabin\n"
    "         Yang and Ye Ding and Yifan Li and Davide Tisi and Qiyu Zeng and "
    "Han \n"
    "         Bao and Yu Xia and Jiameng Huang and Koki Muraoka and Yibo Wang "
    "and \n"
    "         Junhan Chang and Fengbo Yuan and Sigbj{\\o}rn L{\\o}land Bore "
    "and "
    "Chun\n"
    "         Cai and Yinnian Lin and Bo Wang and Jiayan Xu and Jia-Xin Zhu "
    "and \n"
    "         Chenxing Luo and Yuzhi Zhang and Rhys E A Goodall and Wenshuo "
    "Liang\n"
    "         and Anurag Kumar Singh and Sikai Yao and Jingchao Zhang and "
    "Renata\n"
    "         Wentzcovitch and Jiequn Han and Jie Liu and Weile Jia and Darrin "
    "M\n"
    "         York and Weinan E and Roberto Car and Linfeng Zhang and Han "
    "Wang},\n"
    "  journal =  {J. Chem. Phys.},\n"
    "  volume =   159,\n"
    "  issue =    5,  \n"
    "  year =    2023,\n"
    "  pages  =   054801,\n"
    "  doi =      {10.1063/5.0155600},\n"
    "}\n"
    "@Article{Zeng_arXiv_2025_p2502.19161,\n"
    "  annote       = {general purpose},\n"
    "    author =   {Jinzhe Zeng and Duo Zhang and Anyang Peng and Xiangyu "
    "Zhang and Sensen\n"
    "             He and Yan Wang and Xinzijian Liu and Hangrui Bi and Yifan "
    "Li and Chun\n"
    "             Cai and Chengqian Zhang and Yiming Du and Jia-Xin Zhu and "
    "Pinghui Mo\n"
    "             and Zhengtao Huang and Qiyu Zeng and Shaochen Shi and "
    "Xuejian Qin and\n"
    "             Zhaoxi Yu and Chenxing Luo and Ye Ding and Yun-Pei Liu and "
    "Ruosong Shi\n"
    "             and Zhenyu Wang and Sigbj{\\o}rn L{\\o}land Bore and Junhan "
    "Chang and\n"
    "             Zhe Deng and Zhaohan Ding and Siyuan Han and Wanrun Jiang "
    "and Guolin\n"
    "             Ke and Zhaoqing Liu and Denghui Lu and Koki Muraoka and "
    "Hananeh Oliaei\n"
    "             and Anurag Kumar Singh and Haohui Que and Weihong Xu and "
    "Zhangmancang\n"
    "             Xu and Yong-Bin Zhuang and Jiayu Dai and Timothy J. Giese "
    "and Weile\n"
    "             Jia and Ben Xu and Darrin M. York and Linfeng Zhang and Han "
    "Wang},\n"
    "    title =    {{DeePMD-kit v3: A Multiple-Backend Framework for Machine "
    "Learning\n"
    "             Potentials}},\n"
    "    journal =  {arXiv},\n"
    "    year =     2025,\n"
    "    pages =    {2502.19161},\n"
    "    doi =      {10.48550/arXiv.2502.19161},\n"
    "}\n\n";


PairDeepMD::PairDeepMD(CAC *cac) : PairDeepBaseModel(
          cac, cite_user_deepmd_package, &deep_pot, &deep_pot_model_devi) {

}

/* ---------------------------------------------------------------------- */

PairDeepMD::~PairDeepMD() {
  // Ensure base class destructor is called
}

/* ---------------------------------------------------------------------- */

void PairDeepMD::compute(int eflag, int vflag) 
{
  if (numb_models == 0) {
    return;
  }
  // See
  // https://docs.lammps.org/Developer_updating.html#use-ev-init-to-initialize-variables-derived-from-eflag-and-vflag
  ev_init(eflag, vflag);
  if (vflag_atom) {
    error->all(FLERR,
               "6-element atomic virial is not supported. Use compute "
               "centroid/stress/atom command for 9-element atomic virial.");
  }

  // do_ghost is always true currently (not sure why it is that way)
  // could be changed in future versions of deepmd-kit

  bool do_ghost = true;

  //  dpa2 communication
  //  For PyTorch backend only that require passing of ghost data
  //  during NN calculation
  //  Currently comm tiled does not work yet and will not work for PyTorch,
  //    for other backend, just passing dummy variables for communication data 
  //    and thus it does not matter which comm style is used
  if (backend_type == deepmd::DPBackend::PyTorch) {
    if (comm->style == 1) 
      error->all(FLERR, "Comm Tiled does not work with PyTorch backend yet");
    else {
      commbrickdata_ = (CommBrickDeepMD *) comm;
      error->all(FLERR, "Need to reconstruct comm_data");
      // add reconstructed communication data here
    }
  } 

  double **f = atom->f;
  double ****gaussf = element->gaussf;
  int *etype = element->etype;
  int *apc = element->apc;
  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  double *nodal_weight = element->nodal_weight;
  int naghost = 0;
  int neghost = 0;
  if (do_ghost) {
    naghost = atom->nghost;
    neghost = element->nghost;
  }
  int naall = nalocal + naghost;
  int neall = nelocal + neghost;
  int ntotal_all = naall + neall * element->max_nucell * element->maxapc;

  // reestimate the total number of atoms/virtual atoms involved if neighbor list was rebuilt
  // pre-allocate memory to nmax

  if (neighbor->ago == 0) {
    int ntotal_local = nalocal;
    ntotal_local += list->va2nvalist.size() + nelocal * element->max_ngcell * element->maxapc;
    if (ntotal_local > nmax) {
      nmax = ntotal_local;
      memory->destroy(ilist);
      memory->destroy(numneigh);
      memory->sfree(firstneigh);
      memory->create(ilist, nmax, "deepmd:ilist");
      memory->create(numneigh, nmax, "deepmd:numneigh");
      firstneigh = (int **) memory->smalloc(nmax * sizeof(int *), "deepmd:firstneigh");
    }
  }

  // preallocate memory for dcoord and the index_maps with nmax
  vector<double> dcoord;
  vector<int> dtype;
  unordered_map<vector<int>, int, VectorHasher<int>> index_map_cac2deep;
  unordered_map<int, vector<int>> index_map_deep2cac;

  dcoord.reserve(ntotal_all * 3);
  dtype.reserve(ntotal_all);
  index_map_cac2deep.reserve(ntotal_all);
  index_map_deep2cac.reserve(ntotal_all);
  reconstruct_input_nlist(dcoord, dtype, index_map_cac2deep, index_map_deep2cac);
  int nlocal = inum;
  int nghost = nall - nlocal;
  int newton_pair = force->newton_pair;

  //if (atom->sp_flag) {
  //  error->all(
  //      FLERR,
  //      "Pair style 'deepmd' does not support spin atoms, please use pair "
  //      "style 'deepspin' instead.");
  //}

  double dener(0);
  vector<double> dforce(nall * 3);
  vector<double> dvirial(9, 0);
  vector<double> dbox(9, 0);
  vector<double> daparam;

  // get box
  dbox[0] = domain->h[0] / dist_unit_cvt_factor;  // xx
  dbox[4] = domain->h[1] / dist_unit_cvt_factor;  // yy
  dbox[8] = domain->h[2] / dist_unit_cvt_factor;  // zz
  dbox[7] = domain->h[3] / dist_unit_cvt_factor;  // zy
  dbox[6] = domain->h[4] / dist_unit_cvt_factor;  // zx
  dbox[3] = domain->h[5] / dist_unit_cvt_factor;  // yx

  // mapping (for DPA-2 JAX)
  //vector<int> mapping_vec(nall, -1);

  // This mapping seems incomplete (only calculated in serial and not parallel)
  // Thus it is disabled here and if it should be fixed when 
  //   it is fully implemented in future version of deepmd-kit

  //if (comm->nprocs == 1 && atom->map_style != Atom::MAP_NONE) {
  //  for (size_t ii = 0; ii < nall; ++ii) {
  //    mapping_vec[ii] = atom->map(atom->tag[ii]);
  //  }
  //}
  
  if (do_compute_aparam) {
    make_aparam_from_compute(daparam);
  } else if (aparam.size() > 0) {
    // uniform aparam
    make_uniform_aparam(daparam, aparam, nlocal);
  }

  if (do_compute_fparam) {
    make_fparam_from_compute(fparam);
  }
  // neighbor is rebuilt every step
  // need to do this so that the neighbor list internally within deepmd is copied from this new list every step
  // this one bug took 3 days to find...! :(
  
  int ago = 0;

  //int ago = neighbor->ago;
  //if (numb_models > 1) {
  //  if (multi_models_no_mod_devi &&
  //      (out_freq > 0 && update->ntimestep % out_freq == 0)) {
  //    ago = 0;
  //  } else if (multi_models_mod_devi &&
  //      (out_freq == 0 || update->ntimestep % out_freq != 0)) {
  //    ago = 0;
  //  }
  //}

  // compute

  //if (update->ntimestep == 11) return;
  single_model = (numb_models == 1);
  multi_models_no_mod_devi =
    (numb_models > 1 && (out_freq == 0 || update->ntimestep % out_freq != 0));
  multi_models_mod_devi =
    (numb_models > 1 && (out_freq > 0 && update->ntimestep % out_freq == 0));
  if (do_ghost) {

    deepmd_compat::InputNlist cac_list(
        inum, ilist, numneigh, firstneigh,
        nswap, sendnum, recvnum, firstrecv, 
        sendlist, sendproc, recvproc, &world);


    //cac_list.set_mask(NEIGHMASK);
    if (comm->nprocs == 1 && atom->map_style != Atom::MAP_NONE) {
      //cac_list.set_mapping(mapping_vec.data());
      error->all(FLERR,"TEST"); 
    }
    deepmd_compat::InputNlist extend_cac_list;
    if (single_model || multi_models_no_mod_devi) {
      // cvflag_atom is the right flag for the cvatom matrix
      if (!(eflag_atom || cvflag_atom)) {
        try {
          deep_pot.compute(dener, dforce, dvirial, dcoord, dtype, dbox, nghost,
              cac_list, ago, fparam, daparam);
        } catch (deepmd_compat::deepmd_exception &e) {
          error->one(FLERR, e.what());
        }
      }
      // do atomic energy and virial
      else {
        vector<double> deatom(nall * 1, 0);
        vector<double> dvatom(nall * 9, 0);

        try {
          deep_pot.compute(dener, dforce, dvirial, deatom, dvatom, dcoord,
              dtype, dbox, nghost, cac_list, ago, fparam, daparam);
        } catch (deepmd_compat::deepmd_exception &e) {
          error->one(FLERR, e.what());
        }

        // only store virial for atoms and nodes
        // also re-tally dener and dvirial

        if (evflag) {
          dener = 0.0;
          std::fill(dvirial.begin(), dvirial.end(), 0.0);
          //          if (comm->me == 1) printf("step = %d before virial calc dvirial = %g\n",update->ntimestep,dvirial[0]);
          for (int ideep = 0; ideep < nall; ideep++) {
            if (index_map_deep2cac.find(ideep)==index_map_deep2cac.end()) error->one(FLERR,"TEST"); // sanity check
            int i = index_map_deep2cac[ideep][0];
            int iindex = index_map_deep2cac[ideep][1];
            double *eptr, *cvptr;
            double iescale, ivscale;

            if (iindex == -1) {
              if (eflag_atom) eptr = &eatom[i];
              if (cvflag_atom) cvptr = cvatom[i];
              iescale = ivscale = 1.0;
            } else if ((iindex < element->maxapc * element->maxucell) && (iindex >= 0)) {
              int ietype = etype[i];
              int iapc = apc[ietype];
              int ibasis = iindex % iapc;
              int iucell = iindex / iapc;
              int igcell = u2g[ietype][iucell];
              int inode = -1;
              if (igcell >= 0) 
                inode = g2n[ietype][igcell];
              // skip tallying if not a node
              if (inode >= 0) {
                if (eflag_atom) eptr = &enode[i][ibasis][inode];
                if (cvflag_atom) cvptr = cvnode[i][ibasis][inode];
                if (nodal_energy_weight < 0) { 
                  iescale = nodal_weight[ietype] / 2.0;
                } else iescale = nodal_energy_weight / 2.0;
                ivscale = nodal_weight[ietype];
              } else continue;
            } else continue;

            // energy
            if (eflag) {
              double evdwl = scale[1][1] * deatom[ideep] * ener_unit_cvt_factor;
              if (eflag_atom) *eptr +=  evdwl;
              if (eflag_global) dener += evdwl * iescale;
            }
            // virials
            // Added by Davide Tisi 2020
            // interface the atomic virial computed by DeepMD
            // with the one used in centroid atoms

            if (vflag) {
              if (cvflag_atom) {
                cvptr[0] += scale[1][1] * dvatom[9 * ideep + 0] * ener_unit_cvt_factor;  // xx
                cvptr[1] += scale[1][1] * dvatom[9 * ideep + 4] * ener_unit_cvt_factor;  // yy
                cvptr[2] += scale[1][1] * dvatom[9 * ideep + 8] * ener_unit_cvt_factor;  // zz
                cvptr[3] += scale[1][1] * dvatom[9 * ideep + 3] * ener_unit_cvt_factor;  // xy
                cvptr[4] += scale[1][1] * dvatom[9 * ideep + 6] * ener_unit_cvt_factor;  // xz
                cvptr[5] += scale[1][1] * dvatom[9 * ideep + 7] * ener_unit_cvt_factor;  // yz
                cvptr[6] += scale[1][1] * dvatom[9 * ideep + 1] * ener_unit_cvt_factor;  // yx
                cvptr[7] += scale[1][1] * dvatom[9 * ideep + 2] * ener_unit_cvt_factor;  // zx
                cvptr[8] += scale[1][1] * dvatom[9 * ideep + 5] * ener_unit_cvt_factor;  // zy
              } 

              if (vflag_global) {
                for (int d = 0; d < 9; d++)
                  dvirial[d] += scale[1][1] * dvatom[9 * ideep + d] * ener_unit_cvt_factor * ivscale; 
              }
            }
          }
        }
      }
    } else if (multi_models_mod_devi) {
      error->all(FLERR,"TEST"); 
      vector<double> deatom(nall * 1, 0);
      vector<double> dvatom(nall * 9, 0);
      vector<vector<double>> all_virial;
      vector<double> all_energy;
      vector<vector<double>> all_atom_energy;
      vector<vector<double>> all_atom_virial;
      if (!(eflag_atom || cvflag_atom)) {
        try {
          deep_pot_model_devi.compute(all_energy, all_force, all_virial, dcoord,
              dtype, dbox, nghost, cac_list, ago,
              fparam, daparam);
        } catch (deepmd_compat::deepmd_exception &e) {
          error->one(FLERR, e.what());
        }
      } else {
        try {
          deep_pot_model_devi.compute(all_energy, all_force, all_virial,
              all_atom_energy, all_atom_virial, dcoord,
              dtype, dbox, nghost, cac_list, ago,
              fparam, daparam);
        } catch (deepmd_compat::deepmd_exception &e) {
          error->one(FLERR, e.what());
        }
      }

      // only store virial for atoms and nodes
      // also re-tally dener and dvirial

      dener = all_energy[0];
      dforce = all_force[0];
      dvirial = all_virial[0];
      if (evflag) {
        deatom = all_atom_energy[0];
        dener = 0.0;
        dvatom = all_atom_virial[0];
        std::fill(dvirial.begin(), dvirial.end(), 0.0);
        for (int ideep = 0; ideep < nall; ideep++) {
          int i = index_map_deep2cac[ideep][0];
          int iindex = index_map_deep2cac[ideep][1];
          double *eptr, *cvptr;
          double iescale, ivscale;
          if (iindex == -1) {
            if (eflag_atom) eptr = &eatom[i];
            if (cvflag_atom) cvptr = cvatom[i];
            iescale = ivscale = 1;
          } else if (iindex >= 0 && iindex < element->maxapc * element->maxucell) {
            int ietype = etype[i];
            int iapc = apc[ietype];
            int ibasis = iindex % iapc;
            int iucell = iindex / iapc;
            int igcell = u2g[ietype][iucell];
            int inode = -1;
            if (igcell >= 0) 
              inode = g2n[ietype][igcell];
            if (inode >= 0) {
              if (eflag_atom) eptr = &enode[i][ibasis][inode];
              if (cvflag_atom) cvptr = cvnode[i][ibasis][inode];
              if (nodal_energy_weight < 0) iescale = nodal_weight[ietype] / 2.0;
              else iescale = nodal_energy_weight / 2.0;
              ivscale = nodal_weight[ietype];
            } else continue;
          } else continue;
          // energy
          if (eflag) {
            double evdwl = scale[1][1] * deatom[ideep] * ener_unit_cvt_factor;
            if (eflag_atom) *eptr +=  evdwl;
            if (eflag_global) dener += evdwl * iescale;
          }
          // virials
          // Added by Davide Tisi 2020
          // interface the atomic virial computed by DeepMD
          // with the one used in centroid atoms

          if (vflag) {
            if (cvflag_atom) {
              cvptr[0] += scale[1][1] * dvatom[9 * ideep + 0] * ener_unit_cvt_factor;  // xx
              cvptr[1] += scale[1][1] * dvatom[9 * ideep + 4] * ener_unit_cvt_factor;  // yy
              cvptr[2] += scale[1][1] * dvatom[9 * ideep + 8] * ener_unit_cvt_factor;  // zz
              cvptr[3] += scale[1][1] * dvatom[9 * ideep + 3] * ener_unit_cvt_factor;  // xy
              cvptr[4] += scale[1][1] * dvatom[9 * ideep + 6] * ener_unit_cvt_factor;  // xz
              cvptr[5] += scale[1][1] * dvatom[9 * ideep + 7] * ener_unit_cvt_factor;  // yz
              cvptr[6] += scale[1][1] * dvatom[9 * ideep + 1] * ener_unit_cvt_factor;  // yx
              cvptr[7] += scale[1][1] * dvatom[9 * ideep + 2] * ener_unit_cvt_factor;  // zx
              cvptr[8] += scale[1][1] * dvatom[9 * ideep + 5] * ener_unit_cvt_factor;  // zy
            } 

            if (vflag_global) {
              for (int d = 9 * ideep; d < 9 * ideep + 9; d++)
                dvirial[d] += scale[1][1] * dvatom[d] * ener_unit_cvt_factor * ivscale; 
            }
          }
        }
      }

      // Exporting model deviation currently not working yet
      /*
         if (out_freq > 0 && update->ntimestep % out_freq == 0) {
         int rank = comm->me;
      // std force
      if (newton_pair) {
      //#if CAC_VERSION_NUMBER >= 20220324
      //          comm->reverse_comm(this);
      //#else
      comm->reverse_comm_pair(this);
      //#endif
      }
      vector<double> std_f;
      vector<double> tmp_avg_f;
      deep_pot_model_devi.compute_avg(tmp_avg_f, all_force);
      deep_pot_model_devi.compute_std_f(std_f, tmp_avg_f, all_force);
      if (out_rel == 1) {
      deep_pot_model_devi.compute_relative_std_f(std_f, tmp_avg_f, eps);
      }
      double min = numeric_limits<double>::max(), max = 0, avg = 0;
      ana_st(max, min, avg, std_f, nlocal);
      double all_f_min = 0, all_f_max = 0, all_f_avg = 0;
      MPI_Reduce(&min, &all_f_min, 1, MPI_DOUBLE, MPI_MIN, 0, world);
      MPI_Reduce(&max, &all_f_max, 1, MPI_DOUBLE, MPI_MAX, 0, world);
      MPI_Reduce(&avg, &all_f_avg, 1, MPI_DOUBLE, MPI_SUM, 0, world);
      all_f_avg /= double(atom->natoms);
      // std v
      vector<double> send_v(9 * numb_models);
      vector<double> recv_v(9 * numb_models);
      for (int kk = 0; kk < numb_models; ++kk) {
      for (int ii = 0; ii < 9; ++ii) {
      send_v[kk * 9 + ii] = all_virial[kk][ii] / double(atom->natoms);
      }
      }
      MPI_Reduce(&send_v[0], &recv_v[0], 9 * numb_models, MPI_DOUBLE, MPI_SUM,
      0, world);
      vector<vector<double>> all_virial_1(numb_models);
      vector<double> avg_virial, std_virial;
      for (int kk = 0; kk < numb_models; ++kk) {
      all_virial_1[kk].resize(9);
      for (int ii = 0; ii < 9; ++ii) {
      all_virial_1[kk][ii] = recv_v[kk * 9 + ii];
      }
      }
      double all_v_min = numeric_limits<double>::max(), all_v_max = 0,
      all_v_avg = 0;
      if (rank == 0) {
      deep_pot_model_devi.compute_avg(avg_virial, all_virial_1);
      deep_pot_model_devi.compute_std(std_virial, avg_virial, all_virial_1,
      1);
      if (out_rel_v == 1) {
      deep_pot_model_devi.compute_relative_std(std_virial, avg_virial,
      eps_v, 1);
      }
      for (int ii = 0; ii < 9; ++ii) {
      if (std_virial[ii] > all_v_max) {
      all_v_max = std_virial[ii];
      }
      if (std_virial[ii] < all_v_min) {
      all_v_min = std_virial[ii];
      }
      all_v_avg += std_virial[ii] * std_virial[ii];
      }
      all_v_avg = sqrt(all_v_avg / 9);
      }
      if (rank == 0) {
      all_v_max *= ener_unit_cvt_factor;
      all_v_min *= ener_unit_cvt_factor;
      all_v_avg *= ener_unit_cvt_factor;
      all_f_max *= force_unit_cvt_factor;
      all_f_min *= force_unit_cvt_factor;
      all_f_avg *= force_unit_cvt_factor;
      fp << setw(12) << update->ntimestep << " " << setw(18) << all_v_max
        << " " << setw(18) << all_v_min << " " << setw(18) << all_v_avg
        << " " << setw(18) << all_f_max << " " << setw(18) << all_f_min
        << " " << setw(18) << all_f_avg;
    }
    if (out_each == 1) {
      vector<double> std_f_all(atom->natoms);
      // Gather std_f and tags
      tagint *tag = atom->tag;
      int nprocs = comm->nprocs;
      // Grow arrays if necessary
      if (atom->natoms > stdf_comm_buff_size) {
        stdf_comm_buff_size = atom->natoms;
        memory->destroy(stdfsend);
        memory->destroy(stdfrecv);
        memory->destroy(tagsend);
        memory->destroy(tagrecv);
        memory->create(stdfsend, stdf_comm_buff_size, "deepmd:stdfsendall");
        memory->create(stdfrecv, stdf_comm_buff_size, "deepmd:stdfrecvall");
        memory->create(tagsend, stdf_comm_buff_size, "deepmd:tagsendall");
        memory->create(tagrecv, stdf_comm_buff_size, "deepmd:tagrecvall");
      }
      for (int ii = 0; ii < nlocal; ii++) {
        tagsend[ii] = tag[ii];
        stdfsend[ii] = std_f[ii];
      }
      MPI_Gather(&nlocal, 1, MPI_INT, counts, 1, MPI_INT, 0, world);
      displacements[0] = 0;
      for (int ii = 0; ii < nprocs - 1; ii++) {
        displacements[ii + 1] = displacements[ii] + counts[ii];
      }
      MPI_Gatherv(tagsend, nlocal, MPI_CAC_TAGINT, tagrecv, counts,
          displacements, MPI_CAC_TAGINT, 0, world);
      MPI_Gatherv(stdfsend, nlocal, MPI_DOUBLE, stdfrecv, counts,
          displacements, MPI_DOUBLE, 0, world);
      if (rank == 0) {
        for (int dd = 0; dd < atom->natoms; ++dd) {
          std_f_all[tagrecv[dd] - 1] = stdfrecv[dd] * force_unit_cvt_factor;
        }
        for (int dd = 0; dd < atom->natoms; ++dd) {
          fp << " " << setw(18) << std_f_all[dd];
        }
      }
    }
    if (rank == 0) {
      fp << endl;
    }
    }
    */
    } else {
      error->all(FLERR, "unknown computational branch");
    }
  } else {
    if (numb_models == 1) {
      try {
        deep_pot.compute(dener, dforce, dvirial, dcoord, dtype, dbox);
      } catch (deepmd_compat::deepmd_exception &e) {
        error->one(FLERR, e.what());
      }
    } else {
      error->all(FLERR, "Serial version does not support model devi");
    }
  }

  // get force
  // only extract force for atoms/gauss points
  
  for (int ideep = 0; ideep < nall; ideep++) {
    int i = index_map_deep2cac[ideep][0];
    int iindex = index_map_deep2cac[ideep][1];
    double *fptr;

    if (iindex == -1) {
      fptr = f[i];
    } else if (iindex >= 0 && iindex < element->maxapc * element->maxucell) {
      int ietype = etype[i];
      int iapc = apc[ietype];
      int ibasis = iindex % iapc;
      int iucell = iindex / iapc;
      int igcell = u2g[ietype][iucell];
      if (igcell >= 0) {
        fptr = gaussf[i][ibasis][igcell];
      } else continue;
    } else continue;
    for (int dd = 0; dd < 3; ++dd) {
      fptr[dd] += scale[1][1] * dforce[3 * ideep + dd] * force_unit_cvt_factor;
    }

  }

  // accumulate energy and virial
  
  if (eflag_global) {
    eng_vdwl += scale[1][1] * dener * ener_unit_cvt_factor;
  }
  if (vflag_global) {
    virial[0] += 1.0 * dvirial[0] * scale[1][1] * ener_unit_cvt_factor;
    virial[1] += 1.0 * dvirial[4] * scale[1][1] * ener_unit_cvt_factor;
    virial[2] += 1.0 * dvirial[8] * scale[1][1] * ener_unit_cvt_factor;
    virial[3] += 1.0 * dvirial[3] * scale[1][1] * ener_unit_cvt_factor;
    virial[4] += 1.0 * dvirial[6] * scale[1][1] * ener_unit_cvt_factor;
    virial[5] += 1.0 * dvirial[7] * scale[1][1] * ener_unit_cvt_factor;
  }
}

/* ---------------------------------------------------------------------- */

static bool is_key(const string &input) {
  vector<string> keys;
  //keys.push_back("out_freq");
  //keys.push_back("out_file");
  keys.push_back("fparam");
  keys.push_back("aparam");
  keys.push_back("fparam_from_compute");
  keys.push_back("aparam_from_compute");
  //keys.push_back("ttm");
  //keys.push_back("atomic");
  //keys.push_back("relative");
  //keys.push_back("relative_v");
  //keys.push_back("virtual_len");
  //keys.push_back("spin_norm");

  for (int ii = 0; ii < keys.size(); ++ii) {
    if (input == keys[ii]) {
      return true;
    }
  }
  return false;
}

/* ---------------------------------------------------------------------- */

void PairDeepMD::settings(int narg, char **arg) {

  if (narg <= 0) {
    error->all(FLERR, "Illegal pair_style command");
  }

  vector<string> models;
  int iarg = 0;
  while (iarg < narg) {
    if (is_key(arg[iarg])) {
      break;
    }
    iarg++;
  }
  for (int ii = 0; ii < iarg; ++ii) {
    models.push_back(arg[ii]);
  }

  backend_type = deepmd::get_backend(models[0]); 
  numb_models = models.size();
  if (numb_models == 1) {
    try {
      deep_pot.init(arg[0], get_node_rank(), get_file_content(arg[0]));
    } catch (deepmd_compat::deepmd_exception &e) {
      error->one(FLERR, e.what());
    }
    cutoff = deep_pot.cutoff() * dist_unit_cvt_factor;
    numb_types = deep_pot.numb_types();
    numb_types_spin = deep_pot.numb_types_spin();
    dim_fparam = deep_pot.dim_fparam();
    dim_aparam = deep_pot.dim_aparam();
  } else {
    try {
      deep_pot.init(arg[0], get_node_rank(), get_file_content(arg[0]));
      deep_pot_model_devi.init(models, get_node_rank(),
          get_file_content(models));
    } catch (deepmd_compat::deepmd_exception &e) {
      error->one(FLERR, e.what());
    }
    cutoff = deep_pot_model_devi.cutoff() * dist_unit_cvt_factor;
    numb_types = deep_pot_model_devi.numb_types();
    numb_types_spin = deep_pot_model_devi.numb_types_spin();
    dim_fparam = deep_pot_model_devi.dim_fparam();
    dim_aparam = deep_pot_model_devi.dim_aparam();
    assert(cutoff == deep_pot.cutoff() * dist_unit_cvt_factor);
    assert(numb_types == deep_pot.numb_types());
    assert(numb_types_spin == deep_pot.numb_types_spin());
    assert(dim_fparam == deep_pot.dim_fparam());
    assert(dim_aparam == deep_pot.dim_aparam());
  }

  //out_freq = 100;
  out_freq = 0; // temporarily turn off exporting model deviation
  out_file = "model_devi.out";
  out_each = 0;
  out_rel = 0;
  eps = 0.;
  fparam.clear();
  aparam.clear();
  while (iarg < narg) {
    if (!is_key(arg[iarg])) {
      error->all(FLERR,
          "Illegal pair_style command\nwrong number of parameters\n");
    }
    if (string(arg[iarg]) == string("out_freq")) {
      if (iarg + 1 >= narg) 
        error->all(FLERR, "Illegal out_freq, not provided");
      out_freq = universe->inumeric(FLERR, (arg[iarg + 1]));
      if (out_freq < 0) 
        error->all(FLERR, "Illegal out_freq, should be >= 0");
      iarg += 2;
    } else if (string(arg[iarg]) == string("out_file")) {
      if (iarg + 1 >= narg) 
        error->all(FLERR, "Illegal out_file, not provided");
      out_file = string(arg[iarg + 1]);
      iarg += 2;
    } else if (string(arg[iarg]) == string("fparam")) {
      for (int ii = 0; ii < dim_fparam; ++ii) {
        if (iarg + 1 + ii >= narg || is_key(arg[iarg + 1 + ii])) {
          char tmp[1024];
          sprintf(tmp, "Illegal fparam, the dimension should be %d",
              dim_fparam);
          error->all(FLERR, tmp);
        }
        fparam.push_back(atof(arg[iarg + 1 + ii]));
      }
      iarg += 1 + dim_fparam;
    } else if (string(arg[iarg]) == string("aparam")) {
      for (int ii = 0; ii < dim_aparam; ++ii) {
        if (iarg + 1 + ii >= narg || is_key(arg[iarg + 1 + ii])) {
          char tmp[1024];
          sprintf(tmp, "Illegal aparam, the dimension should be %d",
              dim_aparam);
          error->all(FLERR, tmp);
        }
        aparam.push_back(atof(arg[iarg + 1 + ii]));
      }
      iarg += 1 + dim_aparam;

      //    } else if (string(arg[iarg]) == string("ttm")) {
      //#ifdef USE_TTM
      //      for (int ii = 0; ii < 1; ++ii) {
      //        if (iarg + 1 + ii >= narg || is_key(arg[iarg + 1 + ii])) {
      //          error->all(FLERR, "invalid ttm key: should be ttm ttm_fix_id(str)");
      //        }
      //      }
      //      do_ttm = true;
      //      ttm_fix_id = arg[iarg + 1];
      //      iarg += 1 + 1;
      //#else
      //      error->all(FLERR,
      //          "The deepmd-kit was compiled without support for TTM, please "
      //          "rebuild it with CAC version >=20210831");
      //#endif

  }

  ///////////////////////////////////////////////
  // pair_style     deepmd cp.pb fparam_from_compute TEMP
  // compute        TEMP all temp
  //////////////////////////////////////////////
    else if (string(arg[iarg]) == string("fparam_from_compute")) {
      for (int ii = 0; ii < 1; ++ii) {
        if (iarg + 1 + ii >= narg || is_key(arg[iarg + 1 + ii])) {
          error->all(FLERR,
              "invalid fparam_from_compute key: should be "
              "fparam_from_compute compute_fparam_id(str)");
        }
      }
      do_compute_fparam = true;
      compute_fparam_id = arg[iarg + 1];
      iarg += 1 + 1;
    } else if (string(arg[iarg]) == string("aparam_from_compute")) {
      for (int ii = 0; ii < 1; ++ii) {
        if (iarg + 1 + ii >= narg || is_key(arg[iarg + 1 + ii])) {
          error->all(FLERR,
              "invalid aparam_from_compute key: should be "
              "aparam_from_compute compute_aparam_id(str)");
        }
      }
      do_compute_aparam = true;
      compute_aparam_id = arg[iarg + 1];
      iarg += 1 + 1;
      //} else if (string(arg[iarg]) == string("atomic")) {
      //  out_each = 1;
      //  iarg += 1;
      //} else if (string(arg[iarg]) == string("relative")) {
      //  out_rel = 1;
      //  eps = atof(arg[iarg + 1]) / ener_unit_cvt_factor;
      //  iarg += 2;
      //} else if (string(arg[iarg]) == string("relative_v")) {
      //  out_rel_v = 1;
      //  eps_v = atof(arg[iarg + 1]) / ener_unit_cvt_factor;
      //  iarg += 2;
      //} else if (string(arg[iarg]) == string("virtual_len")) {
      //  virtual_len.resize(numb_types_spin);
      //  for (int ii = 0; ii < numb_types_spin; ++ii) {
      //    virtual_len[ii] = atof(arg[iarg + ii + 1]);
      //  }
      //  iarg += numb_types_spin + 1;
      //} else if (string(arg[iarg]) == string("spin_norm")) {
      //  spin_norm.resize(numb_types_spin);
      //  for (int ii = 0; ii < numb_types_spin; ++ii) {
      //    spin_norm[ii] = atof(arg[iarg + ii + 1]);
      //  }
      //  iarg += numb_types_spin + 1;
  }
  }

  if ((int)do_compute_aparam + (int)(aparam.size() > 0) > 1) {
    error->all(FLERR,
        "aparam, aparam_from_compute should NOT be set "
        "simultaneously");
  }
  if (do_compute_fparam && fparam.size() > 0) {
    error->all(
        FLERR,
        "fparam and fparam_from_compute should NOT be set simultaneously");
  }

  if (comm->me == 0) {
    if (numb_models > 1 && out_freq > 0) {
      if (!is_restart) {
        fp.open(out_file);
        fp << scientific;
        fp << "#" << setw(12 - 1) << "step" << setw(18 + 1) << "max_devi_v"
          << setw(18 + 1) << "min_devi_v" << setw(18 + 1) << "avg_devi_v"
          << setw(18 + 1) << "max_devi_f" << setw(18 + 1) << "min_devi_f"
          << setw(18 + 1) << "avg_devi_f";
        if (out_each) {
          // at this time, we don't know how many atoms
          fp << setw(18 + 1) << "atm_devi_f(N)";
        }
        fp << endl;
      } else {
        fp.open(out_file, ofstream::out | ofstream::app);
        fp << scientific;
      }
    }
    string pre = "  ";
    cout << pre << ">>> Info of model(s):" << endl
      << pre << "using " << setw(3) << numb_models << " model(s): ";
    if (narg == 1) {
      cout << arg[0] << " ";
    } else {
      for (int ii = 0; ii < models.size(); ++ii) {
        cout << models[ii] << " ";
      }
    }
    cout << endl
      << pre << "rcut in model:      " << cutoff << endl
      << pre << "ntypes in model:    " << numb_types << endl;
    if (fparam.size() > 0) {
      cout << pre << "using fparam(s):    ";
      for (int ii = 0; ii < dim_fparam; ++ii) {
        cout << fparam[ii] << "  ";
      }
      cout << endl;
    }
    if (do_compute_fparam) {
      cout << pre << "using compute id (fparam):      ";
      cout << compute_fparam_id << "  " << endl;
    }
    if (do_compute_aparam) {
      cout << pre << "using compute id (aparam):      ";
      cout << compute_aparam_id << "  " << endl;
    }
    if (aparam.size() > 0) {
      cout << pre << "using aparam(s):    ";
      for (int ii = 0; ii < aparam.size(); ++ii) {
        cout << aparam[ii] << "  ";
      }
      cout << endl;
    }
    /*
       if (do_ttm) {
       cout << pre << "using ttm fix:      ";
       cout << ttm_fix_id << "  ";
       if (dim_fparam > 0) {
       cout << "(fparam)" << endl;
       } else if (dim_aparam > 0) {
       cout << "(aparam)" << endl;
       }
       }
       */
  }

  comm_atom_reverse = numb_models * 3;
  comm_elem_reverse = element->maxnpe * element->maxapc * comm_atom_reverse;
  all_force.resize(numb_models);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairDeepMD::coeff(int narg, char **arg) {
  if (!allocated) {
    allocate();
  }
  int n = atom->ntypes;
  int ilo, ihi, jlo, jhi;
  ilo = 0;
  jlo = 0;
  ihi = n;
  jhi = n;
  if (narg >= 2) {
    universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
    universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);
    if (ilo != 1 || jlo != 1 || ihi != n || jhi != n) {
      error->all(FLERR,
          "deepmd requires that the scale should be set to all atom "
          "types, i.e. pair_coeff * *.");
    }
  }
  if (narg <= 2) {
    type_idx_map.resize(n);
    for (int ii = 0; ii < n; ++ii) {
      type_idx_map[ii] = ii;
    }
  } else {
    int iarg = 2;

    // type_map is a list of strings with undetermined length
    // note: although we have numb_types from the model, we do not require
    // the number of types in the system matches that in the model
    vector<string> type_map;
    string type_map_str;
    deep_pot.get_type_map(type_map_str);
    // convert the string to a vector of strings
    istringstream iss(type_map_str);
    string type_name;
    while (iss >> type_name) {
      type_map.push_back(type_name);
    }

    type_idx_map.clear();
    type_names.clear();
    while (iarg < narg) {
      string type_name = arg[iarg];
      type_names.push_back(type_name);
      bool found_element = false;
      for (int ii = 0; ii < type_map.size(); ++ii) {
        if (type_map[ii] == type_name) {
          type_idx_map.push_back(ii);
          found_element = true;
          break;
        }
      }
      if (!found_element && "NULL" == type_name) {
        type_idx_map.push_back(type_map.size());  // ghost type
        found_element = true;
      }
      if (!found_element) {
        char error_msg[1024];
        sprintf(error_msg, "pair_coeff: element %s not found in the model", type_name);
        error->all(FLERR, error_msg);
      }
      iarg += 1;
    }
    numb_types = type_idx_map.size();
    if (numb_types < n) {
      type_idx_map.resize(n);
      for (int ii = numb_types; ii < n; ++ii) {
        type_idx_map[ii] = -1;
      }
    }
  }
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      scale[i][j] = 1.0;
      if (i > numb_types || j > numb_types) {
        char warning_msg[1024];
        sprintf(warning_msg,
            "Interaction between types %d and %d is set with deepmd, but "
            "will be ignored.\n Deepmd model has only %d types, it only "
            "computes the mulitbody interaction of types: 1-%d.",
            i, j, numb_types, numb_types);
        error->warning(FLERR, warning_msg);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairDeepMD::pack_atom_reverse_comm(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  //if (atom->sp_flag) {
  //  error->all(
  //      FLERR,
  //      "Pair style 'deepmd' does not support spin atoms, please use pair "
  //      "style 'deepspin' instead.");
  //} else {
  for (i = first; i < last; i++) {
    for (int dd = 0; dd < numb_models; ++dd) {
      buf[m++] = all_force[dd][3 * i + 0];
      buf[m++] = all_force[dd][3 * i + 1];
      buf[m++] = all_force[dd][3 * i + 2];
    }
  }
  //}
  return m;
}

/* ---------------------------------------------------------------------- */

void PairDeepMD::unpack_atom_reverse_comm(int n, int *list, double *buf) {
  int i, j, m;

  m = 0;
  //if (atom->sp_flag) {
  //  error->all(
  //      FLERR,
  //      "Pair style 'deepmd' does not support spin atoms, please use pair "
  //      "style 'deepspin' instead.");
  //} else {
  for (i = 0; i < n; i++) {
    j = list[i];
    for (int dd = 0; dd < numb_models; ++dd) {
      all_force[dd][3 * j + 0] += buf[m++];
      all_force[dd][3 * j + 1] += buf[m++];
      all_force[dd][3 * j + 2] += buf[m++];
    }
  }
  //}
}

void PairDeepMD::print_summary(const string pre) const {
  if (comm->me == 0) {
    // capture cout to a string, then call CAC's utils::logmesg
    // https://stackoverflow.com/a/4043813/9567349
    stringstream buffer;
    streambuf *sbuf = cout.rdbuf();
    cout.rdbuf(buffer.rdbuf());

    cout << "Summary of lammps deepmd module ..." << endl;
    cout << pre << ">>> Info of deepmd-kit:" << endl;
    deep_pot.print_summary(pre);
    cout << pre << ">>> Info of lammps module:" << endl;
    cout << pre << "use deepmd-kit at:  " << STR_DEEPMD_ROOT << endl;
    cout << pre << "source:             " << STR_GIT_SUMM << endl;
    cout << pre << "source branch:      " << STR_GIT_BRANCH << endl;
    cout << pre << "source commit:      " << STR_GIT_HASH << endl;
    cout << pre << "source commit at:   " << STR_GIT_DATE << endl;
    cout << pre << "build with inc:     " << STR_BACKEND_INCLUDE_DIRS << endl;
    cout << pre << "build with lib:     " << STR_BACKEND_LIBRARY_PATH << endl;

    cout.rdbuf(sbuf);
    //utils::logmesg(cac, buffer.str());
  }
}

/* --------------------------------------------------------------------
   Reconstruct input list for deepMD in LAMMPS structure
   Include all atoms and virtual atoms involving in the force calculation
   (neighbors of neighbors)
   The neighbor list are extracted from CAC neighbor list and not built from scratch
   This is done every timestep
   --------------------------------------------------------------------*/
void PairDeepMD::reconstruct_input_nlist(
    vector<double> &dcoord, vector<int> &dtype,
    unordered_map<vector<int>, int, VectorHasher<int>> &index_map_cac2deep,
    unordered_map<int, vector<int>> &index_map_deep2cac
    )
{
  int i, ii, iindex;
  int jj, j, jindex, jnum, *jlist, *jindexlist;
  int ictype, ietype, iapc, iucell, ibasis, igcell, inode;
  int jctype, jetype, japc, jucell, jbasis, jgcell, jnode;

  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;
  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int neall = nelocal + neghost;

  int *atype = atom->type;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *apc = element->apc;
  int **u2g = element->u2g;
  int **g2u = element->g2u;
  double **x = atom->x;
  double ****nodex = element->nodex;

  // scan through the neighbor list for the list of all atoms/virtual atoms involved
  // they must all be included in the input list due to restriction from deepmd backend

  // Add neighbor list of real atoms and gauss points and also their neighbors. 
  // First part: construct DeepMD arrays (dcoord, dtype) and the 2-way mapping between atoms 
  //   and virtual atoms indices in CAC to the index in dcoord array
  // nall after the loop will be the total number of atoms (local + ghost) to pass to DeepMD backend
  // The definition of "ghost" and "local" atoms are tweaked a little here
  // It's the same for atoms/gauss points ("local" if owned and "ghost" if not)
  // But for nva ("1st" neighbor) atoms, they are all treated as: "local", 
  //   all (and only) nva's "ghost" atom/gauss point neighbors are duplicated to "proxy ghost" to avoid force tally on those neighbors
  //   "proxy ghost" has their new index as:
  //   - For gauss points: iindex = iindex + element->maxucell * element->maxapc 
  //   - For atoms:  iindex = -2
  // All "2nd" neighbor atoms (not atom, not gauss, not nva) are treated as "ghost" 
  //   since it is just to fill the deepmd list to be neighbor of "local" atoms
  //   even if it is owned
  // Keep coordinates in dcoord in original values to check rsq <= cutoff in second part
  // update dcoord accordingly after neighbor list is built

  nall = 0;
  unordered_map<vector<int>, int, VectorHasher<int>> gauss_index_map;
  unordered_set<vector<int>, VectorHasher<int>> short_nva_list;
  unordered_set<vector<int>, VectorHasher<int>> ghost_list;
  gauss_index_map.reserve(list->ginum);
  short_nva_list.reserve(list->nvainum);
  ghost_list.reserve(nmax);
  double rsq, delx, dely, delz;
  double cutsq = cutoff * cutoff;

  for (ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    iindex = list->iindexlist[ii]; 

    double icoord[3] = {0.0};
    if (iindex < 0) {
      icoord[0] = x[i][0];     
      icoord[1] = x[i][1];     
      icoord[2] = x[i][2];     
      ictype = type_idx_map[atype[i] - 1];
    } else {

      // note: this index is for gauss point, NOT unit cell 
      // i.e. iindex = igauss * iapc + ibasis

      ietype = etype[i];
      iapc = apc[ietype];
      ibasis = iindex % iapc;
      igcell = iindex / iapc;
      iucell = g2u[ietype][igcell];
      ictype = type_idx_map[ctype[i][ibasis] - 1];
      element->evec->interpolate(icoord, nodex, i, ibasis, iucell, 3);

      // change iindex to unitcell-based index to be consistent with the rest of the code

      iindex = iucell * iapc + ibasis;

      // build index mapping for gauss point

      gauss_index_map[{i, iindex}] = ii;
    }

    // add i to the dcoord vector and stored the index mapping
    // duplicate check is handled in insert_to_deepmd_arrays()
    // all atoms/gauss points from list->ilist are already "local"

    nall += insert_to_deepmd_arrays(i, iindex, icoord, ictype, nall, 
        index_map_cac2deep, index_map_deep2cac, dcoord, dtype);

    // loop through i's neighbors and add them to the reconstructed list 

    jnum = list->numneigh[ii];
    jlist = list->firstneigh[ii];
    jindexlist = list->firstneighindex[ii];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];

      double jcoord[3] = {0.0};
      int nloc;
      jgcell = -1;
      if (jindex < 0) {
        jcoord[0] = x[j][0];     
        jcoord[1] = x[j][1];     
        jcoord[2] = x[j][2];     
        jctype = type_idx_map[atype[j] - 1];
        nloc = nalocal;
      } else {
        jetype = etype[j];
        japc = apc[jetype];
        jbasis = jindex % japc;
        jucell = jindex / japc;
        jgcell = u2g[jetype][jucell];
        jctype = type_idx_map[ctype[j][jbasis] - 1];
        element->evec->interpolate(jcoord, nodex, j, jbasis, jucell, 3);
        nloc = nelocal;
      }

      // skip neighbor outside force cutoff

      delx = jcoord[0] - icoord[0];
      dely = jcoord[1] - icoord[1];
      delz = jcoord[2] - icoord[2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq > cutsq) continue;

      // if j is atom or gauss point,

      if (jindex < 0 || jgcell >= 0) {

        //   check if j is ghost and add to deepmd array accordingly
        if (j >= nloc) {
          // check for duplicate before adding
          if (ghost_list.count({j, jindex}) == 0)
            ghost_list.insert({j, jindex});
        }
      } 

      // j is nva -> add to short_nva_list

      else {
        if (list->va2nvalist.count({j, jindex}) == 0) error->one(FLERR,"TEST"); // sanity check
        int jnva = list->va2nvalist[{j, jindex}];

        //        if (list->va2nvalist[j][jbasis][jucel] < 0) error->one(FLERR,"TEST"); // sanity check
        //        int jnva = list->va2nvalist[j][jbasis][jucel];
        // check for duplicate (skip if it was already added)
        if (short_nva_list.count({jnva, j, jindex}) == 0)
          short_nva_list.insert({jnva, j, jindex}); 
      }
    }
  }

  // loop through short_nva_list and loop through its neighbor list, 

  for (const vector<int> &short_nva : short_nva_list) {
    int inva = short_nva[0];
    i = short_nva[1];
    iindex = short_nva[2];

    ietype = etype[i];
    iapc = apc[ietype];
    ibasis = iindex % iapc;
    iucell = iindex / iapc;
    ictype = type_idx_map[ctype[i][ibasis] - 1];

    double icoord[3] = {0.0};
    element->evec->interpolate(icoord, nodex, i, ibasis, iucell, 3);

    nall += insert_to_deepmd_arrays(i, iindex, icoord, ictype, nall, 
        index_map_cac2deep, index_map_deep2cac, dcoord, dtype);

    jnum = list->nvanumneigh[inva];
    jlist = list->nvafirstneigh[inva];
    jindexlist = list->nvafirstneighindex[inva];

    //int nn = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      double jcoord[3] = {0.0};
      int nloc;
      if (jindex < 0) {
        jcoord[0] = x[j][0];     
        jcoord[1] = x[j][1];     
        jcoord[2] = x[j][2];     
        jgcell = -1;
        nloc = nalocal;
      } else {
        jetype = etype[j];
        japc = apc[jetype];
        jbasis = jindex % japc;
        jucell = jindex / japc;
        jgcell = u2g[jetype][jucell];
        element->evec->interpolate(jcoord, nodex, j, jbasis, jucell, 3);
        nloc = nelocal;
      }

      // skip neighbor outside force cutoff

      delx = icoord[0] - jcoord[0];
      dely = icoord[1] - jcoord[1];
      delz = icoord[2] - jcoord[2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq) continue;
      //nn++;

      int jnva = list->va2nvalist[{j, jindex}];

      // if j is atom/gauss point/short_nva 
      if (jindex < 0 || jgcell >= 0 || short_nva_list.count({jnva, j, jindex}) == 1) {

        // if j is ghost atom/gauss point:
        if (j >= nloc && (jindex < 0 || jgcell >= 0)) {

          //  - add j to "proxy ghost" list (short_nva is handled in the outer i loop)
          if (jindex >= 0) jindex += element->maxucell * element->maxapc;
          else jindex = -2;
          if (ghost_list.count({j, jindex}) == 0)
            ghost_list.insert({j, jindex});
        }
      } 

      // if not, then j is "2nd" neighbor and add it to ghost list
      else {
        if (ghost_list.count({j, jindex}) == 0)
          ghost_list.insert({j, jindex});
      }
    }
  }

  int nlocal = nall;

  // insert ghost atoms to DeepMD arrays
  // do this after to ensure ghost are stored after local
  for (const vector<int> &ghost : ghost_list) {
    i = ghost[0];
    iindex = ghost[1];
    double icoord[3] = {0.0};
    if (iindex < 0) {
      icoord[0] = x[i][0];     
      icoord[1] = x[i][1];     
      icoord[2] = x[i][2];     
      ictype = type_idx_map[atype[i] - 1];
    } else {
      ietype = etype[i];
      iapc = apc[ietype];
      int tmp = iindex % (element->maxapc * element->maxucell);
      ibasis = tmp % iapc;
      iucell = tmp / iapc;
      ictype = type_idx_map[ctype[i][ibasis] - 1];
      element->evec->interpolate(icoord, nodex, i, ibasis, iucell, 3);
    }

    // add i to the dcoord vector and stored the index mapping

    nall += insert_to_deepmd_arrays(i, iindex, icoord, ictype, nall, 
        index_map_cac2deep, index_map_deep2cac, dcoord, dtype);
  }

  // Second part: construct the neighbor list only for local atoms in deepmd arrays
  // only store neighbors within cutoff to save computational cost in DeepMD

  ipage->reset(); 
  inum = 0;
  int n, ideep, jdeep;
  int *neighptr;
  for (ideep = 0; ideep < nlocal; ideep++) {

    //----------- Sanity check, can be removed later -----
    if (index_map_deep2cac.find(ideep) == index_map_deep2cac.end())
      error->one(FLERR,"Incorrect index map");
    //----------------------------------------------------
    i = index_map_deep2cac[ideep][0];
    iindex = index_map_deep2cac[ideep][1];
    int inva = -1;

    // extract neighbor list data: jnum, jlist, jindexlist

    igcell = -1;
    if (iindex < 0) {
      jnum = list->numneigh[i];
      jlist = list->firstneigh[i];
      jindexlist = list->firstneighindex[i];
    } else {
      //check if i is gauss point or nva
      ietype = etype[i];
      iapc = apc[ietype];
      iucell = iindex / iapc;
      igcell = u2g[ietype][iucell];
      if (igcell >= 0) {
        if (gauss_index_map.find({i, iindex}) == gauss_index_map.end()) error->one(FLERR, "TEST"); // sanity check
        ii = gauss_index_map[{i, iindex}];
        jnum = list->numneigh[ii];
        jlist = list->firstneigh[ii];
        jindexlist = list->firstneighindex[ii];
      } else {
        if (list->va2nvalist.count({i, iindex}) == 0) error->one(FLERR,"TEST"); // sanity check
        inva = list->va2nvalist[{i, iindex}];
        if (short_nva_list.count({inva, i, iindex}) == 1) {
          jnum = list->nvanumneigh[inva]; 
          jlist = list->nvafirstneigh[inva];
          jindexlist = list->nvafirstneighindex[inva];
        } 
        // sanity check, can be removed lated
        else {
          error->one(FLERR,"TEST");
        }
      }
    }

    n = 0;
    neighptr = ipage->vget(); 

    for (jj = 0; jj < jnum ; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];

      int nloc;
      if (jindex < 0) {
        jgcell = -1;
        nloc = nalocal;
      } else {
        jetype = etype[j];
        japc = apc[jetype];
        jbasis = jindex % japc;
        jucell = jindex / japc;
        jgcell = u2g[jetype][jucell];
        nloc = nelocal;
      }

      // if i is nva and j is ghost atom/gauss point change jindex to select the "proxy ghost"

      if (inva >= 0 && j >= nloc) {
        if (jindex < 0) {
          jindex = -2;
        } else if (jgcell >= 0) {
          jindex += element->maxucell * element->maxapc;
        }
      }
      auto it = index_map_cac2deep.find({j, jindex});
      if (it != index_map_cac2deep.end()) {
        jdeep = it->second;
        if (jdeep >= nall) error->one(FLERR, "TEST"); // sanity check
        delx = dcoord[ideep * 3 + 0] - dcoord[jdeep * 3 + 0];
        dely = dcoord[ideep * 3 + 1] - dcoord[jdeep * 3 + 1];
        delz = dcoord[ideep * 3 + 2] - dcoord[jdeep * 3 + 2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq <= cutsq) 
          neighptr[n++] = jdeep;
      }
    }

    if (inum >= nmax) error->one(FLERR, "TEST"); // sanity check
    firstneigh[inum] = neighptr;
    numneigh[inum] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    ilist[inum] = ideep;
    inum++;
  }

  //dump_debug(nall, dcoord, index_map_deep2cac, short_nva_list);
  // convert dcoord

  for (int ii = 0; ii < nall; ++ii) {
    for (int dd = 0; dd < 3; ++dd) {
      dcoord[ii * 3 + dd] =
        (dcoord[ii * 3 + dd] - domain->boxlo[dd]) / dist_unit_cvt_factor;
    }
  }
}

/* --------------------------------------------------------------------
   insert atom i into the deepmd_arrays
   return 0 if i already exists in the list
   return 1 if inserted successfully
   --------------------------------------------------------------------*/

int PairDeepMD::insert_to_deepmd_arrays(int i, int iindex, double *coord, int itype, int nall,
    unordered_map<vector<int>, int, VectorHasher<int>> &index_map_cac2deep,
    unordered_map<int, vector<int>> &index_map_deep2cac,
    vector<double> &dcoord, vector<int> &dtype)
{
  vector<int> index_pair = {i, iindex};

  // return if this atom is already added

  if (index_map_cac2deep.find(index_pair) != index_map_cac2deep.end()) {
    error->all(FLERR,"TEST"); //sanity check
    return 0;
  }

  index_map_deep2cac[nall] = index_pair;
  index_map_cac2deep[index_pair] = nall;
  dcoord.push_back(coord[0]);
  dcoord.push_back(coord[1]);
  dcoord.push_back(coord[2]);
  dtype.push_back(itype);

  return 1;
}

/*  ----------------------------------------------------------------------
    same as in pair.cpp 
    the only difference is: 
    eflag_atom = eflag
    vflag_atom = vflag
    -------------------------------------------------------------------------  */

void PairDeepMD::ev_setup(int eflag, int vflag, int alloc)
{
  int maxnpe = element->maxnpe;
  int maxapc = element->maxapc;
  int i, j, k, n;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;

  // need eatom to recaclulate global energy
  eflag_atom = eflag;

  vflag_global = vflag & VIRIAL_PAIR;

  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 1) vflag_global = 1;
  vflag_fdotr = 0;
  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 0) vflag_fdotr = 1;
  vflag_atom = vflag & VIRIAL_ATOM;
  if (vflag & VIRIAL_CENTROID && centroidstressflag != CENTROID_AVAIL) vflag_atom = 1;
  cvflag_atom = vflag;
  vflag_either = vflag_global || vflag_atom || cvflag_atom;

  evflag = eflag_either || vflag_either;

  // reallocate per-atom/per-node arrays if necessary

  if (eflag_atom) {
    if (atom->nmax > maxeatom) {
      maxeatom = atom->nmax;
      if (alloc) {
        memory->destroy(eatom);
        memory->create(eatom, comm->nthreads * maxeatom, "pair:eatom");
      }
    }
    if (element->nmax > maxeelem) {
      maxeelem = element->nmax;
      if (alloc) {
        memory->destroy(enode);
        memory->create(enode, comm->nthreads * maxeelem, maxapc, maxnpe, "pair:enode");
      }
    }
  }

  if (vflag_atom) {
    if (atom->nmax > maxvatom) {
      maxvatom = atom->nmax;
      if (alloc) {
        memory->destroy(vatom);
        memory->create(vatom, comm->nthreads * maxvatom, 6, "pair:vatom");
      }
    }
    if (element->nmax > maxvelem) {
      maxvelem = element->nmax;
      if (alloc) {
        memory->destroy(vnode);
        memory->create(vnode, comm->nthreads * maxvelem, maxapc, maxnpe, 6, "pair:vnode");
      }
    }
  }

  if (cvflag_atom) {
    if (atom->nmax > maxcvatom) {
      maxcvatom = atom->nmax;
      if (alloc) {
        memory->destroy(cvatom);
        memory->create(cvatom, comm->nthreads * maxcvatom, 9, "pair:cvatom");
      }
    }
    if (element->nmax > maxcvelem) {
      maxcvelem = element->nmax;
      if (alloc) {
        memory->destroy(cvnode);
        memory->create(cvnode, comm->nthreads * maxcvelem, maxapc, maxnpe, 9, "pair:cvnode");
      }
    }
  }

  // zero accumulators
  // use force->newton instead of newton_pair
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info

  if (eflag_global) eng_vdwl = eng_coul = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) 
          enode[i][j][k] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) {
          vnode[i][j][k][0] = 0.0;
          vnode[i][j][k][1] = 0.0;
          vnode[i][j][k][2] = 0.0;
          vnode[i][j][k][3] = 0.0;
          vnode[i][j][k][4] = 0.0;
          vnode[i][j][k][5] = 0.0;
        }
  }
  if (cvflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton) n += atom->nghost;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
    }
    n = element->nlocal;
    if (force->newton) n += element->nghost;
    for (i = 0; i < n; i++) 
      for (j = 0; j < maxapc; j++) 
        for (k = 0; k < maxnpe; k++) {
          cvnode[i][j][k][0] = 0.0;
          cvnode[i][j][k][1] = 0.0;
          cvnode[i][j][k][2] = 0.0;
          cvnode[i][j][k][3] = 0.0;
          cvnode[i][j][k][4] = 0.0;
          cvnode[i][j][k][5] = 0.0;
          cvnode[i][j][k][6] = 0.0;
          cvnode[i][j][k][7] = 0.0;
          cvnode[i][j][k][8] = 0.0;
        }
  }

  // if vflag_global = 2 and no element in simulation and pair::compute() calls virial_fdotr_compute() 
  // compute global virial via (F dot r) instead of via pairwise summation
  // unset other flags as appropriate

  if (vflag_global == 2 && no_virial_fdotr_compute == 0 && element->nelements == 0) {
    vflag_fdotr = 1;
  } else vflag_fdotr = 0;
}

void PairDeepMD::dump_debug(int natoms, vector<double> x
    , unordered_map<int, vector<int>> index_map_deep2cac, unordered_set<vector<int>, VectorHasher<int>> short_nva_list) 
{
  char file[128];
  sprintf(file,"dump_deep_atoms_proc_%d.%d.atom",comm->me,update->ntimestep);
  FILE *fp = fopen(file, "w");
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, BIGINT_FORMAT "\n", natoms);
  char boundstr[9];
  domain->boundary_string(boundstr);
  fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
  fprintf(fp, "ITEM: ATOMS id type x y z ghost owned eid aid i iindex iucell igauss ibasis numneigh inva \n");
  for (int ideep = 0; ideep < natoms; ideep++) {
    int i = index_map_deep2cac[ideep][0]; 
    int iindex = index_map_deep2cac[ideep][1]; 
    int ghost;
    int nlocal;
    int itype,iapc;
    int ibasis,igauss,iucell;
    int inva = 0;
    int eid, aid;
    int num;
    eid = aid = -1;
    if (iindex < 0) {
      nlocal = atom->nlocal;
      itype = atom->type[i];
      iucell = ibasis = igauss = -1;
      aid = atom->tag[i];
    } else {
      eid = element->tag[i];
      nlocal = element->nlocal;
      iapc = element->apc[element->etype[i]];
      int tmp = iindex % (element->maxapc * element->maxucell);
      ibasis = tmp % iapc;
      iucell = tmp / iapc;
      //if (iucell >= 1000 || ibasis >= 2) error->all(FLERR,"TEST"); 
      igauss = element->u2g[element->etype[i]][iucell];
      itype = element->ctype[i][ibasis];
      if (list->va2nvalist.count({i, tmp}) == 0) 
        inva = -1;
      else 
        inva = list->va2nvalist[{i, tmp}];
      if (short_nva_list.count({inva, i, tmp}) == 0)
        inva = -1;
    }
    if (ideep < inum) {
      ghost = 0;
      num = numneigh[ideep];
    } else {
      ghost = 1;
      num = 0;
    }
    fprintf(fp, "%d %d %g %g %g %d %d %d %d %d %d %d %d %d %d %d\n",ideep+1,itype
        ,x[3 * ideep]
        ,x[3 * ideep + 1]
        ,x[3 * ideep + 2]
        ,ghost,i < nlocal,eid,aid,i,iindex,iucell,igauss,ibasis,num,inva);
  }
  fclose(fp);
}

void PairDeepMD::dump_debug_neighbor(int ideep, vector<double> x
    , unordered_map<int, vector<int>> index_map_deep2cac, unordered_set<vector<int>, VectorHasher<int>> short_nva_list) 
{
  char file[128];
  sprintf(file,"dump_deep_atoms_neighbor_%d_proc_%d.%d.atom",ideep,comm->me,update->ntimestep);
  FILE *fp = fopen(file, "w");
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, BIGINT_FORMAT "\n", update->ntimestep);
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, BIGINT_FORMAT "\n", numneigh[ideep]+1);
  char boundstr[9];
  domain->boundary_string(boundstr);
  fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
  fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
  fprintf(fp, "ITEM: ATOMS id type x y z ideep numneigh i iindex iucell ibasis\n");

  int jnum = numneigh[ideep];
  int *jlist = firstneigh[ideep];
  int i = index_map_deep2cac[ideep][0]; 
  int iindex = index_map_deep2cac[ideep][1]; 

  fprintf(fp, "%d %d %g %g %g %d %d %d %d %d %d\n",+1,1
      ,x[3 * ideep]
      ,x[3 * ideep + 1]
      ,x[3 * ideep + 2],ideep,jnum,i,iindex,iindex/2,iindex%2);

  for (int jj = 0; jj < jnum; jj++) {
    int jdeep = jlist[jj];
    int nn = 0;
    int j = index_map_deep2cac[jdeep][0]; 
    int jindex = index_map_deep2cac[jdeep][1]; 

    if (jdeep < inum) nn = numneigh[jdeep];
    fprintf(fp, "%d %d %g %g %g %d %d %d %d %d %d\n"
        ,jj+2,2
        ,x[3 * jdeep]
        ,x[3 * jdeep + 1]
        ,x[3 * jdeep + 2]
        ,jdeep,nn,j,jindex,jindex/2,jindex%2
        );
  }
  fclose(fp);
}
