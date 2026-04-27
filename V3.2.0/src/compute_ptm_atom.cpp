// clang-format off
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <algorithm>
#include "compute_ptm_atom.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "update.h"
#include "universe.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"

// PTM library headers
#include "ptm_constants.h"
#include "ptm_functions.h"
#include "ptm_initialize_data.h"

using namespace CAC_NS;

// -------------------------------------------------------------------
// Data passed to the get_neighbours callback
// -------------------------------------------------------------------
// Neighbour positions are stored as absolute coordinates.
// Central particle: d.coordx/coordy/coordz
// Sorted neighbours (distance-ordered): d.delx/dely/delz[k] relative to central.
// For 2-shell structures the callback also needs each first-shell neighbour's
// own neighbourhood; stored in d.neigh_abs[k][m][0..2] (absolute positions)
// and d.neigh_c{x,y,z}[k] (the neighbour's own centre position).
//
// Particle IDs used inside a single ptm_index call:
//   central:          atom_id == central_id
//   neighbour k (0-based): atom_id == k+1   (1 .. ncollect)
//   never overlap because central_id != k+1 for reasonable particle counts.
// -------------------------------------------------------------------

#define CPTM_MAXNEAR (PTM_MAX_INPUT_POINTS)   // 19

struct PTMData {
  // central particle info
  size_t central_id;
  double coordx, coordy, coordz;
  int ncollect;
  double delx[CPTM_MAXNEAR];
  double dely[CPTM_MAXNEAR];
  double delz[CPTM_MAXNEAR];

  // 2-shell: each first-shell neighbour's neighbourhood
  bool need_2shell;
  int    neigh_ncollect[PTM_MAX_NBRS];
  double neigh_coordx[PTM_MAX_NBRS];
  double neigh_coordy[PTM_MAX_NBRS];
  double neigh_coordz[PTM_MAX_NBRS];
  double neigh_abs[PTM_MAX_NBRS][CPTM_MAXNEAR][3]; // absolute positions
};

// -------------------------------------------------------------------
// get_neighbours callback for ptm_index
// -------------------------------------------------------------------
static int get_neighbours(void *vdata,
    size_t central_id, size_t atom_id, int num,
    size_t *neigh_indices, int32_t *numbers, double (*neigh_pos)[3])
{
  PTMData *d = (PTMData *)vdata;

  if (atom_id == central_id) {
    // return central atom's sorted neighbours
    int n = (d->ncollect < num - 1) ? d->ncollect : num - 1;
    neigh_pos[0][0] = neigh_pos[0][1] = neigh_pos[0][2] = 0.0;
    neigh_indices[0] = central_id;
    numbers[0] = 0;
    for (int k = 0; k < n; k++) {
      neigh_pos[k+1][0] = d->delx[k];
      neigh_pos[k+1][1] = d->dely[k];
      neigh_pos[k+1][2] = d->delz[k];
      neigh_indices[k+1] = (size_t)(k + 1);
      numbers[k+1] = 0;
    }
    return n + 1;
  } else {
    // 2-shell: return neighbours of first-shell neighbour atom_id
    int k = (int)atom_id - 1;
    if (!d->need_2shell || k < 0 || k >= d->ncollect || k >= PTM_MAX_NBRS
        || d->neigh_ncollect[k] <= 0)
      return 0;

    double coordx = d->neigh_coordx[k], coordy = d->neigh_coordy[k], coordz = d->neigh_coordz[k];
    int n = (d->neigh_ncollect[k] < num - 1) ? d->neigh_ncollect[k] : num - 1;

    neigh_pos[0][0] = neigh_pos[0][1] = neigh_pos[0][2] = 0.0;
    neigh_indices[0] = atom_id;
    numbers[0] = 0;
    for (int m = 0; m < n; m++) {
      neigh_pos[m+1][0] = d->neigh_abs[k][m][0] - coordx;
      neigh_pos[m+1][1] = d->neigh_abs[k][m][1] - coordy;
      neigh_pos[m+1][2] = d->neigh_abs[k][m][2] - coordz;
      // give each sub-neighbour a unique ID (won't be queried further)
      neigh_indices[m+1] = (size_t)(CPTM_MAXNEAR + k * CPTM_MAXNEAR + m + 1);
      numbers[m+1] = 0;
    }
    return n + 1;
  }
}


// -------------------------------------------------------------------
// ComputePTMAtom constructor
// Syntax: compute ID group ptm/atom <struct-string> <rmsd_threshold> [inner yes/no]
//
// <struct-string>: hyphen-separated list of structure keywords, or
//                  the preset "default" (= fcc-hcp-bcc-ico) or "all".
// <rmsd_threshold>: RMSD cutoff; 0 = INFINITY (accept any match).
// inner: yes = compute for all virtual atoms including interior ucells;
//        no (default) = outer (surface) ucells only.
// -------------------------------------------------------------------
ComputePTMAtom::ComputePTMAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  list(nullptr), atompattern(nullptr), vatompattern(nullptr)
{
  if (narg < 5)
    error->all(FLERR,
        "Illegal compute ptm/atom command: "
        "syntax is  compute ID group ptm/atom <structures> <rmsd_threshold> [inner yes/no]");

  // --- parse structure string (arg[3]) ---
  // accepted formats: "default", "all", or hyphen-separated keywords
  // e.g. "fcc-hcp-bcc" or "fcc" or "fcc-bcc-dcub"

  const char *snames[] = {
    "fcc", "hcp", "bcc", "ico", "sc", "dcub", "dhex", "graphene"};
  const int32_t sflags[] = {
    PTM_CHECK_FCC, PTM_CHECK_HCP, PTM_CHECK_BCC, PTM_CHECK_ICO,
    PTM_CHECK_SC, PTM_CHECK_DCUB, PTM_CHECK_DHEX, PTM_CHECK_GRAPHENE};
  const int nstructs = 8;

  input_flags = 0;

  if (strcmp(arg[3], "default") == 0) {
    input_flags = PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_BCC | PTM_CHECK_ICO;
  } else if (strcmp(arg[3], "all") == 0) {
    input_flags = PTM_CHECK_ALL;
  } else {
    // walk through hyphen-separated tokens
    char *buf = new char[strlen(arg[3]) + 1];
    strcpy(buf, arg[3]);
    char *tok = strtok(buf, "-");
    while (tok != nullptr) {
      bool found = false;
      for (int i = 0; i < nstructs; i++) {
        if (strcmp(tok, snames[i]) == 0) {
          input_flags |= sflags[i];
          found = true;
          break;
        }
      }
      if (!found) {
        delete [] buf;
        error->all(FLERR,
            "Illegal compute ptm/atom command: unknown structure keyword");
      }
      tok = strtok(nullptr, "-");
    }
    delete [] buf;
  }

  if (input_flags == 0)
    error->all(FLERR,
        "Compute ptm/atom: must specify at least one structure type");

  // --- parse RMSD threshold (arg[4]) ---
  rmsd_threshold = universe->numeric(FLERR, arg[4]);
  if (rmsd_threshold < 0.0)
    error->all(FLERR,
        "Compute ptm/atom: rmsd_threshold must be >= 0 (0 = no threshold)");
  if (rmsd_threshold == 0.0)
    rmsd_threshold = INFINITY;

  // --- parse optional keyword arguments ---
  inner_flag = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "inner") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal compute ptm/atom command: missing value for 'inner'");
      if (strcmp(arg[iarg+1], "yes") == 0) inner_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) inner_flag = 0;
      else error->all(FLERR, "Illegal compute ptm/atom command: 'inner' must be yes or no");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal compute ptm/atom command: unknown keyword");
    }
  }

  // --- compute sortnum: max neighbours needed, with a small buffer ---
  sortnum = 4;
  if (input_flags & (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO))
    sortnum = MAX(sortnum, PTM_NUM_NBRS_FCC + 2);      // 14
  if (input_flags & PTM_CHECK_BCC)
    sortnum = MAX(sortnum, PTM_NUM_NBRS_BCC + 2);      // 16
  if (input_flags & PTM_CHECK_SC)
    sortnum = MAX(sortnum, PTM_NUM_NBRS_SC + 2);       // 8
  if (input_flags & (PTM_CHECK_DCUB | PTM_CHECK_DHEX))
    sortnum = MAX(sortnum, PTM_MAX_INPUT_POINTS);      // 19
  if (input_flags & PTM_CHECK_GRAPHENE)
    sortnum = MAX(sortnum, PTM_NUM_NBRS_GRAPHENE + 4); // 13

  peratom_flag = 1;
  size_peratom_cols = PTM_NCOLS;
  element_data_style = VIRTUAL_SURFACE;

  namax = nemax = 0;

  comm_atom_forward = PTM_NCOLS;
  comm_elem_forward = element->maxucell * element->maxapc * PTM_NCOLS;
}

/* ---------------------------------------------------------------------- */

ComputePTMAtom::~ComputePTMAtom()
{
  memory->destroy(atompattern);
  memory->destroy(vatompattern);
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute ptm/atom requires a pair style be defined");

  // PTM global initialisation (idempotent after first call)
  if (ptm_initialize_global() != PTM_NO_ERROR)
    error->all(FLERR, "PTM global initialization failed");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style, "ptm/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute ptm/atom defined");

  // request a full, occasional, vatom-sorted neighbour list
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair      = 0;
  neighbor->requests[irequest]->compute   = 1;
  neighbor->requests[irequest]->half      = 0;
  neighbor->requests[irequest]->full      = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->vatomlist = 1;
  neighbor->requests[irequest]->gausslist = 0;
  neighbor->requests[irequest]->sort      = 1;
  neighbor->requests[irequest]->sortnum   = sortnum;
  neighbor->requests[irequest]->group     = 1;
  neighbor->requests[irequest]->groupbit  = groupbit;
  neighbor->requests[irequest]->noinner   = inner_flag ? 0 : 1;
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow output arrays if needed
  if (atom->nmax > namax) {
    namax = atom->nmax;
    grow_atom(namax);
    array_atom = atompattern;
  }
  if (element->nmax > nemax) {
    nemax = element->nmax;
    grow_vatom(nemax, element->maxapc, element->maxucell);
    array_vatom = vatompattern;
  }

  // zero output arrays
  bigint nbytes = (bigint)namax * PTM_NCOLS * sizeof(double);
  if (nbytes) memset(&atompattern[0][0], 0, nbytes);
  nbytes = (bigint)nemax * element->maxapc * element->maxucell
    * PTM_NCOLS * sizeof(double);
  if (nbytes > 0) memset(&vatompattern[0][0][0][0], 0, nbytes);

  neighbor->build_one(list, preneighflag);

  double **ax    = atom->x;
  double ****nodex = element->nodex;
  int *amask     = atom->mask;
  int *emask     = element->mask;
  int *etype     = element->etype;
  int *apc       = element->apc;
  ElementVec *evec = element->evec;

  int inum  = list->inum;
  int *ilist        = list->ilist;
  int *iindexlist   = list->iindexlist;
  int *numneigh     = list->numneigh;
  int **firstneigh  = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  // --- optional: build a (particle → ilist position) map for 2-shell ---
  bool need_2shell = (input_flags &
      (PTM_CHECK_DCUB | PTM_CHECK_DHEX | PTM_CHECK_GRAPHENE)) != 0;

  std::unordered_map<int64_t, int> particle_map;
  if (need_2shell) {
    particle_map.reserve(inum * 2);
    for (int ii = 0; ii < inum; ii++) {
      int i      = ilist[ii];
      int iindex = iindexlist[ii];
      int64_t key = (int64_t)i * 100000LL + (iindex + 1);
      particle_map[key] = ii;
    }
  }

  // --- PTM local handle (per-call allocation/free) ---
  ptm_initialize_global(); // idempotent
  ptm_local_handle_t lh = ptm_initialize_local();

  // --- main loop over all particles in the neighbour list ---
  for (int ii = 0; ii < inum; ii++) {
    int i      = ilist[ii];
    int iindex = iindexlist[ii];

    int iucell, ibasis, ietype, iapc;
    int imask;
    double *outptr;

    if (iindex < 0) {
      outptr = atompattern[i];
      imask  = amask[i];
    } else {
      ietype = etype[i];
      iapc   = apc[ietype];
      iucell = iindex / iapc;
      ibasis = iindex % iapc;
      outptr = vatompattern[i][ibasis][iucell];
      imask  = emask[i];
    }

    if (!(imask & groupbit)) continue;

    // default output: OTHER, zero rmsd/idist, identity quaternion
    outptr[PTM_COL_TYPE]  = 0.0;
    outptr[PTM_COL_RMSD]  = 0.0;
    outptr[PTM_COL_IDIST] = 0.0;
    outptr[PTM_COL_QUATW] = 1.0;
    outptr[PTM_COL_QUATX] = 0.0;
    outptr[PTM_COL_QUATY] = 0.0;
    outptr[PTM_COL_QUATZ] = 0.0;

    if (numneigh[ii] < 3) continue;

    // --- collect sorted neighbour positions ---
    double coordx, coordy, coordz;
    universe->get_particle_pos(i, iindex, coordx, coordy, coordz);

    int *jlist      = firstneigh[ii];
    int *jindexlist = firstneighindex[ii];
    int nneigh      = numneigh[ii];
    int ncollect    = (nneigh < CPTM_MAXNEAR) ? nneigh : CPTM_MAXNEAR;

    PTMData d;
    d.central_id = (size_t)ii;
    d.coordx = coordx; d.coordy = coordy; d.coordz = coordz;
    d.ncollect  = ncollect;
    d.need_2shell = need_2shell;

    for (int jj = 0; jj < ncollect; jj++) {
      int j      = jlist[jj];
      int jindex = jindexlist[jj];
      double jx, jy, jz;
      universe->get_particle_pos(j, jindex, jx, jy, jz);
      d.delx[jj] = jx - coordx;
      d.dely[jj] = jy - coordy;
      d.delz[jj] = jz - coordz;
    }

    // --- fill 2-shell data (neighbours of neighbours) ---
    if (need_2shell) {
      for (int k = 0; k < ncollect && k < PTM_MAX_NBRS; k++) {
        int j      = jlist[k];
        int jindex = jindexlist[k];
        double jx  = coordx + d.delx[k];
        double jy  = coordy + d.dely[k];
        double jz  = coordz + d.delz[k];
        d.neigh_coordx[k] = jx;
        d.neigh_coordy[k] = jy;
        d.neigh_coordz[k] = jz;
        d.neigh_ncollect[k] = 0;

        int64_t key = (int64_t)j * 100000LL + (jindex + 1);
        auto it = particle_map.find(key);
        if (it == particle_map.end()) continue;

        int kk = it->second;
        int *kjlist      = firstneigh[kk];
        int *kjindexlist = firstneighindex[kk];
        int kjnum        = numneigh[kk];
        int nm = (kjnum < CPTM_MAXNEAR) ? kjnum : CPTM_MAXNEAR;
        for (int m = 0; m < nm; m++) {
          int n2      = kjlist[m];
          int n2index = kjindexlist[m];
          double nx, ny, nz;
          universe->get_particle_pos(n2, n2index, nx, ny, nz);
          d.neigh_abs[k][m][0] = nx;
          d.neigh_abs[k][m][1] = ny;
          d.neigh_abs[k][m][2] = nz;
        }
        d.neigh_ncollect[k] = nm;
      }
    }

    // --- call ptm_index ---
    int32_t type = PTM_MATCH_NONE, alloy_type;
    double  scale = 0.0, rmsd = INFINITY, interatomic_distance = 0.0;
    double  q[4] = {1.0, 0.0, 0.0, 0.0};

    ptm_index(lh, (size_t)ii,
        get_neighbours, (void *)&d,
        input_flags, false,
        &type, &alloy_type, &scale, &rmsd, q,
        nullptr, nullptr, nullptr, nullptr,
        &interatomic_distance, nullptr, nullptr);

    // apply RMSD threshold
    if (rmsd > rmsd_threshold) type = PTM_MATCH_NONE;

    if (type != PTM_MATCH_NONE) {
      outptr[PTM_COL_TYPE]  = (double)type;
      outptr[PTM_COL_RMSD]  = rmsd;
      outptr[PTM_COL_IDIST] = interatomic_distance;
      outptr[PTM_COL_QUATW] = q[0];
      outptr[PTM_COL_QUATX] = q[1];
      outptr[PTM_COL_QUATY] = q[2];
      outptr[PTM_COL_QUATZ] = q[3];
    }
  }

  ptm_uninitialize_local(lh);
}

/* ---------------------------------------------------------------------- */

double ComputePTMAtom::memory_usage()
{
  double bytes = (double)namax * PTM_NCOLS * sizeof(double);
  bytes += (double)nemax * element->maxapc * element->maxucell
    * PTM_NCOLS * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int ComputePTMAtom::pack_atom_forward_comm(int n, int *list, double *buf,
    int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (int k = 0; k < PTM_NCOLS; k++)
      buf[m++] = atompattern[j][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int m = 0, last = first + n;
  for (int i = first; i < last; i++)
    for (int k = 0; k < PTM_NCOLS; k++)
      atompattern[i][k] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int ComputePTMAtom::pack_elem_forward_comm(int n, int *list, double *buf,
    int /*pbc_flag*/, int * /*pbc*/)
{
  int *nucell = element->nucell;
  int *apc    = element->apc;
  int *etype  = element->etype;
  int m = 0;
  for (int ii = 0; ii < n; ii++) {
    int i = list[ii];
    int ietype = etype[i];
    for (int j = 0; j < apc[ietype]; j++)
      for (int k = 0; k < nucell[ietype]; k++)
        for (int c = 0; c < PTM_NCOLS; c++)
          buf[m++] = vatompattern[i][j][k][c];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int *nucell = element->nucell;
  int *apc    = element->apc;
  int *etype  = element->etype;
  int m = 0, last = first + n;
  for (int i = first; i < last; i++) {
    int ietype = etype[i];
    for (int j = 0; j < apc[ietype]; j++)
      for (int k = 0; k < nucell[ietype]; k++)
        for (int c = 0; c < PTM_NCOLS; c++)
          vatompattern[i][j][k][c] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::grow_atom(int nmax)
{
  memory->destroy(atompattern);
  memory->create(atompattern, nmax, PTM_NCOLS, "ptm:atompattern");
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::grow_vatom(int nmax, int maxapc, int maxucell)
{
  memory->destroy(vatompattern);
  memory->create(vatompattern, nmax, maxapc, maxucell, PTM_NCOLS,
      "ptm:vatompattern");
}
