/* Original code from LAMMPS CAC comm_cac.cpp by Adrian Diaz - University of
 * Florida Adopted and optimized by Thanh Phan - Iowa State University
 * */

#include "atom.h"
#include "atom_vec.h"
#include "comm_brick.h"
#include "comm_tiled.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "element.h"
#include "element_vec.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include <string.h>

using namespace CAC_NS;
using namespace MathExtra;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define EPSILON 1.0e-6
#define BOX_EPSILON 1.0e-10
#define OVERLAP_SIZE 10 // size of 1 overlap bounding box in communication
#define DELTA_PROCS 16
#define BIG 1e30

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(CAC *cac) : Comm(cac) {

  style = 1;
  layout = Comm::LAYOUT_UNIFORM;
  init_buffers();
}

/* ---------------------------------------------------------------------- */
// IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommTiled::CommTiled(CAC *cac, Comm *oldcomm) : Comm(*oldcomm) {
  style = 1;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommTiled::~CommTiled() {
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
  deallocate_buffers();
  memory->sfree(rcbinfo);
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
   ------------------------------------------------------------------------- */

void CommTiled::init_buffers() {
  // bufextra = max size of one exchanged atom
  //          = allowed overflow of sendbuf in exchange()
  // atomvec, fix reset these 2 maxexchange values if needed
  // only necessary if their size > BUFEXTRA

  dimension = domain->dimension;

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send, maxsend + bufextra, "comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv, maxrecv, "comm:buf_recv");

  atom_sendnum = atom_recvnum = NULL;
  atom_size_forward_recv = atom_firstrecv = NULL;
  atom_size_reverse_send = atom_size_reverse_recv = NULL;
  atom_forward_recv_offset = atom_reverse_recv_offset = NULL;
  elem_sendnum = elem_recvnum = NULL;
  elem_size_forward_recv = elem_firstrecv = NULL;
  elem_size_reverse_send = elem_size_reverse_recv = NULL;
  elem_forward_recv_offset = elem_reverse_recv_offset = NULL;

  pbc_flag = NULL;
  pbc = NULL;
  minsendbox = NULL;
  minsendbox_flag = NULL;
  maxsendbox = NULL;
  atom_maxsendlist = NULL;
  atom_sendlist = NULL;
  elem_maxsendlist = NULL;
  elem_sendlist = NULL;

  maxreqstat = 0;
  requests = NULL;

  maxoverlap = 0;
  overlap = NULL;
  overlap_pbc = NULL;
  overlap2box = NULL;

  maxforeign = nforeign = 0;
  foreign_overlap = NULL;
  foreign_tag = NULL;
  foreign_boxes = NULL;

  nstencil = 0;
  maxstencil = 0;
  stencil = NULL;
  fboxbins = NULL;
  fboxbinhead = NULL;
  maxbin = 0;

  nprocmax = DELTA_PROCS;
  grow_overlap(DELTA_PROCS, 0);
  overlap_repeat = new int[nprocs];

  rcbinfo = NULL;
  nexchproc = new int[dimension];
  nexchprocmax = new int[dimension];
  exchproc = new int *[dimension];
  exchnum = new int *[dimension];

  for (int i = 0; i < dimension; i++) {
    nexchprocmax[i] = DELTA_PROCS;
    exchproc[i] = new int[DELTA_PROCS];
    exchnum[i] = new int[DELTA_PROCS];
  }
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
   ------------------------------------------------------------------------- */

void CommTiled::deallocate_buffers() {
  delete[] atom_sendnum;
  delete[] atom_recvnum;
  delete[] atom_size_forward_recv;
  delete[] atom_firstrecv;
  delete[] atom_size_reverse_send;
  delete[] atom_size_reverse_recv;
  delete[] atom_forward_recv_offset;
  delete[] atom_reverse_recv_offset;

  delete[] elem_sendnum;
  delete[] elem_recvnum;
  delete[] elem_size_forward_recv;
  delete[] elem_firstrecv;
  delete[] elem_size_reverse_send;
  delete[] elem_size_reverse_recv;
  delete[] elem_forward_recv_offset;
  delete[] elem_reverse_recv_offset;

  delete[] pbc_flag;
  memory->destroy(pbc);
  memory->destroy(minsendbox);
  memory->destroy(minsendbox_flag);
  memory->destroy(maxsendbox);

  delete[] atom_maxsendlist;
  delete[] elem_maxsendlist;
  for (int i = 0; i < nprocmax; i++) {
    memory->destroy(atom_sendlist[i]);
    memory->destroy(elem_sendlist[i]);
  }
  delete[] atom_sendlist;
  delete[] elem_sendlist;

  memory->destroy(foreign_overlap);
  memory->destroy(foreign_tag);
  memory->destroy(foreign_boxes);
  memory->destroy(fboxbins);
  memory->destroy(fboxbinhead);
  memory->destroy(stencil);

  memory->destroy(overlap);
  memory->destroy(overlap_pbc);
  memory->destroy(overlap2box);
  delete[] overlap_repeat;

  delete[] requests;

  delete[] nexchproc;
  delete[] nexchprocmax;

  for (int i = 0; i < dimension; i++) {
    delete[] exchproc[i];
    delete[] exchnum[i];
  }

  delete[] exchproc;
  delete[] exchnum;
}

/* ---------------------------------------------------------------------- */

void CommTiled::init() {

  Comm::init();

  // temporary restrictions

  if (triclinic)
    error->all(FLERR, "Cannot yet use comm_style tiled with triclinic box");
  if (mode == Comm::MULTI)
    error->all(FLERR, "Cannot yet use comm_style tiled with multi-mode comm");
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size and tiling
   1. Find full list of overlap procs with me from all face, edge, corner boxes
      (a proc be in list multiple times with different pbc)
   2. Determine minsendbox for each overlap: atoms/elements overlap this box is
      sent as ghost without checking with foreign box
   3. Determine maxsendbox for each overlap: atoms/elements outside of this box
      wont be send as neighbor (if out of minsendbox and within maxsendbox,
   check with foreign box to see if need to send as ghost)
   4. Determine element box to send to other as foreign_box,
      swap foreign_boxes with other procs
   5. Setup bins, stencil for foreign boxes, bin foreign_boxes
   ------------------------------------------------------------------------- */

void CommTiled::setup_borders() {
  if (debug)
    debug_write_domain_box();
  int i, j, m, n, dim;
  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;

  int *periodicity = domain->periodicity;

  // set function pointers

  if (layout != Comm::LAYOUT_TILED) {
    box_drop = &CommTiled::box_drop_brick;
    box_other = &CommTiled::box_other_brick;
    box_touch = &CommTiled::box_touch_brick;
    point_drop = &CommTiled::point_drop_brick;
  } else {
    box_drop = &CommTiled::box_drop_tiled;
    box_other = &CommTiled::box_other_tiled;
    box_touch = &CommTiled::box_touch_tiled;
    point_drop = &CommTiled::point_drop_tiled;
  }

  // set cutoff for comm forward and comm reverse
  // check that cutoff < any periodic box length

  double cut = MAX(neighbor->cutneighmax, cutghostuser);

  // additional skin from computes

  for (i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->ghostskinflag)
      cut = MAX(modify->compute[i]->ghostskin, cut);

  // if pair exists, check for threebody pair
  // double cut since neighbor of ghost is required

  if (force->pair)
    if (force->pair->threebody_flag)
      cut *= 2.0;

  cutghost[0] = cutghost[1] = cutghost[2] = cut;

  if ((periodicity[0] && cutghost[0] > prd[0]) ||
      (periodicity[1] && cutghost[1] > prd[1]) ||
      (periodicity[2] && cutghost[2] > prd[2] && dimension == 3)) {
    error->all(FLERR, "Communication cutoff for comm_style tiled "
                      "cannot exceed periodic box length");
  }

  // update element bounding box
  // max_cutghost = maximum search range for ghost
  // NOTE: NEED TO OPTIMIZE max_cutghost LATER

  element->update_element_bound_box();

  double *box_size = element->max_element_bound_box_size;
  double max = MAX(box_size[0], box_size[1]);
  if (dimension == 3)
    max = MAX(box_size[2], max);
  MPI_Allreduce(&max, max_cutghost, 1, MPI_DOUBLE, MPI_MAX, world);
  max_cutghost[0] /= 2.0;
  max_cutghost[0] += cut;
  max_cutghost[1] = max_cutghost[2] = max_cutghost[0];

  // limit max_cutghost to prd

  int limit_flag[3];
  for (dim = 0; dim < dimension; dim++)
    if (max_cutghost[dim] >= prd[dim]) {
      max_cutghost[dim] = MIN(prd[dim], max_cutghost[dim]);
      limit_flag[dim] = 1;
    } else {
      limit_flag[dim] = 0;
    }

  //***************************************************************************
  // 1. Setup borders for full corner and edge swaps
  // determine which procs I will send to and receive from in each overlap
  // done by intersecting ghost box with all proc sub-boxes it overlaps
  // sets noverlap, overlap
  // sets overlapother, overlelf, pbc_flag, pbc
  // resets nprocmax
  //***************************************************************************

  // setup the ghost box grid
  // +/- BOX_EPSILON to account for rounding off error when max_cutghost is
  // limit to prd

  double box_grid[3][4];
  for (dim = 0; dim < dimension; dim++) {
    box_grid[dim][0] =
        sublo[dim] - max_cutghost[dim] + BOX_EPSILON * limit_flag[dim];
    box_grid[dim][1] = sublo[dim];
    box_grid[dim][2] = subhi[dim];
    box_grid[dim][3] =
        subhi[dim] + max_cutghost[dim] - BOX_EPSILON * limit_flag[dim];
  }

  if (dimension == 2)
    box_grid[2][0] = box_grid[2][1] = box_grid[2][2] = box_grid[2][3] =
        sublo[2];

  // one = first ghost box in same periodic image
  // two = second ghost box wrapped across periodic boundary
  // either may not exist
  // there can be multiple second ghost box

  int noverlap1, indexme;
  double lo1[3], hi1[3], lo2[8][3], hi2[8][3];
  int current_pbc[3], current_pbc_flag[3];
  int pbc2[8][3]; // pbc flag for second ghost box set
  int one, two;

  noverlap = 0;

  // initialize repeat check array

  for (i = 0; i < nprocs; i++)
    overlap_repeat[i] = 0;

  // loop over ghost box grids
  // NOTE: this does not works with 2D yet
  if (domain->dimension == 2)
    error->all(FLERR, "Comm_tiled does not work with 2D yet");

  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      for (int iz = 0; iz < 3; iz++) {

        // skip if it's my subbox

        if (ix == 1 && iy == 1 && iz == 1)
          continue;

        lo1[0] = box_grid[0][ix];
        lo1[1] = box_grid[1][iy];
        lo1[2] = box_grid[2][iz];
        hi1[0] = box_grid[0][ix + 1];
        hi1[1] = box_grid[1][iy + 1];
        hi1[2] = box_grid[2][iz + 1];

        one = 1;
        two = 0;

        current_pbc[0] = current_pbc[1] = current_pbc[2] = 0;
        current_pbc_flag[0] = current_pbc_flag[1] = current_pbc_flag[2] = 0;
        for (dim = 0; dim < dimension; dim++) {
          if (periodicity[dim]) {
            if (lo1[dim] < boxlo[dim])
              current_pbc[dim] = 1;
            if (hi1[dim] > boxhi[dim])
              current_pbc[dim] = -1;
            if (current_pbc[dim] != 0) {
              two = 1;
              current_pbc_flag[dim] = 1;
            }
          }
        }

        if (nprocs == 1)
          one = 0;

        int n2box = 0; // image box counter

        // compute how many and which images of lo1 to test for overlap

        if (two) {
          for (int pbx = 0; pbx <= current_pbc_flag[0]; pbx++)
            for (int pby = 0; pby <= current_pbc_flag[1]; pby++)
              for (int pbz = 0; pbz <= current_pbc_flag[2]; pbz++) {
                if (!pbx && !pby && !pbz)
                  continue;
                copy3(lo1, lo2[n2box]);
                copy3(hi1, hi2[n2box]);
                pbc2[n2box][0] = 0;
                pbc2[n2box][1] = 0;
                pbc2[n2box][2] = 0;

                if (pbx) {
                  lo2[n2box][0] += prd[0] * current_pbc[0];
                  hi2[n2box][0] += prd[0] * current_pbc[0];
                  pbc2[n2box][0] = current_pbc[0];
                }

                if (pby) {
                  lo2[n2box][1] += prd[1] * current_pbc[1];
                  hi2[n2box][1] += prd[1] * current_pbc[1];
                  pbc2[n2box][1] = current_pbc[1];
                }

                if (pbz) {
                  lo2[n2box][2] += prd[2] * current_pbc[2];
                  hi2[n2box][2] += prd[2] * current_pbc[2];
                  pbc2[n2box][2] = current_pbc[2];
                }

                for (dim = 0; dim < dimension; dim++) {
                  lo2[n2box][dim] = MAX(lo2[n2box][dim], boxlo[dim]);
                  hi2[n2box][dim] = MIN(hi2[n2box][dim], boxhi[dim]);
                }

                if (lo2[n2box][0] == hi2[n2box][0] ||
                    lo2[n2box][1] == hi2[n2box][1] ||
                    lo2[n2box][2] == hi2[n2box][2])
                  continue;

                n2box++;
              }
        }

        if (one) {
          for (dim = 0; dim < dimension; dim++) {
            lo1[dim] = MAX(lo1[dim], boxlo[dim]);
            hi1[dim] = MIN(hi1[dim], boxhi[dim]);
            if (hi1[dim] == lo1[dim])
              one = 0;
          }
        }

        // find overlap from ghost box one and two

        if (one)
          (this->*box_drop)(0, lo1, hi1, NULL, indexme, 1);
        noverlap1 = noverlap;

        if (two)
          for (i = 0; i < n2box; i++)
            (this->*box_drop)(0, lo2[i], hi2[i], pbc2[i], indexme, 2);
      }
    }
  }

  overlapself = 0;
  int tmp[6];
  int swapindex = noverlap - 1;

  // move all overlaps with me to the end of the list with
  // their respective pbc images and other overlap arrays

  if (nprocs != 1) {
    for (i = 0; i < noverlap; i++) {
      if (overlap[i] == me) {
        while (overlap[swapindex] == me) {
          swapindex--;
          overlapself++;
        }
        if (i >= swapindex)
          break;
        tmp[0] = overlap[swapindex];
        tmp[1] = overlap_pbc[swapindex][0];
        tmp[2] = overlap_pbc[swapindex][1];
        tmp[3] = overlap_pbc[swapindex][2];
        overlap[swapindex] = overlap[i];
        overlap_pbc[swapindex][0] = overlap_pbc[i][0];
        overlap_pbc[swapindex][1] = overlap_pbc[i][1];
        overlap_pbc[swapindex][2] = overlap_pbc[i][2];
        overlap[i] = tmp[0];
        overlap_pbc[i][0] = tmp[1];
        overlap_pbc[i][1] = tmp[2];
        overlap_pbc[i][2] = tmp[3];

        if (layout != Comm::LAYOUT_TILED) {
          tmp[0] = overlap2box[swapindex][0];
          tmp[1] = overlap2box[swapindex][1];
          tmp[2] = overlap2box[swapindex][2];
          tmp[3] = overlap2box[swapindex][3];
          tmp[4] = overlap2box[swapindex][4];
          tmp[5] = overlap2box[swapindex][5];
          overlap2box[swapindex][0] = overlap2box[i][0];
          overlap2box[swapindex][1] = overlap2box[i][1];
          overlap2box[swapindex][2] = overlap2box[i][2];
          overlap2box[swapindex][3] = overlap2box[i][3];
          overlap2box[swapindex][4] = overlap2box[i][4];
          overlap2box[swapindex][5] = overlap2box[i][5];
          overlap2box[i][0] = tmp[0];
          overlap2box[i][1] = tmp[1];
          overlap2box[i][2] = tmp[2];
          overlap2box[i][3] = tmp[3];
          overlap2box[i][4] = tmp[4];
          overlap2box[i][5] = tmp[5];
        }
        swapindex--;
        overlapself++;
      }
    }
  } else
    overlapself = noverlap;

  // reallocate 2nd dimensions of all send/recv arrays, based on noverlap
  // grow overlap if needed

  if (noverlap > nprocmax) {
    int oldmax = nprocmax;
    while (nprocmax < noverlap)
      nprocmax += DELTA_PROCS;
    grow_overlap(nprocmax, oldmax);
  }

  // overlap how has list of noverlap procs
  // includes pbc effects

  if (noverlap && noverlap - overlapself)
    overlapother = 1;
  else
    overlapother = 0;

  //*************************************************************************************
  // 2. Compute minsendbox for each of my overlap,
  // this is the absolute send box (elements/atoms overlap this minsendbox is
  // sent to other proc) obox = other proc subdomain box sbox = minimum box I
  // need to send to other proc
  //*************************************************************************************

  double oboxlo[3], oboxhi[3], sbox[6], sbox_multi[6];

  if (mode == Comm::SINGLE) {
    for (i = 0; i < noverlap; i++) {
      minsendbox_flag[i] = 1;
      pbc_flag[i] = 0;
      pbc[i][0] = pbc[i][1] = pbc[i][2] = pbc[i][3] = pbc[i][4] = pbc[i][5] = 0;
      if (overlap_pbc[i][0] || overlap_pbc[i][1] || overlap_pbc[i][2]) {
        pbc_flag[i] = 1;
        pbc[i][0] = overlap_pbc[i][0];
        pbc[i][1] = overlap_pbc[i][1];
        pbc[i][2] = overlap_pbc[i][2];
      }
      (this->*box_other)(0, 0, overlap[i], oboxlo, oboxhi, i);
      oboxlo[0] -= pbc[i][0] * prd[0] + cutghost[0];
      oboxlo[1] -= pbc[i][1] * prd[1] + cutghost[1];
      oboxlo[2] -= pbc[i][2] * prd[2] + cutghost[2];
      oboxhi[0] -= pbc[i][0] * prd[0] - cutghost[0];
      oboxhi[1] -= pbc[i][1] * prd[1] - cutghost[1];
      oboxhi[2] -= pbc[i][2] * prd[2] - cutghost[2];
      for (dim = 0; dim < dimension; dim++) {
        sbox[3 + dim] = MIN(oboxhi[dim], subhi[dim]);
        if (oboxlo[dim] <= sublo[dim]) {
          sbox[dim] = sublo[dim];
          if (sbox[3 + dim] <= sbox[dim])
            sbox[3 + dim] = sbox[dim];
        } else {
          sbox[dim] = oboxlo[dim];
          if (sbox[dim] >= sbox[3 + dim])
            sbox[dim] = sbox[3 + dim];
        }
      }

      // determine if minsendbox has zero thickness due to a lack of overlap
      // with force cutoff radius

      for (dim = 0; dim < dimension; dim++) {
        if (sbox[dim] == sbox[dim + 3])
          minsendbox_flag[i] = 0;
      }

      memcpy(minsendbox[i], sbox, 6 * sizeof(double));
    }
  } else {
  }

  int nmax = noverlap;
  if (nmax > maxreqstat) {
    maxreqstat = nmax;
    delete[] requests;
    requests = new MPI_Request[maxreqstat];
  }

  //*************************************************************************************
  // 3. Determine maxsendbox for each overlap
  // accquire expanded subbox due to element bounding boxes from other procs
  // send my expanded subbox and extend it by cutghost to other procs as my
  // maxsendbox
  //*************************************************************************************

  expand_subbox[0] = MIN(element->local_element_bound_box[0], sublo[0]);
  expand_subbox[1] = MIN(element->local_element_bound_box[1], sublo[1]);
  expand_subbox[2] = MIN(element->local_element_bound_box[2], sublo[2]);
  expand_subbox[3] = MAX(element->local_element_bound_box[3], subhi[0]);
  expand_subbox[4] = MAX(element->local_element_bound_box[4], subhi[1]);
  expand_subbox[5] = MAX(element->local_element_bound_box[5], subhi[2]);

  int nother = noverlap - overlapself;

  for (m = 0; m < nother; m++)
    MPI_Irecv(&maxsendbox[m][0], 6, MPI_DOUBLE, overlap[m], 0, world,
              &requests[m]);
  for (m = 0; m < nother; m++)
    MPI_Send(&expand_subbox[0], 6, MPI_DOUBLE, overlap[m], 0, world);
  if (overlapother)
    MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);

  for (i = 0; i < nother; i++) {
    maxsendbox[i][0] -= pbc[i][0] * prd[0] + cutghost[0];
    maxsendbox[i][1] -= pbc[i][1] * prd[1] + cutghost[1];
    maxsendbox[i][2] -= pbc[i][2] * prd[2] + cutghost[2];
    maxsendbox[i][3] -= pbc[i][0] * prd[0] - cutghost[0];
    maxsendbox[i][4] -= pbc[i][1] * prd[1] - cutghost[1];
    maxsendbox[i][5] -= pbc[i][2] * prd[2] - cutghost[2];
  }

  // copy to overlap with self

  for (i = nother; i < noverlap; i++) {
    maxsendbox[i][0] = expand_subbox[0] - pbc[i][0] * prd[0] - cutghost[0];
    maxsendbox[i][1] = expand_subbox[1] - pbc[i][1] * prd[1] - cutghost[1];
    maxsendbox[i][2] = expand_subbox[2] - pbc[i][2] * prd[2] - cutghost[2];
    maxsendbox[i][3] = expand_subbox[3] - pbc[i][0] * prd[0] + cutghost[0];
    maxsendbox[i][4] = expand_subbox[4] - pbc[i][1] * prd[1] + cutghost[1];
    maxsendbox[i][5] = expand_subbox[5] - pbc[i][2] * prd[2] + cutghost[2];
  }

  //*************************************************************************************
  // 4. Communicating foreign_boxes for each overlap proc
  //    obox = expanded subbox by cutghost and element bound box of other box
  //    element bound box overlap with obox will be sent to other as foreign box
  //    send/recv foreign_box
  //*************************************************************************************

  int nelocal = element->nlocal;
  double **element_bound_box = element->element_bound_box;
  int max_send_size = 0;
  int total_recv_size = 0;
  nforeign = 0;
  for (m = 0; m < noverlap; m++) {
    overlap_sendnum[m] = overlap_recvnum[m] = 0;

    oboxlo[0] = maxsendbox[m][0];
    oboxlo[1] = maxsendbox[m][1];
    oboxlo[2] = maxsendbox[m][2];
    oboxhi[0] = maxsendbox[m][3];
    oboxhi[1] = maxsendbox[m][4];
    oboxhi[2] = maxsendbox[m][5];

    // skip if no overlap between oboxlo and my expanded subbox

    if (expand_subbox[0] > oboxhi[0] || expand_subbox[3] < oboxlo[0])
      continue;
    if (expand_subbox[1] > oboxhi[1] || expand_subbox[4] < oboxlo[1])
      continue;
    if (expand_subbox[2] > oboxhi[2] || expand_subbox[5] < oboxlo[2])
      continue;

    // gather list of elements with bounding box overlapping obox
    // grow overlap_sendlist if needed

    n = 0;
    for (i = 0; i < nelocal; i++) {
      if (element_bound_box[i][0] > oboxhi[0] ||
          element_bound_box[i][3] < oboxlo[0])
        continue;
      if (element_bound_box[i][1] > oboxhi[1] ||
          element_bound_box[i][4] < oboxlo[1])
        continue;
      if (element_bound_box[i][2] > oboxhi[2] ||
          element_bound_box[i][5] < oboxlo[2])
        continue;
      if (n == overlap_maxsendlist[m])
        grow_overlap_list(m, n);
      overlap_sendlist[m][n++] = i;
    }
    overlap_sendnum[m] = n;
    max_send_size = MAX(n * OVERLAP_SIZE, max_send_size);
  }

  // communicate number of overlap element bound boxes to recv from others
  // copy if send to me

  for (m = 0; m < nother; m++)
    MPI_Irecv(&overlap_recvnum[m], 1, MPI_INT, overlap[m], 1, world,
              &requests[m]);
  for (m = 0; m < nother; m++)
    MPI_Send(&overlap_sendnum[m], 1, MPI_INT, overlap[m], 1, world);

  // copy from self

  for (m = nother; m < noverlap; m++)
    overlap_recvnum[m] = overlap_sendnum[m];

  if (overlapother) {
    MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);
    for (m = 0; m < nother; m++)
      total_recv_size += overlap_recvnum[m] * OVERLAP_SIZE;
  }

  if (noverlap) {
    overlap_firstrecv[0] = 0;
    for (m = 1; m < noverlap; m++)
      overlap_firstrecv[m] = overlap_firstrecv[m - 1] + overlap_recvnum[m - 1];
  }

  // ensure size of buf_send/recv is enough

  if (max_send_size > maxsend)
    grow_send(max_send_size, 0);
  if (total_recv_size > maxrecv)
    grow_recv(total_recv_size);

  // setup recv requests

  for (m = 0; m < nother; m++)
    if (overlap_recvnum[m])
      MPI_Irecv(&buf_recv[overlap_firstrecv[m] * OVERLAP_SIZE],
                overlap_recvnum[m] * OVERLAP_SIZE, MPI_DOUBLE, overlap[m], 2,
                world, &requests[m]);

  // pack element bound boxes to buf_send and send
  // only send/recv if sendnum/recvnum for this overlap is > 0

  for (m = 0; m < nother; m++)
    if (overlap_sendnum[m]) {
      pack_element_bound_box(overlap_sendnum[m], m, overlap_sendlist[m],
                             buf_send);
      MPI_Send(buf_send, overlap_sendnum[m] * OVERLAP_SIZE, MPI_DOUBLE,
               overlap[m], 2, world);
    }

  if (overlapother)
    MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);

  // unpack element bound boxes from buf_recv

  for (m = 0; m < nother; m++)
    if (overlap_recvnum[m])
      unpack_element_bound_box(overlap_recvnum[m], overlap_firstrecv[m],
                               &buf_recv[overlap_firstrecv[m] * OVERLAP_SIZE],
                               m);

  // pack and unpack to self

  for (m = nother; m < noverlap; m++)
    if (overlap_sendnum[m]) {
      pack_element_bound_box(overlap_sendnum[m], m, overlap_sendlist[m],
                             buf_send);
      unpack_element_bound_box(overlap_recvnum[m], overlap_firstrecv[m],
                               buf_send, m);
    }

  if (debug)
    debug_write_foreign_box(1);

  //*************************************************************************************
  // 5. Setup bins/stencils for foreign boxes
  //    Similar to nbin.cpp and nstencil.cppj
  //*************************************************************************************

  setup_bins();
  if (nforeign) {
    bin_foreign_boxes();
    create_stencil();
  }
}

/* ---------------------------------------------------------------------- */

void CommTiled::create_stencil() {
  double cutx =
      (maxfbox[0] + element->max_element_bound_box_size[0]) / 2 + cutghost[0];
  double cuty =
      (maxfbox[1] + element->max_element_bound_box_size[1]) / 2 + cutghost[1];
  double cutz =
      (maxfbox[2] + element->max_element_bound_box_size[2]) / 2 + cutghost[2];

  int sx = static_cast<int>(cutx * bininvx);
  if (sx * binsizex < cutx)
    sx++;
  int sy = static_cast<int>(cuty * bininvy);
  if (sy * binsizey < cuty)
    sy++;
  int sz;
  if (dimension == 3) {
    sz = static_cast<int>(cutz * bininvz);
    if (sz * binsizez < cutz)
      sz++;
  } else
    sz = 0;

  int smax = (2 * sx + 1) * (2 * sy + 1) * (2 * sz + 1);

  if (smax > maxstencil) {
    maxstencil = smax;
    memory->destroy(stencil);
    memory->create(stencil, maxstencil, "comm:stencil");
  }

  nstencil = 0;
  if (dimension == 3) {
    for (int k = -sz; k <= sz; k++)
      for (int j = -sy; j <= sy; j++)
        for (int i = -sx; i <= sx; i++)
          if (bin_distance(i, 0) < cutx && bin_distance(j, 1) < cuty &&
              bin_distance(k, 2) < cutz)
            stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
  } else {
    for (int j = -sy; j <= sy; j++)
      for (int i = -sx; i <= sx; i++)
        if (bin_distance(i, 0) < cutx && bin_distance(j, 1) < cuty)
          stencil[nstencil++] = j * mbinx + i;
  }
}

/* ---------------------------------------------------------------------- */

double CommTiled::bin_distance(int i, int dim) {
  double distance;
  double binsize;
  if (dim == 0)
    binsize = binsizex;
  if (dim == 1)
    binsize = binsizey;
  if (dim == 2)
    binsize = binsizez;

  if (i > 0)
    distance = (i - 1) * binsize;
  else if (i == 0)
    distance = 0.0;
  else
    distance = -(i + 1) * binsize;

  return distance;
}

/* ---------------------------------------------------------------------- */

void CommTiled::setup_bins() {
  int i, dim;

  // bsubbox = bounding box of all my foreign boxes

  double bsubbox[6], bbox[3];
  double binsize_optimal, tmp;

  bsubbox[0] = bsubbox[1] = bsubbox[2] = BIG;
  bsubbox[3] = bsubbox[4] = bsubbox[5] = -BIG;
  maxfbox[0] = maxfbox[1] = maxfbox[2] = 0;
  if (nforeign)
    for (i = 0; i < nforeign; i++) {
      for (dim = 0; dim < 3; dim++) {
        bsubbox[dim] = MIN(bsubbox[dim], foreign_boxes[i][dim]);
        bsubbox[3 + dim] = MAX(bsubbox[3 + dim], foreign_boxes[i][3 + dim]);
        tmp = foreign_boxes[i][3 + dim] - foreign_boxes[i][dim];
        maxfbox[dim] = MAX(maxfbox[dim], tmp);
      }
    }
  else {
    bsubbox[0] = sublo[0];
    bsubbox[1] = sublo[1];
    bsubbox[2] = sublo[2];
    bsubbox[3] = subhi[0];
    bsubbox[4] = subhi[1];
    bsubbox[5] = subhi[2];
  }
  tmp = MAX((maxfbox[0] + element->max_element_bound_box_size[0]) / 2.0 +
                cutghost[0],
            (maxfbox[1] + element->max_element_bound_box_size[1]) / 2.0 +
                cutghost[1]);
  tmp = MAX((maxfbox[2] + element->max_element_bound_box_size[2]) / 2.0 +
                cutghost[2],
            tmp);

  MPI_Allreduce(&tmp, &binsize_optimal, 1, MPI_DOUBLE, MPI_MAX, world);

  if (nforeign) {
    bsubbox[0] = MIN(bsubbox[0], sublo[0] - binsize_optimal);
    bsubbox[1] = MIN(bsubbox[1], sublo[1] - binsize_optimal);
    bsubbox[2] = MIN(bsubbox[2], sublo[2] - binsize_optimal);
    bsubbox[3] = MAX(bsubbox[0], subhi[0] + binsize_optimal);
    bsubbox[4] = MAX(bsubbox[4], subhi[1] + binsize_optimal);
    bsubbox[5] = MAX(bsubbox[5], subhi[2] + binsize_optimal);

    double binsizeinv =
        2.0 / binsize_optimal; // binsize_optimal is doubled here hence the 2.0

    bbox[0] = boxhi[0] - boxlo[0];
    bbox[1] = boxhi[1] - boxlo[1];
    bbox[2] = boxhi[2] - boxlo[2];

    if (bbox[0] * binsizeinv > MAXSMALLINT ||
        bbox[1] * binsizeinv > MAXSMALLINT ||
        bbox[2] * binsizeinv > MAXSMALLINT)
      error->all(FLERR, "Domain too large for foreign bins");

    nbinx = static_cast<int>(bbox[0] * binsizeinv);
    nbiny = static_cast<int>(bbox[1] * binsizeinv);
    if (dimension == 3)
      nbinz = static_cast<int>(bbox[2] * binsizeinv);
    else
      nbinz = 1;

    if (nbinx == 0)
      nbinx = 1;
    if (nbiny == 0)
      nbiny = 1;
    if (nbinz == 0)
      nbinz = 1;

    binsizex = bbox[0] / nbinx;
    binsizey = bbox[1] / nbiny;
    binsizez = bbox[2] / nbinz;

    bininvx = 1.0 / binsizex;
    bininvy = 1.0 / binsizey;
    bininvz = 1.0 / binsizez;

    int mbinxhi, mbinyhi, mbinzhi;
    double coord;

    coord = bsubbox[0] - EPSILON * bbox[0];
    mbinxlo = static_cast<int>((coord - boxlo[0]) * bininvx);
    if (coord < boxlo[0])
      mbinxlo--;
    coord = bsubbox[3] + EPSILON * bbox[0];
    mbinxhi = static_cast<int>((coord - boxlo[0]) * bininvx);

    coord = bsubbox[1] - EPSILON * bbox[1];
    mbinylo = static_cast<int>((coord - boxlo[1]) * bininvy);
    if (coord < boxlo[1])
      mbinylo--;
    coord = bsubbox[4] + EPSILON * bbox[1];
    mbinyhi = static_cast<int>((coord - boxlo[1]) * bininvy);

    if (dimension == 3) {
      coord = bsubbox[2] - EPSILON * bbox[2];
      mbinzlo = static_cast<int>((coord - boxlo[2]) * bininvz);
      if (coord < boxlo[2])
        mbinzlo--;
      coord = bsubbox[5] + EPSILON * bbox[2];
      mbinzhi = static_cast<int>((coord - boxlo[2]) * bininvz);
    }

    mbinxlo--;
    mbinxhi++;
    mbinx = mbinxhi - mbinxlo + 1;

    mbinylo--;
    mbinyhi++;
    mbiny = mbinyhi - mbinylo + 1;

    if (dimension == 3) {
      mbinzlo--;
      mbinzhi++;
    } else
      mbinzlo = mbinzhi = 0;
    mbinz = mbinzhi - mbinzlo + 1;

    bigint bbin = ((bigint)mbinx) * ((bigint)mbiny) * ((bigint)mbinz) + 1;
    if (bbin > MAXSMALLINT)
      error->one(FLERR, "Too many foreign bins");
    mbins = bbin;

    if (mbins > maxbin) {
      maxbin = mbins;
      memory->destroy(fboxbinhead);
      memory->create(fboxbinhead, maxbin, "comm:fboxbinhead");
    }
  }
}

/* ---------------------------------------------------------------------- */

void CommTiled::bin_foreign_boxes() {
  int i, dim, ibin;
  for (i = 0; i < mbins; i++)
    fboxbinhead[i] = -1;

  double coord[3];
  for (i = nforeign - 1; i >= 0; i--) {
    for (dim = 0; dim < 3; dim++)
      coord[dim] = (foreign_boxes[i][dim] + foreign_boxes[i][dim + 3]) / 2.0;
    ibin = coord2bin(coord[0], coord[1], coord[2]);
    fboxbins[i] = fboxbinhead[ibin];
    fboxbinhead[ibin] = i;
  }
}

/* ---------------------------------------------------------------------- */

int CommTiled::coord2bin(double x, double y, double z) {
  int ix, iy, iz;

  if (!ISFINITE(x) || !ISFINITE(y) || !ISFINITE(z))
    error->one(FLERR, "Non-numeric positions - simulation unstable");

  if (x >= boxhi[0])
    ix = static_cast<int>((x - boxhi[0]) * bininvx) + nbinx;
  else if (x >= boxlo[0]) {
    ix = static_cast<int>((x - boxlo[0]) * bininvx);
    ix = MIN(ix, nbinx - 1);
  } else
    ix = static_cast<int>((x - boxlo[0]) * bininvx) - 1;

  if (y >= boxhi[1])
    iy = static_cast<int>((y - boxhi[1]) * bininvy) + nbiny;
  else if (y >= boxlo[1]) {
    iy = static_cast<int>((y - boxlo[1]) * bininvy);
    iy = MIN(iy, nbiny - 1);
  } else
    iy = static_cast<int>((y - boxlo[1]) * bininvy) - 1;

  if (z >= boxhi[2])
    iz = static_cast<int>((z - boxhi[2]) * bininvz) + nbinz;
  else if (z >= boxlo[2]) {
    iz = static_cast<int>((z - boxlo[2]) * bininvz);
    iz = MIN(iz, nbinz - 1);
  } else
    iz = static_cast<int>((z - boxlo[2]) * bininvz) - 1;
  return (iz - mbinzlo) * mbiny * mbinx + (iy - mbinylo) * mbinx +
         (ix - mbinxlo);
}

/* ---------------------------------------------------------------------- */

void CommTiled::grow_foreign_box() {
  maxforeign += BUFEXTRA;
  memory->grow(foreign_overlap, maxforeign, "comm:foreign_overlap");
  memory->grow(foreign_tag, maxforeign, "comm:foreign_tag");
  memory->grow(foreign_boxes, maxforeign, 6, "comm:foreign_boxes");
  memory->grow(fboxbins, maxforeign, "comm:fboxbins");
}

/* ---------------------------------------------------------------------- */

int CommTiled::pack_element_bound_box(int n, int overlap_index, int *list,
                                      double *buf) {
  int i, j, m;
  double dx, dy, dz;
  double **element_bound_box = element->element_bound_box;

  m = 0;
  dx = dy = dz = 0.0;

  if (pbc_flag[overlap_index]) {
    if (domain->triclinic == 0) {
      dx = pbc[overlap_index][0] * domain->xprd;
      dy = pbc[overlap_index][1] * domain->yprd;
      dz = pbc[overlap_index][2] * domain->zprd;
    } else {
      dx = pbc[overlap_index][0];
      dy = pbc[overlap_index][1];
      dz = pbc[overlap_index][2];
    }
  }
  for (j = 0; j < n; j++) {
    i = list[j];
    // buf[m++] = element->tag[i];    // for debugging
    buf[m++] = pbc[overlap_index][0];
    buf[m++] = pbc[overlap_index][1];
    buf[m++] = pbc[overlap_index][2];
    buf[m++] = element_bound_box[i][0] + dx;
    buf[m++] = element_bound_box[i][1] + dy;
    buf[m++] = element_bound_box[i][2] + dz;
    buf[m++] = element_bound_box[i][3] + dx;
    buf[m++] = element_bound_box[i][4] + dy;
    buf[m++] = element_bound_box[i][5] + dz;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int CommTiled::unpack_element_bound_box(int n, int first, double *buf,
                                        int overlap_index) {
  int i, j, m, last;
  int foreign_proc;   // which proc this foreign box comes from
  int foreign_pbc[3]; // pbc image of this foreign box
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {

    if (nforeign == maxforeign)
      grow_foreign_box();

    // foreign_tag[i] = (int) buf[m++];        // for debugging
    foreign_proc = overlap[overlap_index];
    foreign_pbc[0] = (int)buf[m++];
    foreign_pbc[1] = (int)buf[m++];
    foreign_pbc[2] = (int)buf[m++];
    foreign_overlap[i] = -1;

    // determine which overlap I need to send atoms/elements overlaping with
    // this foreign box

    for (j = 0; j < noverlap; j++)
      if (foreign_proc == overlap[j] && foreign_pbc[0] == -pbc[j][0] &&
          foreign_pbc[1] == -pbc[j][1] && foreign_pbc[2] == -pbc[j][2]) {
        foreign_overlap[i] = j;
      }
    if (foreign_overlap[i] < 0) {
      error->one(FLERR,
                 "Cannot find which overlap this foreign box is relevant to");
    }

    // foreign box = element bound box extended by cutghost

    foreign_boxes[i][0] = buf[m++] - cutghost[0];
    foreign_boxes[i][1] = buf[m++] - cutghost[1];
    foreign_boxes[i][2] = buf[m++] - cutghost[2];
    foreign_boxes[i][3] = buf[m++] + cutghost[0];
    foreign_boxes[i][4] = buf[m++] + cutghost[1];
    foreign_boxes[i][5] = buf[m++] + cutghost[2];
    nforeign++;
  }

  return m;
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition exchange communication patterns
   function of current box size and tiling
   ------------------------------------------------------------------------- */

void CommTiled::setup_exchange() {
  int i, j;

  // domain properties used in setup method and methods it calls

  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  sublo = domain->sublo;
  subhi = domain->subhi;

  int *periodicity = domain->periodicity;

  // set function pointers

  if (layout != Comm::LAYOUT_TILED) {
    box_drop = &CommTiled::box_drop_brick;
    box_other = &CommTiled::box_other_brick;
    box_touch = &CommTiled::box_touch_brick;
    point_drop = &CommTiled::point_drop_brick;
  } else {
    box_drop = &CommTiled::box_drop_tiled;
    box_other = &CommTiled::box_other_tiled;
    box_touch = &CommTiled::box_touch_tiled;
    point_drop = &CommTiled::point_drop_tiled;
  }

  // if RCB decomp exists and just changed, gather needed global RCB info

  if (layout == Comm::LAYOUT_TILED)
    coord2proc_setup();

  // setup exchange communication
  // find nearest neighbor comm to infer touching exchange procs
  // determine which procs I will send to and receive from in each exchange
  // direction done by intersecting epsilon box on each face with all proc
  // sub-boxes it overlaps sets nexchproc & exchproc, resets nexchprocmax

  double cut = MIN(prd[0], prd[1]);
  if (dimension == 3)
    cut = MIN(cut, prd[2]);
  cut *= EPSILON * EPSILON;

  int indexme, noverlap1;
  double lo[3], hi[3];

  for (int idim = 0; idim < dimension; idim++) {

    // noverlap = # of overlaps of box with procs via box_drop()
    // overlap = list of overlapping procs
    // if overlap with self, indexme = index of me in list

    noverlap = 0;

    // west/south/down directions

    if (periodicity[idim] || sublo[idim] != boxlo[idim]) {

      copy3(sublo, lo);
      copy3(subhi, hi);
      lo[idim] = MAX(sublo[idim] - cut, boxlo[idim]);
      hi[idim] = sublo[idim];

      if (lo[idim] == hi[idim]) {
        lo[idim] = boxhi[idim] - cut;
        hi[idim] = boxhi[idim];
      }

      indexme = -1;
      (this->*box_drop)(idim, lo, hi, NULL, indexme, 0);

      // if self is in overlap list remove it from list

      if (indexme >= 0) {
        int tmp = overlap[noverlap - 1];
        overlap[noverlap - 1] = overlap[indexme];
        overlap[indexme] = tmp;
        noverlap--;
      }

      // loop through overlap list, if procs does not touch me, remove from list

      i = 0;
      while (i < noverlap) {
        if (!(this->*box_touch)(overlap[i], idim, 0)) {
          overlap[i] = overlap[noverlap - 1];
          noverlap--;
        } else
          i++;
      }
    }
    noverlap1 = noverlap;

    // east/north/up directions

    if (periodicity[idim] || subhi[idim] != boxhi[idim]) {

      copy3(sublo, lo);
      copy3(subhi, hi);
      lo[idim] = subhi[idim];
      hi[idim] = MIN(subhi[idim] + cut, boxhi[idim]);

      if (lo[idim] == hi[idim]) {
        lo[idim] = boxlo[idim];
        hi[idim] = boxlo[idim] + cut;
      }

      indexme = -1;
      (this->*box_drop)(idim, lo, hi, NULL, indexme, 0);

      // if self is in overlap list remove it from list

      if (indexme >= 0) {
        int tmp = overlap[noverlap - 1];
        overlap[noverlap - 1] = overlap[indexme];
        overlap[indexme] = tmp;
        noverlap--;
      }

      // loop through second half of overlap list, if procs does not touch me,
      // remove from list compare with first half of overlap list to insure each
      // proc appears exactly once remove repeats

      i = noverlap1;
      while (i < noverlap) {
        if (!(this->*box_touch)(overlap[i], idim, 1)) {
          overlap[i] = overlap[noverlap - 1];
          noverlap--;
        } else {
          for (j = 0; j < noverlap1; j++)
            if (overlap[i] == overlap[j]) {
              overlap[i] = overlap[noverlap - 1];
              noverlap--;
              break;
            }
          if (j == noverlap1)
            i++;
        }
      }
    }

    // reallocate exchproc and exchnum if needed based on noverlap
    // copy overlap list to exchproc list

    if (noverlap > nexchprocmax[idim]) {
      while (nexchprocmax[idim] < noverlap)
        nexchprocmax[idim] += DELTA_PROCS;
      delete[] exchproc[idim];
      exchproc[idim] = new int[nexchprocmax[idim]];
      delete[] exchnum[idim];
      exchnum[idim] = new int[nexchprocmax[idim]];
    }

    nexchproc[idim] = noverlap;
    for (i = 0; i < noverlap; i++)
      exchproc[idim][i] = overlap[i];
  }

  // reallocate MPI Requests and Statuses as needed

  int nmax = 0;
  for (i = 0; i < dimension; i++)
    nmax = MAX(nmax, nexchprocmax[i]);
  if (nmax > maxreqstat) {
    maxreqstat = nmax;
    delete[] requests;
    requests = new MPI_Request[maxreqstat];
  }
}

/* ----------------------------------------------------------------------
   forward communication of atom/element coords every timestep
   other per-atom/element attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm(int dummy) {
  int i, irecv, n, nother;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **x = atom->x;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if overlapself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  nother = noverlap - overlapself;

  // exchange atoms

  if (atom->natoms) {

    for (i = 0; i < nother; i++)
      MPI_Irecv(x[atom_firstrecv[i]], atom_size_forward_recv[i], MPI_DOUBLE,
                overlap[i], 0, world, &requests[i]);
    for (i = 0; i < nother; i++) {
      n = avec->pack_comm(atom_sendnum[i], atom_sendlist[i], buf_send,
                          pbc_flag[i], pbc[i]);
      MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
    }

    for (i = nother; i < noverlap; i++)
      avec->pack_comm(atom_sendnum[i], atom_sendlist[i], x[atom_firstrecv[i]],
                      pbc_flag[i], pbc[i]);

    if (overlapother)
      MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);
  }

  // exchange elements

  if (element->nelements) {
    for (i = 0; i < nother; i++)
      MPI_Irecv(&buf_recv[elem_size_forward * elem_forward_recv_offset[i]],
                elem_size_forward_recv[i], MPI_DOUBLE, overlap[i], 0, world,
                &requests[i]);

    for (i = 0; i < nother; i++) {
      n = evec->pack_comm(elem_sendnum[i], elem_sendlist[i], buf_send,
                          pbc_flag[i], pbc[i]);
      MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
    }

    for (i = nother; i < noverlap; i++) {
      evec->pack_comm(elem_sendnum[i], elem_sendlist[i], buf_send, pbc_flag[i],
                      pbc[i]);
      evec->unpack_comm(elem_recvnum[i], elem_firstrecv[i], buf_send);
    }

    for (i = 0; i < nother; i++) {
      MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
      evec->unpack_comm(
          elem_recvnum[irecv], elem_firstrecv[irecv],
          &buf_recv[elem_size_forward * elem_forward_recv_offset[irecv]]);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms (only) every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm() {

  int i, irecv, n, nother;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;
  double **f = atom->f;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if overlapself is set
  // wait on all procs except self and unpack received data
  // if comm_f_only set, exchange or copy directly from f, don't pack

  nother = noverlap - overlapself;

  if (atom->natoms) {
    if (comm_f_only) {
      for (i = 0; i < nother; i++) {
        MPI_Irecv(&buf_recv[atom_size_reverse * atom_reverse_recv_offset[i]],
                  atom_size_reverse_recv[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      }
      for (i = 0; i < nother; i++)
        MPI_Send(f[atom_firstrecv[i]], atom_size_reverse_send[i], MPI_DOUBLE,
                 overlap[i], 0, world);

      for (i = nother; i < noverlap; i++)
        avec->unpack_reverse(atom_sendnum[i], atom_sendlist[i],
                             f[atom_firstrecv[i]]);

      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        avec->unpack_reverse(
            atom_sendnum[irecv], atom_sendlist[irecv],
            &buf_recv[atom_size_reverse * atom_reverse_recv_offset[irecv]]);
      }
    } else {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[atom_size_reverse * atom_reverse_recv_offset[i]],
                  atom_size_reverse_recv[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = avec->pack_reverse(atom_recvnum[i], atom_firstrecv[i], buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }
      for (i = nother; i < noverlap; i++) {
        avec->pack_reverse(atom_recvnum[i], atom_firstrecv[i], buf_send);
        avec->unpack_reverse(atom_sendnum[i], atom_sendlist[i], buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        avec->unpack_reverse(
            atom_sendnum[irecv], atom_sendlist[irecv],
            &buf_recv[atom_size_reverse * atom_reverse_recv_offset[irecv]]);
      }
    }
  }

  // exchange elements if threebody potential

  if (element->nelements && force->pair->threebody_flag) {

    for (i = 0; i < nother; i++)
      MPI_Irecv(&buf_recv[elem_size_reverse * elem_reverse_recv_offset[i]],
                elem_size_reverse_recv[i], MPI_DOUBLE, overlap[i], 0, world,
                &requests[i]);
    for (i = 0; i < nother; i++) {
      n = evec->pack_reverse(elem_recvnum[i], elem_firstrecv[i], buf_send);
      MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
    }

    for (i = nother; i < noverlap; i++) {
      evec->pack_reverse(elem_recvnum[nother], elem_firstrecv[nother],
                         buf_send);
      evec->unpack_reverse(elem_sendnum[nother], elem_sendlist[nother],
                           buf_send);
    }
    for (i = 0; i < nother; i++) {
      MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
      evec->unpack_reverse(
          elem_sendnum[irecv], elem_sendlist[irecv],
          &buf_recv[elem_size_reverse * elem_reverse_recv_offset[irecv]]);
    }
  }
}

/* ----------------------------------------------------------------------
   exchange move atoms/elements to correct processors
   atoms/elements exchanged with procs that touch sub-box in each of 3 dims
   send out atoms/elements that have left my box, receive ones entering my box
   atoms/elements will be lost if not inside a touching proc's box
   can happen if atom/element moves outside of non-periodic bounary
   or if atom/element moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms/elements must be in lamda coords (0-1) before exchange
   is called
   ------------------------------------------------------------------------- */

void CommTiled::exchange() {
  bigint n, sum;
  int i, m, nexch, nsend, nrecv, nlocal, proc, offset;
  double lo, hi, value;
  double **x;
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (atom_map_style)
    atom->map_clear();
  if (elem_map_style)
    element->map_clear();
  atom->nghost = 0;
  element->nghost = 0;

  // insure send buf is large enough for single atom
  // bufextra = max size of one atom = allowed overflow of sendbuf
  // fixes can change per-atom size requirement on-the-fly

  int bufextra_old = bufextra;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
  if (bufextra > bufextra_old)
    memory->grow(buf_send, maxsend + bufextra, "comm:buf_send");

  // domain properties used in exchange method and methods it calls
  // subbox bounds for orthogonal or triclinic

  prd = domain->prd;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  // exchange atoms

  if (atom->natoms) {

    for (int dim = 0; dim < dimension; dim++) {

      // fill buffer with atoms leaving my box, using < and >=
      // when atom is deleted, fill it in with last atom

      x = atom->x;
      nlocal = atom->nlocal;
      lo = sublo[dim];
      hi = subhi[dim];

      i = nsend = 0;

      while (i < nlocal) {
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend)
            grow_send(nsend, 1);
          proc = (this->*point_drop)(dim, x[i]);
          if (proc != me) {
            buf_send[nsend++] = proc;
            nsend += avec->pack_exchange(i, &buf_send[nsend]);
          }
          avec->copy(nlocal - 1, i, 1);
          nlocal--;
        } else
          i++;
      }
      atom->nlocal = nlocal;

      // send and recv atoms from neighbor procs that touch my sub-box in dim
      // no send/recv with self
      // send size of message first
      // receiver may receive multiple messages, realloc buf_recv if needed

      nexch = nexchproc[dim];
      if (!nexch)
        continue;

      for (m = 0; m < nexch; m++)
        MPI_Irecv(&exchnum[dim][m], 1, MPI_INT, exchproc[dim][m], 0, world,
                  &requests[m]);
      for (m = 0; m < nexch; m++)
        MPI_Send(&nsend, 1, MPI_INT, exchproc[dim][m], 0, world);
      MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);

      nrecv = 0;
      for (m = 0; m < nexch; m++)
        nrecv += exchnum[dim][m];
      if (nrecv > maxrecv)
        grow_recv(nrecv);

      offset = 0;
      for (m = 0; m < nexch; m++) {
        MPI_Irecv(&buf_recv[offset], exchnum[dim][m], MPI_DOUBLE,
                  exchproc[dim][m], 0, world, &requests[m]);
        offset += exchnum[dim][m];
      }
      for (m = 0; m < nexch; m++)
        MPI_Send(buf_send, nsend, MPI_DOUBLE, exchproc[dim][m], 0, world);
      MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);

      // check incoming atoms to see if I own it and they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      //   atom may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        proc = static_cast<int>(buf_recv[m++]);
        if (proc == me) {
          value = buf_recv[m + dim + 1];
          if (value >= lo && value < hi) {
            m += avec->unpack_exchange(&buf_recv[m]);
            continue;
          }
        }
        m += static_cast<int>(buf_recv[m]);
      }
    }
  }

  // exchange elements

  if (element->nelements) {

    for (int dim = 0; dim < dimension; dim++) {

      // fill buffer with elements leaving my box, using < and >=
      // when element is deleted, fill it in with last element

      x = element->x;
      nlocal = element->nlocal;

      lo = sublo[dim];
      hi = subhi[dim];
      i = nsend = 0;

      while (i < nlocal) {
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend)
            grow_send(nsend, 1);
          proc = (this->*point_drop)(dim, x[i]);
          if (proc != me) {
            buf_send[nsend++] = proc;
            nsend += evec->pack_exchange(i, &buf_send[nsend]);
          }
          evec->copy(nlocal - 1, i, 1);
          nlocal--;
        } else
          i++;
      }

      element->nlocal = nlocal;

      // send and recv elements from neighbor procs that touch my sub-box in dim
      // no send/recv with self
      // send size of message first
      // receiver may receive multiple messages, realloc buf_recv if needed

      nexch = nexchproc[dim];
      if (!nexch)
        continue;

      for (m = 0; m < nexch; m++)
        MPI_Irecv(&exchnum[dim][m], 1, MPI_INT, exchproc[dim][m], 0, world,
                  &requests[m]);
      for (m = 0; m < nexch; m++)
        MPI_Send(&nsend, 1, MPI_INT, exchproc[dim][m], 0, world);
      MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);
      nrecv = 0;
      for (m = 0; m < nexch; m++)
        nrecv += exchnum[dim][m];
      if (nrecv > maxrecv)
        grow_recv(nrecv);

      offset = 0;
      for (m = 0; m < nexch; m++) {
        MPI_Irecv(&buf_recv[offset], exchnum[dim][m], MPI_DOUBLE,
                  exchproc[dim][m], 0, world, &requests[m]);
        offset += exchnum[dim][m];
      }
      for (m = 0; m < nexch; m++)
        MPI_Send(buf_send, nsend, MPI_DOUBLE, exchproc[dim][m], 0, world);
      MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);

      // check incoming elements to see if I own it and they are in my box
      // if so, add to my list
      // box check is only for this dimension,
      //   element may be passed to another proc in later dims

      m = 0;
      while (m < nrecv) {
        proc = static_cast<int>(buf_recv[m++]);
        if (proc == me) {
          value = buf_recv[m + dim + 1];
          if (value >= lo && value < hi) {
            m += evec->unpack_exchange(&buf_recv[m]);
            continue;
          }
        }
        m += static_cast<int>(buf_recv[m]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   borders list nearby atoms/elements to send to neighboring procs at every
   timestep one list is created per overlap proc that will be made as list is
   made, actually do communication this does equivalent of a forward_comm(), so
   don't need to explicitly call forward_comm() on reneighboring timestep this
   routine is called before every reneighboring for triclinic, atoms must be in
   lamda coords (0-1) before borders is called
   ------------------------------------------------------------------------- */

void CommTiled::borders() {

  int i, j, k, m, n, nother, ibin;
  int checkflag; // flag if atoms/elements need to check with foreign boxes
  double *lobox; // min sendbox I need to send
  double
      *hibox; // max sendbox out of which atoms/elements are not needed as ghost
  double *coord;
  double **x;
  int nlocal, nscountall, nrcountall, nscount, nrcount;
  int *sendflag =
      new int[noverlap]; // flag per overlap if I will be send to that overlap
  AtomVec *avec = atom->avec;
  ElementVec *evec = element->evec;

  // send/recv max one = max # of atoms/elements in single send/recv
  // send/recv max all = max # of atoms/elements in all sends/recvs

  smaxone = smaxall = 0;
  rmaxone = rmaxall = 0;

  nother = noverlap - overlapself;

  for (m = 0; m < noverlap; m++) {
    atom_sendnum[m] = 0;
    elem_sendnum[m] = 0;
  }

  // Loop through all my local atoms/elements
  // 1. loop throught noverlap and check with minsendbox and maxsendbox
  //    set sendflag and checkflag
  // 2. If checkflag = 1, loop through foreign box in stencil to check
  //    if foreign box belong to overlap with sendflag = 1, skip and don't check
  // Grow buf_send/recv if necessary

  // swap atoms

  if (atom->natoms) {

    // find atoms within minsendboxes using >= and <
    // hi test with ">" is important b/c don't want to send an atom
    //   in lower dim (on boundary) that a proc will recv again in higher dim
    // store sent atom indices in atom_sendlist for use in future timesteps

    x = atom->x;
    nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      coord = x[i];
      checkflag = 0;
      for (m = 0; m < noverlap; m++) {
        sendflag[m] = -1;
        lobox = minsendbox[m];
        hibox = maxsendbox[m];
        if (minsendbox_flag[m])
          if (coord[0] >= lobox[0] && coord[0] < lobox[3] &&
              coord[1] >= lobox[1] && coord[1] < lobox[4] &&
              coord[2] >= lobox[2] && coord[2] < lobox[5]) {
            if (atom_sendnum[m] == atom_maxsendlist[m])
              grow_atom_list(m, atom_sendnum[m]);

            atom_sendlist[m][atom_sendnum[m]++] = i;
            sendflag[m] = 1;
            continue;
          }

        if (coord[0] >= hibox[0] && coord[0] < hibox[3] &&
            coord[1] >= hibox[1] && coord[1] < hibox[4] &&
            coord[2] >= hibox[2] && coord[2] < hibox[5]) {

          sendflag[m] = 0;
          checkflag = 1;
        }
      }

      if (checkflag && nforeign) {
        ibin = coord2bin(coord[0], coord[1], coord[2]);
        for (k = 0; k < nstencil; k++) {
          if (ibin + stencil[k] < 0 || ibin + stencil[k] >= mbins)
            error->one(FLERR, "Incorrect bin & stencil setup"); // sanity check
          for (j = fboxbinhead[ibin + stencil[k]]; j >= 0; j = fboxbins[j]) {
            m = foreign_overlap[j];
            if (sendflag[m] == 0) {
              if (coord[0] >= foreign_boxes[j][0] &&
                  coord[0] <= foreign_boxes[j][3] &&
                  coord[1] >= foreign_boxes[j][1] &&
                  coord[1] <= foreign_boxes[j][4] &&
                  coord[2] >= foreign_boxes[j][2] &&
                  coord[2] <= foreign_boxes[j][5]) {
                if (atom_sendnum[m] == atom_maxsendlist[m])
                  grow_atom_list(m, atom_sendnum[m]);
                atom_sendlist[m][atom_sendnum[m]++] = i;
                sendflag[m] = 1;
              }
            }
          }
        }
      }
    }

    // send atom_sendnum counts to procs who recv from me except self
    // copy data to self if overlapself is set

    if (overlapother) {
      for (m = 0; m < nother; m++)
        MPI_Irecv(&atom_recvnum[m], 1, MPI_INT, overlap[m], 0, world,
                  &requests[m]);
      for (m = 0; m < nother; m++)
        MPI_Send(&atom_sendnum[m], 1, MPI_INT, overlap[m], 0, world);
    }
    if (overlapself)
      for (m = nother; m < noverlap; m++)
        atom_recvnum[m] = atom_sendnum[m];
    if (overlapother)
      MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);

    // setup other per proc values from atom_sendnum and atom_recvnum

    nscountall = nrcountall = 0;
    for (m = 0; m < noverlap; m++) {
      nscount = atom_sendnum[m];
      nrcount = atom_recvnum[m];
      smaxone = MAX(smaxone, nscount);
      rmaxone = MAX(rmaxone, nrcount);
      nscountall += nscount;
      nrcountall += nrcount;

      atom_size_reverse_recv[m] = nscount * atom_size_reverse;
      if (m == 0)
        atom_reverse_recv_offset[0] = 0;
      else
        atom_reverse_recv_offset[m] =
            atom_reverse_recv_offset[m - 1] + atom_sendnum[m - 1];

      atom_size_forward_recv[m] = nrcount * atom_size_forward;
      atom_size_reverse_send[m] = nrcount * atom_size_reverse;
      if (m == 0) {
        atom_firstrecv[0] = nlocal;
        atom_forward_recv_offset[0] = 0;
      } else {
        atom_firstrecv[m] = atom_firstrecv[m - 1] + atom_recvnum[m - 1];
        atom_forward_recv_offset[m] =
            atom_forward_recv_offset[m - 1] + atom_recvnum[m - 1];
      }
    }
    smaxall = MAX(smaxall, nscountall);
    rmaxall = MAX(rmaxall, nrcountall);

    // insure send/recv buffers are large enough for this border comm swap

    if (smaxone * atom_size_border > maxsend)
      grow_send(smaxone * atom_size_border, 0);
    if (rmaxall * atom_size_border > maxrecv)
      grow_recv(rmaxall * atom_size_border);

    // swap atoms with other procs using pack_border(), unpack_border()
    // use Waitall() instead of Waitany() because calls to unpack_border()
    //   must increment per-atom arrays in ascending order

    if (overlapother) {
      for (m = 0; m < nother; m++)
        MPI_Irecv(&buf_recv[atom_size_border * atom_forward_recv_offset[m]],
                  atom_recvnum[m] * atom_size_border, MPI_DOUBLE, overlap[m], 0,
                  world, &requests[m]);
      for (m = 0; m < nother; m++) {
        n = avec->pack_border(atom_sendnum[m], atom_sendlist[m], buf_send,
                              pbc_flag[m], pbc[m]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[m], 0, world);
      }
    }
    if (overlapself) {
      for (m = nother; m < noverlap; m++) {
        avec->pack_border(atom_sendnum[m], atom_sendlist[m], buf_send,
                          pbc_flag[m], pbc[m]);
        avec->unpack_border(atom_recvnum[m], atom_firstrecv[m], buf_send);
      }
    }

    if (overlapother) {
      MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);
      for (m = 0; m < nother; m++)
        avec->unpack_border(
            atom_recvnum[m], atom_firstrecv[m],
            &buf_recv[atom_size_border * atom_forward_recv_offset[m]]);
    }

    // increment ghost atoms

    if (noverlap)
      atom->nghost +=
          atom_forward_recv_offset[noverlap - 1] + atom_recvnum[noverlap - 1];

    // insure send/recv buffers are long enough for all forward & reverse comm
    // send buf is for one forward or reverse sends to one proc
    // recv buf is for all forward or reverse recvs in one swap

    int max = MAX(atom_maxforward * smaxone, atom_maxreverse * rmaxone);
    if (max > maxsend)
      grow_send(max, 0);
    max = MAX(atom_maxforward * rmaxall, atom_maxreverse * smaxall);
    if (max > maxrecv)
      grow_recv(max);

    // reset global->local map

    if (atom_map_style)
      atom->map_set();
  }

  // swap elements

  if (element->nelements) {

    // find atoms within minotherboxes using >= and <
    // hi test with ">" is important b/c don't want to send an atom
    //   in lower dim (on boundary) that a proc will recv again in higher dim
    // store sent atom indices in atom_sendlist for use in future timesteps
    // 2 bound boxes overlap if hi1 >= lo2 and lo1 <= hi2 in all dims

    x = element->element_bound_box;
    nlocal = element->nlocal;

    for (i = 0; i < nlocal; i++) {
      coord = x[i];
      checkflag = 0;

      for (m = 0; m < noverlap; m++) {
        sendflag[m] = -1;
        lobox = minsendbox[m];
        hibox = maxsendbox[m];
        if (minsendbox_flag[m])
          if (coord[3] >= lobox[0] && coord[0] <= lobox[3] &&
              coord[4] >= lobox[1] && coord[1] <= lobox[4] &&
              coord[5] >= lobox[2] && coord[2] <= lobox[5]) {
            if (elem_sendnum[m] == elem_maxsendlist[m])
              grow_elem_list(m, elem_sendnum[m]);
            elem_sendlist[m][elem_sendnum[m]++] = i;
            sendflag[m] = 1;
            continue;
          }
        if (coord[3] >= hibox[0] && coord[0] <= hibox[3] &&
            coord[4] >= hibox[1] && coord[1] <= hibox[4] &&
            coord[5] >= hibox[2] && coord[2] <= hibox[5]) {
          sendflag[m] = 0;
          checkflag = 1;
        }
      }

      if (checkflag && nforeign) {
        ibin =
            coord2bin((coord[0] + coord[3]) / 2.0, (coord[1] + coord[4]) / 2.0,
                      (coord[2] + coord[5]) / 2.0);
        for (k = 0; k < nstencil; k++) {
          for (j = fboxbinhead[ibin + stencil[k]]; j >= 0; j = fboxbins[j]) {
            m = foreign_overlap[j];
            if (sendflag[m] == 0) {
              if (coord[3] >= foreign_boxes[j][0] &&
                  coord[0] <= foreign_boxes[j][3] &&
                  coord[4] >= foreign_boxes[j][1] &&
                  coord[1] <= foreign_boxes[j][4] &&
                  coord[5] >= foreign_boxes[j][2] &&
                  coord[2] <= foreign_boxes[j][5]) {
                if (elem_sendnum[m] == elem_maxsendlist[m])
                  grow_elem_list(m, elem_sendnum[m]);
                elem_sendlist[m][elem_sendnum[m]++] = i;
                sendflag[m] = 1;
              }
            }
          }
        }
      }
    }

    // send elem_sendnum counts to procs who recv from me except self
    // copy data to self if overlapself is set

    if (overlapother) {
      for (m = 0; m < nother; m++)
        MPI_Irecv(&elem_recvnum[m], 1, MPI_INT, overlap[m], 0, world,
                  &requests[m]);
      for (m = 0; m < nother; m++)
        MPI_Send(&elem_sendnum[m], 1, MPI_INT, overlap[m], 0, world);
    }
    if (overlapself)
      for (m = nother; m < noverlap; m++)
        elem_recvnum[m] = elem_sendnum[m];
    if (overlapother)
      MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);

    // setup other per proc values from elem_sendnum and elem_recvnum

    nscountall = nrcountall = 0;
    for (m = 0; m < noverlap; m++) {
      nscount = elem_sendnum[m];
      nrcount = elem_recvnum[m];
      smaxone = MAX(smaxone, nscount);
      rmaxone = MAX(rmaxone, nrcount);
      nscountall += nscount;
      nrcountall += nrcount;

      elem_size_reverse_recv[m] = nscount * elem_size_reverse;
      if (m == 0)
        elem_reverse_recv_offset[0] = 0;
      else
        elem_reverse_recv_offset[m] =
            elem_reverse_recv_offset[m - 1] + elem_sendnum[m - 1];

      elem_size_forward_recv[m] = nrcount * elem_size_forward;
      elem_size_reverse_send[m] = nrcount * elem_size_reverse;
      if (m == 0) {
        elem_firstrecv[0] = nlocal;
        elem_forward_recv_offset[0] = 0;
      } else {
        elem_firstrecv[m] = elem_firstrecv[m - 1] + elem_recvnum[m - 1];
        elem_forward_recv_offset[m] =
            elem_forward_recv_offset[m - 1] + elem_recvnum[m - 1];
      }
    }

    smaxall = MAX(smaxall, nscountall);
    rmaxall = MAX(rmaxall, nrcountall);

    // insure send/recv buffers are large enough for this border comm swap

    if (smaxone * elem_size_border > maxsend)
      grow_send(smaxone * elem_size_border, 0);
    if (rmaxall * elem_size_border > maxrecv)
      grow_recv(rmaxall * elem_size_border);

    // swap elems with other procs using pack_border(), unpack_border()
    // use Waitall() instead of Waitany() because calls to unpack_border()
    //   must increment per-elem arrays in ascending order

    if (overlapother) {
      for (m = 0; m < nother; m++)
        MPI_Irecv(&buf_recv[elem_size_border * elem_forward_recv_offset[m]],
                  elem_recvnum[m] * elem_size_border, MPI_DOUBLE, overlap[m], 0,
                  world, &requests[m]);
      for (m = 0; m < nother; m++) {
        n = evec->pack_border(elem_sendnum[m], elem_sendlist[m], buf_send,
                              pbc_flag[m], pbc[m]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[m], 0, world);
      }
    }
    if (overlapself) {
      for (m = nother; m < noverlap; m++) {
        evec->pack_border(elem_sendnum[m], elem_sendlist[m], buf_send,
                          pbc_flag[m], pbc[m]);
        evec->unpack_border(elem_recvnum[m], elem_firstrecv[m], buf_send);
      }
    }

    if (overlapother) {
      MPI_Waitall(nother, requests, MPI_STATUS_IGNORE);
      for (m = 0; m < nother; m++)
        evec->unpack_border(
            elem_recvnum[m], elem_firstrecv[m],
            &buf_recv[elem_size_border * elem_forward_recv_offset[m]]);
    }

    // increment ghost elems

    if (noverlap)
      element->nghost =
          elem_forward_recv_offset[noverlap - 1] + elem_recvnum[noverlap - 1];

    // insure send/recv buffers are long enough for all forward & reverse comm
    // send buf is for one forward or reverse sends to one proc
    // recv buf is for all forward or reverse recvs in one swap

    int max = MAX(elem_maxforward * smaxone, elem_maxreverse * rmaxone);
    if (max > maxsend)
      grow_send(max, 0);
    max = MAX(elem_maxforward * rmaxall, elem_maxreverse * smaxall);
    if (max > maxrecv)
      grow_recv(max);

    // reset global->local map

    if (elem_map_style)
      element->map_set();
  }

  delete[] sendflag;
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_pair(Pair *pair) {
  int i, irecv, n, nother;

  int nsize;

  nother = noverlap - overlapself;

  // forward communicate atoms from pair

  if (atom->natoms) {
    nsize = pair->comm_atom_forward;
    for (i = 0; i < nother; i++)
      MPI_Irecv(&buf_recv[nsize * atom_forward_recv_offset[i]],
                nsize * atom_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                &requests[i]);
    for (i = 0; i < nother; i++) {
      n = pair->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i],
                                       buf_send, pbc_flag[i], pbc[i]);
      MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
    }

    for (i = nother; i < noverlap; i++) {
      pair->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i], buf_send,
                                   pbc_flag[i], pbc[i]);
      pair->unpack_atom_forward_comm(atom_recvnum[i], atom_firstrecv[i],
                                     buf_send);
    }

    for (i = 0; i < nother; i++) {
      MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
      pair->unpack_atom_forward_comm(
          atom_recvnum[irecv], atom_firstrecv[irecv],
          &buf_recv[nsize * atom_forward_recv_offset[irecv]]);
    }
  }

  // forward communicate elements from pair

  if (element->nelements) {
    nsize = pair->comm_elem_forward;
    for (i = 0; i < nother; i++)
      MPI_Irecv(&buf_recv[nsize * elem_forward_recv_offset[i]],
                nsize * elem_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                &requests[i]);
    for (i = 0; i < nother; i++) {
      n = pair->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i],
                                       buf_send, pbc_flag[i], pbc[i]);
      MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
    }

    for (i = nother; i < noverlap; i++) {
      pair->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i], buf_send,
                                   pbc_flag[i], pbc[i]);
      pair->unpack_elem_forward_comm(elem_recvnum[i], elem_firstrecv[i],
                                     buf_send);
    }

    for (i = 0; i < nother; i++) {
      MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
      pair->unpack_elem_forward_comm(
          elem_recvnum[irecv], elem_firstrecv[irecv],
          &buf_recv[nsize * elem_forward_recv_offset[irecv]]);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_pair(Pair *pair) {
  int i, irecv, n, nother;

  int nsize;

  // reverse communicate atoms from pair

  nother = noverlap - overlapself;

  if (atom->natoms) {
    nsize = MAX(pair->comm_atom_reverse, pair->comm_atom_reverse_off);

    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * atom_reverse_recv_offset[i]],
                  nsize * atom_sendnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = pair->pack_atom_reverse_comm(atom_recvnum[i], atom_firstrecv[i],
                                         buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        pair->pack_atom_reverse_comm(atom_recvnum[i], atom_firstrecv[i],
                                     buf_send);
        pair->unpack_atom_reverse_comm(atom_sendnum[i], atom_sendlist[i],
                                       buf_send);
      }

      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        pair->unpack_atom_reverse_comm(
            atom_sendnum[irecv], atom_sendlist[irecv],
            &buf_recv[nsize * atom_reverse_recv_offset[irecv]]);
      }
    }
  }

  // reverse communicate elements from pair

  if (element->nelements) {
    nsize = MAX(pair->comm_elem_reverse, pair->comm_elem_reverse_off);

    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * elem_reverse_recv_offset[i]],
                  nsize * elem_sendnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = pair->pack_elem_reverse_comm(elem_recvnum[i], elem_firstrecv[i],
                                         buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        pair->pack_elem_reverse_comm(elem_recvnum[i], elem_firstrecv[i],
                                     buf_send);
        pair->unpack_elem_reverse_comm(elem_sendnum[i], elem_sendlist[i],
                                       buf_send);
      }

      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        pair->unpack_elem_reverse_comm(
            elem_sendnum[irecv], elem_sendlist[irecv],
            &buf_recv[nsize * elem_reverse_recv_offset[irecv]]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_atom_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_atom_forward
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_fix(Fix *fix, int size) {
  error->all(FLERR, "Forward comm fix not yet supported by CommTiled");
  //  int i,irecv,n,nsize,nsend,nrecv;
  //
  //  if (size) nsize = size;
  //  else nsize = fix->comm_atom_forward;
  //
  //  for (int iswap = 0; iswap < nswap; iswap++) {
  //    nsend = noverlap - overlapself;
  //    nrecv = noverlap - overlapself;
  //
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[i]],
  //            nsize*atom_recvnum[i],
  //            MPI_DOUBLE,overlap[i],0,world,&requests[i]);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++) {
  //        n = fix->pack_atom_forward_comm(atom_sendnum[i],atom_sendlist[i],
  //            buf_send,pbc_flag[i],pbc[i]);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,overlap[i],0,world);
  //      }
  //    }
  //    if (overlapself) {
  //      fix->pack_atom_forward_comm(atom_sendnum[nsend],atom_sendlist[nsend],
  //          buf_send,pbc_flag[nsend],
  //          pbc[nsend]);
  //      fix->unpack_atom_forward_comm(atom_recvnum[nrecv],atom_firstrecv[nrecv],
  //          buf_send);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++) {
  //        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
  //        fix->unpack_atom_forward_comm(atom_recvnum[irecv],atom_firstrecv[irecv],
  //            &buf_recv[nsize*
  //            atom_forward_recv_offset[irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_atom_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
   some are smaller than max stored in its comm_atom_forward
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix(Fix *fix, int size) {
  error->all(FLERR, "Reverse comm fix not yet supported by CommTiled");
  //  int i,irecv,n,nsize,nsend,nrecv;
  //
  //  if (size) nsize = size;
  //  else nsize = fix->comm_atom_reverse;
  //
  //  for (int iswap = nswap-1; iswap >= 0; iswap--) {
  //    nsend = noverlap - overlapself;
  //    nrecv = noverlap - overlapself;
  //
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[i]],
  //            nsize*atom_sendnum[i],MPI_DOUBLE,
  //            overlap[i],0,world,&requests[i]);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++) {
  //        n = fix->pack_atom_reverse_comm(atom_recvnum[i],atom_firstrecv[i],
  //            buf_send);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,overlap[i],0,world);
  //      }
  //    }
  //    if (overlapself) {
  //      fix->pack_atom_reverse_comm(atom_recvnum[nrecv],atom_firstrecv[nrecv],
  //          buf_send);
  //      fix->unpack_atom_reverse_comm(atom_sendnum[nsend],atom_sendlist[nsend],
  //          buf_send);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++) {
  //        MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
  //        fix->unpack_atom_reverse_comm(atom_sendnum[irecv],atom_sendlist[irecv],
  //            &buf_recv[nsize*
  //            atom_reverse_recv_offset[irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for all pack sizes to insure buf_send is big enough
   handshake sizes before irregular comm to insure buf_recv is big enough
NOTE: how to setup one big buf recv with correct offsets ??
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix_variable(Fix *fix) {
  error->all(FLERR, "Reverse comm fix variable not yet supported by CommTiled");
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_compute(Compute *compute) {
  int i, irecv, n, nother;

  int nsize;

  nother = noverlap - overlapself;
  if (atom->natoms) {
    nsize = compute->comm_atom_forward;
    if (nsize) {

      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * atom_forward_recv_offset[i]],
                  nsize * atom_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = compute->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i],
                                            buf_send, pbc_flag[i], pbc[i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        compute->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i],
                                        buf_send, pbc_flag[i], pbc[i]);
        compute->unpack_atom_forward_comm(atom_recvnum[i], atom_firstrecv[i],
                                          buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_atom_forward_comm(
            atom_recvnum[irecv], atom_firstrecv[irecv],
            &buf_recv[nsize * atom_forward_recv_offset[irecv]]);
      }
    }
  }

  if (element->nelements) {
    nsize = compute->comm_elem_forward;
    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * elem_forward_recv_offset[i]],
                  nsize * elem_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = compute->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i],
                                            buf_send, pbc_flag[i], pbc[i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        compute->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i],
                                        buf_send, pbc_flag[i], pbc[i]);
        compute->unpack_elem_forward_comm(elem_recvnum[i], elem_firstrecv[i],
                                          buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_elem_forward_comm(
            elem_recvnum[irecv], elem_firstrecv[irecv],
            &buf_recv[nsize * elem_forward_recv_offset[irecv]]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_compute(Compute *compute) {
  int i, irecv, n, nother;

  int nsize;

  // reverse communicate atoms from compute

  nother = noverlap - overlapself;
  if (atom->natoms) {
    nsize = compute->comm_atom_reverse;
    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * atom_reverse_recv_offset[i]],
                  nsize * atom_sendnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = compute->pack_atom_reverse_comm(atom_recvnum[i], atom_firstrecv[i],
                                            buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        compute->pack_atom_reverse_comm(atom_recvnum[i], atom_firstrecv[i],
                                        buf_send);
        compute->unpack_atom_reverse_comm(atom_sendnum[i], atom_sendlist[i],
                                          buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_atom_reverse_comm(
            atom_sendnum[irecv], atom_sendlist[irecv],
            &buf_recv[nsize * atom_reverse_recv_offset[irecv]]);
      }
    }
  }

  // reverse communicate elements from compute

  if (element->nelements) {
    nsize = compute->comm_elem_reverse;
    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * elem_reverse_recv_offset[i]],
                  nsize * elem_sendnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = compute->pack_elem_reverse_comm(elem_recvnum[i], elem_firstrecv[i],
                                            buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        compute->pack_elem_reverse_comm(elem_recvnum[i], elem_firstrecv[i],
                                        buf_send);
        compute->unpack_elem_reverse_comm(elem_sendnum[i], elem_sendlist[i],
                                          buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_elem_reverse_comm(
            elem_sendnum[irecv], elem_sendlist[irecv],
            &buf_recv[nsize * elem_reverse_recv_offset[irecv]]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_dump(Dump *dump) {
  int i, irecv, n, nother;
  int nsize;

  nother = noverlap - overlapself;
  if (atom->natoms) {
    nsize = dump->comm_atom_forward;
    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * atom_forward_recv_offset[i]],
                  nsize * atom_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = dump->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i],
                                         buf_send, pbc_flag[i], pbc[i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        dump->pack_atom_forward_comm(atom_sendnum[i], atom_sendlist[i],
                                     buf_send, pbc_flag[i], pbc[i]);
        dump->unpack_atom_forward_comm(atom_recvnum[i], atom_firstrecv[i],
                                       buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        dump->unpack_atom_forward_comm(
            atom_recvnum[irecv], atom_firstrecv[irecv],
            &buf_recv[nsize * atom_forward_recv_offset[irecv]]);
      }
    }
  }

  if (element->nelements) {
    nsize = dump->comm_elem_forward;
    if (nsize) {
      for (i = 0; i < nother; i++)
        MPI_Irecv(&buf_recv[nsize * elem_forward_recv_offset[i]],
                  nsize * elem_recvnum[i], MPI_DOUBLE, overlap[i], 0, world,
                  &requests[i]);
      for (i = 0; i < nother; i++) {
        n = dump->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i],
                                         buf_send, pbc_flag[i], pbc[i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, overlap[i], 0, world);
      }

      for (i = nother; i < noverlap; i++) {
        dump->pack_elem_forward_comm(elem_sendnum[i], elem_sendlist[i],
                                     buf_send, pbc_flag[i], pbc[i]);
        dump->unpack_elem_forward_comm(elem_recvnum[i], elem_firstrecv[i],
                                       buf_send);
      }
      for (i = 0; i < nother; i++) {
        MPI_Waitany(nother, requests, &irecv, MPI_STATUS_IGNORE);
        dump->unpack_elem_forward_comm(
            elem_recvnum[irecv], elem_firstrecv[irecv],
            &buf_recv[nsize * elem_forward_recv_offset[irecv]]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer limit
   ------------------------------------------------------------------------- */

void CommTiled::reverse_comm_dump(Dump *dump) {
  error->all(FLERR, "Reverse comm dump not yet supported by CommTiled");
  //  int i,irecv,n,nother,nrecv;
  //
  //  int nsize = dump->comm_atom_reverse;
  //
  //  for (int iswap = nswap-1; iswap >= 0; iswap--) {
  //    nsend = noverlap - overlapself;
  //    nrecv = noverlap - overlapself;
  //
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_reverse_recv_offset[i]],
  //            nsize*atom_sendnum[i],MPI_DOUBLE,
  //            overlap[i],0,world,&requests[i]);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++) {
  //        n = dump->pack_atom_reverse_comm(atom_recvnum[i],atom_firstrecv[i],
  //            buf_send);
  //        MPI_Send(buf_send,n,MPI_DOUBLE,overlap[i],0,world);
  //      }
  //    }
  //    if (overlapself) {
  //      dump->pack_atom_reverse_comm(atom_recvnum[nrecv],atom_firstrecv[nrecv],
  //          buf_send);
  //      dump->unpack_atom_reverse_comm(atom_sendnum[nsend],atom_sendlist[nsend],
  //          buf_send);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++) {
  //        MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
  //        dump->unpack_atom_reverse_comm(atom_sendnum[irecv],atom_sendlist[irecv],
  //            &buf_recv[nsize*
  //            atom_reverse_recv_offset[irecv]]);
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------
   forward communication of Nsize values in per-atom array
   ------------------------------------------------------------------------- */

void CommTiled::forward_comm_array(int nsize, double **array) {
  error->all(FLERR, "Forward comm array not yet supported by CommTiled");
  //  int i,j,k,m,iatom,last,irecv,nsend,nrecv;
  //
  //  // insure send/recv bufs are big enough for nsize
  //  // based on smaxone/rmaxall from most recent borders() invocation
  //
  //  if (nsize > atom_maxforward) {
  //    atom_maxforward = nsize;
  //    if (atom_maxforward*smaxone > maxsend)
  //    grow_send(atom_maxforward*smaxone,0); if (atom_maxforward*rmaxall >
  //    maxrecv) grow_recv(atom_maxforward*rmaxall);
  //  }
  //
  //  for (int iswap = 0; iswap < nswap; iswap++) {
  //    nsend = noverlap - overlapself;
  //    nrecv = noverlap - overlapself;
  //
  //    MPI_Barrier(world);
  //
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++)
  //        MPI_Irecv(&buf_recv[nsize*atom_forward_recv_offset[i]],
  //            nsize*atom_recvnum[i],
  //            MPI_DOUBLE,overlap[i],0,world,&requests[i]);
  //    }
  //    if (overlapother) {
  //      for (i = 0; i < nsend; i++) {
  //        m = 0;
  //        for (iatom = 0; iatom < atom_sendnum[i]; iatom++) {
  //          j = atom_sendlist[i][iatom];
  //          for (k = 0; k < nsize; k++)
  //            buf_send[m++] = array[j][k];
  //        }
  //        MPI_Send(buf_send,nsize*atom_sendnum[i],
  //            MPI_DOUBLE,overlap[i],0,world);
  //      }
  //    }
  //    if (overlapself) {
  //      m = 0;
  //      for (iatom = 0; iatom < atom_sendnum[nsend]; iatom++) {
  //        j = atom_sendlist[nsend][iatom];
  //        for (k = 0; k < nsize; k++)
  //          buf_send[m++] = array[j][k];
  //      }
  //      m = 0;
  //      last = atom_firstrecv[nrecv] + atom_recvnum[nrecv];
  //      for (iatom = atom_firstrecv[nrecv]; iatom < last; iatom++)
  //        for (k = 0; k < nsize; k++)
  //          array[iatom][k] = buf_send[m++];
  //    }
  //
  //    if (overlapother) {
  //      for (i = 0; i < nrecv; i++) {
  //        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
  //        m = nsize*atom_forward_recv_offset[irecv];
  //        last = atom_firstrecv[irecv] + atom_recvnum[irecv];
  //        for (iatom = atom_firstrecv[irecv]; iatom < last; iatom++)
  //          for (k = 0; k < nsize; k++)
  //            array[iatom][k] = buf_recv[m++];
  //      }
  //    }
  //  }
}

/* ----------------------------------------------------------------------*/

void CommTiled::grow_overlap_array() {
  maxoverlap += DELTA_PROCS;
  memory->grow(overlap, maxoverlap, "comm:overlap");
  memory->grow(overlap_pbc, maxoverlap, 3, "comm:overlap_pbc");
  memory->grow(overlap2box, maxoverlap, 6, "comm:overlap2box");
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   box is owned by me and extends in dim
   flag =
   + 0: normal drop for setup_exchange
   + 1: box 1 drop
   + 2: box 2 drop
   ------------------------------------------------------------------------- */

void CommTiled::box_drop_brick(int idim, double *lo, double *hi,
                               int *current_pbc, int &indexme, int flag) {
  // NOTE: this is not triclinic compatible
  // NOTE: these error messages are internal sanity checks
  //       should not occur, can be removed at some point

  if (flag) {
    int pbc_bit;
    double subbox_size[3];
    double procbox_lo[3];
    double procbox_hi[3];
    int proc;
    double *split_array[3];
    subbox_size[0] = subhi[0] - sublo[0];
    subbox_size[1] = subhi[1] - sublo[1];
    subbox_size[2] = subhi[2] - sublo[2];
    int xi[3], xf[3];
    split_array[0] = xsplit;
    split_array[1] = ysplit;
    split_array[2] = zsplit;
    for (int i = 0; i < dimension; i++) {
      xi[i] = static_cast<int>((lo[i] - boxlo[i]) / subbox_size[i]);
      xf[i] = static_cast<int>((hi[i] - boxlo[i]) / subbox_size[i]);

      if (xf[i] == procgrid[i])
        xf[i] = procgrid[i] - 1;
      if (xi[i] < 0)
        xi[i] = 0;
    }

    for (int x = xi[0]; x <= xf[0]; x++) {
      for (int y = xi[1]; y <= xf[1]; y++) {
        for (int z = xi[2]; z <= xf[2]; z++) {
          for (int dim = 0; dim < dimension; dim++) {
            int split_index;
            if (dim == 0)
              split_index = x;
            if (dim == 1)
              split_index = y;
            if (dim == 2)
              split_index = z;
            procbox_lo[dim] =
                boxlo[dim] + prd[dim] * split_array[dim][split_index];
            if (split_index < procgrid[dim] - 1)
              procbox_hi[dim] =
                  boxlo[dim] + prd[dim] * split_array[dim][split_index + 1];
            else
              procbox_hi[dim] = boxhi[dim];
          }
          proc = grid2proc[x][y][z];

          // skip repeats from overlap

          if (flag == 1) {
            if (proc == me)
              continue;
            pbc_bit = 13;
          } else {
            pbc_bit = ((current_pbc[0] + 1) + 3 * (current_pbc[1] + 1) +
                       9 * (current_pbc[2] + 1));
          }
          if (!(overlap_repeat[proc] & (1 << pbc_bit)))
            overlap_repeat[proc] += 1 << pbc_bit;
          else
            continue;

          if (noverlap >= maxoverlap)
            grow_overlap_array();
          overlap2box[noverlap][0] = procbox_lo[0];
          overlap2box[noverlap][1] = procbox_lo[1];
          overlap2box[noverlap][2] = procbox_lo[2];
          overlap2box[noverlap][3] = procbox_hi[0];
          overlap2box[noverlap][4] = procbox_hi[1];
          overlap2box[noverlap][5] = procbox_hi[2];
          if (flag == 2) {
            overlap_pbc[noverlap][0] = current_pbc[0];
            overlap_pbc[noverlap][1] = current_pbc[1];
            overlap_pbc[noverlap][2] = current_pbc[2];
          } else if (flag == 1) {
            overlap_pbc[noverlap][0] = 0;
            overlap_pbc[noverlap][1] = 0;
            overlap_pbc[noverlap][2] = 0;
          }
          overlap[noverlap++] = proc;
        }
      }
    }

  } else {
    int index = -1;
    int dir;
    if (hi[idim] == sublo[idim]) {
      index = myloc[idim] - 1;
      dir = -1;
    } else if (lo[idim] == subhi[idim]) {
      index = myloc[idim] + 1;
      dir = 1;
    } else if (hi[idim] == boxhi[idim]) {
      index = procgrid[idim] - 1;
      dir = -1;
    } else if (lo[idim] == boxlo[idim]) {
      index = 0;
      dir = 1;
    } else
      error->one(FLERR, "Comm tiled mis-match in box drop brick");

    int other1, other2, proc;
    double lower, upper;
    double *split;

    if (idim == 0) {
      other1 = myloc[1];
      other2 = myloc[2];
      split = xsplit;
    } else if (idim == 1) {
      other1 = myloc[0];
      other2 = myloc[2];
      split = ysplit;
    } else {
      other1 = myloc[0];
      other2 = myloc[1];
      split = zsplit;
    }

    if (index < 0 || index > procgrid[idim])
      error->one(FLERR, "Comm tiled invalid index in box drop brick");

    while (1) {
      lower = boxlo[idim] + prd[idim] * split[index];
      if (index < procgrid[idim] - 1)
        upper = boxlo[idim] + prd[idim] * split[index + 1];
      else
        upper = boxhi[idim];
      if (lower >= hi[idim] || upper <= lo[idim])
        break;

      if (idim == 0)
        proc = grid2proc[index][other1][other2];
      else if (idim == 1)
        proc = grid2proc[other1][index][other2];
      else
        proc = grid2proc[other1][other2][index];

      if (noverlap == maxoverlap)
        grow_overlap_array();
      if (proc == me)
        indexme = noverlap;
      overlap[noverlap++] = proc;
      index += dir;
      if (index < 0 || index >= procgrid[idim])
        break;
    }
  }
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
   flag =
   + 0: normal drop for setup_exchange
   + 1: box 1 drop
   + 2: box 2 drop
   ------------------------------------------------------------------------- */

void CommTiled::box_drop_tiled(int idim, double *lo, double *hi,
                               int *current_pbc, int &indexme, int flag) {
  box_drop_tiled_recurse(lo, hi, 0, nprocs - 1, current_pbc, indexme, flag);
}

void CommTiled::box_drop_tiled_recurse(double *lo, double *hi, int proclower,
                                       int procupper, int *current_pbc,
                                       int &indexme, int flag) {
  // end recursion when partition is a single proc
  // add proc to overlap list

  if (proclower == procupper) {
    if (noverlap == maxoverlap)
      grow_overlap_array();

    if (flag) {

      // skip repeat overlap

      int pbc_bit;
      if (flag == 1) {
        if (proclower == me)
          return;
        pbc_bit = 13;
      } else {
        pbc_bit = ((current_pbc[0] + 1) + 3 * (current_pbc[1] + 1) +
                   9 * (current_pbc[2] + 1));
      }
      if (!(overlap_repeat[proclower] & (1 << pbc_bit)))
        overlap_repeat[proclower] += 1 << pbc_bit;
      else
        return;

      // set pbc image for this overlap
      // no need to check overlap with self here
      // will be sorted later

      if (flag == 2) {
        overlap_pbc[noverlap][0] = current_pbc[0];
        overlap_pbc[noverlap][1] = current_pbc[1];
        overlap_pbc[noverlap][2] = current_pbc[2];
      } else if (flag == 1) {
        overlap_pbc[noverlap][0] = 0;
        overlap_pbc[noverlap][1] = 0;
        overlap_pbc[noverlap][2] = 0;
      }

    } else {
      if (proclower == me)
        indexme = noverlap;
    }
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use > and < criteria and add +/- EPSILON so does not include a box it only
  // touches procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim] * rcbinfo[procmid].cutfrac;

  if (lo[idim] < cut)
    box_drop_tiled_recurse(lo, hi, proclower, procmid - 1, current_pbc, indexme,
                           flag);
  if (hi[idim] > cut)
    box_drop_tiled_recurse(lo, hi, procmid, procupper, current_pbc, indexme,
                           flag);
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
   ------------------------------------------------------------------------- */

void CommTiled::box_other_brick(int idim, int idir, int proc, double *lo,
                                double *hi, int overlap_counter) {
  lo[0] = overlap2box[overlap_counter][0];
  lo[1] = overlap2box[overlap_counter][1];
  lo[2] = overlap2box[overlap_counter][2];
  hi[0] = overlap2box[overlap_counter][3];
  hi[1] = overlap2box[overlap_counter][4];
  hi[2] = overlap2box[overlap_counter][5];
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
   ------------------------------------------------------------------------- */

void CommTiled::box_other_tiled(int idim, int idir, int proc, double *lo,
                                double *hi, int overlap_counter) {
  double(*split)[2] = rcbinfo[proc].mysplit;

  lo[0] = boxlo[0] + prd[0] * split[0][0];
  if (split[0][1] < 1.0)
    hi[0] = boxlo[0] + prd[0] * split[0][1];
  else
    hi[0] = boxhi[0];

  lo[1] = boxlo[1] + prd[1] * split[1][0];
  if (split[1][1] < 1.0)
    hi[1] = boxlo[1] + prd[1] * split[1][1];
  else
    hi[1] = boxhi[1];

  lo[2] = boxlo[2] + prd[2] * split[2][0];
  if (split[2][1] < 1.0)
    hi[2] = boxlo[2] + prd[2] * split[2][1];
  else
    hi[2] = boxhi[2];
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
   procneigh stores 6 procs that touch me
   ------------------------------------------------------------------------- */

int CommTiled::box_touch_brick(int proc, int idim, int idir) {
  if (procneigh[idim][idir] == proc)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if proc's box touches me, else 0
   ------------------------------------------------------------------------- */

int CommTiled::box_touch_tiled(int proc, int idim, int idir) {
  // sending to left
  // only touches if proc hi = my lo, or if proc hi = boxhi and my lo = boxlo

  if (idir == 0) {
    if (rcbinfo[proc].mysplit[idim][1] == rcbinfo[me].mysplit[idim][0])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][1] == 1.0 &&
             rcbinfo[me].mysplit[idim][0] == 0.0)
      return 1;

    // sending to right
    // only touches if proc lo = my hi, or if proc lo = boxlo and my hi = boxhi

  } else {
    if (rcbinfo[proc].mysplit[idim][0] == rcbinfo[me].mysplit[idim][1])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][0] == 0.0 &&
             rcbinfo[me].mysplit[idim][1] == 1.0)
      return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_brick(int idim, double *x) {
  if (closer_subbox_edge(idim, x))
    return procneigh[idim][1];
  return procneigh[idim][0];
}

/* ----------------------------------------------------------------------
   determine which proc owns point x via recursion thru RCB tree
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_tiled(int idim, double *x) {
  double xnew[3];
  xnew[0] = x[0];
  xnew[1] = x[1];
  xnew[2] = x[2];

  if (idim == 0) {
    if (xnew[1] < sublo[1] || xnew[1] > subhi[1]) {
      if (closer_subbox_edge(1, x))
        xnew[1] = subhi[1];
      else
        xnew[1] = sublo[1];
    }
  }
  if (idim <= 1) {
    if (xnew[2] < sublo[2] || xnew[2] > subhi[2]) {
      if (closer_subbox_edge(2, x))
        xnew[2] = subhi[2];
      else
        xnew[2] = sublo[2];
    }
  }

  int proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
  if (proc == me)
    return me;

  if (idim == 0) {
    int done = 1;
    if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
      xnew[1] -= EPSILON * (subhi[1] - sublo[1]);
      done = 0;
    }
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
      done = 0;
    }
    if (!done) {
      proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
      done = 1;
      if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
        xnew[1] -= EPSILON * (subhi[1] - sublo[1]);
        done = 0;
      }
      if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
        xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
        done = 0;
      }
      if (!done)
        proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
    }
  } else if (idim == 1) {
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
      proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
    }
  }

  return proc;
}

/* ----------------------------------------------------------------------
   recursive point drop thru RCB tree
   ------------------------------------------------------------------------- */

int CommTiled::point_drop_tiled_recurse(double *x, int proclower,
                                        int procupper) {
  // end recursion when partition is a single proc
  // return proc

  if (proclower == procupper)
    return proclower;

  // drop point on side of cut it is on
  // use < criterion so point is not on high edge of proc sub-domain
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim] * rcbinfo[procmid].cutfrac;

  if (x[idim] < cut)
    return point_drop_tiled_recurse(x, proclower, procmid - 1);
  else
    return point_drop_tiled_recurse(x, procmid, procupper);
}

/* ----------------------------------------------------------------------
   assume x[idim] is outside subbox bounds in same dim
   ------------------------------------------------------------------------- */

int CommTiled::closer_subbox_edge(int idim, double *x) {
  double deltalo, deltahi;

  if (sublo[idim] == boxlo[idim])
    deltalo = fabs(x[idim] - prd[idim] - sublo[idim]);
  else
    deltalo = fabs(x[idim] - sublo[idim]);

  if (subhi[idim] == boxhi[idim])
    deltahi = fabs(x[idim] + prd[idim] - subhi[idim]);
  else
    deltahi = fabs(x[idim] - subhi[idim]);

  if (deltalo < deltahi)
    return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   if RCB decomp exists and just changed, gather needed global RCB info
   ------------------------------------------------------------------------- */

void CommTiled::coord2proc_setup() {
  if (!rcbnew)
    return;

  if (!rcbinfo)
    rcbinfo =
        (RCBinfo *)memory->smalloc(nprocs * sizeof(RCBinfo), "comm:rcbinfo");
  rcbnew = 0;
  RCBinfo rcbone;
  memcpy(&rcbone.mysplit[0][0], &mysplit[0][0], 6 * sizeof(double));
  rcbone.cutfrac = rcbcutfrac;
  rcbone.dim = rcbcutdim;
  MPI_Allgather(&rcbone, sizeof(RCBinfo), MPI_CHAR, rcbinfo, sizeof(RCBinfo),
                MPI_CHAR, world);
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with coord x[3] based on current decomp
   x will be in box (orthogonal) or lamda coords (triclinic)
   if layout = UNIFORM or NONUNIFORM, invoke parent method
   if layout = TILED, use point_drop_recurse()
   return owning proc ID, ignore igx,igy,igz
   ------------------------------------------------------------------------- */

int CommTiled::coord2proc(double *x, int &igx, int &igy, int &igz) {
  if (layout != Comm::LAYOUT_TILED)
    return Comm::coord2proc(x, igx, igy, igz);
  return point_drop_tiled_recurse(x, 0, nprocs - 1);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
   ------------------------------------------------------------------------- */

void CommTiled::grow_send(int n, int flag) {
  maxsend = static_cast<int>(BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send, maxsend + bufextra, "comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send, maxsend + bufextra, "comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_recv(int n) {
  maxrecv = static_cast<int>(BUFFACTOR * n);
  if (maxrecv < 0 || maxrecv >= MAXBIGINT)
    error->one(FLERR, "Data size too large for comm:buf_recv");
  memory->destroy(buf_recv);
  memory->create(buf_recv, maxrecv, "comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the atom_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_atom_list(int iwhich, int n) {
  atom_maxsendlist[iwhich] = static_cast<int>(BUFFACTOR * n);
  memory->grow(atom_sendlist[iwhich], atom_maxsendlist[iwhich],
               "comm:atom_sendlist[i]");
}

/* ----------------------------------------------------------------------
   realloc the size of the elem_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_elem_list(int iwhich, int n) {
  elem_maxsendlist[iwhich] = static_cast<int>(BUFFACTOR * n);
  memory->grow(elem_sendlist[iwhich], elem_maxsendlist[iwhich],
               "comm:elem_sendlist[i]");
}

/* ----------------------------------------------------------------------
   realloc the size of the overlap_sendlist as needed with BUFFACTOR
   ------------------------------------------------------------------------- */

void CommTiled::grow_overlap_list(int iwhich, int n) {
  overlap_maxsendlist[iwhich] = static_cast<int>(BUFFACTOR * n);
  memory->grow(overlap_sendlist[iwhich], overlap_maxsendlist[iwhich],
               "comm:overlap_sendlist[i]");
}

/* ----------------------------------------------------------------------
   grow info for overlap, to allow for N procs to communicate with
   ------------------------------------------------------------------------- */

void CommTiled::grow_overlap(int n, int nold) {
  delete[] atom_sendnum;
  delete[] atom_recvnum;
  delete[] atom_size_forward_recv;
  delete[] atom_firstrecv;
  delete[] atom_size_reverse_send;
  delete[] atom_size_reverse_recv;
  delete[] atom_forward_recv_offset;
  delete[] atom_reverse_recv_offset;
  atom_sendnum = new int[n];
  atom_recvnum = new int[n];
  atom_size_forward_recv = new int[n];
  atom_firstrecv = new int[n];
  atom_size_reverse_send = new int[n];
  atom_size_reverse_recv = new int[n];
  atom_forward_recv_offset = new int[n];
  atom_reverse_recv_offset = new int[n];

  delete[] elem_sendnum;
  delete[] elem_recvnum;
  delete[] elem_size_forward_recv;
  delete[] elem_firstrecv;
  delete[] elem_size_reverse_send;
  delete[] elem_size_reverse_recv;
  delete[] elem_forward_recv_offset;
  delete[] elem_reverse_recv_offset;
  elem_sendnum = new int[n];
  elem_recvnum = new int[n];
  elem_size_forward_recv = new int[n];
  elem_firstrecv = new int[n];
  elem_size_reverse_send = new int[n];
  elem_size_reverse_recv = new int[n];
  elem_forward_recv_offset = new int[n];
  elem_reverse_recv_offset = new int[n];

  delete[] pbc_flag;
  pbc_flag = new int[n];
  memory->destroy(pbc);
  memory->create(pbc, n, 6, "comm:pbc_flag");
  memory->destroy(minsendbox);
  memory->create(minsendbox, n, 6, "comm:minsendbox");
  memory->destroy(minsendbox_flag);
  memory->create(minsendbox_flag, n, "comm:minsendbox_flag");
  memory->destroy(maxsendbox);
  memory->create(maxsendbox, n, 6, "comm:maxsendbox");

  delete[] atom_maxsendlist;
  delete[] elem_maxsendlist;
  atom_maxsendlist = new int[n];
  elem_maxsendlist = new int[n];

  for (int i = 0; i < nold; i++) {
    memory->destroy(atom_sendlist[i]);
    memory->destroy(elem_sendlist[i]);
  }

  delete[] atom_sendlist;
  delete[] elem_sendlist;
  atom_sendlist = new int *[n];
  elem_sendlist = new int *[n];

  for (int i = 0; i < n; i++) {
    atom_maxsendlist[i] = BUFMIN;
    elem_maxsendlist[i] = BUFMIN;
    memory->create(atom_sendlist[i], BUFMIN, "comm:atom_sendlist[i]");
    memory->create(elem_sendlist[i], BUFMIN, "comm:elem_sendlist[i]");
  }

  delete[] overlap_sendnum;
  delete[] overlap_recvnum;
  delete[] overlap_maxsendlist;
  delete[] overlap_firstrecv;
  for (int i = 0; i < nold; i++) {
    memory->destroy(overlap_sendlist[i]);
  }
  delete[] overlap_sendlist;

  overlap_sendnum = new int[n];
  overlap_recvnum = new int[n];
  overlap_maxsendlist = new int[n];
  overlap_sendlist = new int *[n];
  overlap_firstrecv = new int[n];

  for (int i = 0; i < n; i++) {
    overlap_maxsendlist[i] = BUFMIN;
    memory->create(overlap_sendlist[i], BUFMIN, "comm:overlap_sendlist[i]");
  }
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
NOTE: this method is currently not used
------------------------------------------------------------------------- */

int CommTiled::exchange_variable(int n, double * /*inbuf*/,
                                 double *& /*outbuf*/) {
  int nrecv = n;
  return nrecv;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint CommTiled::memory_usage() {
  bigint bytes = 0;
  for (int i = 0; i < nprocmax; i++) {
    bytes += memory->usage(atom_sendlist[i], atom_maxsendlist[i]);
    bytes += memory->usage(elem_sendlist[i], elem_maxsendlist[i]);
    bytes += memory->usage(overlap_sendlist[i], overlap_maxsendlist[i]);
  }
  bytes += memory->usage(buf_send, maxsend + bufextra);
  bytes += memory->usage(buf_recv, maxrecv);

  return bytes;
}

/* ----------------------------------------------------------------------
   debug function to output foreign boxes
   flag = 0: dump file
          1: tecplot file
---------------------------------------------------------------------- */

void CommTiled::debug_write_foreign_box(int flag) {
  FILE *fp;

  char file[50];
  if (flag) {
    sprintf(file, "foreign_box_%d.dat", me);
    fp = fopen(file, "w");
    fprintf(fp, "title = \"%d\"\n", me);
    fprintf(fp, "variables = \"x\", \"y\", \"z\"\n");
    fprintf(fp, "zone t = \"Coarse Element\",");
    fprintf(fp, " n = %d e = %d", nforeign * 8, nforeign);
    fprintf(fp, " datapacking = point, zonetype = febrick\n");
    for (int i = 0; i < nforeign; i++) {
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][0], foreign_boxes[i][1],
              foreign_boxes[i][2]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][3], foreign_boxes[i][1],
              foreign_boxes[i][2]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][3], foreign_boxes[i][4],
              foreign_boxes[i][2]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][0], foreign_boxes[i][4],
              foreign_boxes[i][2]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][0], foreign_boxes[i][1],
              foreign_boxes[i][5]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][3], foreign_boxes[i][1],
              foreign_boxes[i][5]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][3], foreign_boxes[i][4],
              foreign_boxes[i][5]);
      fprintf(fp, "%g %g %g\n", foreign_boxes[i][0], foreign_boxes[i][4],
              foreign_boxes[i][5]);
    }
    int count = 1;
    for (int i = 0; i < nforeign; i++) {
      for (int j = 0; j < 8; j++)
        fprintf(fp, "%d ", count++);
      fprintf(fp, "\n");
    }
    fclose(fp);

  } else {
    sprintf(file, "foreign_box_%d.atom", me);
    fp = fopen(file, "w");

    fprintf(fp, "ITEM: TIMESTEP\n");

    fprintf(fp,
            "%d"
            "\n",
            me);

    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,
            "%d"
            "\n",
            nforeign * 9);
    char boundstr[9];
    domain->boundary_string(boundstr);
    if (triclinic == 0) {
      fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
      fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
      fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
      fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
    } else {
      fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0],
              domain->boxhi_bound[0], domain->xy);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1],
              domain->boxhi_bound[1], domain->xz);
      fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2],
              domain->boxhi_bound[2], domain->yz);
    }
    fprintf(fp, "ITEM: ATOMS id type x y z tag overlap\n");

    int n = 1;
    for (int i = 0; i < nforeign; i++) {
      fprintf(fp, "%d 4 %g %g %g %d %d\n", n++,
              (foreign_boxes[i][0] + foreign_boxes[i][3]) / 2,
              (foreign_boxes[i][1] + foreign_boxes[i][4]) / 2,
              (foreign_boxes[i][2] + foreign_boxes[i][5]) / 2, foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][0],
              foreign_boxes[i][1], foreign_boxes[i][2], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][3],
              foreign_boxes[i][1], foreign_boxes[i][2], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][0],
              foreign_boxes[i][4], foreign_boxes[i][2], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][3],
              foreign_boxes[i][4], foreign_boxes[i][2], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][0],
              foreign_boxes[i][1], foreign_boxes[i][5], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][3],
              foreign_boxes[i][1], foreign_boxes[i][5], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][0],
              foreign_boxes[i][4], foreign_boxes[i][5], foreign_tag[i],
              foreign_overlap[i]);
      fprintf(fp, "%d 5 %g %g %g %d %d\n", n++, foreign_boxes[i][3],
              foreign_boxes[i][4], foreign_boxes[i][5], foreign_tag[i],
              foreign_overlap[i]);
    }

    fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   debug function to output domain boxes
   ---------------------------------------------------------------------- */

void CommTiled::debug_write_domain_box() {
  double *buf = new double[6];
  int tmp;

  buf[0] = domain->sublo[0];
  buf[1] = domain->sublo[1];
  buf[2] = domain->sublo[2];
  buf[3] = domain->subhi[0];
  buf[4] = domain->subhi[1];
  buf[5] = domain->subhi[2];

  if (me == 0) {
    FILE *fp;
    fp = fopen("domain_box_partition.dat", "w");
    fprintf(fp, "title = domain_box_partition\n");
    fprintf(fp, "variables = \"x\", \"y\", \"z\", \"procID\"\n");
    fprintf(fp, "zone t = \"Coarse Element\",");
    fprintf(fp, " n = %d e = %d", nprocs * 8, nprocs);
    fprintf(fp, " datapacking = point, zonetype = febrick\n");

    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf, 6, MPI_DOUBLE, iproc, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, iproc, 0, world);
        MPI_Wait(&request, &status);
      }
      fprintf(fp, "%g %g %g %d\n", buf[0], buf[1], buf[2], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[3], buf[1], buf[2], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[3], buf[4], buf[2], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[0], buf[4], buf[2], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[0], buf[1], buf[5], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[3], buf[1], buf[5], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[3], buf[4], buf[5], iproc);
      fprintf(fp, "%g %g %g %d\n", buf[0], buf[4], buf[5], iproc);
    }

    int count = 1;
    for (int iproc = 0; iproc < nprocs; iproc++) {
      fprintf(fp, "%d %d %d %d %d %d %d %d\n", count, count + 1, count + 2,
              count + 3, count + 4, count + 5, count + 6, count + 7);
      count += 8;
    }

    fclose(fp);
  } else {
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(buf, 6, MPI_DOUBLE, 0, 0, world);
  }

  delete[] buf;
}
