#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "compute_cna_atom.h"
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

using namespace CAC_NS;

#define MAXNEAR 14
#define MAXCOMMON 6

enum{UNKNOWN, OTHER, FCC, HCP, BCC, ICOS};
enum{NCOMMON, NBOND, MAXBOND, MINBOND};

/*  ----------------------------------------------------------------------  */

ComputeCNAAtom::ComputeCNAAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg), 
  list(NULL), atompattern(NULL), vatompattern(NULL)
{
  if (narg != 4) error->all(FLERR, "Illegal compute cna/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff;
  int adapt_cut_flag = 0;
  
  if (strcmp(arg[3], "adaptive") == 0) {
    cutoff = 0.0;
  } else if (strcmp(arg[3], "Ac/fcc") == 0) 
    cutoff = 4.53237;
  else if (strcmp(arg[3], "Ag/fcc") == 0) 
    cutoff = 3.49103;
  else if (strcmp(arg[3], "Al/fcc") == 0) 
    cutoff = 3.45689;
  else if (strcmp(arg[3], "Ar/fcc") == 0) 
    cutoff = 4.48969;
  else if (strcmp(arg[3], "Au/fcc") == 0) 
    cutoff = 3.48250;
  else if (strcmp(arg[3], "Ba/bcc") == 0) 
    cutoff = 6.05968;
  else if (strcmp(arg[3], "C/dia") == 0) 
    cutoff = 2.74223;
  else if (strcmp(arg[3], "Ca/fcc") == 0) 
    cutoff = 4.76283;
  else if (strcmp(arg[3], "Ce/fcc") == 0) 
    cutoff = 4.40434;
  else if (strcmp(arg[3], "Cr/bcc") == 0) 
    cutoff = 3.47647;
  else if (strcmp(arg[3], "Cs/bcc") == 0) 
    cutoff = 7.30300;
  else if (strcmp(arg[3], "Cu/fcc") == 0) 
    cutoff = 3.08133;
  else if (strcmp(arg[3], "Eu/bcc") == 0) 
    cutoff = 5.56476;
  else if (strcmp(arg[3], "Fe/bcc") == 0) 
    cutoff = 3.46440;
  else if (strcmp(arg[3], "Ge/dia") == 0) 
    cutoff = 4.34762;
  else if (strcmp(arg[3], "Ir/fcc") == 0) 
    cutoff = 3.27764;
  else if (strcmp(arg[3], "K/bcc") == 0) 
    cutoff = 6.31317;
  else if (strcmp(arg[3], "Kr/fcc") == 0) 
    cutoff = 4.88233;
  else if (strcmp(arg[3], "Li/bcc") == 0) 
    cutoff = 4.21280;
  else if (strcmp(arg[3], "Mo/bcc") == 0) 
    cutoff = 3.80239;
  else if (strcmp(arg[3], "Na/bcc") == 0) 
    cutoff = 5.10606;
  else if (strcmp(arg[3], "Nb/bcc") == 0) 
    cutoff = 3.98345;
  else if (strcmp(arg[3], "Ne/fcc") == 0) 
    cutoff = 3.78124;
  else if (strcmp(arg[3], "Ni/fcc") == 0) 
    cutoff = 3.00451;
  else if (strcmp(arg[3], "Pb/fcc") == 0) 
    cutoff = 4.22509;
  else if (strcmp(arg[3], "Pd/fcc") == 0) 
    cutoff = 3.32032;
  else if (strcmp(arg[3], "Pt/fcc") == 0) 
    cutoff = 3.34593;
  else if (strcmp(arg[3], "Rb/bcc") == 0) 
    cutoff = 6.74773;
  else if (strcmp(arg[3], "Rh/fcc") == 0) 
    cutoff = 3.24350;
  else if (strcmp(arg[3], "Si/dia") == 0) 
    cutoff = 4.19095;
  else if (strcmp(arg[3], "Sr/fcc") == 0) 
    cutoff = 5.18960;
  else if (strcmp(arg[3], "Ta/bcc") == 0) 
    cutoff = 3.99552;
  else if (strcmp(arg[3], "Th/fcc") == 0) 
    cutoff = 4.33605;
  else if (strcmp(arg[3], "V/bcc") == 0) 
    cutoff = 3.64546;
  else if (strcmp(arg[3], "W/bcc") == 0) 
    cutoff = 3.81446;
  else if (strcmp(arg[3], "Xe/fcc") == 0) 
    cutoff = 5.29203;
  else if (strcmp(arg[3], "Yb/fcc") == 0) 
    cutoff = 4.68601;
  else {
    cutoff = universe->numeric(FLERR, arg[3]);
    if (cutoff < 0.0) error->all(FLERR, "Illegal compute cna/atom command");
  }
  if (cutoff == 0.0) adapt_cut_flag = 1;
  cutsq = cutoff * cutoff;

  namax = nemax = 0;
  nalocalmax = nelocalmax = 0;

  // set comm size needed by this compute

  comm_atom_forward = 1;
  comm_elem_forward = element->maxucell * element->maxapc;
}

/*  ----------------------------------------------------------------------  */

ComputeCNAAtom::~ComputeCNAAtom()
{
  memory->destroy(atompattern);
  memory->destroy(vatompattern);
}

/*  ----------------------------------------------------------------------  */

void ComputeCNAAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR, "Compute cna/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute cna/atom cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if (2.0 * sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR, "Compute cna/atom cutoff may be too large to find "
        "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style, "cna/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute cna/atom defined");

  // need an occasional full vatom neighbor list

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->vatomlist = 1;
  neighbor->requests[irequest]->gausslist = 0;
  //neighbor->requests[irequest]->neighneigh = 1;
  if (cutsq > 0) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = 2.0 * sqrt(cutsq);
  }
  neighbor->requests[irequest]->sort = 1;
  neighbor->requests[irequest]->sortnum = 14;
  neighbor->requests[irequest]->group = 1;
  neighbor->requests[irequest]->groupbit = groupbit;
  neighbor->requests[irequest]->noinner = 1;


}

/*  ----------------------------------------------------------------------  */

void ComputeCNAAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/*  ----------------------------------------------------------------------  */

void ComputeCNAAtom::compute_peratom()
{
  int nfcc, nhcp, nbcc4, nbcc6, nico;
  int cj, ck, cl, cm;
  int pattern;
  int i, j, k, l;
  int ii, jj, kk, ll, n, inear;
  int num;
  int ietype, jetype;
  int iucell, jucell, iindex, jindex, iapc, japc;
  int inode, jnode;
  int ibasis, jbasis;
  int firstflag, nnearest;
  int fccflag, bccflag;
  int ncommonfcc[MAXNEAR], commonfcc[MAXNEAR][MAXCOMMON];
  int ncommonbcc[MAXNEAR], commonbcc[MAXNEAR][MAXCOMMON];

  double rcutfccsq, rcutbccsq, rcutbcc1, rcutbcc2;
  double rscalefcc = (3.0 + sqrt(8.0)) / 576.0;
  double rscalebcc1 = (1.0 + sqrt(2.0)) / (8.0 * sqrt(3));
  double rscalebcc2 = (1.0 + sqrt(2.0)) / 12.0;

  int cnafcc[MAXNEAR][MAXCOMMON];
  int cnabcc[MAXNEAR][MAXCOMMON];
  int bonds[MAXCOMMON], nbonds, maxbonds, minbonds;

  //int onenearest[MAXNEAR], onenearesttype[MAXNEAR];
  double mynearestx[14][3];

  bigint common[MAXCOMMON];
  int commontype[MAXCOMMON];
  double ix, iy, iz, jx, jy, jz, kx, ky, kz;
  double delx, dely, delz, rsq, r;
  int *jlist, *jindexlist;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary
  
  if (atom->nmax > namax) {
    namax = atom->nmax;
    grow_atom(namax);
    vector_atom = atompattern;
  }

  if (element->nmax > nemax) {
    nemax = element->nmax;
    grow_vatom(nemax, element->maxapc, element->maxucell);
    vector_vatom = vatompattern;
  }

  int nbytes = namax * sizeof(double);
  if (nbytes)
    memset(&atompattern[0], 0, nbytes);
  nbytes = nemax * element->maxapc * element->maxucell * sizeof(double); 
  if (nbytes)
    memset(&vatompattern[0][0][0], 0, nbytes);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list, preneighflag);

  double **ax = atom->x;
  double **ex = element->x;
  double ****nodex = element->nodex;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  int *npe = element->npe;
  int *apc = element->apc;
  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int *nucell = element->nucell;
  ElementVec *evec = element->evec;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  //int *element_firstindex = list->element_firstindex;

  int imask;
  double *patternptr;

  int nerror = 0;

  if (adapt_cut_flag) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      iindex = iindexlist[ii];
      // i is atom

      if (iindex < 0) {
        patternptr = &atompattern[i];
        imask = amask[i];
      } 

      // i is virtual atom

      else {
        ietype = etype[i];
        iapc = apc[ietype];
        iucell = iindex / iapc;
        ibasis = iindex % iapc;
        patternptr = &vatompattern[i][ibasis][iucell];
        imask = emask[i];
      }

      // skip if atom I has less than 16 neighbors
      // or not in group

      if (!(imask & groupbit)) continue;
      
      *patternptr = OTHER;
      if (numneigh[ii] < 12) continue;
      fccflag = 1;
      if (numneigh[ii] < 14) {
        num = 12;
        bccflag = 0;
      } else {
        num = 14;
        bccflag = 1;
      } 

      if (iindex < 0) {
        ix = ax[i][0];
        iy = ax[i][1];
        iz = ax[i][2];
      } else {
        ix = evec->interpolate(nodex, i, ibasis, iucell, 0);
        iy = evec->interpolate(nodex, i, ibasis, iucell, 1);
        iz = evec->interpolate(nodex, i, ibasis, iucell, 2);
      }

      // average nearest neighbor distance
      // to find rcut for CNA calculation (see OVITO doc)
      // NOTE: neighbor list is already sorted by distance

      rcutfccsq = 0.0;
      rcutbcc1 = 0.0;
      rcutbcc2 = 0.0;
      jlist = firstneigh[ii];
      jindexlist = firstneighindex[ii];
      for (jj = 0; jj < num; jj++) {
        j = jlist[jj];
        jindex = jindexlist[jj];
        if (jindex < 0) {
          mynearestx[jj][0] = ax[j][0];
          mynearestx[jj][1] = ax[j][1];
          mynearestx[jj][2] = ax[j][2];
        } else {
          jetype = etype[j];
          japc = apc[jetype];
          jucell = jindex / japc;
          jbasis = jindex % japc;
          mynearestx[jj][0] = evec->interpolate(nodex, j, jbasis, jucell, 0);
          mynearestx[jj][1] = evec->interpolate(nodex, j, jbasis, jucell, 1);
          mynearestx[jj][2] = evec->interpolate(nodex, j, jbasis, jucell, 2);
        }
        delx = mynearestx[jj][0] - ix;
        dely = mynearestx[jj][1] - iy;
        delz = mynearestx[jj][2] - iz;
        r = sqrt(delx * delx + dely * dely + delz * delz);
        if (jj < 12) rcutfccsq += r;
        if (bccflag) {
          if (jj < 8) rcutbcc1 += r;
          else rcutbcc2 += r;
        }
        ncommonfcc[jj] = 0;
        ncommonbcc[jj] = 0;
      }
      rcutfccsq = rcutfccsq * rcutfccsq * rscalefcc;
      if (bccflag) {
        rcutbcc1 *= rscalebcc1;
        rcutbcc2 *= rscalebcc2;
        rcutbccsq = (rcutbcc1 + rcutbcc2) * (rcutbcc1 + rcutbcc2);
      }

      // loop over near neighbors of I to build cna data structure
      // build both data structure for FCC and BCC at the same time
      // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
      // cna[k][NBONDS] = # of bonds between those common neighbors
      // cna[k][MAXBOND] = max # of bonds of any common neighbor

      for (jj = 0; jj < num; jj++) {

        // common = list of neighbors common to atom I and atom J
        // use near neighbor list of I to find common
        // also update common list for K 
        // since if K is common neighbor of I and J then
        // J is common neighbor if I and K

        for (kk = jj + 1; kk < num; kk++) {
          delx = mynearestx[jj][0] - mynearestx[kk][0];
          dely = mynearestx[jj][1] - mynearestx[kk][1];
          delz = mynearestx[jj][2] - mynearestx[kk][2];
          rsq = delx * delx + dely * dely + delz * delz;

          // ensure that this is done for only the 12 nearest neighbors
          if (fccflag && jj < 12 && kk < 12) 
            if (rsq < rcutfccsq) {
              if (ncommonfcc[jj] == 5 || 
                  ncommonfcc[kk] == 5) {
                fccflag = 0;
              } else {
                commonfcc[jj][ncommonfcc[jj]++] = kk;
                commonfcc[kk][ncommonfcc[kk]++] = jj;
              }
            }

          if (bccflag) 
            if (rsq < rcutbccsq) {
              if (ncommonbcc[jj] == 6 || 
                  ncommonbcc[kk] == 6) {
                bccflag = 0;
              } else { 
                commonbcc[jj][ncommonbcc[jj]++] = kk;
                commonbcc[kk][ncommonbcc[kk]++] = jj;
              }
            }

          // if fail to detect fcc or bcc, skip to next I atom;

          if (!fccflag && !bccflag) goto NEXT;
        }

        if (fccflag && ncommonfcc[jj] != 4 && ncommonfcc[jj] != 5) 
          fccflag = 0; 
        if (bccflag && ncommonbcc[jj] != 4 && ncommonbcc[jj] != 6) 
          bccflag = 0;

        // if fail to detect fcc or bcc, skip to next I atom;

        if (!fccflag && !bccflag) goto NEXT;

        // calculate total # of bonds between common neighbor atoms
        // also max and min # of common atoms any common atom is bonded to
        // bond = pair of atoms within cutoff

        if (fccflag) {
          for (n = 0; n < MAXCOMMON; n++) bonds[n] = 0;
          nbonds = 0;
          for (kk = 0; kk < ncommonfcc[jj]-1; kk++) {
            k = commonfcc[jj][kk];
            for (ll = kk + 1; ll < ncommonfcc[jj]; ll++) {
              l = commonfcc[jj][ll];
              delx = mynearestx[k][0] - mynearestx[l][0];
              dely = mynearestx[k][1] - mynearestx[l][1];
              delz = mynearestx[k][2] - mynearestx[l][2];
              rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < rcutfccsq) {
                nbonds++;
                bonds[kk]++;
                bonds[ll]++;
              }
            }
          }
          if (nbonds == 2 || nbonds == 5) {
            maxbonds = 0;
            minbonds = MAXCOMMON;
            for (n = 0; n < ncommonfcc[jj]; n++) {
              maxbonds = MAX(bonds[n], maxbonds);
              minbonds = MIN(bonds[n], minbonds);
            }
            cnafcc[jj][NCOMMON] = ncommonfcc[jj];
            cnafcc[jj][NBOND] = nbonds;
            cnafcc[jj][MAXBOND] = maxbonds;
            cnafcc[jj][MINBOND] = minbonds;
          } else fccflag = 0;
        }

        if (bccflag) {
          for (n = 0; n < MAXCOMMON; n++) bonds[n] = 0;
          nbonds = 0;
          for (kk = 0; kk < ncommonbcc[jj]-1; kk++) {
            k = commonbcc[jj][kk];
            for (ll = kk + 1; ll < ncommonbcc[jj]; ll++) {
              l = commonbcc[jj][ll];
              delx = mynearestx[k][0] - mynearestx[l][0];
              dely = mynearestx[k][1] - mynearestx[l][1];
              delz = mynearestx[k][2] - mynearestx[l][2];
              rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < rcutbccsq) {
                nbonds++;
                bonds[kk]++;
                bonds[ll]++;
              }
            }
          }
          if (nbonds == 4 || nbonds == 6) {
            maxbonds = 0;
            minbonds = MAXCOMMON;
            for (n = 0; n < ncommonbcc[jj]; n++) {
              maxbonds = MAX(bonds[n], maxbonds);
              minbonds = MIN(bonds[n], minbonds);
            }
            cnabcc[jj][NCOMMON] = ncommonbcc[jj];
            cnabcc[jj][NBOND] = nbonds;
            cnabcc[jj][MAXBOND] = maxbonds;
            cnabcc[jj][MINBOND] = minbonds;
          } else bccflag = 0;
        }

        // if fail to detect fcc or bcc, skip to next I atom;
        if (!fccflag && !bccflag) goto NEXT;
      }
      // detect CNA pattern of the atom


      if (fccflag) {
        nfcc = nhcp = nico = 0;
        for (inear = 0; inear < 12; inear++) {
          cj = cnafcc[inear][NCOMMON];
          ck = cnafcc[inear][NBOND];
          cl = cnafcc[inear][MAXBOND];
          cm = cnafcc[inear][MINBOND];
          if (cj == 4 && ck == 2 && cl == 1 && cm == 1)
            nfcc++;
          else if (cj == 4 && ck == 2 && cl == 2 && cm == 0)
            nhcp++;
          else if (cj == 5 && ck == 5 && cl == 2 && cm == 2)
            nico++;
        }
        if (nfcc == 12)
          *patternptr = FCC;
        else if (nfcc == 6 && nhcp == 6)
          *patternptr = HCP;
        else if (nico == 12)
          *patternptr = ICOS;
        else fccflag = 0;
      }

      if (!fccflag && bccflag) {
        nbcc4 = nbcc6 = 0;
        for (inear = 0; inear < 14; inear++) {
          cj = cnabcc[inear][NCOMMON];
          ck = cnabcc[inear][NBOND];
          cl = cnabcc[inear][MAXBOND];
          cm = cnabcc[inear][MINBOND];
          if (cj == 4 && ck == 4 && cl == 2 && cm == 2)
            nbcc4++;
          else if (cj == 6 && ck == 6 && cl == 2 && cm == 2)
            nbcc6++;
        }
        if (nbcc4 == 6 && nbcc6 == 8) *patternptr = BCC;
      }

NEXT: continue;
    } 
  } else {

  }
}

/*  ----------------------------------------------------------------------
   memory usage of local atom-based array
   -------------------------------------------------------------------------  */

double ComputeCNAAtom::memory_usage()
{
  double bytes = nalocalmax * sizeof(int);
  bytes += 2 * nalocalmax * MAXNEAR * sizeof(int);
  bytes += namax * sizeof(double);
  bytes += 2 * nelocalmax * element->maxucell * sizeof(int);
  bytes += nelocalmax * element->maxucell * MAXNEAR * sizeof(int);
  bytes += nemax * element->maxucell * sizeof(double);

  return bytes;
}


/* ----------------------------------------------------------------------------------------------- */

int ComputeCNAAtom::pack_atom_forward_comm(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = atompattern[j];
  }
  return m;
}
/* ------------------------------------------------------------------------------------------ */

void ComputeCNAAtom::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) atompattern[i] = buf[m++];
}


/* --------------------------------------------------------------------------------------------- */

int ComputeCNAAtom::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, k, l, m;
  int *nucell = element->nucell;
  int *apc = element->apc;
  int *etype = element->etype;
  int ietype;
  m = 0;
  for (j = 0; j < n; j++) {
    i = list[j];
    ietype = etype[i];
    for (l = 0; l < apc[ietype]; l++)
      for (k = 0; k < nucell[ietype]; k++) 
        buf[m++] = vatompattern[i][k][l];

  }
  return m;
}

/* ---------------------------------------------------------------------------------------------- */

void ComputeCNAAtom::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i, j, k, m, last;
  int *nucell = element->nucell;
  int *apc = element->apc;
  int *etype = element->etype;
  int ietype;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    ietype = etype[i];
    for (j = 0; j < apc[ietype]; j++)
      for (k = 0; k < nucell[ietype]; k++) 
        vatompattern[i][j][k] = buf[m++];

  }
}

/*  ----------------------------------------------------------------------
   grow atom-based arrays
   -------------------------------------------------------------------------  */

void ComputeCNAAtom::grow_atom(int nmax)  
{
  memory->destroy(atompattern);
  memory->create(atompattern, nmax, "cna:atompattern");
}

/*  ----------------------------------------------------------------------
   grow virtual atom-based arrays
   -------------------------------------------------------------------------  */

void ComputeCNAAtom::grow_vatom(int nmax, int maxapc, int nucell) 
{
  memory->destroy(vatompattern);
  memory->create(vatompattern, nmax, maxapc, nucell, "cna:vatompattern");
}
