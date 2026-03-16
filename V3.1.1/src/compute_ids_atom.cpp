#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "compute_ids_atom.h"
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

#define MAXNEAR 20

enum{UNKNOWN, OTHER, CUBIC, HEX};
enum{MAXBOND, MINBOND};

/*  ----------------------------------------------------------------------  */

ComputeIDSAtom::ComputeIDSAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg), 
  list(nullptr), atompattern(nullptr), vatompattern(nullptr), 
  ibucket(nullptr), indexbucket(nullptr), rbucket(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute ids/atom command");

  peratom_flag = 1;
  ghostskinflag = 1;
  size_peratom_cols = 0;

  // rscale = ((1+sqrt(2))/2 * (1/12))^2
  // see OVITO
  
  rscale = (3.0 + sqrt(8))/576.0; 

  namax = nemax = 0;

  // allocate buckets for bubble sort
  
  maxbucket = 50;
  memory->create(ibucket, maxbucket, "npair:ibucket");
  memory->create(indexbucket, maxbucket, "npair:indexbucket");
  memory->create(rbucket, maxbucket, "npair:rbucket");
}

/*  ----------------------------------------------------------------------  */

ComputeIDSAtom::~ComputeIDSAtom()
{
  memory->destroy(ibucket);
  memory->destroy(indexbucket);
  memory->destroy(rbucket);

  memory->destroy(atompattern);
  memory->destroy(vatompattern);
}

/*  ----------------------------------------------------------------------  */

void ComputeIDSAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute ids/atom requires a pair style be defined");

  int cutoff = force->pair->cutforce + neighbor->skin;
  cutsq = cutoff * cutoff;
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style, "ids/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute ids/atom defined");

  // need an occasional full ucell neighbor list with sorting upto 4 nearest neighbors

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->vatomlist = 1;
  neighbor->requests[irequest]->gausslist = 0;
  neighbor->requests[irequest]->neighneigh = 1;
  neighbor->requests[irequest]->sort = 1;
  neighbor->requests[irequest]->sortnum = 4;
  neighbor->requests[irequest]->group = 1;
  neighbor->requests[irequest]->groupbit = groupbit;
  //neighbor->requests[irequest]->noinner = 1;

}

/*  ----------------------------------------------------------------------  */

void ComputeIDSAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/*  ----------------------------------------------------------------------  */

void ComputeIDSAtom::compute_peratom()
{

  int nfcc, nhcp;
  int cl, cm;
  int i, j, k, l;
  int ii, jj, kk, ll, n, inear;
  int jnum, knum;
  int ietype, jetype, ketype;
  int iucell, jucell, kucell;
  int iindex, jindex, kindex;
  int ibasis, jbasis, kbasis;
  int inode;
  int iapc, japc, kapc;
  int jlistindex;
  double rcutsq; 

  // cna[k][MAXBOND] = max # of bonds of any common neighbor
  // cna[k][MINBOND] = min # of bonds of any common neighbor

  int cna[12][2], bonds[4], nbonds, maxbonds, minbonds;

  int mynearest[12], mynearestindex[12];
  double mynearestx[12][3];

  int ncommon[12];
  int common[12][4];
  double ix, iy, iz, jx, jy, jz;
  double delx, dely, delz, rsq;
  int *jlist, *jindexlist, *klist, *kindexlist;

  invoked_peratom = update->ntimestep;


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
  int naghost = atom->nghost;
  ElementVec *evec = element->evec;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int *nucell = element->nucell;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  int *element_firstindex = list->element_firstindex;

  int imask;
  double *patternptr;


  /* -------------------------------------------
  *Calculate IDS structure for atoms
  *---------------------------------------- */

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

    if (!(imask & groupbit)) {
      *patternptr = UNKNOWN;
      continue;
    } 
    *patternptr = OTHER;
    if (numneigh[ii] < 16) continue;

    if (iindex < 0) {
      ix = ax[i][0];
      iy = ax[i][1];
      iz = ax[i][2];
    } else {
      ix = evec->interpolate(nodex, i, ibasis, iucell, 0);
      iy = evec->interpolate(nodex, i, ibasis, iucell, 1);
      iz = evec->interpolate(nodex, i, ibasis, iucell, 2);
    }

    // loop over 4 nearest neighbors J of I
    // take 4 nearest neighbors of each J (exclude I) as 2nd nearest neighbors of atom I

    n = 0;
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];
    for (jj = 0; jj < 4; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];      

      // if J is owned

      //if ((j < nalocal && jindex == -1) || (j < nelocal && jindex >= 0)) {

      if (0) {
        // loop over 4 nearest neighbors of J

        if (jindex < 0) 
          jlistindex = j;
        else 
          jlistindex = element_firstindex[j] + jindex;

        knum = numneigh[jlistindex];
        klist = firstneigh[jlistindex];
        kindexlist = firstneighindex[jlistindex];

        if (knum < 4) goto NEXT;
        
        for (kk = 0; kk < 4; kk++) {
          k = klist[kk];
          kindex = kindexlist[kk];

          // skip if K is I

          if (k == i && kindex == iindex) continue;
          if (n == 12) goto NEXT;
          mynearest[n] = k;
          mynearestindex[n++] = kindex;
        }
      } 

      // if J is ghost, use full neighbor list of I (including I) to find 4 nearest neighbor of J (exclude J)
      // do include I since I might not be one of the 4 nearest neighbors of J after sorting

      else {
        if (jindex < 0) {
          jx = ax[j][0];
          jy = ax[j][1];
          jz = ax[j][2];
        } else {
          jetype = etype[j];
          japc = apc[jetype];
          jucell = jindex / japc;
          jbasis = jindex % japc;
          jx = evec->interpolate(nodex, j, jbasis, jucell, 0);
          jy = evec->interpolate(nodex, j, jbasis, jucell, 1);
          jz = evec->interpolate(nodex, j, jbasis, jucell, 2);
        }

        int nbuck = 1;   
        ibucket[0] = i;
        indexbucket[0] = iindex;
        delx = ix - jx;
        dely = iy - jy;
        delz = iz - jz;
        rbucket[0] = delx * delx + dely * dely + delz * delz;

        knum = numneigh[ii];
        klist = jlist;
        kindexlist = jindexlist;
        for (kk = 0; kk < knum; kk++) {
          k = klist[kk];
          kindex = kindexlist[kk];
          delx = jx; dely = jy; delz = jz;

          // skip if K == J

          if (k == j && kindex == jindex) continue;
          if (kindex < 0) {
            delx -= ax[k][0];
            dely -= ax[k][1];
            delz -= ax[k][2];
          } else {
            ketype = etype[k];
            kapc = apc[ketype];
            kucell = kindex / kapc;
            kbasis = kindex % kapc;
            delx -= evec->interpolate(nodex, k, kbasis, kucell, 0);
            dely -= evec->interpolate(nodex, k, kbasis, kucell, 1);
            delz -= evec->interpolate(nodex, k, kbasis, kucell, 2);
          }
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            if (nbuck == maxbucket) grow_bucket();
            ibucket[nbuck] = k;
            indexbucket[nbuck] = kindex;
            rbucket[nbuck] = rsq;
            nbuck++;
          }
        }

        // if number of neighbors < 4, just add them to list
        // start from 1 to skip I

        if (nbuck < 4) {
          if (n + nbuck - 1 > 12) goto NEXT;
          for (l = 1; l < nbuck; l++) {
            mynearest[n] = ibucket[l];
            mynearestindex[n++] = indexbucket[l];
          }
        } 

        // if more than 4 neighbors, perform bubble sort to find 4 nearest neighbors of J
          
        else {
          int permute[nbuck];
          for (l = 0; l < nbuck; l++)
            permute[l] = l;
          for (l = 0; l < 4; l++) {
            for (k = nbuck-1; k > l; k--) {
              int current = permute[k];
              int next = permute[k-1];
              if (rbucket[current] < rbucket[next]) {
                permute[k] = next;
                permute[k-1] = current;
              }
            }
            k = ibucket[permute[l]];
            kindex = indexbucket[permute[l]];

            // skip if K == I

            if (k == i && kucell == -1) continue;
            if (n == 12) goto NEXT;
            mynearest[n] = k;
            mynearestindex[n++] = kindex;
          }
        }
      }
    }

    // if number of 2nd nearest != 12, continue to next atom

    if (n != 12) continue;

    // average second nearest neighbor distance
    // to find rcut for CNA calculation (see OVITO doc)

    rcutsq = 0;
    for (jj = 0; jj < 12; jj++) {
      j = mynearest[jj];
      jindex = mynearestindex[jj];
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
      rcutsq += sqrt(delx * delx + dely * dely + delz * delz);
      ncommon[jj] = 0;
    }
    rcutsq = rcutsq * rcutsq * rscale;

    // loop over second nearest neighbors of I to build cna data structure

    for (jj = 0; jj < 12; jj++) {

      // common = list of neighbors common to atom I and atom J
      // loop through 2nd nearest neighbor of atom I and check if
      // distance to atom J < local cutoff (--> common)
      // if ncommon == 4 goto next I atom

      for (kk = jj+1; kk < 12; kk++) {
        delx = mynearestx[jj][0] - mynearestx[kk][0];
        dely = mynearestx[jj][1] - mynearestx[kk][1];
        delz = mynearestx[jj][2] - mynearestx[kk][2];
        rsq = delx * delx + dely * dely + delz * delz;

        // if within local cutoff add to common list

        if (rsq < rcutsq) {
          if (ncommon[jj] == 4 || ncommon[kk] == 4) 
            goto NEXT;
          common[jj][ncommon[jj]++] = kk;
          common[kk][ncommon[kk]++] = jj;
        }
      }

      if (ncommon[jj] != 4) goto NEXT;

      bonds[0] = bonds[1] = bonds[2] = bonds[3] = 0;

      nbonds = 0;
      for (kk = 0; kk < 3; kk++) {
        k = common[jj][kk];
        for (ll = kk + 1; ll < 4; ll++) {
          l = common[jj][ll];
          delx = mynearestx[k][0] - mynearestx[l][0];
          dely = mynearestx[k][1] - mynearestx[l][1];
          delz = mynearestx[k][2] - mynearestx[l][2];
          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < rcutsq) {
            if (nbonds == 2) goto NEXT;
            nbonds++;
            bonds[kk]++;
            bonds[ll]++;
          }
        }
      }

      if (nbonds != 2) goto NEXT;

      maxbonds = 0;
      maxbonds = MAX(bonds[0], maxbonds);
      maxbonds = MAX(bonds[1], maxbonds);
      maxbonds = MAX(bonds[2], maxbonds);
      cna[jj][MAXBOND] = MAX(bonds[3], maxbonds);
      minbonds = 4;
      minbonds = MIN(bonds[0], minbonds);
      minbonds = MIN(bonds[1], minbonds);
      minbonds = MIN(bonds[2], minbonds);
      cna[jj][MINBOND] = MIN(bonds[3], minbonds);
    }

    // detect IDS pattern of atom I

    nfcc = nhcp = 0;
    for (inear = 0; inear < 12; inear++) {
      cl = cna[inear][MAXBOND];
      cm = cna[inear][MINBOND];
      if (cl == 1 && cm == 1) nfcc++;
      else if (cl == 2 && cm == 0) nhcp++;
    }
    if (nfcc == 12) *patternptr = CUBIC;
    else if (nfcc == 6 && nhcp == 6) *patternptr = HEX;
NEXT: continue;
  }

  if (update->ntimestep == 100) error->all(FLERR,"TEST"); 
}

/*  ----------------------------------------------------------------------
   memory usage of local atom-based array
   -------------------------------------------------------------------------  */

double ComputeIDSAtom::memory_usage()
{
  double bytes = namax * sizeof(double);
  bytes += nemax * element->maxapc * element->maxucell * sizeof(double);

  return bytes;
}

/*  ----------------------------------------------------------------------
   grow atom-based arrays
   -------------------------------------------------------------------------  */

void ComputeIDSAtom::grow_bucket()  
{
  maxbucket += maxbucket/2;
  memory->grow(ibucket, maxbucket, "ids:ibucket");
  memory->grow(indexbucket, maxbucket, "ids:indexbucket");
  memory->grow(rbucket, maxbucket, "ids:rbucket");
}
/*  ----------------------------------------------------------------------
   grow atom-based arrays
   -------------------------------------------------------------------------  */

void ComputeIDSAtom::grow_atom(int nmax)  
{
  memory->destroy(atompattern);
  memory->create(atompattern, nmax, "ids:atompattern");
}

/*  ----------------------------------------------------------------------
   grow virtual atom-based arrays
   -------------------------------------------------------------------------  */

void ComputeIDSAtom::grow_vatom(int nmax, int maxapc, int nucell) 
{
  memory->destroy(vatompattern);
  memory->create(vatompattern, nmax, maxapc, nucell, "ids:vatompattern");
}
