#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "compute_ids_atom.h"
#include "atom.h"
#include "element.h"
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

enum{OTHER,CUBIC,HEX,UNKNOWN};
enum{MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeIDSAtom::ComputeIDSAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  list(NULL), apattern(NULL), iapattern(NULL), 
  ibucket(NULL), indexbucket(NULL), rbucket(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ids/atom command");

  peratom_flag = 1;
  ghostskinflag = 1;
  size_peratom_cols = 0;

  // rscale = ((1+sqrt(2))/2 * (1/12))^2
  // see OVITO
  
  rscale = (3.0 + sqrt(8))/576.0; 

  namax = nemax = 0;

  // allocate buckets for bubble sort
  
  maxbucket = 50;
  memory->create(ibucket,maxbucket,"npair:ibucket");
  memory->create(indexbucket,maxbucket,"npair:indexbucket");
  memory->create(rbucket,maxbucket,"npair:rbucket");
}

/* ---------------------------------------------------------------------- */

ComputeIDSAtom::~ComputeIDSAtom()
{
  memory->destroy(ibucket);
  memory->destroy(indexbucket);
  memory->destroy(rbucket);

  memory->destroy(apattern);
  memory->destroy(iapattern);
}

/* ---------------------------------------------------------------------- */

void ComputeIDSAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute ids/atom requires a pair style be defined");

  int cutoff = force->pair->cutforce + neighbor->skin;
  cutsq = cutoff*cutoff;
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ids/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ids/atom defined");

  // need an occasional full intpl neighbor list with sorting upto 4 nearest neighbors

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->intpllist = 1;
  neighbor->requests[irequest]->intglist = 0;
  neighbor->requests[irequest]->sort = 1;
  neighbor->requests[irequest]->sortnum = 4;
 
}

/* ---------------------------------------------------------------------- */

void ComputeIDSAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeIDSAtom::compute_peratom()
{

  int nfcc,nhcp;
  int cl,cm;
  int i,j,k,l;
  int ii,jj,kk,ll,n,m,p,inear;
  int jnum,knum;
  int ietype,jetype,ketype,letype;
  int iintpl,jintpl,kintpl,lintpl;
  int inode,jnode,knode,lnode;
  bigint iintpl_local,jintpl_local,kintpl_local;
  double rcutsq; 

  // cna[k][MAXBOND] = max # of bonds of any common neighbor
  // cna[k][MINBOND] = min # of bonds of any common neighbor

  int cna[12][2],bonds[4],nbonds,maxbonds,minbonds;

  int mynearest[12],mynearestintpl[12];

  int ncommon[12];
  int common[12][4];
  int commonintpl[12][4];
  double ix,iy,iz,jx,jy,jz,kx,ky,kz;
  double delx,dely,delz,rsq;
  int *jlist,*jintpllist,*klist,*kintpllist;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary (max size should be nlocal)

  if (atom->nmax > namax) {
    namax = atom->nmax;
    grow_atom(namax);
    vector_atom = apattern;
  }

  if (element->nmax > nemax) {
    nemax = element->nmax;
    grow_intpl(nemax,element->maxintpl);
    vector_intpl_atom = iapattern;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list,preneighflag);

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  double ***shape_array = element->shape_array;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  int *npe = element->npe;
  int nalocal = atom->nlocal;
  int naghost = atom->nghost;

  int nelocal = element->nlocal;
  int neghost = element->nghost;
  int *nintpl = element->nintpl;

  int ainum = list->ainum;
  int einum = list->einum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;
  bigint *e2ialist = list->e2ialist;
  int *numneigha = list->numneigha;
  int *numneighia = list->numneighia;
  int **firstneigha = list->firstneigha;
  int **firstneigha_index = list->firstneigha_index;
  int **firstneighia = list->firstneighia;
  int **firstneighia_index = list->firstneighia_index;

  // compute cna for each atom/interpolated atom in group
  // only performed if # of nearest neighbors = 12

  /*-------------------------------------------
   * Calculate IDS structure for atoms
   * ----------------------------------------*/

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];

    // skip if atom I has less than 16 neighbors
    // or not in group

    if (!(amask[i] & groupbit)) {
      apattern[i] = UNKNOWN;
      continue;
    } 
    
    apattern[i] = OTHER;

    if (numneigha[i] < 16) continue;

    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];
    
    // loop over 4 nearest neighbors atom J of atom I
    // take 4 nearest neighbor of each atom J (exclude I) as 2nd nearest neighbors of atom I
    
    n = 0;
    jlist = firstneigha[i];
    jintpllist = firstneigha_index[i];
    for (jj = 0; jj < 4; jj++) {
      j = jlist[jj];
      jintpl = jintpllist[jj];      

      // if J is owned

      if ((j < nalocal && jintpl == -1) || (j < nelocal && jintpl >= 0)) {

        // loop over 4 nearest neighbors of J

        if (jintpl == -1) {
          knum = numneigha[j];
          klist = firstneigha[j];
          kintpllist = firstneigha_index[j];
        } else {
          jintpl_local = e2ialist[j] + jintpl;
          knum = numneighia[jintpl_local];
          klist = firstneighia[jintpl_local];
          kintpllist = firstneighia_index[jintpl_local];
        }

        knum = MIN(knum,4);
        for (kk = 0; kk < 4; kk++) {
          k = klist[kk];
          kintpl = kintpllist[kk];

          // skip if K is I

          if (k == i && kintpl == -1) continue;
          if (n == 12) goto NEXTATOM;
          mynearest[n] = k;
          mynearestintpl[n++] = kintpl;
        }

      // if J is ghost, use full neighbor list of I (include I) to find 4 nearest neighbor of J (exclude J)
      // do include I since I might not be one of the 4 nearest neighbors of J after sorting

      } else {
        if (jintpl == -1) {
          jx = ax[j][0];
          jy = ax[j][1];
          jz = ax[j][2];
        } else {
          jx = jy = jz = 0.0;
          jetype = etype[j];
          for (jnode = 0; jnode < npe[jetype]; jnode++) {
            jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
            jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
            jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
          }
        }

        int nbuck = 1;   
        ibucket[0] = i;
        indexbucket[0] = -1;
        delx = ix - jx;
        dely = iy - jy;
        delz = iz - jz;
        rbucket[0] = delx*delx + dely*dely + delz*delz;

        knum = numneigha[i];
        klist = jlist;
        kintpllist = jintpllist;
        for (kk = 0; kk < knum; kk++) {
          k = klist[kk];
          kintpl = kintpllist[kk];
          delx = jx; dely = jy; delz = jz;

          // skip if K == J

          if (k == j && kintpl == jintpl) continue;
          if (kintpl == -1) {
            delx -= ax[k][0];
            dely -= ax[k][1];
            delz -= ax[k][2];
          } else {
            ketype = etype[k];
            for (knode = 0; knode < npe[ketype]; knode++) {
              delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
              dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
              delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
            }
          }
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            if (nbuck == maxbucket) grow_bucket();
            ibucket[nbuck] = k;
            indexbucket[nbuck] = kintpl;
            rbucket[nbuck] = rsq;
            nbuck++;
          }
        }

        // if number of neighbors < 4, just add them to list
        // start from 1 to skip I

        if (nbuck < 4) {
          if (n + nbuck - 1 > 12) goto NEXTATOM;
          for (l = 1; l < nbuck; l++) {
            mynearest[n] = ibucket[l];
            mynearestintpl[n++] = indexbucket[l];
          }

          // if more than 4 neighbors, perform bubble sort to find 4 nearest neighbors of J

        } else {
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
            kintpl = indexbucket[permute[l]];

            // skip if K == I

            if (k == i && kintpl == -1) continue;
            if (n == 12) goto NEXTATOM;
            mynearest[n] = k;
            mynearestintpl[n++] = kintpl;
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
      jintpl = mynearestintpl[jj];
      delx = ix; dely = iy; delz = iz;
      if (jintpl == -1) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
      } else {
        jetype = etype[j]; 
        for (jnode = 0; jnode < npe[jetype]; jnode++) {
          delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
          dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
          delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
        }
      }
      rcutsq += sqrt(delx*delx + dely*dely + delz*delz);
      ncommon[jj] = 0;
    }
    rcutsq = rcutsq*rcutsq*rscale;

    // loop over second nearest neighbors of I to build cna data structure

    for (jj = 0; jj < 12; jj++) {
      j = mynearest[jj];
      jintpl = mynearestintpl[jj];

      // common = list of neighbors common to atom I and atom J
      // loop through 2nd nearest neighbor of atom I and check if
      // distance to atom J < local cutoff (--> common)
      // if ncommon != 4 goto next I atom

      if (jj < 11) {
        if (jintpl == -1) {
          jx = ax[j][0];
          jy = ax[j][1];
          jz = ax[j][2];
        } else {
          jetype = etype[j];
          jx = jy = jz = 0.0;
          for (jnode = 0; jnode < npe[jetype]; jnode++) {
            jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0]; 
            jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1]; 
            jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2]; 
          }
        }

        for (kk = jj+1; kk < 12; kk++) {
          k = mynearest[kk];
          kintpl = mynearestintpl[kk];
          delx = jx; dely = jy; delz = jz;
          if (kintpl == -1) {
            delx -= ax[k][0];
            dely -= ax[k][1];
            delz -= ax[k][2];
          } else {
            ketype = etype[k];
            for (knode = 0; knode < npe[ketype]; knode++) {
              delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
              dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
              delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
            }
          }
          rsq = delx*delx + dely*dely + delz*delz;

          // if within local cutoff add to common list

          if (rsq < rcutsq) {
            if (ncommon[jj] == 4 || ncommon[kk] == 4) 
              goto NEXTATOM;
            common[jj][ncommon[jj]] = k;
            commonintpl[jj][ncommon[jj]++] = kintpl;
            common[kk][ncommon[kk]] = j;
            commonintpl[kk][ncommon[kk]++] = jintpl;
          }
        }
      }
      if (ncommon[jj] != 4) goto NEXTATOM;

      bonds[0] = bonds[1] = bonds[2] = bonds[3] = 0;

      nbonds = 0;
      for (kk = 0; kk < 3; kk++) {
        k = common[jj][kk];
        kintpl = commonintpl[jj][kk];

        if (kintpl == -1) {
          kx = ax[k][0];
          ky = ax[k][1];
          kz = ax[k][2];
        } else {
          ketype = etype[k];
          kx = ky = kz = 0.0;
          for (knode = 0; knode < npe[ketype]; knode++) {
            kx += shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
            ky += shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
            kz += shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
          }
        }

        for (ll = kk + 1; ll < 4; ll++) {
          l = common[jj][ll];
          lintpl = commonintpl[jj][ll];
          delx = kx; dely = ky; delz = kz;

          if (lintpl == -1) {
            delx -= ax[l][0];
            dely -= ax[l][1];
            delz -= ax[l][2];
          } else {
            letype = etype[l];
            for (lnode = 0; lnode < npe[letype]; lnode++) {
              delx -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][0];
              dely -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][1];
              delz -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][2];
            }
          }
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rcutsq) {
            if (nbonds == 2) goto NEXTATOM;
            nbonds++;
            bonds[kk]++;
            bonds[ll]++;
          }
        }
      }

      if (nbonds != 2) goto NEXTATOM;

      maxbonds = 0;
      maxbonds = MAX(bonds[0],maxbonds);
      maxbonds = MAX(bonds[1],maxbonds);
      maxbonds = MAX(bonds[2],maxbonds);
      cna[jj][MAXBOND] = MAX(bonds[3],maxbonds);
      minbonds = 4;
      minbonds = MIN(bonds[0],minbonds);
      minbonds = MIN(bonds[1],minbonds);
      minbonds = MIN(bonds[2],minbonds);
      cna[jj][MINBOND] = MIN(bonds[3],minbonds);
    }

    // detect IDS pattern of atom I

    nfcc = nhcp = 0;
    for (inear = 0; inear < 12; inear++) {
      cl = cna[inear][MAXBOND];
      cm = cna[inear][MINBOND];
      if (cl == 1 && cm == 1) nfcc++;
      else if (cl == 2 && cm == 0) nhcp++;
    }
    if (nfcc == 12) apattern[i] = CUBIC;
    else if (nfcc == 6 && nhcp == 6) apattern[i] = HEX;
NEXTATOM: continue;
  }

  /*-------------------------------------------
   * Calculate IDS structure for interpolated atoms
   * ----------------------------------------*/

  // loop over elements in compute group

  for (ii = 0; ii < einum; ii++) {

    i = eilist[ii];
    ietype = etype[i];

    if (!(emask[i] & groupbit)) {
      for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) 
        iapattern[i][iintpl] = UNKNOWN;
      continue;
    }

    for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) { 
      iintpl_local = e2ialist[i] + iintpl;

      iapattern[i][iintpl] = OTHER;
      // skip if atom I has less than 16 neighbors
      // or not in group

      if (numneighia[iintpl_local] < 16) continue;

      ix = iy = iz = 0.0;
      for (inode = 0; inode < npe[ietype]; inode++) {
        ix += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        iy += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        iz += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
      }

      // loop over 4 nearest neighbors J of I
      // take 4 nearest neighbor of each J (exclude I) as 2nd nearest neighbors of I

      n = 0;
      jlist = firstneighia[iintpl_local];
      jintpllist = firstneighia_index[iintpl_local];
      for (jj = 0; jj < 4; jj++) {
        j = jlist[jj];
        jintpl = jintpllist[jj];      

        // if J is owned

        if ((j < nalocal && jintpl == -1) || (j < nelocal && jintpl >= 0)) {

          // loop over 4 nearest neighbors of J

          if (jintpl == -1) {
            knum = numneigha[j];
            klist = firstneigha[j];
            kintpllist = firstneigha_index[j];
          } else {
            jintpl_local = e2ialist[j] + jintpl;
            knum = numneighia[jintpl_local];
            klist = firstneighia[jintpl_local];
            kintpllist = firstneighia_index[jintpl_local];
          }

          knum = MIN(knum,4);
          for (kk = 0; kk < 4; kk++) {
            k = klist[kk];
            kintpl = kintpllist[kk];

            // skip if K is I

            if (k == i && kintpl == iintpl) continue;
            if (n == 12) goto NEXTINTPL;
            mynearest[n] = k;
            mynearestintpl[n++] = kintpl;
          }

          // if J is ghost, use full neighbor list of I (include I) to find 4 nearest neighbor of J (exclude J)
          // do include I since I might not be one of the 4 nearest neighbors of J after sorting

        } else {
          if (jintpl == -1) {
            jx = ax[j][0];
            jy = ax[j][1];
            jz = ax[j][2];
          } else {
            jx = jy = jz = 0.0;
            jetype = etype[j];
            for (jnode = 0; jnode < npe[jetype]; jnode++) {
              jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
              jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
              jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
            }
          }

          int nbuck = 1;   
          ibucket[0] = i;
          indexbucket[0] = iintpl;
          delx = ix - jx;
          dely = iy - jy;
          delz = iz - jz;
          rbucket[0] = delx*delx + dely*dely + delz*delz;

          knum = numneighia[iintpl_local];
          klist = jlist;
          kintpllist = jintpllist;
          for (kk = 0; kk < knum; kk++) {
            k = klist[kk];
            kintpl = kintpllist[kk];
            delx = jx; dely = jy; delz = jz;

            // skip if K == J

            if (k == j && kintpl == jintpl) continue;
            if (kintpl == -1) {
              delx -= ax[k][0];
              dely -= ax[k][1];
              delz -= ax[k][2];
            } else {
              ketype = etype[k];
              for (knode = 0; knode < npe[ketype]; knode++) {
                delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
              }
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutsq) {
              if (nbuck == maxbucket) grow_bucket();
              ibucket[nbuck] = k;
              indexbucket[nbuck] = kintpl;
              rbucket[nbuck] = rsq;
              nbuck++;
            }
          }

          // if number of neighbors < 4, just add them to list
          // start from 1 to skip I

          if (nbuck < 4) {
            if (n  + nbuck - 1 > 12) goto NEXTINTPL;
            for (l = 1; l < nbuck; l++) {
              mynearest[n] = ibucket[l];
              mynearestintpl[n++] = indexbucket[l];
            }

            // if more than 4 neighbors, perform bubble sort to find 4 nearest neighbors of J

          } else {
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
              kintpl = indexbucket[permute[l]];

              // skip if K == I

              if (k == i && kintpl == iintpl) continue;
              if (n == 12) goto NEXTINTPL;
              mynearest[n] = k;
              mynearestintpl[n++] = kintpl;
            }
          }
        }
      }

      // if number of 2nd nearest != 12, continue to next intpl 

      if (n != 12) continue;

      // average second nearest neighbor distance
      // to find rcut for CNA calculation (see OVITO doc)

      rcutsq = 0;
      for (jj = 0; jj < 12; jj++) {
        j = mynearest[jj];
        jintpl = mynearestintpl[jj];
        delx = ix; dely = iy; delz = iz;
        if (jintpl == -1) {
          delx -= ax[j][0];
          dely -= ax[j][1];
          delz -= ax[j][2];
        } else {
          jetype = etype[j]; 
          for (jnode = 0; jnode < npe[jetype]; jnode++) {
            delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
            dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
            delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
          }
        }
        rcutsq += sqrt(delx*delx + dely*dely + delz*delz);
        ncommon[jj] = 0;
      }
      rcutsq = rcutsq*rcutsq*rscale;

      // loop over second nearest neighbors of I to build cna data structure

      for (jj = 0; jj < 12; jj++) {
        j = mynearest[jj];
        jintpl = mynearestintpl[jj];

        // common = list of neighbors common to I and J
        // loop through 2nd nearest neighbor of I and check if
        // distance to J < local cutoff (--> common)
        // if ncommon != 4 goto next intpl

        if (jj < 11) {
          if (jintpl == -1) {
            jx = ax[j][0];
            jy = ax[j][1];
            jz = ax[j][2];
          } else {
            jetype = etype[j];
            jx = jy = jz = 0.0;
            for (jnode = 0; jnode < npe[jetype]; jnode++) {
              jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0]; 
              jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1]; 
              jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2]; 
            }
          }

          for (kk = jj+1; kk < 12; kk++) {
            k = mynearest[kk];
            kintpl = mynearestintpl[kk];
            delx = jx; dely = jy; delz = jz;
            if (kintpl == -1) {
              delx -= ax[k][0];
              dely -= ax[k][1];
              delz -= ax[k][2];
            } else {
              ketype = etype[k];
              for (knode = 0; knode < npe[ketype]; knode++) {
                delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
              }
            }
            rsq = delx*delx + dely*dely + delz*delz;

            // if within local cutoff add to common list

            if (rsq < rcutsq) {
              if (ncommon[jj] == 4 || ncommon[kk] == 4) 
                goto NEXTINTPL;
              common[jj][ncommon[jj]] = k;
              commonintpl[jj][ncommon[jj]++] = kintpl;
              common[kk][ncommon[kk]] = j;
              commonintpl[kk][ncommon[kk]++] = jintpl;
            }
          }
        }

        if (ncommon[jj] != 4) goto NEXTINTPL;

        bonds[0] = bonds[1] = bonds[2] = bonds[3] = 0;

        nbonds = 0;
        for (kk = 0; kk < 3; kk++) {
          k = common[jj][kk];
          kintpl = commonintpl[jj][kk];

          if (kintpl == -1) {
            kx = ax[k][0];
            ky = ax[k][1];
            kz = ax[k][2];
          } else {
            ketype = etype[k];
            kx = ky = kz = 0.0;
            for (knode = 0; knode < npe[ketype]; knode++) {
              kx += shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
              ky += shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
              kz += shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
            }
          }

          for (ll = kk + 1; ll < 4; ll++) {
            l = common[jj][ll];
            lintpl = commonintpl[jj][ll];
            delx = kx; dely = ky; delz = kz;

            if (lintpl == -1) {
              delx -= ax[l][0];
              dely -= ax[l][1];
              delz -= ax[l][2];
            } else {
              letype = etype[l];
              for (lnode = 0; lnode < npe[letype]; lnode++) {
                delx -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][0];
                dely -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][1];
                delz -= shape_array[letype][lintpl][lnode]*nodex[l][lnode][2];
              }
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < rcutsq) {
              if (nbonds == 2) goto NEXTINTPL;
              nbonds++;
              bonds[kk]++;
              bonds[ll]++;
            }
          }
        }

        if (nbonds != 2) goto NEXTINTPL;

        maxbonds = 0;
        maxbonds = MAX(bonds[0],maxbonds);
        maxbonds = MAX(bonds[1],maxbonds);
        maxbonds = MAX(bonds[2],maxbonds);
        cna[jj][MAXBOND] = MAX(bonds[3],maxbonds);
        minbonds = 4;
        minbonds = MIN(bonds[0],minbonds);
        minbonds = MIN(bonds[1],minbonds);
        minbonds = MIN(bonds[2],minbonds);
        cna[jj][MINBOND] = MIN(bonds[3],minbonds);

      }

      // detect IDS pattern of the interpolated atom I

      nfcc = nhcp = 0;
      for (inear = 0; inear < 12; inear++) {
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
        if (cl == 1 && cm == 1) nfcc++;
        else if (cl == 2 && cm == 0) nhcp++;
      }
      if (nfcc == 12) iapattern[i][iintpl] = CUBIC;
      else if (nfcc == 6 && nhcp == 6) iapattern[i][iintpl] = HEX;
NEXTINTPL: continue;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeIDSAtom::memory_usage()
{
  double bytes = namax * sizeof(double);
  bytes += nemax * element->maxintpl * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   grow atom-based arrays
   ------------------------------------------------------------------------- */

void ComputeIDSAtom::grow_bucket()  
{
  maxbucket += maxbucket/2;
  memory->grow(ibucket,maxbucket,"ids:ibucket");
  memory->grow(indexbucket,maxbucket,"ids:indexbucket");
  memory->grow(rbucket,maxbucket,"ids:rbucket");
}
/* ----------------------------------------------------------------------
   grow atom-based arrays
   ------------------------------------------------------------------------- */

void ComputeIDSAtom::grow_atom(int nmax)  
{
  memory->destroy(apattern);
  memory->create(apattern,nmax,"ids:apattern");
}

/* ----------------------------------------------------------------------
   grow interpolated atom-based arrays
   ------------------------------------------------------------------------- */

void ComputeIDSAtom::grow_intpl(int nmax, int nintpl) 
{
  memory->destroy(iapattern);
  memory->create(iapattern,nmax,nintpl,"ids:iapattern");
}
