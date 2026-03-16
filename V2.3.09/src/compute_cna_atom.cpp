#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "compute_cna_atom.h"
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

#define MAXNEAR 16
#define MAXCOMMON 8

enum{UNKNOWN,FCC,HCP,BCC,ICOS,OTHER};
enum{NCOMMON,NBOND,MAXBOND,MINBOND};
enum{ATOM,INTPL};

/* ---------------------------------------------------------------------- */

ComputeCNAAtom::ComputeCNAAtom(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  list(NULL), anearest(NULL), anearesttype(NULL), ianearest(NULL), ianearesttype(NULL), annearest(NULL),iannearest(NULL), apattern(NULL), iapattern(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal compute cna/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff = universe->numeric(FLERR,arg[3]);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute cna/atom command");
  cutsq = cutoff*cutoff;

  namax = nemax = 0;

  // set comm size needed by this compute

  comm_atom_forward = 1;
  comm_elem_forward = element->maxintpl;
}

/* ---------------------------------------------------------------------- */

ComputeCNAAtom::~ComputeCNAAtom()
{
  memory->destroy(anearest);
  memory->destroy(anearesttype);
  memory->destroy(annearest);
  memory->destroy(apattern);
  memory->destroy(ianearest);
  memory->destroy(ianearesttype);
  memory->destroy(iannearest);
  memory->destroy(iapattern);

}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute cna/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute cna/atom cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute cna/atom cutoff may be too large to find "
        "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cna/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cna/atom defined");

  // need an occasional full intpl neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->intpllist = 1;
  neighbor->requests[irequest]->intglist = 0;

  // use custom cutoff to save computational time
  
  neighbor->requests[irequest]->cut = 1;
  neighbor->requests[irequest]->cutoff = sqrt(cutsq) + neighbor->skin;
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCNAAtom::compute_peratom()
{
  int nfcc,nhcp,nbcc4,nbcc6,nico;
  int cj,ck,cl,cm;
  int pattern;
  bigint i,j,k,l;
  int ii,jj,kk,ll,n,inear,jnear;
  int jnum,knum,ei,ej,ek,el;
  int ietype,jetype,ketype,letype;
  int iintpl,jintpl,kintpl,lintpl;
  int inode,jnode,knode,lnode;
  bigint iintpl_local,jintpl_local,kintpl_local;
  int firstflag,ncommon,nnearest;

  int cna[MAXNEAR][4],bonds[MAXCOMMON],nbonds,maxbonds,minbonds;
  int onenearest[MAXNEAR],onenearesttype[MAXNEAR];

  bigint common[MAXCOMMON];
  int commontype[MAXCOMMON];
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
  int npe = element->npe;
  int nalocal = atom->nlocal;

  int nelocal = element->nlocal;
  int *nintpl = element->nintpl;

  int ainum = list->ainum;
  int einum = list->einum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;
  bigint *e2ialist = list->e2ialist;
  int *ia2elist = list->ia2elist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighia2a = list->numneighia2a;
  int *numneighe2e = list->numneighe2e;
  int *numneighia2ia = list->numneighia2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighia2a = list->firstneighia2a;
  int **firstneighe2e = list->firstneighe2e;
  int **firstneighia2ia = list->firstneighia2ia;
  int **firstneighia2ia_index = list->firstneighia2ia_index;

  int nerror = 0;

  // find the neigbours of each atoms/interpolated atoms within cutoff using full neighbor list
  // anearest[] = atom indices of nearest neighbors, up to MAXNEAR
  // ianearest[] = interpolated atom indices of nearest neighbors, up to MAXNEAR
  // do this for all atoms/interpolated atoms, not just compute group
  // since CNA calculation requires neighbors of neighbors

  /*---------------------------------------------
   * find atoms' nearest neighbors 
   * ------------------------------------------*/

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];
    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];

    n = 0;

    // loop over atom neighbors
    
    jnum = numneigha2a[i];
    if (jnum) {
      jlist = firstneigha2a[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        //j &= NEIGHMASK;
        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          if (n < MAXNEAR) {
            anearest[i][n] = j;
            anearesttype[i][n++] = ATOM;
          } else {
            nerror++;
            break;
          }
        }
      }
    }

    // loop over interpolated atom neighbors
    
    jnum = numneigha2ia[i];
    if (jnum) {
      jlist = firstneigha2ia[i];
      jintpllist = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jetype = etype[j];
        //j &= NEIGHMASK;
        jintpl = jintpllist[jj];
        jintpl_local = e2ialist[j] + jintpl; 
        delx = ix; dely = iy; delz = iz;
        for (jnode = 0; jnode < npe; jnode++) {
          delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
          dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
          delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
        }
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) 
          if (n < MAXNEAR) {
            anearest[i][n] = jintpl_local;
            anearesttype[i][n++] = INTPL;
          } else {
            nerror++;
            break;
          }
      }
    }
    annearest[i] = n;
  }

  /*---------------------------------------------
   * find interpolated atoms' nearest neighbors 
   * ------------------------------------------*/

  // loop over elements then loop over interpolated atoms inside each element

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    for (iintpl = 0; iintpl < nintpl[ietype]; iintpl++) {
      n = 0;
      ix = iy = iz = 0.0;
      for (inode = 0; inode < npe; inode++) {
        ix += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        iy += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        iz += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
      }

      iintpl_local = e2ialist[i] + iintpl;

      // loop over neigboring atoms

      jnum = numneighia2a[iintpl_local];
      if (jnum) {
        jlist = firstneighia2a[iintpl_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          //j &= NEIGHMASK;
          delx = ix - ax[j][0];
          dely = iy - ax[j][1];
          delz = iz - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            if (n < MAXNEAR) { 
              ianearest[i][iintpl][n] = j;
              ianearesttype[i][iintpl][n++] = ATOM;
            } else {
              nerror++;
              break;
            }
          }
        }
      }

      // loop over neighboring interpolated atoms
      
      jnum = numneighia2ia[iintpl_local];
      if (jnum) {
        jlist = firstneighia2ia[iintpl_local]; 
        jintpllist = firstneighia2ia_index[iintpl_local];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          jintpl = jintpllist[jj];
          jintpl_local = e2ialist[j] + jintpl;
          jetype = etype[j];
          delx = ix; dely = iy; delz = iz;
          for (jnode = 0; jnode < npe; jnode++) {
            delx -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
            dely -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
            delz -= shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
          }
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < cutsq) 
            if (n < MAXNEAR) {
              ianearest[i][iintpl][n] = jintpl_local;
              ianearesttype[i][iintpl][n++] = INTPL;
            } else {
              nerror++;
              break;
            }
        }
      }
      iannearest[i][iintpl] = n;
    }
  }

  // warning message

  int nerrorall;
  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0) {
    char str[128];
    sprintf(str,"Too many neighbors in CNA for %d atoms",nerrorall);
    error->warning(FLERR,str,0);
  }

  // compute CNA for each atom/interpolated atom in group
  // only performed if # of nearest neighbors = 12 or 14 (fcc,hcp)

  nerror = 0;

  /*-------------------------------------------
   * Calculate CNA structure for atoms
   * ----------------------------------------*/

  int flag = 0;

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];

    if (!(amask[i] & groupbit)) {
      apattern[i] = UNKNOWN;
      continue;
    }

    apattern[i] = OTHER;
    if (annearest[i] != 12 && annearest[i] != 14) continue;

    // loop over near neighbors of I to build cna data structure
    // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
    // cna[k][NBONDS] = # of bonds between those common neighbors
    // cna[k][MAXBOND] = max # of bonds of any common neighbor
    // cna[k][MINBOND] = min # of bonds of any common neighbor

    for (jj = 0; jj < annearest[i]; jj++) {
      j = anearest[i][jj];

      // common = list of neighbors common to atom I and atom J
      // J neighbor is atom

      if (anearesttype[i][jj] == ATOM) {
        if (j < nalocal) {

          // if J is an owned atom, use its near neighbor list to find them

          firstflag = 1;
          ncommon = 0;
          for (inear = 0; inear < annearest[i]; inear++)
            for (jnear = 0; jnear < annearest[j]; jnear++)

              // neighbor index and type must match to be common

              if (anearest[i][inear] == anearest[j][jnear] &&
                  anearesttype[i][inear] == anearesttype[j][jnear] ) 
                if (ncommon < MAXCOMMON) {
                  common[ncommon] = anearest[i][inear];
                  commontype[ncommon++] = anearesttype[i][inear];
                } else if (firstflag) {
                  nerror++;
                  firstflag = 0;
                }

        } else {

          // if J is a ghost atom, use full neighbor list of I to find them
          // must exclude J from I's neighbor list

          jx = ax[j][0];
          jy = ax[j][1];
          jz = ax[j][2];

          n = 0;

          // loop over atom neighbor K of I

          knum = numneigha2a[i];
          if (knum) {
            klist = firstneigha2a[i];
            for (kk = 0; kk < knum; kk++) {
              k = klist[kk];
              if (k == j) continue;
              delx = jx - ax[k][0];
              dely = jy - ax[k][1];
              delz = jz - ax[k][2];
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) 
                if (n < MAXNEAR) {
                  onenearest[n] = k;
                  onenearesttype[n++] = ATOM;
                } else break;
            }
          }

          // loop over interpolated atom neighbor K of I

          knum = numneigha2ia[i];
          if (knum) {
            klist = firstneigha2ia[i];
            kintpllist = firstneigha2ia_index[i];
            for (kk = 0; kk < knum; kk++) {
              k = klist[kk];
              ketype = etype[k];
              kintpl = kintpllist[kk];
              kintpl_local = e2ialist[k] + kintpl;
              //k &= NEIGHMASK;
              delx = jx; dely = jy; delz = jz;
              for (knode = 0; knode < npe; knode++) {
                delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) 
                if (n < MAXNEAR) {
                  onenearest[n] = kintpl_local;
                  onenearesttype[n++] = INTPL;
                } else break;
            }
          }
          firstflag = 1;
          ncommon = 0;
          for (inear = 0; inear < annearest[i]; inear++)
            for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
              if (anearest[i][inear] == onenearest[jnear] && 
                  anearesttype[i][inear] == onenearesttype[jnear]) 
                if (ncommon < MAXCOMMON) {
                  common[ncommon] = anearest[i][inear];
                  commontype[ncommon++] = anearesttype[i][inear];
                } else if (firstflag) {
                  nerror++;
                  firstflag = 0;
                }
        }

        cna[jj][NCOMMON] = ncommon;

        for (n = 0; n < ncommon; n++) bonds[n] = 0;

        nbonds = 0;
        for (kk = 0; kk < ncommon - 1; kk++) {
          k = common[kk];

          // get the coords of L based on its type (ATOM or INTPL)

          if (commontype[kk] == ATOM) {
            kx = ax[k][0];
            ky = ax[k][1];
            kz = ax[k][2];
          } else {
            ek = ia2elist[k];
            kintpl = k - e2ialist[ek];
            ketype = etype[ek];
            kx = ky = kz = 0.0;
            for (knode = 0; knode < npe; knode++) {
              kx += shape_array[ketype][kintpl][knode]*nodex[ek][knode][0];
              ky += shape_array[ketype][kintpl][knode]*nodex[ek][knode][1];
              kz += shape_array[ketype][kintpl][knode]*nodex[ek][knode][2];
            }
          }

          for (ll = kk + 1; ll < ncommon; ll++) {
            l = common[ll];

            // get the coords of L based on its type (ATOM or INTPL)


            delx = kx; dely = ky; delz = kz;
            if (commontype[ll] == ATOM) {
              delx -= ax[l][0];
              dely -= ax[l][1];
              delz -= ax[l][2];
            } else {
              el = ia2elist[l];
              lintpl = l - e2ialist[el];
              letype = etype[el];
              for (lnode = 0; lnode < npe; lnode++) {
                delx -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][0];
                dely -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][1];
                delz -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][2];
              }
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutsq) {
              nbonds++;
              bonds[kk]++;
              bonds[ll]++;
            }
          }
        }

        cna[jj][NBOND] = nbonds;

        maxbonds = 0;
        minbonds = MAXCOMMON;
        for (n = 0; n < ncommon; n++) {
          maxbonds = MAX(bonds[n],maxbonds);
          minbonds = MIN(bonds[n],minbonds);
        }
        cna[jj][MAXBOND] = maxbonds;
        cna[jj][MINBOND] = minbonds;

      } else {

        // J neighbor is interpolated atom

        // index of parent element of J

        ej = ia2elist[j];

        // index of J in its parent element

        jintpl = j - e2ialist[ej]; 

        if (ej < nelocal) {

          // if J is an owned interpolated atom, use its near neighbor list to find common

          firstflag = 1;
          ncommon = 0;
          for (inear = 0; inear < annearest[i]; inear++)
            for (jnear = 0; jnear < iannearest[ej][jintpl]; jnear++)

              // neighbor index and type must match to be common

              if (anearest[i][inear] == ianearest[ej][jintpl][jnear] &&
                  anearesttype[i][inear] == ianearesttype[ej][jintpl][jnear] ) {
                if (ncommon < MAXCOMMON) {
                  common[ncommon] = anearest[i][inear];
                  commontype[ncommon++] = anearesttype[i][inear];
                } else if (firstflag) {
                  nerror++;
                  firstflag = 0;
                }
              }
        } else {

          // if J is a ghost interpolated atom, use full neighbor list of I to find them

          jetype = etype[ej]; 
          jx = jy = jz = 0.0;
          for (jnode = 0; jnode < npe; jnode++) {
            jx += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][0];
            jy += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][1];
            jz += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][2];
          }

          n = 0;

          // loop over atom neighbor K of I

          knum = numneigha2a[i];
          if (knum) {
            klist = firstneigha2a[i];
            for (kk = 0; kk < knum; kk++) {
              k = klist[kk];
              //k &= NEIGHMASK;
              delx = jx - ax[k][0];
              dely = jy - ax[k][1];
              delz = jz - ax[k][2];
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) 
                if (n < MAXNEAR) {
                  onenearest[n] = k;
                  onenearesttype[n++] = ATOM;
                } else break;
            }
          }

          // loop over interpolated atom neighbor K of I
          // must exclude J from full neighbor list of I

          knum = numneigha2ia[i];
          if (knum) {
            klist = firstneigha2ia[i];
            kintpllist = firstneigha2ia_index[i];
            for (kk = 0; kk < knum; kk++) {
              k = klist[kk];
              ketype = etype[k];
              kintpl = kintpllist[kk];
              kintpl_local = e2ialist[k] + kintpl;
              if (k == ej && kintpl == jintpl) continue;
              //k &= NEIGHMASK;
              delx = jx; dely = jy; delz = jz;
              for (knode = 0; knode < npe; knode++) {
                delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) 
                if (n < MAXNEAR) {
                  onenearest[n] = kintpl_local;
                  onenearesttype[n++] = INTPL;
                } else break;
            }
          }
          firstflag = 1;
          ncommon = 0;
          for (inear = 0; inear < annearest[i]; inear++)
            for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
              if (anearest[i][inear] == onenearest[jnear] && 
                  anearesttype[i][inear] == onenearesttype[jnear]) {
                if (ncommon < MAXCOMMON) {
                  common[ncommon] = anearest[i][inear];
                  commontype[ncommon++] = anearesttype[i][inear];
                } else if (firstflag) {
                  nerror++;
                  firstflag = 0;
                }
              }
        }
        cna[jj][NCOMMON] = ncommon;

        for (n = 0; n < ncommon; n++) bonds[n] = 0;

        nbonds = 0;
        for (kk = 0; kk < ncommon - 1; kk++) {
          k = common[kk];

          // get the coords of L based on its type (ATOM or INTPL)

          if (commontype[kk] == ATOM) {
            kx = ax[k][0];
            ky = ax[k][1];
            kz = ax[k][2];
          } else {
            ek = ia2elist[k];
            kintpl = k - e2ialist[ek];
            ketype = etype[ek];
            kx = ky = kz = 0.0;
            for (knode = 0; knode < npe; knode++) {
              kx += shape_array[ketype][kintpl][knode]*nodex[ek][knode][0];
              ky += shape_array[ketype][kintpl][knode]*nodex[ek][knode][1];
              kz += shape_array[ketype][kintpl][knode]*nodex[ek][knode][2];
            }
          }

          for (ll = kk + 1; ll < ncommon; ll++) {
            l = common[ll];
            // get the coords of L based on its type (ATOM or INTPL)

            delx = kx; dely = ky; delz = kz; 
            if (commontype[ll] == ATOM) {
              delx -= ax[l][0];
              dely -= ax[l][1];
              delz -= ax[l][2];
            } else {
              el = ia2elist[l];
              lintpl = l - e2ialist[el];
              letype = etype[el];
              for (lnode = 0; lnode < npe; lnode++) {
                delx -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][0];
                dely -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][1];
                delz -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][2];
              }
            }
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < cutsq) {
              nbonds++;
              bonds[kk]++;
              bonds[ll]++;
            }
          }
        }

        cna[jj][NBOND] = nbonds;

        maxbonds = 0;
        minbonds = MAXCOMMON;
        for (n = 0; n < ncommon; n++) {
          maxbonds = MAX(bonds[n],maxbonds);
          minbonds = MIN(bonds[n],minbonds);
        }
        cna[jj][MAXBOND] = maxbonds;
        cna[jj][MINBOND] = minbonds;

      }
    }

    // detect CNA pattern of atom I
    nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
    if (annearest[i] == 12) {
      for (inear = 0; inear < 12; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 2 && cl == 1 && cm == 1) nfcc++;
        else if (cj == 4 && ck == 2 && cl == 2 && cm == 0) nhcp++;
        else if (cj == 5 && ck == 5 && cl == 2 && cm == 2) nico++;
      }
      if (nfcc == 12) apattern[i] = FCC;
      else if (nfcc == 6 && nhcp == 6) apattern[i] =  HCP;
      else if (nico == 12) apattern[i] = ICOS;
    } else if (annearest[i] == 14) {
      for (inear = 0; inear < 14; inear++) {
        cj = cna[inear][NCOMMON];
        ck = cna[inear][NBOND];
        cl = cna[inear][MAXBOND];
        cm = cna[inear][MINBOND];
        if (cj == 4 && ck == 4 && cl == 2 && cm == 2) nbcc4++;
        else if (cj == 6 && ck == 6 && cl == 2 && cm == 2) nbcc6++;
      }
      if (nbcc4 == 6 && nbcc6 == 8) apattern[i] =  BCC;
    }

  }

  /*-------------------------------------------
   * find cna structure for interpolated atoms
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
      iapattern[i][iintpl] = OTHER;
      if (iannearest[i][iintpl] != 12 && iannearest[i][iintpl] != 14)
        continue;

      // loop over near neighbors of I to build cna data structure
      // cna[k][NCOMMON] = # of common neighbors of I with each of its neighs
      // cna[k][NBONDS] = # of bonds between those common neighbors
      // cna[k][MAXBOND] = max # of bonds of any common neighbor
      // cna[k][MINBOND] = min # of bonds of any common neighbor

      for (jj = 0; jj < iannearest[i][iintpl]; jj++) {

        j = ianearest[i][iintpl][jj];

        // common = list of neighbors common to I and J

        // J neighbor is atom
        if (ianearesttype[i][iintpl][jj] == ATOM) {

          if (j < nalocal) {

            // if J is an owned atom, use its near neighbor list to find common

            firstflag = 1;
            ncommon = 0;
            for (inear = 0; inear < iannearest[i][iintpl]; inear++)
              for (jnear = 0; jnear < annearest[j]; jnear++)

                // neighbor index and type must match to be common
                if (ianearest[i][iintpl][inear] == anearest[j][jnear] &&
                    ianearesttype[i][iintpl][inear] == anearesttype[j][jnear] ) {
                  if (ncommon < MAXCOMMON) {
                    common[ncommon] = anearest[j][jnear];
                    commontype[ncommon++] = anearesttype[j][jnear];
                  } else if (firstflag) {
                    nerror++;
                    firstflag = 0;
                  }
                }
          } else {

            // if J is a ghost atom, use full neighbor list of I to find them
            // must exclude J from I's neighbor list

            jx = ax[j][0];
            jy = ax[j][1];
            jz = ax[j][2];

            n = 0;

            iintpl_local = e2ialist[i] + iintpl;
            // loop over atom neighbor K of I
            knum = numneighia2a[iintpl_local];
            if (knum) {
              klist = firstneighia2a[iintpl_local];
              for (kk = 0; kk < knum; kk++) {
                k = klist[kk];
                //k &= NEIGHMASK;
                if (k == j) continue;
                delx = jx - ax[k][0];
                dely = jy - ax[k][1];
                delz = jz - ax[k][2];
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutsq) 
                  if (n < MAXNEAR) {
                    onenearest[n] = k;
                    onenearesttype[n++] = ATOM;
                  } else break;

              }
            }
            // loop over interpolated atom neighbor K of I
            knum = numneighia2ia[iintpl_local];
            if (knum) {
              klist = firstneighia2ia[iintpl_local];
              kintpllist = firstneighia2ia_index[iintpl_local];
              for (kk = 0; kk < knum; kk++) {
                k = klist[kk];
                ketype = etype[k];
                kintpl = kintpllist[kk];
                kintpl_local = e2ialist[k] + kintpl;
                //k &= NEIGHMASK;
                delx = jx; dely = jy; delz = jz;
                for (knode = 0; knode < npe; knode++) {
                  delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                  dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                  delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
                }
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutsq) {
                  if (n < MAXNEAR) {
                    onenearest[n] = kintpl_local;
                    onenearesttype[n++] = INTPL;
                  } else break;
                }
              }
            }
            firstflag = 1;
            ncommon = 0;

            for (inear = 0; inear < iannearest[i][iintpl]; inear++)
              for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
                if (ianearest[i][iintpl][inear] == onenearest[jnear] && 
                    ianearesttype[i][iintpl][inear] == onenearesttype[jnear]) 
                  if (ncommon < MAXCOMMON) {
                    common[ncommon] = onenearest[jnear];
                    commontype[ncommon++] = onenearesttype[jnear];
                  } else if (firstflag) {
                    nerror++;
                    firstflag = 0;
                  }
          }
          cna[jj][NCOMMON] = ncommon;

          for (n = 0; n < ncommon; n++) bonds[n] = 0;

          nbonds = 0;
          for (kk = 0; kk < ncommon - 1; kk++) {
            k = common[kk];

            // get the coords of K based on its type (ATOM or INTPL)

            if (commontype[kk] == ATOM) {
              kx = ax[k][0];
              ky = ax[k][1];
              kz = ax[k][2];
            } else {
              ek = ia2elist[k];
              kintpl = k - e2ialist[ek];
              ketype = etype[ek];
              kx = ky = kz = 0.0;
              for (knode = 0; knode < npe; knode++) {
                kx += shape_array[ketype][kintpl][knode]*nodex[ek][knode][0];
                ky += shape_array[ketype][kintpl][knode]*nodex[ek][knode][1];
                kz += shape_array[ketype][kintpl][knode]*nodex[ek][knode][2];
              }
            }

            for (ll = kk + 1; ll < ncommon; ll++) {
              l = common[ll];

              // get the coords of L based on its type (ATOM or INTPL)

              delx = kx; dely = ky; delz = kz; 
              if (commontype[ll] == ATOM) {
                delx -= ax[l][0];
                dely -= ax[l][1];
                delz -= ax[l][2];
              } else {
                el = ia2elist[l];
                lintpl = l - e2ialist[el];
                letype = etype[el];
                for (lnode = 0; lnode < npe; lnode++) {
                  delx -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][0];
                  dely -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][1];
                  delz -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][2];
                }
              }
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) {
                nbonds++;
                bonds[kk]++;
                bonds[ll]++;
              }
            }
          }

          cna[jj][NBOND] = nbonds;

          maxbonds = 0;
          minbonds = MAXCOMMON;
          for (n = 0; n < ncommon; n++) {
            maxbonds = MAX(bonds[n],maxbonds);
            minbonds = MIN(bonds[n],minbonds);
          }
          cna[jj][MAXBOND] = maxbonds;
          cna[jj][MINBOND] = minbonds;
        } else {

          // J neighbor is interpolated atom

          // index of parent element of J

          ej = ia2elist[j];

          // index of J in its parent element

          jintpl = j - e2ialist[ej]; 

          if (ej < nelocal) {

            // if J is an owned interpolated atom, use its near neighbor list to find common

            firstflag = 1;
            ncommon = 0;
            for (inear = 0; inear < iannearest[i][iintpl]; inear++)
              for (jnear = 0; jnear < iannearest[ej][jintpl]; jnear++)

                // neighbor index and type must match to be common

                if (ianearest[i][iintpl][inear] == ianearest[ej][jintpl][jnear] &&
                    ianearesttype[i][iintpl][inear] == ianearesttype[ej][jintpl][jnear] ) {
                  if (ncommon < MAXCOMMON) {
                    common[ncommon] = ianearest[i][iintpl][inear];
                    commontype[ncommon++] = ianearesttype[i][iintpl][inear];
                  } else if (firstflag) {
                    nerror++;
                    firstflag = 0;
                  }
                }
          } else {

            // if J is a ghost interpolated atom, use full neighbor list of I to find them
            // must exclude J from I's neighbor list

            jetype = etype[ej]; 

            jx = jy = jz = 0.0;
            for (jnode = 0; jnode < npe; jnode++) {
              jx += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][0];
              jy += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][1];
              jz += shape_array[jetype][jintpl][jnode]*nodex[ej][jnode][2];
            }

            n = 0;

            iintpl_local = e2ialist[i] + iintpl;

            // loop over atom neighbor K of I

            knum = numneighia2a[iintpl_local];
            if (knum) {
              klist = firstneighia2a[iintpl_local];
              for (kk = 0; kk < knum; kk++) {
                k = klist[kk];
                //k &= NEIGHMASK;
                delx = jx - ax[k][0];
                dely = jy - ax[k][1];
                delz = jz - ax[k][2];
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutsq) 
                  if (n < MAXNEAR) {
                    onenearest[n] = k;
                    onenearesttype[n++] = ATOM;
                  } else break;
              }
            }

            // loop over interpolated atom neighbor K of I

            knum = numneighia2ia[iintpl_local];
            if (knum) {
              klist = firstneighia2ia[iintpl_local];
              kintpllist = firstneighia2ia_index[iintpl_local];
              for (kk = 0; kk < knum; kk++) {
                k = klist[kk];
                kintpl = kintpllist[kk];

                if (ej == k && jintpl == kintpl) continue;
                ketype = etype[k];
                kintpl_local = e2ialist[k] + kintpl;
                //k &= NEIGHMASK;
                delx = jx; dely = jy; delz = jz;
                for (knode = 0; knode < npe; knode++) {
                  delx -= shape_array[ketype][kintpl][knode]*nodex[k][knode][0];
                  dely -= shape_array[ketype][kintpl][knode]*nodex[k][knode][1];
                  delz -= shape_array[ketype][kintpl][knode]*nodex[k][knode][2];
                }
                rsq = delx*delx + dely*dely + delz*delz;
                if (rsq < cutsq) 
                  if (n < MAXNEAR) {
                    onenearest[n] = kintpl_local;
                    onenearesttype[n++] = INTPL;
                  } else break;
              }
            }

            firstflag = 1;
            ncommon = 0;
            for (inear = 0; inear < iannearest[i][iintpl]; inear++)
              for (jnear = 0; (jnear < n) && (n < MAXNEAR); jnear++)
                if (ianearest[i][iintpl][inear] == onenearest[jnear] && 
                    ianearesttype[i][iintpl][inear] == onenearesttype[jnear]) {
                  if (ncommon < MAXCOMMON) {
                    common[ncommon] = onenearest[jnear];
                    commontype[ncommon++] = onenearesttype[jnear];
                  } else if (firstflag) {
                    nerror++;
                    firstflag = 0;
                  }
                }
          }

          cna[jj][NCOMMON] = ncommon;

          for (n = 0; n < ncommon; n++) bonds[n] = 0;
          nbonds = 0;
          for (kk = 0; kk < ncommon - 1; kk++) {
            k = common[kk];

            // get the coords of K based on its type (ATOM or INTPL)

            if (commontype[kk] == ATOM) {
              kx = ax[k][0];
              ky = ax[k][1];
              kz = ax[k][2];
            } else {
              ek = ia2elist[k];
              kintpl = k - e2ialist[ek];
              ketype = etype[ek];
              kx = ky = kz = 0.0;
              for (knode = 0; knode < npe; knode++) {
                kx += shape_array[ketype][kintpl][knode]*nodex[ek][knode][0];
                ky += shape_array[ketype][kintpl][knode]*nodex[ek][knode][1];
                kz += shape_array[ketype][kintpl][knode]*nodex[ek][knode][2];
              }
            }

            for (ll = kk + 1; ll < ncommon; ll++) {
              l = common[ll];

              // get the coords of L based on its type (ATOM or INTPL)

              delx = kx; dely = ky; delz = kz; 
              if (commontype[ll] == ATOM) {
                delx -= ax[l][0];
                dely -= ax[l][1];
                delz -= ax[l][2];
              } else {
                el = ia2elist[l];
                lintpl = l - e2ialist[el];
                letype = etype[el];
                for (lnode = 0; lnode < npe; lnode++) {
                  delx -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][0];
                  dely -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][1];
                  delz -= shape_array[letype][lintpl][lnode]*nodex[el][lnode][2];
                }
              }

              // check if L is bonded with K

              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) {
                nbonds++;
                bonds[kk]++;
                bonds[ll]++;
              }
            }
          }
          cna[jj][NBOND] = nbonds;
          maxbonds = 0;
          minbonds = MAXCOMMON;
          for (n = 0; n < ncommon; n++) {
            maxbonds = MAX(bonds[n],maxbonds);
            minbonds = MIN(bonds[n],minbonds);
          }
          cna[jj][MAXBOND] = maxbonds;
          cna[jj][MINBOND] = minbonds;
        }
      }


      // detect CNA pattern of the interpolated atom I

      nfcc = nhcp = nbcc4 = nbcc6 = nico = 0;
      if (iannearest[i][iintpl] == 12) {
        for (inear = 0; inear < 12; inear++) {
          cj = cna[inear][NCOMMON];
          ck = cna[inear][NBOND];
          cl = cna[inear][MAXBOND];
          cm = cna[inear][MINBOND];
          if (cj == 4 && ck == 2 && cl == 1 && cm == 1) nfcc++;
          else if (cj == 4 && ck == 2 && cl == 2 && cm == 0) nhcp++;
          else if (cj == 5 && ck == 5 && cl == 2 && cm == 2) nico++;
        }
        if (nfcc == 12) iapattern[i][iintpl] = FCC;
        else if (nfcc == 6 && nhcp == 6) iapattern[i][iintpl] =  HCP;
        else if (nico == 12) iapattern[i][iintpl] = ICOS;

      } else if (iannearest[i][iintpl] == 14) {
        for (inear = 0; inear < 14; inear++) {
          cj = cna[inear][NCOMMON];
          ck = cna[inear][NBOND];
          cl = cna[inear][MAXBOND];
          cm = cna[inear][MINBOND];
          if (cj == 4 && ck == 4 && cl == 2 && cm == 2) nbcc4++;
          else if (cj == 6 && ck == 6 && cl == 2 && cm == 2) nbcc6++;
        }
        if (nbcc4 == 6 && nbcc6 == 8) iapattern[i][iintpl] =  BCC;
      }
    }
  }

  // warning message

  MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
  if (nerrorall && comm->me == 0) {
    char str[128];
    sprintf(str,"Too many common neighbors in CNA %d times",nerrorall);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeCNAAtom::memory_usage()
{
  double maxintpl = element->maxintpl;
  double bytes = namax * sizeof(int);
  bytes += 2 * namax * MAXNEAR * sizeof(int);
  bytes += namax * sizeof(double);
  bytes += 2 * nemax * maxintpl * sizeof(int);
  bytes += nemax * maxintpl * MAXNEAR * sizeof(int);
  bytes += nemax * maxintpl * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   grow atom-based arrays
   ------------------------------------------------------------------------- */

void ComputeCNAAtom::grow_atom(int nmax) {
  memory->destroy(anearest);
  memory->destroy(anearesttype);
  memory->destroy(annearest);
  memory->destroy(apattern);
  memory->create(anearest,nmax,MAXNEAR,"cna:anearest");
  memory->create(anearesttype,nmax,MAXNEAR,"cna:anearesttype");
  memory->create(annearest,nmax,"cna:annearest");
  memory->create(apattern,nmax,"cna:apattern");
}

/* ----------------------------------------------------------------------
   grow interpolated atom-based arrays
   ------------------------------------------------------------------------- */

void ComputeCNAAtom::grow_intpl(int nmax,int nintpl) {
  memory->destroy(ianearest);
  memory->destroy(ianearesttype);
  memory->destroy(iannearest);
  memory->destroy(iapattern);
  memory->create(ianearest,nmax,nintpl,MAXNEAR,"cna:ianearest");
  memory->create(ianearesttype,nmax,nintpl,MAXNEAR,"cna:ianearesttype");
  memory->create(iannearest,nmax,nintpl,"cna:iannearest");
  memory->create(iapattern,nmax,nintpl,"cna:iapattern");
}

/*-----------------------------------------------------------------------------------------------*/

int ComputeCNAAtom::pack_atom_forward_comm(int n, int *list, double *buf,
    int pbc_flag, int *pbc)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = apattern[j];
  }
  return m;
}
/*------------------------------------------------------------------------------------------*/

void ComputeCNAAtom::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) apattern[i] = buf[m++];
}


/*---------------------------------------------------------------------------------------------*/

int ComputeCNAAtom::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k,m;
  int nintpl = element->maxintpl;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nintpl; k++) buf[m++] = iapattern[j][k];
  }
  return m;
}

/*----------------------------------------------------------------------------------------------*/

void ComputeCNAAtom::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i,j,m,last;
  int nintpl = element->maxintpl;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (j = 0; j < nintpl; j++) iapattern[i][j] = buf[m++];
}


