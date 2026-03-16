#include <stdlib.h>
#include <string.h>
#include "compute_stress_mech.h"
#include "atom.h"
#include "element.h"
#include "comm.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"
#include "my_page.h"
#include "universe.h"

using namespace CAC_NS;

#define PGDELTA 1
#define DELTA 1024
#define EPSILON 1e-6

/* ---------------------------------------------------------------------- */

ComputeStressMech::ComputeStressMech(CAC *cac, int narg, char **arg) :
  Compute(cac, narg, arg),
  surface(NULL)
{
  // surface neighbor list initialization

  a_ab_page = ia_ab_page = ia_index_ab_page = NULL;
  a_bl_page = ia_bl_page = ia_index_bl_page = NULL;
  numneigh_a_ab = numneigh_ia_ab = NULL;
  numneigh_a_bl = numneigh_ia_bl = NULL;
  firstneigh_a_ab = firstneigh_ia_ab = firstneigh_ia_index_ab = NULL;
  firstneigh_a_bl = firstneigh_ia_bl = firstneigh_ia_index_bl = NULL;
  pgsize = 0;
  oneatom = 0;
  oneatom_user = 0;

  snum = 0;

  // surface at middle of cell by default

  xsurface[0] = xsurface[1] = xsurface[2] = 0.5;

  if (narg < 13) error->all(FLERR,"Illegal compute stress/mech command");

  // compute flag and default values

  scalar_flag = array_flag = 1;
  size_array_cols = 6; // columns 0 1 2 = x y z, columns 3 4 5 = sx sy sz
  size_array_rows = 0;
  size_array_rows_variable = 1;
  ghostskinflag = 1;

  surface_arrange_style = 1;
  triclinic = domain->triclinic;

  // process arguments 
  // list of argument in the following order:
  // dir nx ny nz xlo xhi ylo yhi zlo zhi + optional arguments

  int iarg; 
  if (strcmp(arg[3],"x") == 0)
    surface_direction = 0;
  else if (strcmp(arg[3],"y") == 0)
    surface_direction = 1;
  else if (strcmp(arg[3],"z") == 0)
    surface_direction = 2;
  else error->all(FLERR,"Illegal surface direction in compute stress/mech command");

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int *periodicity = domain->periodicity;

  for (int i = 0; i < 3; i++) {
    nsurface[i] = universe->inumeric(FLERR,arg[i+4]);
    if (nsurface[i] <= 0) error->all(FLERR,"Illegal value for number of surfaces in compute stress/mech command");

    iarg = i*2 + 7;
    if (strcmp(arg[iarg],"EDGE") == 0) {
      if (domain->box_exist == 0)
        error->all(FLERR,"Cannot use EDGE when box does not exist");
      if (triclinic) 
        error->all(FLERR,"Cannot use EDGE for triclinic box");
      lo[i] = boxlo[i];
    } else {
      lo[i] = universe->numeric(FLERR,arg[iarg]);
      if (periodicity[i] && lo[i] < boxlo[i])
        error->all(FLERR,"Surface bound must not exceed box bound for periodicity in compute stress/mech command");
    }

    iarg++;
    if (strcmp(arg[iarg],"EDGE") == 0) {
      if (domain->box_exist == 0)
        error->all(FLERR,"Cannot use EDGE when box does not exist");
      if (triclinic)
        error->all(FLERR,"Cannot use EDGE for triclinic box");
      hi[i] = boxhi[i];
    } else {
      hi[i] = universe->numeric(FLERR,arg[iarg]);
      if (periodicity[i] && hi[i] > boxhi[i])
        error->all(FLERR,"Surface bound must not exceed box bound for periodicity in compute stress/mech command");
    }

    if (hi[i] < lo[i]) error->all(FLERR,"Illegal bound values in compute stress/mech command");
  }

  if (triclinic) {
    error->all(FLERR,"Compute stress/mech does not work with triclinic box yet");
    // check surface box for triclinic case
  }

  surface_area = 1;
  for (int i = 0; i < 3; i++) {
    delta[i] = (hi[i]-lo[i])/nsurface[i];
    if (surface_direction != i) {
      ghostskin = MAX(ghostskin,delta[i]/2.0);
      surface_area *= delta[i];
    }
  } 

  // process optional args

  iarg = 13;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"style") == 0) {
      if (strcmp(arg[iarg+1],"standard") == 0)
        surface_arrange_style = 1;
      else if (0) {

        // add to include more styles

      } else
        error->all(FLERR,"Illegal compute stress/mech command");
    } else if (strcmp(arg[iarg],"position") == 0) {
      double x = universe->numeric(FLERR,arg[iarg+1]);
      if (x < 0.0 || x > 1.0) 
        error->all(FLERR,"Surface position must be between 0 and 1 of the cell");
      xsurface[surface_direction] = x;
    } else if (strcmp(arg[iarg],"one") == 0) {
      oneatom_user = universe->inumeric(FLERR,arg[iarg+1]);
    } else 
      error->all(FLERR,"Illegal compute stress/mech command");
    iarg += 2;
  }

  snum = nsurface[0]*nsurface[1]*nsurface[2];
  size_array_rows = snum;
  allocate();

}

/* ---------------------------------------------------------------------- */

ComputeStressMech::~ComputeStressMech()
{
  memory->destroy(surface);

  memory->destroy(numneigh_a_ab);
  memory->destroy(numneigh_ia_ab);
  memory->sfree(firstneigh_a_ab);
  memory->sfree(firstneigh_ia_ab);
  memory->sfree(firstneigh_ia_index_ab);

  memory->destroy(numneigh_a_bl);
  memory->destroy(numneigh_ia_bl);
  memory->sfree(firstneigh_a_bl);
  memory->sfree(firstneigh_ia_bl);
  memory->sfree(firstneigh_ia_index_bl);

  delete [] a_ab_page;
  delete [] ia_ab_page;
  delete [] ia_index_ab_page;
  delete [] a_bl_page;
  delete [] ia_bl_page;
  delete [] ia_index_bl_page;
}

/* ---------------------------------------------------------------------- */

void ComputeStressMech::init()
{

  if (force->pair == NULL)
    error->all(FLERR,"Compute stress/mech requires a pair style be defined");
  if (force->pair->threebody_flag) 
    error->all(FLERR,"Compute stress/mech requires doesn't work with 3-body pair potential yet");
  if (domain->triclinic) error->all(FLERR,"Compute stress/mech doesn't work with triclinic box yet");

  triclinic = domain->triclinic;

  for (int i = 0; i < 3; i++) {
    cutneighsq[i] = force->pair->cutforce;
    if (i != surface_direction) cutneighsq[i] += delta[i];
    cutneighsq[i] = cutneighsq[i]*cutneighsq[i];
  }

  // local compute stress/mech surface neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;

  // check one page == NULL is enough

  if (a_ab_page == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {

    delete [] a_ab_page;
    delete [] a_bl_page;
    delete [] ia_ab_page;
    delete [] ia_bl_page;
    delete [] ia_index_ab_page;
    delete [] ia_index_bl_page;

    pgsize = neighbor->pgsize;
    if (oneatom_user) oneatom = oneatom_user;
    else oneatom = neighbor->oneatom;

    int nmypage = comm->nthreads;
    a_ab_page = new MyPage<int>[nmypage];
    ia_ab_page = new MyPage<int>[nmypage];
    ia_index_ab_page = new MyPage<int>[nmypage];
    a_bl_page = new MyPage<int>[nmypage];
    ia_bl_page = new MyPage<int>[nmypage];
    ia_index_bl_page = new MyPage<int>[nmypage];

    for (int i = 0; i < nmypage; i++) {
      a_ab_page[i].init(oneatom,pgsize);
      ia_ab_page[i].init(oneatom,pgsize);
      ia_index_ab_page[i].init(oneatom,pgsize);
      a_bl_page[i].init(oneatom,pgsize);
      ia_bl_page[i].init(oneatom,pgsize);
      ia_index_bl_page[i].init(oneatom,pgsize);
    }
  }
}

/* -------------------------------------
   Check if force between 2 atoms (x1 & x2) crosses surface at x0 
   if force vector is parallel to surface, 
   crossing if atoms are on the surface
   the crossing point must be within box to avoid double count from other procs
   return 0 if not cross surface or distance > force cutoff
          1 if cross surface
   ------------------------------------- */

int ComputeStressMech::is_crossing(double *x0, double *x1, double *x2) 
{
  int d = surface_direction;
  double a[3];
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  for (int i = 0; i < 3; i++) a[i] = x2[i] - x1[i];

  // check if distance < force cutoff

  double cutforcesq = cutneighsq[d];
  double rsq = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
  if (rsq >= cutforcesq) return 0;

  // check if parallel

  if (a[d] == 0.0) {
    if (fabs(x1[d]-x0[d]) < EPSILON) return 1;
    else return 0;
  }


  double t = (x0[d]-x1[d])/a[d];
  double x_cross[3];
  for (int i = 0; i < 3; i++) {
    x_cross[i] = a[i]*t + x1[i];

    // check if force vector is crossing surface

    if (i != d && fabs(x_cross[i]-x0[i]) > delta[i]/2.0) return 0;

    // check if intersection point is within bound

    if (sublo[i] > x_cross[i] || x_cross[i] >= subhi[i]) return 0;

  }

  return 1;
}

/* -------------------------------------------------------------------------- */

void ComputeStressMech::setup_surfaces()
{
  int s = 0;

  if (surface_arrange_style == 1) {
    double xtemp,ytemp,ztemp;
    for (int i = 0; i < nsurface[0]; i++) {
      xtemp = lo[0] + (i + xsurface[0])*delta[0];
      for (int j = 0; j < nsurface[1]; j++) {
        ytemp = lo[1] + (j + xsurface[1])*delta[1];
        for (int k = 0; k < nsurface[2]; k++) {
          ztemp = lo[2] + (k + xsurface[2])*delta[2];

          surface[s][0] = xtemp;          
          surface[s][1] = ytemp;          
          surface[s][2] = ztemp;          
          s++;

        }
      }
    }
  } else {

    // add more to include more style

  }
}

/* -------------------------------------------------------------------------- */

void ComputeStressMech::allocate()
{

  array = memory->create(surface,snum,size_array_cols,"compute:surface");

  memory->create(numneigh_a_ab,snum,"compute:numneigha_ab");
  memory->create(numneigh_ia_ab,snum,"compute:numneighia_ab");
  memory->create(numneigh_a_bl,snum,"compute:numneigha_bl");
  memory->create(numneigh_ia_bl,snum,"compute:numneighia_bl");
  firstneigh_a_ab = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneigha_ab");
  firstneigh_ia_ab = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneighia_ab");
  firstneigh_ia_index_ab = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneighia_index_ab");
  firstneigh_a_bl = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneigha_bl");
  firstneigh_ia_bl = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneighia_bl");
  firstneigh_ia_index_bl = (int **) memory->smalloc(snum*sizeof(int *),"compute:firstneighia_index_bl");
}

/* ----------------------------------------------------------------------
   find local below neighbor to surface center 
   currently loop through all atoms, elements
   can implement better binning and stencil method to save time if needed
   (rectangular box neighboring instead of spherical)
   ---------------------------------------------------------------------- */

void ComputeStressMech::surface_neighbor()
{
  int atom_ntotal = atom->nlocal + atom->nghost;
  int elem_ntotal = element->nlocal + element->nghost;

  double **ax = atom->x;
  double **ex = element->x;
  int *etype = element->etype;
  int *npe = element->npe;
  double ***nodex = element->nodex;
  double ***shape_array = element->shape_array;
  double ***shape_array_center_subelem = element->shape_array_center_subelem;
  int ***ias2ia = element->ias2ia;
  int **natom_subelem = element->natom_subelem;
  int *nsubelem = element->nsubelem;

  int i,j,k,l,jetype,jsub,jnode,jintpl,ictype,jctype;
  double ix,iy,iz,jx,jy,jz,delx,dely,delz,rsq,delxsq,delysq,delzsq;
  int na_ab,na_bl,nia_ab,nia_bl;
  int *a_ab_neighptr,*a_bl_neighptr,*ia_ab_neighptr,*ia_bl_neighptr;
  int *ia_index_ab_neighptr,*ia_index_bl_neighptr;

  //*** NEED FIX ***
  double cutelemsq = force->pair->cutforce + ghostskin + element->max_diag_size/2.0;
  double cutsubelemsq = force->pair->cutforce + ghostskin + element->max_diag_size/2.0/element->subsplit[1];
  cutelemsq = cutelemsq*cutelemsq;
  cutsubelemsq = cutsubelemsq*cutsubelemsq;

  int d = surface_direction;

  // reset all pages

  a_ab_page->reset();
  ia_ab_page->reset();
  ia_index_ab_page->reset();
  a_bl_page->reset();
  ia_bl_page->reset();
  ia_index_bl_page->reset();

  // loop through all surfaces and find neighbor atoms and interpolated atoms for surface I

  for (i = 0; i < snum; i++) {
    ix = surface[i][0]; 
    iy = surface[i][1]; 
    iz = surface[i][2]; 

    // find neighboring atoms to surface

    if (atom_ntotal) {
      na_ab = na_bl = 0;
      a_ab_neighptr = a_ab_page->vget(); 
      a_bl_neighptr = a_bl_page->vget(); 
      for (j = 0; j < atom_ntotal; j++) {
        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        delxsq = delx*delx;
        delysq = dely*dely;
        delzsq = delz*delz;

        if (delxsq < cutneighsq[0] && delysq < cutneighsq[1] && delzsq < cutneighsq[2]) {
          if (ax[j][d] < surface[i][d]) 
            a_bl_neighptr[na_bl++] = j;
          else
            a_ab_neighptr[na_ab++] = j;
        }
      }
      firstneigh_a_ab[i] = a_ab_neighptr;
      firstneigh_a_bl[i] = a_bl_neighptr;
      numneigh_a_ab[i] = na_ab;
      numneigh_a_bl[i] = na_bl;

      a_ab_page->vgot(na_ab);
      a_bl_page->vgot(na_bl);
      if (a_ab_page->status() || a_bl_page->status()) 
        error->one(FLERR,"Compute stress/mech neighbor list overflow, boost one");
    } else {
      firstneigh_a_ab[i] = NULL;
      firstneigh_a_bl[i] = NULL;
      numneigh_a_ab[i] = 0;
      numneigh_a_bl[i] = 0;
    }

    // find neighboring interpolated atoms to surface

    if (elem_ntotal) {
      nia_ab = nia_bl = 0;
      ia_ab_neighptr = ia_ab_page->vget(); 
      ia_bl_neighptr = ia_bl_page->vget(); 
      ia_index_ab_neighptr = ia_index_ab_page->vget(); 
      ia_index_bl_neighptr = ia_index_bl_page->vget(); 
      for (j = 0; j < elem_ntotal; j++) {
        jetype = etype[j];
        delx = ix - ex[j][0];
        dely = iy - ex[j][1];
        delz = iz - ex[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        // check if element J is a neighbor, then loop through its sub-elements

        if (rsq < cutelemsq) {
          for (jsub = 0; jsub < nsubelem[jetype]; jsub++) {

            // interpolating sub-element center coordinates

            delx = ix; dely = iy; delz = iz;
            for (jnode = 0; jnode < npe[jetype]; jnode++) {
              delx -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][0];
              dely -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][1];
              delz -= shape_array_center_subelem[jetype][jsub][jnode]*nodex[j][jnode][2];
            }
            rsq = delx*delx + dely*dely + delz*delz;

            // check if sub-element J is a neighbor, then loop through its interpolated atoms

            if (rsq < cutsubelemsq) {
              for (l = 0; l < natom_subelem[jetype][jsub]; l++) {

                // interpolating interpolated atom coordinates

                jx = jy = jz = 0.0;
                jintpl = ias2ia[jetype][jsub][l];
                for (jnode = 0; jnode < npe[jetype]; jnode++) {
                  jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
                  jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
                  jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
                }
                delx = ix - jx; 
                dely = iy - jy; 
                delz = iz - jz; 
                delxsq = delx*delx;
                delysq = dely*dely;
                delzsq = delz*delz;

                // check if interpolated atom is a neighbor

                if (delxsq < cutneighsq[0] && delysq < cutneighsq[1] && delzsq < cutneighsq[2]) {
                  if ((d == 0 && jx < ix) ||
                      (d == 1 && jy < iy) ||
                      (d == 2 && jz < iz)) {
                    ia_bl_neighptr[nia_bl] = j;
                    ia_index_bl_neighptr[nia_bl++] = jintpl;
                  } else {
                    ia_ab_neighptr[nia_ab] = j;
                    ia_index_ab_neighptr[nia_ab++] = jintpl;
                  }
                }
              }
            }
          }
        }
      }

      firstneigh_ia_ab[i] = ia_ab_neighptr;
      firstneigh_ia_index_ab[i] = ia_index_ab_neighptr;
      firstneigh_ia_bl[i] = ia_bl_neighptr;
      firstneigh_ia_index_bl[i] = ia_index_bl_neighptr;
      numneigh_ia_ab[i] = nia_ab;
      numneigh_ia_bl[i] = nia_bl;
      ia_ab_page->vgot(nia_ab);
      ia_index_ab_page->vgot(nia_ab);
      ia_bl_page->vgot(nia_bl);
      ia_index_bl_page->vgot(nia_bl);
      if (ia_ab_page->status() ||
          ia_bl_page->status() ||
          ia_index_ab_page->status() ||
          ia_index_bl_page->status())
        error->one(FLERR,"Compute stress/mech neighbor list overflow, boost one");
    } else {
      firstneigh_ia_ab[i] = NULL;
      firstneigh_ia_bl[i] = NULL;
      firstneigh_ia_index_ab[i] = NULL;
      firstneigh_ia_index_bl[i] = NULL;
      numneigh_ia_ab[i] = 0;
      numneigh_ia_bl[i] = 0;
    }
  }
}



/* ---------------------------------------------------------------------- */

double ComputeStressMech::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  scalar = surface_direction;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeStressMech::compute_array()
{
  int i,j,k,jnum,knum,jj,kk,node,jetype,ketype,jintpl,kintpl;
  int *jlist,*jindexlist,*klist,*kindexlist;
  double jx[3],kx[3],f[3],stress[3];
  double **atomx = atom->x;
  double ***nodex = element->nodex;
  double ***shape_array = element->shape_array;
  int *etype = element->etype;
  int *npe = element->npe;

  double nktv2p = force->nktv2p;

  invoked_array = update->ntimestep;

  // setup surfaces

  setup_surfaces();

  // build neighbor

  surface_neighbor();

  int count = 0;
  for (i = 0; i < snum; i++) {

    // clear stress

    stress[0] = stress[1] = stress[2] = 0.0;

    // loop through below neighboring atoms (J)

    jnum = numneigh_a_bl[i];
    if (jnum) {
      jlist = firstneigh_a_bl[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj]; 

        // loop through all above neighboring atoms (K)

        knum = numneigh_a_ab[i];
        if (knum) {
          klist = firstneigh_a_ab[i];
          for (kk = 0; kk < knum; kk++) {
            k = klist[kk];

            // check if force is crossing surface

            if (is_crossing(surface[i],atomx[j],atomx[k])) {

              // calculate force from above atom (K) acting on below atom (J)

              force->pair->pair_a2a(j,k,f);
              for (int l = 0; l < 3; l++)
                stress[l] += f[l];
            }
          }
        }

        // loop through all above neighboring interpolated atoms

        knum = numneigh_ia_ab[i];
        if (knum) {
          klist = firstneigh_ia_ab[i];
          kindexlist = firstneigh_ia_index_ab[i];
          for (kk = 0; kk < knum; kk++) {
            k = klist[kk];
            kintpl = kindexlist[kk];
            ketype = etype[k];

            // interpolate position of interpolated atom K

            kx[0] = kx[1] = kx[2] = 0.0;
            for (node = 0; node < npe[ketype]; node++)
              for (int l = 0; l < 3; l++)
                kx[l] += shape_array[ketype][kintpl][node]*nodex[k][node][l];

            // check if force is crossing surface && distance < force cutoff

            if (is_crossing(surface[i],atomx[j],atomx[k])) {

              // calculate force from above atom (K) acting on below atom (J)

              force->pair->pair_a2ia(j,k,kintpl,kx,f);
              for (int l = 0; l < 3; l++)
                stress[l] += f[l];
            }
          }
        }
      }
    }

    // loop through below neighboring interpolated atoms (J)

    jnum = numneigh_ia_bl[i];
    if (jnum) {
      jlist = firstneigh_ia_bl[i];
      jindexlist = firstneigh_ia_index_bl[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jindexlist[jj];
        jetype = etype[j];

        // interpolate position of interpolated atom J

        jx[0] = jx[1] = jx[2] = 0.0;
        for (node = 0; node < npe[jetype]; node++)
          for (int l = 0; l < 3; l++)
            jx[l] += shape_array[jetype][jintpl][node]*nodex[j][node][l];

        // loop through all above neighboring atoms K

        knum = numneigh_a_ab[i];
        if (knum) {
          klist = firstneigh_a_ab[i];
          for (kk = 0; kk < knum; kk++) {
            k = klist[kk];

            // check if force is crossing surface

            if (is_crossing(surface[i],jx,atomx[k])) {

              // calculate force from above atom (K) acting on below atom (J) 
              // (negative sign since calculated force is from atom J acting on atom K)

              force->pair->pair_a2ia(k,j,jintpl,jx,f);
              for (int l = 0; l < 3; l++)
                stress[l] -= f[l];
            }
          }
        }

        // loop through all above neighboring interpolated atoms

        knum = numneigh_ia_ab[i];
        if (knum) {
          klist = firstneigh_ia_ab[i];
          kindexlist = firstneigh_ia_index_ab[i];
          for (kk = 0; kk < knum; kk++) {
            k = klist[kk];
            kintpl = kindexlist[kk];
            ketype = etype[k];

            // interpolate position of interpolated atom K

            kx[0] = kx[1] = kx[2] = 0.0;
            for (node = 0; node < npe[ketype]; node++)
              for (int l = 0; l < 3; l++)
                kx[l] += shape_array[ketype][kintpl][node]*nodex[k][node][l];

            // check if force is crossing surface && distance < force cutoff

            if (is_crossing(surface[i],jx,kx)) {

              // calculate force from above atom (K) acting on below atom (J)

              force->pair->pair_ia2ia(j,jintpl,jx,k,kintpl,kx,f);
              for (int l = 0; l < 3; l++)
                stress[l] += f[l];
            }
          }
        }
      }
    }

    // convert to stress (pressure) units

    stress[0] *= nktv2p/surface_area;
    stress[1] *= nktv2p/surface_area;
    stress[2] *= nktv2p/surface_area;

    // sum up stress components from all procs

    MPI_Allreduce(&stress[0],&surface[i][3],3,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

double ComputeStressMech::memory_usage()
{
  double bytes = 0;
  bytes += memory->usage(surface,snum,size_array_cols);
  bytes += memory->usage(numneigh_a_ab,snum);
  bytes += memory->usage(numneigh_ia_ab,snum);
  bytes += memory->usage(numneigh_a_bl,snum);
  bytes += memory->usage(numneigh_ia_bl,snum);

  int nmypage = comm->nthreads;

  if (a_ab_page) for (int i = 0; i < nmypage; i++) bytes += a_ab_page[i].size();
  if (ia_ab_page) for (int i = 0; i < nmypage; i++) bytes += ia_ab_page[i].size();
  if (ia_index_ab_page) for (int i = 0; i < nmypage; i++) bytes += ia_index_ab_page[i].size();
  if (a_bl_page) for (int i = 0; i < nmypage; i++) bytes += a_bl_page[i].size();
  if (ia_bl_page) for (int i = 0; i < nmypage; i++) bytes += ia_bl_page[i].size();
  if (ia_index_bl_page) for (int i = 0; i < nmypage; i++) bytes += ia_index_bl_page[i].size();

  return bytes;
}



