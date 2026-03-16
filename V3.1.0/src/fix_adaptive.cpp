#include <cstring>
#include <cstdlib>
#include "fix_adaptive.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "compute.h"
#include "domain.h"
#include "memory.h"
#include "group.h"
#include "variable.h"
#include "error.h"
#include "irregular_comm.h"
#include "universe.h"
#include "math_extra.h"

using namespace CAC_NS;
using namespace FixConst;
using namespace MathExtra;

#define MAXCUT 2
#define MAXSUB  20
#define BIG     1e30

/* ---------------------------------------------------------------------- */

FixAdaptive::FixAdaptive(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg), 
  id_structure(nullptr), splitlist(nullptr)
{
  boxindices[0][0] = 1;
  boxindices[0][1] = 4;
  boxindices[0][2] = 7;
  boxindices[0][3] = 10;
  boxindices[0][4] = 16;
  boxindices[0][5] = 19;
  boxindices[0][6] = 22;
  boxindices[0][7] = 25;
  boxindices[1][0] = 3;
  boxindices[1][1] = 4;
  boxindices[1][2] = 5;
  boxindices[1][3] = 12;
  boxindices[1][4] = 14;
  boxindices[1][5] = 21;
  boxindices[1][6] = 22;
  boxindices[1][7] = 23;
  boxindices[2][0] = 9;
  boxindices[2][1] = 10;
  boxindices[2][2] = 11;
  boxindices[2][3] = 12;
  boxindices[2][4] = 14;
  boxindices[2][5] = 15;
  boxindices[2][6] = 16;
  boxindices[2][7] = 17;

  maxplanes = 10;
  memory->create(planexi, maxplanes, "fix:planexi");
  memory->create(planecount, maxplanes, "fix:planecount");

  debug = 0;
  create_attribute = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  //if (strcmp(element->element_style, "cac") != 0) 
  //  error->all(FLERR, "Fix adaptive only works for rhb element style");
  //if (element->maxapc > 1) error->all(FLERR, "Fix adaptive does not work for multiple atoms per unit cell yet");

  if (narg < 6) error->all(FLERR, "Illegal fix adaptive command");

  nevery = universe->inumeric(FLERR, arg[3]);
  if (nevery <= 0) error->all(FLERR, "Illegal fix adaptive command");

  int n = strlen(arg[4]) + 1; 
  id_structure = new char[n];
  strcpy(id_structure, arg[4]);
  istructure = modify->find_compute(id_structure);
  if (istructure < 0)
    error->all(FLERR, "Compute ID for fix adative does not exist");
  structure = modify->compute[istructure];
  perfect_crystal = universe->inumeric(FLERR, arg[5]);

  /*
     id_stress = new char[n];
     strcpy(id_stress, arg[6]);
     istress = modify->find_compute(id_stress);
     if (istress < 0)
     error->all(FLERR, "Stress Compute ID for fix adative does not exist");
     stress = modify->compute[istress];
     */
  // optional args

  tol = 1e-3;
  eig_frac_thres = 0.5;
  cut_atom_neigh = 10;
  split_width = 0.5;
  plane_width_fraction = 0.2;
  discretize_flag = 0;
  dislocation_flag = 1;
  thres = 10;
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tol") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      tol = universe->numeric(FLERR, arg[iarg+1]);
      if (tol <= 0) error->all(FLERR, "Tolerant for fix adaptive command must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "cut") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      cut_atom_neigh = universe->numeric(FLERR, arg[iarg+1]);
      if (cut_atom_neigh < 0) error->all(FLERR, "Neighbor cut for fix adaptive command must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "swidth") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      split_width = universe->numeric(FLERR, arg[iarg+1]);
      if (split_width <= 0 || split_width > 1) error->all(FLERR, "Split width for fix adaptive command must be between 0 and 1 (unit cell)");
      iarg += 2;
    } else if (strcmp(arg[iarg], "pwidth") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      plane_width_fraction = universe->numeric(FLERR, arg[iarg+1]);
      if (plane_width_fraction <= 0 || plane_width_fraction > 1) error->all(FLERR, "Plane width fraction for fix adaptive command must be between 0 and 1 (unit cell)");
      iarg += 2;

    } else if (strcmp(arg[iarg], "thres") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      thres = universe->inumeric(FLERR, arg[iarg+1]);
      if (thres < 0) error->all(FLERR, "Invalid thres for fix adaptive command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "discretize") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      if (strcmp(arg[iarg+1], "yes") == 0)
        discretize_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0)
        discretize_flag = 0;
      else error->all(FLERR, "Illegal fix adaptive command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dislocation") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      if (strcmp(arg[iarg+1], "yes") == 0)
        dislocation_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0)
        dislocation_flag = 0;
      else error->all(FLERR, "Illegal fix adaptive command");
      iarg += 2;

    } else if (strcmp(arg[iarg], "debug") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix adaptive command");
      if (strcmp(arg[iarg+1], "yes") == 0)
        debug = 1;
      else if (strcmp(arg[iarg+1], "no") == 0)
        debug = 0;
      else error->all(FLERR, "Illegal fix adaptive command");
      iarg += 2;
    } else error->all(FLERR, "Illegal fix adaptive command");
  }

  if (dislocation_flag == 0 && discretize_flag == 0) 
    error->all(FLERR, "Illegal fix adaptive command, either dislocation and/or discretize flag must be set");
  maxsplit = 1024;
  grow_splitlist(); 

  pending = 0;
  nsplits_total = 0;
}

/* ---------------------------------------------------------------------- */

FixAdaptive::~FixAdaptive()
{
  delete [] id_structure;
  //delete [] id_stress;
  delete [] splitlist;
  memory->destroy(planexi);
  memory->destroy(planecount);
}

/* ---------------------------------------------------------------------- */

int FixAdaptive::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdaptive::init()
{
  // need an occasional full ucell neighbor list

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->atomlist = 0;
  neighbor->requests[irequest]->elemlist = 1;
  neighbor->requests[irequest]->elemonlylist = 1;
  neighbor->requests[irequest]->gausslist = 0;
  neighbor->requests[irequest]->atom_filter = 1;
  neighbor->requests[irequest]->atom_filter_icompute = istructure;
  neighbor->requests[irequest]->atom_filter_value = perfect_crystal;
  neighbor->requests[irequest]->atom_filter_preneighflag = 1;
  neighbor->requests[irequest]->group = 1;
  neighbor->requests[irequest]->groupbit = groupbit;
  neighbor->requests[irequest]->force_rebuild = 1;
  if (cut_atom_neigh > 0) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = cut_atom_neigh;
  }

  if (cut_atom_neigh > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR, "Fix adaptive cutoff may be too large to find "
        "ghost atom neighbors");

}

/* ---------------------------------------------------------------------- */

void FixAdaptive::init_list(int id, NeighList *ptr)
{
  list = ptr;
}


/* ---------------------------------------------------------------------- */

void FixAdaptive::setup(int vflag)
{
  //next_reneighbor = (update->ntimestep/nevery)*nevery;
}

void FixAdaptive::setup_pre_exchange()
{

  // need to exchange, acquire ghosts, build neigh list on the 
  // first step in order for compute to build neighbor list
  // same as in Verlet::setup()
  if (domain->triclinic) {
    domain->x2lamda(atom->nlocal, atom->x);
    domain->x2lamda(element->nlocal, element->x);
    domain->nodex2lamda(element->nlocal, element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();
  comm->exchange();
  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 
  //comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal + atom->nghost, atom->x);
    domain->lamda2x(element->nlocal + element->nghost, element->x);
    domain->lamda2nodex(element->nlocal + element->nghost, element->nodex);
  }
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();

  neighbor->build();
  neighbor->ncalls = 0;

  next_reneighbor = (update->ntimestep/nevery)*nevery;
  post_integrate();
  next_reneighbor = (update->ntimestep/nevery)*nevery+nevery;

  // if splitting
  // reset neighbor->lastcall to -1 
  // to force neigh list rebuild in Veret::setup()

  if (pending) {
    neighbor->lastcall = -1; 
    pending = 0;
  }
  
}

/* ----------------------------------------------------------------------
   if pending, dislocation detected at this time step,  
   set new next_reneighbor here after forcing neigh list rebuild
   ------------------------------------------------------------------------- */

void FixAdaptive::pre_exchange()
{
  if (!pending) return;
  next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
  pending = 0;
}

/* ---------------------------------------------------------------------- */
/*
void FixAdaptive::post_integrate()
{
  system_changed = 0;

  // return if not an adaptive timestep

  if (update->ntimestep < next_reneighbor) return;

  // compute geometric properties for all elements

  element->evec->compute_element_geometry(element->nlocal + element->nghost);

  // invoke stress and structure compute peratom

  if (structure->invoked_peratom != update->ntimestep) {
    structure->compute_peratom();
    structure->preneighflag = 1;
  }

  // acquire ghost structure id 

  comm->forward_comm_compute(structure);

  // build neighbor list for elements

  neighbor->build_one(list, 1);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  int i, j, k, l, n, ii, jj, kk, ll, jindex, jbasis, jucell, jetype, japc;
  int jnum, knum, *jlist, *jindexlist, *klist, *kindexlist;
  int break_flag;
  double delx, dely, delz;
  double **ax = atom->x;
  double **ex = element->x;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  int *apc = element->apc;
  double *apattern = structure->vector_atom;
  double ***vapattern = structure->vector_vatom;
  int jouter;
  int *nsurface = element->nsurface;
  int *element_shape_ids = element->element_shape_ids;
  int ***surface_ucell = element->surface_ucell;
  int **nsurface_ucell = element->nsurface_ucell;
  int cut[MAXCUT], cutdim;
  int nsplit = 0;
  int disc_flag;
  int natom;

  double P[3];                    // Point on splip plane (calculated as centroid)
  double normal[3];
  double A[3][3];                 // 3x3 matrix for SVD, will be V after SVD
  double eigen_values[3];


  // loop through neighbors of elements to find dislocation nearby
  // add element to splitlist

  for (ii = 0; ii < inum; ii++) {
    disc_flag = 0;
    i = ilist[ii];

    if (!(emask[i] & groupbit)) continue;

    if (element_shape_ids[etype[i]] != Element::HEXAHEDRON) 
      continue;

    // loop over atom neighbors

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];

    natom = 0;

    if (dislocation_flag) {
      P[0] = P[1] = P[2] = 0.0;
      A[0][0] = A[0][1] = A[0][2] = 0.0;
      A[1][0] = A[1][1] = A[1][2] = 0.0;
      A[2][0] = A[2][1] = A[2][2] = 0.0;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      if (jindex < 0) {
        if (dislocation_flag) {
          P[0] += ax[j][0];
          P[1] += ax[j][1];
          P[2] += ax[j][2];
        }
      } else {
        if (n == 0) natom = jj;
        break;
      }
    }

    if (natom >= thres) {
      disc_flag = 1;
      if (dislocation_flag) {
        P[0] /= natom;
        P[1] /= natom;
        P[2] /= natom;

        for (jj = 0; jj < natom; jj++) {
          j = jlist[jj];
          delx = ax[j][0] - P[0];
          dely = ax[j][1] - P[1];
          delz = ax[j][2] - P[2];
          A[0][0] += delx * delx;
          A[1][1] += dely * dely;
          A[2][2] += delz * delz;
          A[0][1] += delx * dely;
          A[0][2] += delx * delz;
          A[1][2] += dely * delz;
        }

        A[1][0] = A[0][1];
        A[2][0] = A[0][2];
        A[2][1] = A[1][2];

        eigen_decomposition(A, eigen_values);

        // the smallest eigen_values must be much smaller than the next

        if (eigen_values[0]/eigen_values[1] < eig_frac_thres) {
          normal[0] = A[0][0];
          normal[1] = A[1][0];
          normal[2] = A[2][0];
          int ncut = element->evec->check_split_element(i, normal, P, tol, split_width, cut, cutdim);

          if (ncut) disc_flag = 0;
          for (int icut = 0; icut < ncut; icut++) {
            if (nsplit == maxsplit) grow_splitlist();
            splitlist[nsplit][0] = i;
            splitlist[nsplit][1] = cut[icut];
            splitlist[nsplit][2] = cutdim;
            if (debug) {
              write_split_atom(i, jnum, jlist);
              printf("Split #%d me = %d element %d from atomistic region P = %g %g %g n = %g %g %g\n", 
                  nsplits_total, comm->me, element->tag[i], P[0], P[1], P[2], normal[0], normal[1], normal[2]);
            }
            nsplits_total++;
            nsplit++;
          }
        }
      }
    }
    if (disc_flag && discretize_flag) {
      if (nsplit == maxsplit) grow_splitlist();
      splitlist[nsplit][0] = i;
      splitlist[nsplit][1] = 0;
      splitlist[nsplit][2] = -1;
      if (debug) {
        write_split_atom(i, jnum, jlist);
        printf("Discretize #%d me = %d element %d\n", 
            nsplits_total, comm->me, element->tag[i]);
      }
      nsplit++;
      nsplits_total++;
      continue;
    }

    if (!dislocation_flag) continue;

    // continue checking for split from dislocations in CG domain

    for (jj = natom; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      jetype = etype[j];
      japc = apc[jetype];
      if (element_shape_ids[jetype] == Element::HEXAHEDRON) {
        for (kk = 0; kk < nsurface[element_shape_ids[jetype]]; kk++) {
          knum = nsurface_ucell[jetype][kk];
          klist = surface_ucell[jetype][kk];
          break_flag = 0;
          for (k = 0; k < knum; k++) {
            jucell = klist[k];
            for (jbasis = 0; jbasis < japc; jbasis++) {
              if (vapattern[j][jbasis][jucell] == perfect_crystal) continue;
              int ncut = element->evec->check_split_element(i, j, kk, tol, cut, cutdim, MAXCUT);
              for (int icut = 0; icut < ncut; icut++) {
                if (nsplit == maxsplit) grow_splitlist();
                splitlist[nsplit][0] = i;
                splitlist[nsplit][1] = cut[icut];
                splitlist[nsplit][2] = cutdim;
                if (debug) {
                  write_split_element(i, j);
                  printf("Split #%d me = %d element %d from element %d surface %d\n"
                      , nsplits_total, comm->me, element->tag[i], element->tag[j], kk);
                }
                nsplits_total++;
                nsplit++;
              }

              break_flag = 1;
              break;
            }
            if (break_flag) break;
          }
        }
      } else if (element_shape_ids[jetype] == Element::QUADRILATERAL) {
        if (domain->dimension == 3) {
          knum = nsurface_ucell[jetype][0];
          klist = surface_ucell[jetype][0];
          break_flag = 0;
          for (k = 0; k < knum; k++) {
            jucell = klist[k];

            for (jbasis = 0; jbasis < japc; jbasis++) {
              if (vapattern[j][jbasis][jucell] == perfect_crystal) continue;
              int ncut = element->evec->check_split_element(i, j, 0, tol, cut, cutdim, MAXCUT);
              if (ncut) { 
                if (nsplit + ncut > maxsplit) grow_splitlist();
                for (int icut = 0; icut < ncut; icut++) {
                  splitlist[nsplit][0] = i;
                  splitlist[nsplit][1] = cut[icut];
                  splitlist[nsplit][2] = cutdim;
                  nsplits_total++;
                  nsplit++;
                  if (debug) {
                    write_split_element(i, j);
                    printf("Split #%d me = %d element %d from element %d surface %d\n"
                        , nsplits_total, comm->me, element->tag[i], element->tag[j], kk);
                  }

                }
              }
              break_flag = 1;
              break;
            }
            if (break_flag) break;
          }
        } else {
          // 2D Quad element here
        }
      } else {
        // add new element types here
      }
    }
  }

  int nsplit_all;
  MPI_Allreduce(&nsplit, &nsplit_all, 1, MPI_INT, MPI_MAX, world);

  // return and update next reneighbor if no plane need to split
  // this will not force reneighbor at this time step

  if (nsplit_all == 0) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen, "No new dislocation found from fix adaptive timestep = %d\n", update->ntimestep);
    return;
  } else 
    system_changed = 1;

  pending = 1;
  int nelocal = element->nlocal;
  int nalocal = atom->nlocal;

  // split elements

  int nesplit, nesplit_all;
  nesplit = element->evec->split_element(splitlist, nsplit);

  // extend tags to new atoms and elements

  atom->tag_extend();
  element->tag_extend();

  int nnew_atoms = atom->nlocal - nalocal;
  int nnew_elems = element->nlocal - nelocal;
  int nnew_atoms_all, nnew_elems_all;
  MPI_Allreduce(&nnew_atoms, &nnew_atoms_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nnew_elems, &nnew_elems_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nesplit, &nesplit_all, 1, MPI_INT, MPI_SUM, world);

  // update number of elements and number of atoms

  atom->natoms += nnew_atoms_all;
  element->nelements += nnew_elems_all;

  // recount # of nodes, node cells, vatoms, ucells

  element->count_nodes(1);
  element->count_node_cells(1);
  element->count_vatoms();
  element->count_ucells();

  // migrates new elements if there are any

  if (nnew_elems_all) {
    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal, atom->x);
      domain->x2lamda(element->nlocal, element->x);
      domain->nodex2lamda(element->nlocal, element->nodex);
    }
    domain->reset_box();

    IrregularComm *irrcomm = new IrregularComm(cac);
    irrcomm->migrate(1);
    delete irrcomm;

    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal, atom->x);
      domain->lamda2x(element->nlocal, element->x);
      domain->lamda2nodex(element->nlocal, element->nodex);
    }
  }

  if (comm->me == 0) {
    if (screen) fprintf(screen, "Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n", nesplit_all, nnew_atoms_all, nnew_elems_all+nesplit_all);
    if (logfile) fprintf(logfile, "Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n", nesplit_all, nnew_atoms_all, nnew_elems_all+nesplit_all);
  }

}
*/
/* ---------------------------------------------------------------------- */

void FixAdaptive::post_integrate()
{
  system_changed = 0;

  // return if not an adaptive timestep

  if (update->ntimestep < next_reneighbor) return;

  // compute geometric properties for all elements

  element->evec->compute_element_geometry(element->nlocal + element->nghost);

  // invoke stress and structure compute peratom

  if (structure->invoked_peratom != update->ntimestep) {
    structure->compute_peratom();
    structure->preneighflag = 1;
  }

  // acquire ghost structure id 

  comm->forward_comm_compute(structure);

  // build neighbor list for elements

  neighbor->build_one(list, 1);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  int i, j, k, l, n, ii, jj, kk, ll, jindex, jbasis, jucell, jetype, japc;
  int jnum, knum, *jlist, *jindexlist, *klist;;
  int break_flag;
  double delx, dely, delz;
  double **ax = atom->x;
  double **ex = element->x;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  int *apc = element->apc;
  double *apattern = structure->vector_atom;
  double ***vapattern = structure->vector_vatom;
  int *nsurface = element->nsurface;
  int *element_shape_ids = element->element_shape_ids;
  int ***surface_ucell = element->surface_ucell;
  int **nsurface_ucell = element->nsurface_ucell;
  int *cut = new int[element->max_edge_ucell];
  int nsplit = 0;
  int disc_flag;
  int natom;
  int cutdim, dim, dir;

  double P[3];                    // Point on splip plane (calculated as centroid)
  double normal[3];
  double A[3][3];                 // 3x3 matrix for SVD, will be V after SVD
  double eigen_values[3];

  double products[27][6];         // tally of x_i * x_j products for each neighboring region
  double sums[27][3];             // tally of x_i for each neighboring region
  double num[27];                 // number of atoms in each neighboring region

  // loop through neighbors of elements to find dislocation nearby
  // add element to splitlist
  
  for (ii = 0; ii < inum; ii++) {
    disc_flag = 0;
    i = ilist[ii];

    if (!(emask[i] & groupbit)) continue;

    if (element_shape_ids[etype[i]] != Element::HEXAHEDRON) 
      continue;

    // loop over atom neighbors

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];

    natom = 0;

    // initiallize arrays

    if (dislocation_flag) {
      for (j = 0; j < 27; j++) {
        sums[j][0] = sums[j][1] = sums[j][2] = 0.0;
        products[j][0] = products[j][1] = products[j][2] = 0.0;
        products[j][3] = products[j][4] = products[j][5] = 0.0;
        num[j] = 0;
      }
    } 

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];

      if (jindex >= 0) {
        if (dislocation_flag) {
          sums[jindex][0] += ax[j][0];
          sums[jindex][1] += ax[j][1];
          sums[jindex][2] += ax[j][2];
          products[jindex][0] += ax[j][0] * ax[j][0];
          products[jindex][1] += ax[j][1] * ax[j][1];
          products[jindex][2] += ax[j][2] * ax[j][2];
          products[jindex][3] += ax[j][0] * ax[j][1];
          products[jindex][4] += ax[j][0] * ax[j][2];
          products[jindex][5] += ax[j][1] * ax[j][2];
          num[jindex]++;
        }
        natom++;
      } else break;
    }
    //if (element->tag[i] == 17260) {
    //  printf("natom = %d me = %d\n",natom,comm->me);
    //  write_neighbor(i);
    //  error->one(FLERR,"TEST");
    //}

    if (natom >= thres) {
      disc_flag = 1;
      if (dislocation_flag) {

        // check dislocations coming from 6 directions

        for (dim = 0; dim < 3; dim++) {
          P[0] = P[1] = P[2] = 0.0;
          A[0][0] = A[0][1] = A[0][2] = 0.0;
          A[1][0] = A[1][1] = A[1][2] = 0.0;
          A[2][0] = A[2][1] = A[2][2] = 0.0;
          n = 0;          

          // compute A matrix and P 
          for (j = 0; j < 8; j++) {
            k = boxindices[dim][j];
            P[0] += sums[k][0];
            P[1] += sums[k][1];
            P[2] += sums[k][2];
            A[0][0] += products[k][0];
            A[1][1] += products[k][1];
            A[2][2] += products[k][2];
            A[0][1] += products[k][3];
            A[0][2] += products[k][4];
            A[1][2] += products[k][5];
            n += num[k];
          }
          if (n == 0) continue;
          P[0] /= n;
          P[1] /= n;
          P[2] /= n;
          A[0][0] -= P[0] * P[0] * n;
          A[1][1] -= P[1] * P[1] * n;
          A[2][2] -= P[2] * P[2] * n;
          A[0][1] -= P[0] * P[1] * n;
          A[0][2] -= P[0] * P[2] * n;
          A[1][2] -= P[1] * P[2] * n;

          // A is symmetric

          A[1][0] = A[0][1];
          A[2][0] = A[0][2];
          A[2][1] = A[1][2];

          eigen_decomposition(A, eigen_values);

          // the smallest eigen_values must be much smaller than the next

          if (eigen_values[0] / eigen_values[1] < eig_frac_thres) {
            normal[0] = A[0][0];
            normal[1] = A[1][0];
            normal[2] = A[2][0];

            int ncut = check_split_element(i, normal, tol, cut, dim, natom, jlist, jindexlist);

            if (ncut) disc_flag = 0;
            else continue;
            for (int icut = 0; icut < ncut; icut++) {
              if (nsplit == maxsplit) grow_splitlist();
              splitlist[nsplit][0] = i;
              splitlist[nsplit][1] = cut[icut];
              splitlist[nsplit][2] = dim;
              if (debug) {
                write_split_atom(i, jnum, jlist);
                printf("Split #%d me = %d element %d from atomistic region P = %g %g %g n = %g %g %g\n", 
                    nsplits_total, comm->me, element->tag[i], P[0], P[1], P[2], normal[0], normal[1], normal[2]);
              }
              nsplits_total++;
              nsplit++;
            }
            break;
          }
        }
      }
    }
    if (disc_flag && discretize_flag) {
      if (nsplit == maxsplit) grow_splitlist();
      splitlist[nsplit][0] = i;
      splitlist[nsplit][1] = 0;
      splitlist[nsplit][2] = -1;
      if (debug) {
        write_split_atom(i, jnum, jlist);
        printf("Discretize #%d me = %d element %d\n", 
            nsplits_total, comm->me, element->tag[i]);
      }
      nsplit++;
      nsplits_total++;
      continue;
    }

    if (!dislocation_flag) continue;

    // continue checking for split from dislocations in CG domain

    for (jj = natom; jj < jnum; jj++) {
      j = jlist[jj];
      jetype = etype[j];
      japc = apc[jetype];
      if (element_shape_ids[jetype] == Element::HEXAHEDRON) {
        for (kk = 0; kk < nsurface[element_shape_ids[jetype]]; kk++) {
          knum = nsurface_ucell[jetype][kk];
          klist = surface_ucell[jetype][kk];
          break_flag = 0;
          for (k = 0; k < knum; k++) {
            jucell = klist[k];
            for (jbasis = 0; jbasis < japc; jbasis++) {
              if (vapattern[j][jbasis][jucell] == perfect_crystal) continue;
              int ncut = check_split_element(i, j, kk, tol, cut, cutdim);
              for (int icut = 0; icut < ncut; icut++) {
                if (nsplit == maxsplit) grow_splitlist();
                splitlist[nsplit][0] = i;
                splitlist[nsplit][1] = cut[icut];
                splitlist[nsplit][2] = cutdim;
                if (debug) {
                  write_split_element(i, j);
                  printf("Split #%d me = %d element %d from element %d surface %d\n"
                      , nsplits_total, comm->me, element->tag[i], element->tag[j], kk);
                }
                nsplits_total++;
                nsplit++;
              }

              break_flag = 1;
              break;
            }
            if (break_flag) break;
          }
        }
      } else if (element_shape_ids[jetype] == Element::QUADRILATERAL) {
        if (domain->dimension == 3) {
          knum = nsurface_ucell[jetype][0];
          klist = surface_ucell[jetype][0];
          break_flag = 0;
          for (k = 0; k < knum; k++) {
            jucell = klist[k];

            for (jbasis = 0; jbasis < japc; jbasis++) {
              if (vapattern[j][jbasis][jucell] == perfect_crystal) continue;
              int ncut = check_split_element(i, j, 0, tol, cut, cutdim);
              if (ncut) { 
                if (nsplit + ncut > maxsplit) grow_splitlist();
                for (int icut = 0; icut < ncut; icut++) {
                  splitlist[nsplit][0] = i;
                  splitlist[nsplit][1] = cut[icut];
                  splitlist[nsplit][2] = cutdim;
                  nsplits_total++;
                  nsplit++;
                  if (debug) {
                    write_split_element(i, j);
                    printf("Split #%d me = %d element %d from element %d surface %d\n"
                        , nsplits_total, comm->me, element->tag[i], element->tag[j], kk);
                  }

                }
              }
              break_flag = 1;
              break;
            }
            if (break_flag) break;
          }
        } else {
          // 2D Quad element here
        }
      } else {
        // add new element types here
      }
    }
  }
  //error->all(FLERR,"TEST"); 
  int nsplit_all;
  MPI_Allreduce(&nsplit, &nsplit_all, 1, MPI_INT, MPI_MAX, world);

  // return and update next reneighbor if no plane need to split
  // this will not force reneighbor at this time step

  if (nsplit_all == 0) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen, "No new dislocation found from fix adaptive timestep = " BIGINT_FORMAT "\n", update->ntimestep);
    return;
  } else 
    system_changed = 1;

  pending = 1;
  int nelocal = element->nlocal;
  int nalocal = atom->nlocal;

  // split elements

  int nesplit, nesplit_all;
  nesplit = split_element(splitlist, nsplit);

  // extend tags to new atoms and elements

  atom->tag_extend();
  element->tag_extend();

  int nnew_atoms = atom->nlocal - nalocal;
  int nnew_elems = element->nlocal - nelocal;
  int nnew_atoms_all, nnew_elems_all;
  MPI_Allreduce(&nnew_atoms, &nnew_atoms_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nnew_elems, &nnew_elems_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nesplit, &nesplit_all, 1, MPI_INT, MPI_SUM, world);

  // update number of elements and number of atoms

  atom->natoms += nnew_atoms_all;
  element->nelements += nnew_elems_all;

  // recount # of nodes, node cells, vatoms, ucells

  element->count_nodes(1);
  element->count_node_cells(1);
  element->count_vatoms();
  element->count_ucells();

  // migrates new elements if there are any

  if (nnew_elems_all) {
    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal, atom->x);
      domain->x2lamda(element->nlocal, element->x);
      domain->nodex2lamda(element->nlocal, element->nodex);
    }
    domain->reset_box();

    IrregularComm *irrcomm = new IrregularComm(cac);
    irrcomm->migrate(1);
    delete irrcomm;

    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal, atom->x);
      domain->lamda2x(element->nlocal, element->x);
      domain->lamda2nodex(element->nlocal, element->nodex);
    }
  }
  if (comm->me == 0) {
    if (screen) fprintf(screen, "Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n", nesplit_all, nnew_atoms_all, nnew_elems_all+nesplit_all);
    if (logfile) fprintf(logfile, "Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n", nesplit_all, nnew_atoms_all, nnew_elems_all+nesplit_all);
  }
  delete [] cut;
}

/*  ----------------------------------------------------------------------
    check if a plane with normal vector <n1, n2, n3> that
    go through point P intersects with element I
    -------------------------------------------------------------------------  */

int FixAdaptive::check_split_element(int i, double n[3], double tol, 
    int *cut, int dim, int jnum, int *jlist, int *jindexlist)
{
  int ietype = element->etype[i];

  // check if cut plane is parallel to surface plane of i
  // if not, return 0

  double *isurface = element->surface_plane[i][dim];

  //printf("n = %g %g %g isurface = %g %g %g dot3 = %g\n",n[0],n[1],n[2],isurface[0],isurface[1],isurface[2],dot3(isurface,n));
  //return 0;
  if (1.0 - fabs(dot3(isurface, n)) > tol) return 0;

  int ncell = element->evec->ncells[ietype][dim];
  if (ncell <= 2) return 0;
  double a = 2.0 / ncell;

  // convert neighbor coords to xi and group into planes

  int nplanes = 0;
  int j, jj, jindex, newbin_flag;
  double plane_width = a * plane_width_fraction;
  double tmp[3], xi;
  double **ax = atom->x;
  for (jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    jindex = jindexlist[jj];
    for (int ibox = 0; ibox < 8; ibox ++)
      if (jindex == boxindices[dim][ibox]) {
        element->box2natural(ax[j], tmp, i);
        xi = tmp[dim];

        // check if this atom is in any planes 
        // if yes, update that plane center position
        // otherwise, use it to start a new plane

        newbin_flag = 1;
        for (int iplane = 0; iplane < nplanes; iplane++) {
          if (fabs(xi - planexi[iplane]) < plane_width) {
            if (newbin_flag == 0) {
              printf("i = %d j = %d xi = %g planexi[%d] = %g plane_width = %g n = %g %g %g dim = %d\n isurface[0] = %g %g %g\n isurface[1] = %g %g %g\n isurface[2] = %g %g %g\n"
                  ,element->tag[i],atom->tag[j],xi,iplane,planexi[iplane],plane_width,n[0],n[1],n[2],dim,
                  element->surface_plane[i][0][0],
                  element->surface_plane[i][0][1],
                  element->surface_plane[i][0][2],
                  element->surface_plane[i][1][0],
                  element->surface_plane[i][1][1],
                  element->surface_plane[i][1][2],
                  element->surface_plane[i][2][0],
                  element->surface_plane[i][2][1],
                  element->surface_plane[i][2][2]
                  );

              error->warning(FLERR,"plane_width too large, overlapping planes exist");
              return 0;
            }
            newbin_flag = 0;
            //if (element->tag[i] == 832)
            //  printf("j = %d xi = %g planexi[%d] = %g dim = %d\n",
            //      atom->tag[j],xi,iplane,planexi[iplane],dim);
            planexi[iplane] = planexi[iplane] * planecount[iplane] + xi;
            planexi[iplane] /= ++planecount[iplane];
          }
        }

        if (newbin_flag) {
          if (nplanes == maxplanes) {
            maxplanes *= 2;
            memory->grow(planexi, maxplanes, "fix:planexi");
            memory->grow(planecount, maxplanes, "fix:planecount");
          }
          planexi[nplanes] = xi;
          planecount[nplanes] = 1;
          nplanes++;
        }
        break;  
      }
  }

  int *count = new int[ncell];
  int *cutflag = new int[ncell];
  memset(count, 0, ncell * sizeof(int));
  memset(cutflag, 0, ncell * sizeof(int));

  // convert plane coord to unit cell index  
  // and count how many planes are in each unit cell

  for (int iplane = 0; iplane < nplanes; iplane++) {
    planexi[iplane] = floor((planexi[iplane] + 1) / a);
    if (planexi[iplane] >= ncell || planexi[iplane] < 0)
      error->one(FLERR,"TEST");
    count[(int) planexi[iplane]]++;
  }

  for (int icell = 0; icell < ncell; icell++) {

    // if more than 1 plane in a unit cell, split above and below that unit cell
    if (count[icell] > 1) {
      if (icell != ncell - 1)
        cutflag[icell + 1] = 1;
      if (icell != 0)
        cutflag[icell] = 1;
    } 

    // if only 1 plane in a unit cell, split above and/or below that unit cell 
    // if only the neighbor cell also contains 1 plane
    else if (count[icell] == 1) {
      if (icell != ncell - 1)
        if (count[icell + 1] == 1)
          cutflag[icell + 1] = 1;
      if (icell >= 1)
        if (count[icell - 1] == 1)
          cutflag[icell] = 1;
    }
  }

  int ncut = 0;
  for (int icut = 1; icut < ncell; icut++) 
    if (cutflag[icut])
      cut[ncut++] = icut;

  delete [] count;
  delete [] cutflag;

  return ncut;
}

/* ----------------------------------------------------------------------
   check if element J's surface K cut through element I
   -------------------------------------------------------------------------  */

int FixAdaptive::check_split_element(int i, int j, int surface, double tol, int *cut, int &cutdim)
{
  int jetype = element->etype[j];
  int j_shape_id = element->element_shape_ids[jetype];
  // check if cut surface of J parallel to any surface plane of I
  // if not, return 0

  int jcutdim, dir;
  double **isurface = element->surface_plane[i];
  double **jsurface = element->surface_plane[j];

  if (j_shape_id == Element::HEXAHEDRON) {
    jcutdim = surface / 2;
    dir = surface % 2;
  } else if (j_shape_id == Element::QUADRILATERAL) {
    if (domain->dimension == 3) {
      jcutdim = 2;
      dir = 0;
    } else {
      jcutdim = surface / 2;
      dir = surface % 2;
    }
  }

  // check for surface of I that is parallel to the cut surface from J
  // skip if none are parallel (within the tolerant)

  cutdim = -1;
  for (int dim = 0; dim < 3; dim++) {
    double tmp = 1.0 - fabs(dot3(isurface[dim], jsurface[jcutdim]));
    if (tmp < tol) {
      cutdim = dim;
      break;
    }
  }

  if (cutdim < 0) return 0;

  // look for node (center of the unit cell) on the cut surface of  element J 
  // closest to center element I and use it as the reference point

  int jnode, jj, jbasis;
  double delx, dely, delz, rsq;
  double xtmp, ytmp, ztmp;
  double rsqmin = BIG;
  double *ix = element->x[i];
  double *jx = element->x[j];
  double ***jnodex = element->nodex[j];
  double closest_coord[3];
  int ietype = element->etype[i];
  int ncell = element->evec->ncells[ietype][cutdim];
  int japc = element->apc[jetype];
  double cell_size = element->cell_size[i][cutdim];

  if (j_shape_id == Element::HEXAHEDRON) {
    for (jj = 0; jj < 4; jj++) {
      jnode = element->node_set_3D[jcutdim][dir][jj];
      xtmp = ytmp = ztmp = 0.0;
      for (jbasis = 0; jbasis < japc; jbasis++) {
        xtmp += jnodex[jbasis][jnode][0];
        ytmp += jnodex[jbasis][jnode][1];
        ztmp += jnodex[jbasis][jnode][2];
      }
      xtmp /= japc;
      ytmp /= japc;
      ztmp /= japc;
      delx = xtmp - ix[0];
      dely = ytmp - ix[1];
      delz = ztmp - ix[2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < rsqmin) { 
        rsqmin = rsq;
        closest_coord[0] = xtmp;
        closest_coord[1] = ytmp;
        closest_coord[2] = ztmp;
      }
    }

    // project closest node coord of element J onto slip plane normal line of element I
    // if projection falls inside element --> need splitting

    double a = cell_size / (ncell - 1.0);
    double dist = (dot3(isurface[cutdim], closest_coord) -
        dot3(isurface[cutdim], ix));
    dist /= a;
    double orient = dot3(isurface[cutdim], closest_coord) -
      dot3(isurface[cutdim], jx);
    double hi = ncell / 2.0;
    double lo = -hi;

    if (orient > 0)
      hi -= 1.0;
    else
      lo += 1.0;


    if (dist < lo || dist > hi) return 0;

    cut[0] = static_cast<int> (dist - lo);
    cut[0]++;

    if (cut[0] > 0 && cut[0] < ncell) 
      return 1;

  } else if (j_shape_id == Element::QUADRILATERAL) {
    if (domain->dimension == 3) {
      for (jnode = 0; jnode < 4; jnode++) {
        xtmp = ytmp = ztmp = 0.0;
        for (jbasis = 0; jbasis < japc; jbasis++) {
          xtmp += jnodex[jbasis][jnode][0];
          ytmp += jnodex[jbasis][jnode][1];
          ztmp += jnodex[jbasis][jnode][2];
        }
        xtmp /= japc;
        ytmp /= japc;
        ztmp /= japc;
        delx = xtmp - ix[0];
        dely = ytmp - ix[1];
        delz = ztmp - ix[2];
        if (rsq < rsqmin) { 
          rsqmin = rsq;
          closest_coord[0] = xtmp;
          closest_coord[1] = ytmp;
          closest_coord[2] = ztmp;
        }
      }
      double a = cell_size / (ncell - 1);
      double dist = (dot3(isurface[cutdim], closest_coord) -
          dot3(isurface[cutdim], ix));
      dist /= a;
      double hi = ncell / 2.0;
      double lo = -hi;

      if (dist < lo || dist > hi) return 0;

      cut[0] = static_cast<int> (dist - lo);

      if (cut[0] == 0) {
        cut[0] = 1;
        return 1;
      } else if (cut[0] == ncell - 1) {
        return 1;
      } else if (cut[0] > 0 && cut[0] < ncell - 1) {
        cut[1] = cut[0] + 1;
        return 2;
      }
    } else {
      // Add 2d here
    }

  } else {
    // add new element types here
  }
  return 0;
}

/*  ----------------------------------------------------------------------
    called from fix_adaptive to split elements in list
    return number of elements being splitted
    -------------------------------------------------------------------------  */

int FixAdaptive::split_element(int **splitlist, int n)
{
  int ietype, iapc;
  int *ictype = new int[element->maxapc];
  int type_info[8];
  int *nelem, **cellsize;
  int ii = 0;
  int i, icut, icutdim;
  int *isplitlist, **nelem_array, ***cellsize_array;
  int *idisclist,*del_flag_list;
  ElementVec *evec = element->evec;
  int **ncells = evec->ncells;
  int **ngcells = evec->ngcells;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int *apc = element->apc;

  // counter for number of elements being splitted

  int nsplits = 0;

  // counter for number of elements being discretized

  int ndiscs = 0;

  // first pass: check for new etypes
  // loop through all splits in list
  // each element might have several splits and
  // they are all next to each other in splitlist
  // if icutdim = -1, discretize element i to atoms

  int max_del_flag = element->nlocal;
  if (n > 0) {
    memory->create(isplitlist, n, "evec:isplitlist");
    memory->create(idisclist, n, "evec:idisclist");
    memory->create(nelem_array, n, 3, "evec:nelem_array");
    memory->create(cellsize_array, n, 3, MAXSUB, "evec:cellsize_array");

    memory->create(del_flag_list, max_del_flag, "evec:max_del_flag");
    for (i = 0; i < max_del_flag; i++)
      del_flag_list[i] = 0;

    while (ii < n) {
      nelem = nelem_array[nsplits];
      cellsize = cellsize_array[nsplits];
      nelem[0] = nelem[1] = nelem[2] = 1;

      i = splitlist[ii][0];

      if (i >= max_del_flag) error->all(FLERR, "TEST"); 
      del_flag_list[i] = 1;
      if (splitlist[ii][2] == -1) {
        idisclist[ndiscs++] = i;
        ii++;
        continue;
      }
      ietype = etype[i];
      ietype = 1;
      type_info[6] = element->element_shape_ids[ietype];
      type_info[7] = element->apc[ietype];

      // apply first split to element I

      icut = splitlist[ii][1];
      icutdim = splitlist[ii][2];

      isplitlist[nsplits] = i;
      ietype = etype[i];
      cellsize[0][0] = ncells[ietype][0];
      cellsize[1][0] = ncells[ietype][1];
      cellsize[2][0] = ncells[ietype][2];
      cellsize[icutdim][0] = icut;
      cellsize[icutdim][1] = ncells[ietype][icutdim]-icut;
      nelem[icutdim]++;

      // apply each additional split to element cellsize list 
      // if next split still cuts element I

      if (ii < n-1) {
        while (i == splitlist[ii+1][0]) {

          icut = splitlist[ii+1][1];
          icutdim = splitlist[ii+1][2];

          if (nelem[icutdim] == MAXSUB) 
            error->one(FLERR, "Too many split in one direction");
          int sum = 0;

          for (int l = 0; l < nelem[icutdim]; l++) {
            if (icut > sum + cellsize[icutdim][l])
              sum += cellsize[icutdim][l];
            else if (icut < sum + cellsize[icutdim][l]) {
              for (int m = nelem[icutdim]-1; m >= l; m--) 
                cellsize[icutdim][m+1] = cellsize[icutdim][m];
              cellsize[icutdim][l] = icut-sum;
              cellsize[icutdim][l+1] -= icut-sum;
              if (cellsize[icutdim][l] < 1 || cellsize[icutdim][l] > ncells[ietype][icutdim]
                  || cellsize[icutdim][l+1] < 1 || cellsize[icutdim][l+1] > ncells[ietype][icutdim])
              {
                printf("icut = %d sum = %d cellsize[%d] = %d cellsize[%d] = %d icutdim = %d\n", icut, sum, l, cellsize[icutdim][l], l+1, cellsize[icutdim][l+1], icutdim);
                error->one(FLERR, "TEST");
              }
              nelem[icutdim]++;
              break;
            } else break;
          }
          ii++;
          if (ii == n - 1) break;
        }
      }

      // check rearrange splits and split any elements of size 2

      int index;
      for (int idim = 0; idim < domain->dimension; idim++) {
        index = 0; 
        while (index < nelem[idim]) {
          //printf("cellsize[%d][%d] = %d\n", idim, index, cellsize[idim][index]);
          if (cellsize[idim][index] == 2) {
            if (nelem[icutdim] == MAXSUB) 
              error->one(FLERR, "Too many split in one direction");
            for (int m = nelem[idim]; m > index+1; m--) {
              cellsize[idim][m] = cellsize[idim][m-1];
            }
            cellsize[idim][index] = 1;
            cellsize[idim][index+1] = 1;
            nelem[idim]++;
          } 
          index++;
        }
      }     

      // make sure the splits tally up to the original element sizes

      int checksum = 0;
      for (int ix = 0; ix < nelem[0]; ix++) 
        checksum += cellsize[0][ix];
      if (checksum > ncells[ietype][0]) {
        printf("Ncells = %d cellsize x =", ncells[ietype][0]);
        for (int ix = 0; ix < nelem[0]; ix++) 
          printf(" %d", cellsize[0][ix]);
        printf("\n");
        error->one(FLERR, "Incorrect split");
      }

      checksum = 0;
      for (int iy = 0; iy < nelem[1]; iy++) 
        checksum += cellsize[1][iy];
      if (checksum > ncells[ietype][1]) {
        printf("Ncells = %d cellsize y =", ncells[ietype][1]);
        for (int iy = 0; iy < nelem[1]; iy++) 
          printf(" %d", cellsize[1][iy]);
        printf("\n");
        error->one(FLERR, "Incorrect split");
      }

      checksum = 0;
      for (int iz = 0; iz < nelem[2]; iz++) 
        checksum += cellsize[2][iz];
      if (checksum > ncells[ietype][2]) {
        printf("Ncells = %d cellsize z =", ncells[ietype][2]);
        for (int iz = 0; iz < nelem[2]; iz++) 
          printf(" %d", cellsize[2][iz]);
        printf("\n");
        error->one(FLERR, "Incorrect split");
      }

      for (int ix = 0; ix < nelem[0]; ix++) {
        for (int iy = 0; iy < nelem[1]; iy++) {
          for (int iz = 0; iz < nelem[2]; iz++) {
            type_info[0] = cellsize[0][ix];
            type_info[3] = MIN(type_info[0], ngcells[ietype][0]);
            type_info[1] = cellsize[1][iy];
            type_info[4] = MIN(type_info[1], ngcells[ietype][1]);
            type_info[2] = cellsize[2][iz];
            type_info[5] = MIN(type_info[2], ngcells[ietype][2]);
            if (type_info[0] <= 0 || type_info[0] > ncells[ietype][0]) {
              printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
                  , type_info[0]
                  , type_info[1]
                  , type_info[2]);
              error->one(FLERR, "TEST");
            }
            if (type_info[1] <= 0 || type_info[1] > ncells[ietype][1]) {
              printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
                  , type_info[0]
                  , type_info[1]
                  , type_info[2]);
              error->one(FLERR, "TEST");
            }
            if (type_info[2] <= 0 || type_info[2] > ncells[ietype][2]) {
              printf("Error i = %d me = %d type_info %d %d %d\n", element->tag[i], comm->me
                  , type_info[0]
                  , type_info[1]
                  , type_info[2]);
              error->one(FLERR, "TEST");
            }

            if (domain->dimension == 3) {

              // skip if 2 or more directions have size 1
              // will split to atoms

              if (((type_info[0] == 1) + (type_info[1] == 1) + (type_info[2] == 1)) >= 2) 
                continue;

              // if x direction size is 1, yz -> xy
              if (type_info[0] == 1) {
                type_info[0] = type_info[1];
                type_info[1] = type_info[2];
                type_info[2] = 1;
                type_info[3] = type_info[4];
                type_info[4] = type_info[5];
                type_info[5] = 1;
                type_info[6] = Element::QUADRILATERAL;
              }

              // if y direction size is 1, xz -> xy
              if (type_info[1] == 1) {
                type_info[1] = type_info[2];
                type_info[2] = 1;
                type_info[4] = type_info[5];
                type_info[5] = 1;
                type_info[6] = Element::QUADRILATERAL;
              }

              if (type_info[2] == 1) 
                type_info[6] = Element::QUADRILATERAL;

            }
            evec->request_new_etype(type_info);
            type_info[6] = Element::HEXAHEDRON;
          }
        }
      }
      nsplits++;
      ii++;
    }
  }

  evec->add_requested_etype(0);

  // redefine pointers in case they were resized
  
  ncells = evec->ncells;
  ngcells = evec->ngcells;
  etype = element->etype;
  ctype = element->ctype;
  apc = element->apc;

  if (n > 0) {
    // second pass: split elements

    int index_offset[3], index_node0, iucell;
    int myetype, ncellx, ncelly, ncellz;

    // iucell_list = list of ucell index in old element for each node of new element

    int *iucell_list = new int[8]; 
    int *nucell = element->nucell;

    for (ii = 0; ii < nsplits; ii++) {
      int nalocal_previous = atom->nlocal;
      int nelocal_previous = element->nlocal;
      i = isplitlist[ii];
      ietype = etype[i];
      ncellx = ncells[ietype][0];
      ncelly = ncells[ietype][1];
      ncellz = ncells[ietype][2];
      iapc = apc[ietype];
      for (int j = 0; j < iapc; j++)
        ictype[j] = ctype[i][j];
      nelem = nelem_array[ii];
      cellsize = cellsize_array[ii];
      index_offset[0] = index_offset[1] = index_offset[2] = 0;
      index_node0 = 0;
      for (int ix = 0; ix < nelem[0]; ix++) {
        index_offset[0] = ncelly * ncellz * (cellsize[0][ix] - 1);
        type_info[0] = cellsize[0][ix];
        type_info[3] = MIN(type_info[0], ngcells[ietype][0]);
        for (int iy = 0; iy < nelem[1]; iy++) {
          index_offset[1] = ncellz * (cellsize[1][iy] - 1);
          type_info[1] = cellsize[1][iy];
          type_info[4] = MIN(type_info[1], ngcells[ietype][1]);
          if (domain->dimension == 3) {
            for (int iz = 0; iz < nelem[2]; iz++) {
              index_offset[2] = cellsize[2][iz] - 1;
              type_info[2] = cellsize[2][iz];
              type_info[5] = MIN(type_info[2], ngcells[ietype][2]);
              if (((type_info[0] == 1) + (type_info[1] == 1) + (type_info[2] == 1)) > 1) {
                for (int iix = 0; iix < type_info[0]; iix++)
                  for (int iiy = 0; iiy < type_info[1]; iiy++)
                    for (int iiz = 0; iiz < type_info[2]; iiz++) {
                      iucell = index_node0 
                        + iix * ncelly * ncellz 
                        + iiy * ncellz 
                        + iiz;
                      if (iucell > nucell[ietype]) {
                        printf("Split to atom iucell = %d nucell = %d ietype = %d ncells = %d %d %d\n"
                            , iucell, nucell[ietype], ietype, ncells[ietype][0], ncells[ietype][1], ncells[ietype][2]);
                        error->one(FLERR, "Wrong ucell index");
                      }
                      for (int ibasis = 0; ibasis < iapc; ibasis++)
                        evec->create_pass_on_atom(i, ibasis, iucell, ictype[ibasis], 0);
                    }
              } else {
                if (type_info[0] == 1) {
                  type_info[0] = type_info[1];
                  type_info[1] = type_info[2];
                  type_info[2] = 1;
                  type_info[3] = type_info[4];
                  type_info[4] = type_info[5];
                  type_info[5] = 1;
                  type_info[6] = Element::QUADRILATERAL;
                  iucell_list[0] = index_node0; 
                  iucell_list[1] = index_node0 + index_offset[1];
                  iucell_list[2] = index_node0 + index_offset[1] 
                    + index_offset[2];
                  iucell_list[3] = index_node0 + index_offset[2];
                } else if (type_info[1] == 1) {
                  type_info[1] = type_info[2];
                  type_info[2] = 1;
                  type_info[4] = type_info[5];
                  type_info[5] = 1;
                  type_info[6] = Element::QUADRILATERAL;
                  iucell_list[0] = index_node0; 
                  iucell_list[1] = index_node0 + index_offset[0];
                  iucell_list[2] = index_node0 + index_offset[0] 
                    + index_offset[2];
                  iucell_list[3] = index_node0 + index_offset[2];
                } else if (type_info[2] == 1) {
                  type_info[6] = Element::QUADRILATERAL;
                  iucell_list[0] = index_node0; 
                  iucell_list[1] = index_node0 + index_offset[0];
                  iucell_list[2] = index_node0 + index_offset[0] 
                    + index_offset[1];
                  iucell_list[3] = index_node0 + index_offset[1];
                } else {

                  type_info[6] = Element::HEXAHEDRON;
                  // calculate ucell index in element I for each node of new element

                  iucell_list[0] = index_node0; 
                  iucell_list[1] = index_node0 + index_offset[0];
                  iucell_list[2] = index_node0 + index_offset[0] 
                    + index_offset[1];
                  iucell_list[3] = index_node0 + index_offset[1];
                  iucell_list[4] = iucell_list[0] + index_offset[2];
                  iucell_list[5] = iucell_list[1] + index_offset[2];
                  iucell_list[6] = iucell_list[2] + index_offset[2];
                  iucell_list[7] = iucell_list[3] + index_offset[2];

                }
                myetype = evec->find_etype(type_info);
                if (myetype == 0) {
                  printf("type_info = %d %d %d %d %d %d %d %d\n"
                      , type_info[0]
                      , type_info[1]
                      , type_info[2]
                      , type_info[3]
                      , type_info[4]
                      , type_info[5]
                      , type_info[6]
                      , type_info[7]);
                  error->one(FLERR, "Etype not defined yet");
                } 

                for (int aa = 0; aa < element->npe[myetype]; aa++)
                  if (iucell_list[aa] > nucell[ietype]) {
                    printf("Split to element iucell_list[%d] = %d nucell = %d ietype = %d ncells = %d %d %d\n"
                        , aa, iucell_list[aa], nucell[ietype], ietype, ncells[ietype][0], ncells[ietype][1], ncells[ietype][2]);
                    error->one(FLERR, "Wrong ucell_list");
                  }
                evec->create_pass_on_element(i, iucell_list, myetype, ictype, 0);

              } 
              index_node0 += index_offset[2]+1;
            }
            index_node0 += index_offset[1];
          } else {
            type_info[2] = type_info[5] = 1;
            if (type_info[0] == 1) {
              for (int iiy = 0; iiy < type_info[1]; iiy++)
                for (int ibasis = 0; ibasis < iapc; ibasis++)
                  evec->create_pass_on_atom(i, ibasis, iucell, ictype[ibasis], 0);
            } else if (type_info[1] == 1) {
              error->one(FLERR, "TEST"); 
              //            for (int iix = 0; iix < info[0]; iix++) 
              //              create_pass_on_atom(i, index_node0+iix * ncellyz, ictype, 0);
            } else {
              myetype = evec->find_etype(type_info);
              if (myetype == 0) error->one(FLERR, "Etype not defined yet");
              iucell_list[0] = index_node0; 
              iucell_list[1] = index_node0 + index_offset[0];
              iucell_list[2] = index_node0 + index_offset[0] 
                + index_offset[1];
              iucell_list[3] = index_node0 + index_offset[1];
              evec->create_pass_on_element(i, iucell_list, myetype, ictype, 0);
            }
            index_node0 += index_offset[1]+1;
          }
        }
        index_node0 += index_offset[0];
      }

      // check mass conservation 
      // number of unitcells in splitted elements and atoms must be equal 
      // to the number of atoms in the original element
      int sum = atom->nlocal - nalocal_previous;
      for (int j = nelocal_previous; j < element->nlocal; j++)
        sum += nucell[etype[j]];
      if (sum != nucell[ietype]) {
        printf("split element tag i = %d i = %d ietype = %d ncells = %d %d %d sum = %d new atom = %d\n", element->tag[i], i, ietype, ncellx, ncelly, ncellz, sum, atom->nlocal-nalocal_previous);
        for (int j = 0; j < n; j++)
          if (splitlist[j][0] == i)
            printf("split %d %d\n", splitlist[j][1], splitlist[j][2]);
        printf("Ncells = %d cellsize x =", ncells[ietype][0]);
        for (int iz = 0; iz < nelem[0]; iz++) 
          printf(" %d", cellsize[0][iz]);
        printf("\n");
        printf("Ncells = %d cellsize y =", ncells[ietype][1]);
        for (int iz = 0; iz < nelem[1]; iz++) 
          printf(" %d", cellsize[1][iz]);
        printf("\n");
        printf("Ncells = %d cellsize z =", ncells[ietype][2]);
        for (int iz = 0; iz < nelem[2]; iz++) 
          printf(" %d", cellsize[2][iz]);
        printf("\n");

        error->one(FLERR, "Mass not conserved during splitting");
      } 
    }

    for (ii = 0; ii < ndiscs; ii++) {
      i = idisclist[ii];
      evec->element2atom(i);
    }

    // delete all splitted and discretized elements

    i = 0;
    while (i < max_del_flag && i < element->nlocal) {
      if (del_flag_list[i]) {
        element->nlocal--;
        evec->copy(element->nlocal, i, 1);
        if (element->nlocal < max_del_flag)
          del_flag_list[i] = del_flag_list[element->nlocal];
        else del_flag_list[i] = 0;
      } else i++;
    }

    // clean up

    memory->destroy(del_flag_list);
    memory->destroy(idisclist); 
    memory->destroy(isplitlist); 
    memory->destroy(nelem_array); 
    memory->destroy(cellsize_array);
    delete [] iucell_list;
  }
  delete [] ictype; 
  return nsplits + ndiscs;
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

void FixAdaptive::grow_splitlist()
{
  maxsplit += maxsplit/2;
  memory->grow(splitlist, maxsplit, 3, "fix:splitlist");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

double FixAdaptive::memory_usage()
{
  return (atom->nmax + element->nmax) * sizeof(double);
}


// DEBUG FUNCTION
void FixAdaptive::write_neighbor(int ii)
{

  int i = list->ilist[ii];
  int jnum = list->numneigh[ii];

  if (jnum == 0) 
    printf("No neighbor for element %d me = %d\n", element->tag[i], comm->me); 

  int *jlist = list->firstneigh[ii];
  int *jindexlist = list->firstneighindex[ii];

  char *file = new char[100];
  sprintf(file, "dump_fix_adaptive_neighbor_%d_%d_%d.atom", element->tag[i], comm->me, update->ntimestep);
  FILE *fp = fopen(file, "w");


  fprintf(fp, "ITEM: TIMESTEP\n");

  fprintf(fp, "%d\n", update->ntimestep);

  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");

  int count = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    if (jindexlist[jj] >= 0) 
      count++;
  }

  fprintf(fp, "%d\n", count);

  char boundstr[9];
  domain->boundary_string(boundstr);

  if (domain->triclinic == 0) {
    fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
    fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
    fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
    fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
  } else {
    fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
    fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0], domain->boxhi_bound[0], domain->xy);
    fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1], domain->boxhi_bound[1], domain->xz);
    fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2], domain->boxhi_bound[2], domain->yz);
  }

  fprintf(fp, "ITEM: ATOMS id type x y z quad\n");

  int j;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    if (jindexlist[jj] >= 0) 
      fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %d\n", tag[j], type[j], x[j][0], x[j][1], x[j][2],jindexlist[jj]);
  }

  fclose(fp);
  delete [] file;

}

void FixAdaptive::write_surface_ucell(int i)
{
  /*
     char *file = new char[100];
     sprintf(file, "dump_fix_adaptive_surface_ucell_%d_%d_%d.atom", element->tag[i], comm->me, update->ntimestep);
     FILE *fp = fopen(file, "w");

     int *nsurface = element->nsurface;
     int ***surface_ucell = element->surface_ucell;
     int **nsurface_ucell = element->nsurface_ucell;

     int itype = element->etype[i];
     int sum = 0;
     for (int ii = 0; ii < nsurface[itype]; ii++)
     sum += nsurface_ucell[itype][ii];

     fprintf(fp, "ITEM: TIMESTEP\n");

     fprintf(fp, "%d\n", update->ntimestep);

     fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
     fprintf(fp, "%d\n", sum);

     char boundstr[9];
     domain->boundary_string(boundstr);

     if (domain->triclinic == 0) {
     fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
     } else {
     fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0], domain->boxhi_bound[0], domain->xy);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1], domain->boxhi_bound[1], domain->xz);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2], domain->boxhi_bound[2], domain->yz);
     }

     fprintf(fp, "ITEM: ATOMS id type x y z\n");

     int j, *jlist;
     double ***nodex = element->nodex;
     int *type = element->ctype;
     tagint *tag = element->tag;
     double coord[3];
     for (int ii = 0; ii < nsurface[itype]; ii++) {
     jlist = surface_ucell[itype][ii];
     for (int jj = 0; jj < nsurface_ucell[itype][ii]; jj++) {
     j = jlist[jj];
     element->evec->interpolate(coord, nodex, i, j, 3);
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e\n", tag[j], type[j], coord[0], coord[1], coord[2]);
     }
     }
     fclose(fp);
     delete [] file;
     */
}

void FixAdaptive::write_split_element(int i, int j)
{
  /*
     char *file = new char[100];
     sprintf(file, "dump_fix_adaptive_split_%d_me_%d.atom", nsplits_total, comm->me);
     FILE *fp = fopen(file, "w");

     int ietype = element->etype[i];
     int jetype = element->etype[j];
     int inucell = element->nucell[ietype];
     int jnucell = element->nucell[jetype];
     fprintf(fp, "ITEM: TIMESTEP\n");
     fprintf(fp, "%d\n", update->ntimestep);
     fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
     fprintf(fp, "%d\n", inucell+jnucell);
     char boundstr[9];
     domain->boundary_string(boundstr);
     if (domain->triclinic == 0) {
     fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
     } else {
     fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0], domain->boxhi_bound[0], domain->xy);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1], domain->boxhi_bound[1], domain->xz);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2], domain->boxhi_bound[2], domain->yz);
     }
     fprintf(fp, "ITEM: ATOMS id type x y z structure_id eid inode\n");

     int count = 1;
     tagint *tag = element->tag;
     double ***nodex = element->nodex;
     double **vapattern = structure->vector_ucell_atom;
     int **is_outer = element->is_outer;
     int **ia2i = element->ia2i;
     int **i2n = element->i2n;
     int iouter, jouter, iintg, jintg, inode, jnode;
     double coord[3];
     for (int iucell; iucell < inucell; iucell++) {
     iouter = is_outer[ietype][iucell]-1;
     iintg = ia2i[ietype][iucell];
     inode = -1;
     if (iintg >= 0 && iouter >= 0) 
     inode = i2n[ietype][iintg];
     element->evec->interpolate(coord, nodex, i, iucell, 3);
     if (iouter >= 0) 
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 1, coord[0], coord[1], coord[2], vapattern[i][iouter], tag[i], inode);
     else
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 1, coord[0], coord[1], coord[2], -1, tag[i], inode);
     }
     for (int jucell; jucell < jnucell; jucell++) {
     jouter = is_outer[jetype][jucell]-1;
     jintg = ia2i[jetype][jucell];
     jnode = -1;
     if (jintg >= 0 && jouter >= 0) 
     jnode = i2n[jetype][jintg];
     element->evec->interpolate(coord, nodex, j, jucell, 3);
     if (iouter >= 0) 
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 2, coord[0], coord[1], coord[2], vapattern[j][jouter], tag[j], jnode);
     else
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 2, coord[0], coord[1], coord[2], -1, tag[j], jnode);
     }
     fclose(fp);
     delete [] file;
     */
}

void FixAdaptive::write_split_atom(int i, int jnum, int *jlist)
{
  /*
     char *file = new char[100];
     sprintf(file, "dump_fix_adaptive_split_%d_me_%d.atom", nsplits_total, comm->me);
     FILE *fp = fopen(file, "w");

     int ietype = element->etype[i];
     int inucell = element->nucell[ietype];
     fprintf(fp, "ITEM: TIMESTEP\n");
     fprintf(fp, "%d\n", update->ntimestep);
     fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
     fprintf(fp, "%d\n", inucell+2);
     char boundstr[9];
     domain->boundary_string(boundstr);
     if (domain->triclinic == 0) {
     fprintf(fp, "ITEM: BOX BOUNDS %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[0], domain->boxhi[0]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[1], domain->boxhi[1]);
     fprintf(fp, "%-1.16e %-1.16e\n", domain->boxlo[2], domain->boxhi[2]);
     } else {
     fprintf(fp, "ITEM: BOX BOUNDS xy xz yz %s\n", boundstr);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[0], domain->boxhi_bound[0], domain->xy);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[1], domain->boxhi_bound[1], domain->xz);
     fprintf(fp, "%-1.16e %-1.16e %-1.16e\n", domain->boxlo_bound[2], domain->boxhi_bound[2], domain->yz);
     }
     fprintf(fp, "ITEM: ATOMS id type x y z structure_id eid\n");

     int count = 0;
     tagint *etag = element->tag;
     tagint *atag = atom->tag;
     double ***nodex = element->nodex;
     double **x = atom->x;
     double **vapattern = structure->vector_ucell_atom;
     double *apattern = structure->vector_atom;
     int **is_outer = element->is_outer;
     int iouter, jouter, inode, iintg;
     int **ia2i = element->ia2i;
     int **i2n = element->i2n;
     double coord[3];
     for (int jj = 0; jj < jnum; jj++) {
     int j = jlist[jj];
     count = MAX(count, atag[j]);
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n"
     , atag[j], 2, x[j][0], x[j][1], x[j][2], apattern[j], -1, -1);
     }
     count++;
     for (int iucell; iucell < inucell; iucell++) {
     iouter = is_outer[ietype][iucell]-1;
     inode = -1;
     iintg = ia2i[ietype][iucell];
     if (iintg >= 0)
     inode = i2n[ietype][iintg];
     element->evec->interpolate(coord, nodex, i, iucell, 3);
     if (iouter >= 0) 
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 1, coord[0], coord[1], coord[2], :apattern[i][iouter], etag[i], inode);
     else
     fprintf(fp, "%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n", count++, 1, coord[0], coord[1], coord[2], -1, etag[i], inode);
     }
     fclose(fp);
     delete [] file;
     */
}
