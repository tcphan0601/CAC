#include <mpi.h>
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

/* ---------------------------------------------------------------------- */

FixAdaptive::FixAdaptive(CAC *cac, int narg, char **arg) :
  Fix(cac, narg, arg),
  id_structure(NULL), splitlist(NULL)
{
  debug = 0;
  create_attribute = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  //if (strcmp(element->element_style,"cac") != 0) 
  //  error->all(FLERR,"Fix adaptive only works for rhb element style");
  if (element->max_apc > 1) error->all(FLERR,"Fix adaptive does not work for multiple atoms per unit cell yet");

  if (narg < 6) error->all(FLERR,"Illegal fix adaptive command");

  nevery = universe->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix adaptive command");

  int n = strlen(arg[4]) + 1; 
  id_structure = new char[n];
  strcpy(id_structure,arg[4]);
  istructure = modify->find_compute(id_structure);
  if (istructure < 0)
    error->all(FLERR,"Compute ID for fix adative does not exist");
  structure = modify->compute[istructure];
  perfect_crystal = universe->inumeric(FLERR,arg[5]);

  /*
     id_stress = new char[n];
     strcpy(id_stress,arg[6]);
     istress = modify->find_compute(id_stress);
     if (istress < 0)
     error->all(FLERR,"Stress Compute ID for fix adative does not exist");
     stress = modify->compute[istress];
     */
  // optional args

  tol = 1e-3;
  eig_frac_thres = 0.5;
  cut_atom_neigh = 10;
  split_width = 0.5;
  thres = 10;
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adaptive command");
      tol = universe->numeric(FLERR,arg[iarg+1]);
      if (tol <= 0) error->all(FLERR,"Tolerant for fix adaptive command must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cut") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adaptive command");
      cut_atom_neigh = universe->numeric(FLERR,arg[iarg+1]);
      if (cut_atom_neigh < 0) error->all(FLERR,"Neighbor cut for fix adaptive command must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"swidth") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adaptive command");
      split_width = universe->numeric(FLERR,arg[iarg+1]);
      if (split_width <= 0 || split_width > 1) error->all(FLERR,"Split width for fix adaptive command must be between 0 and 1 (unit cell)");
      iarg += 2;
    } else if (strcmp(arg[iarg],"thres") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adaptive command");
      thres = universe->inumeric(FLERR,arg[iarg+1]);
      if (thres < 0) error->all(FLERR,"Invalid thres for fix adaptive command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"debug") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix adaptive command");
      if (strcmp(arg[iarg+1],"yes"))
        debug = 1;
      else if (strcmp(arg[iarg+1],"no"))
        debug = 0;
      else error->all(FLERR,"Illegal fix adaptive command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix adaptive command");
  }

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
  // need an occasional full intpl neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->atomlist = 0;
  neighbor->requests[irequest]->elemonlylist = 1;
  neighbor->requests[irequest]->intglist = 0;
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
    error->warning(FLERR,"Fix adaptive cutoff may be too large to find "
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
    domain->x2lamda(atom->nlocal,atom->x);
    domain->x2lamda(element->nlocal,element->x);
    domain->nodex2lamda(element->nlocal,element->nodex);
  }
  domain->pbc();
  domain->reset_box();
  comm->setup_exchange();
  comm->exchange();
  comm->setup_borders();
  if (neighbor->style) neighbor->setup_bins(); 
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();

  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal,atom->x);
    domain->lamda2x(element->nlocal,element->x);
    domain->lamda2nodex(element->nlocal,element->nodex);
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

void FixAdaptive::post_integrate()
{

  system_changed = 0;

  // return if not an adaptive timestep

  if (update->ntimestep < next_reneighbor) return;

  // update slip plane for all elements

  element->evec->update_surface_plane();

  // invoke stress and structure compute peratom

  if (structure->invoked_peratom != update->ntimestep) {
    structure->compute_peratom();
    structure->preneighflag = 1;
  }

  //if (stress->invoked_peratom != update->ntimestep)
  //  stress->compute_peratom();

  // acquire ghost structure id & stress
  //for (int i = 0; i < element->nlocal+element->nghost; i++)
  //  if (element->tag[i] == 54)
  //    printf("**** me = %d before before comm j = %d jintpl = %d pattern = %g\n",comm->me,i,11,structure->vector_intpl_atom[i][11]);

  comm->forward_comm_compute(structure);
  //comm->forward_comm_compute(stress);
  //
  //for (int i = 0; i < element->nlocal+element->nghost; i++)
  //  if (element->tag[i] == 54)
  //    printf("**** me = %d before after comm j = %d jintpl = %d pattern = %g\n",comm->me,i,11,structure->vector_intpl_atom[i][11]);

  // build neighbor list for elements

  //MPI_Barrier(world);
  //if (comm->me == 0) printf("Step %d Start build neighbor (fix_adaptive)\n",update->ntimestep);
  neighbor->build_one(list,1);
  //MPI_Barrier(world);
  //if (comm->me == 0) printf("Step %d Finish build neighbor (fix_adaptive)\n",update->ntimestep);

  int inum = list->einum;
  int *ilist = list->eilist;
  int *numneighe2a = list->numneighe2a;
  int *numneighe2e = list->numneighe2e;
  int **firstneighe2a = list->firstneighe2a;
  int **firstneighe2e = list->firstneighe2e;
  int i,j,k,l,n,ii,jj,kk,ll,jintpl,jetype;
  int jnum,knum,*jlist,*klist;
  double delx,dely,delz;
  double **ax = atom->x;
  double **ex = element->x;
  int *amask = atom->mask;
  int *emask = element->mask;
  int *etype = element->etype;
  double *apattern = structure->vector_atom;
  double **iapattern = structure->vector_intpl_atom;
  int **is_outer = element->is_outer;
  int jouter;
  //double **atom_stress = stress->array_atom;
  int *nsurface = element->nsurface;
  int *element_shape_ids = element->element_shape_ids;
  int ***surface_intpl = element->surface_intpl;
  int **nsurface_intpl = element->nsurface_intpl;
  int cut[MAXCUT],cutdim;
  int nsplit = 0;

  double P[3];                    // Point on splip plane (calculated as centroid)
  double normal[3];
  double A[3][3];                 // 3x3 matrix for SVD, will be V after SVD
  double eigen_values[3];


  // loop through neighbors of elements to find dislocation nearby
  // add element to splitlist
  //MPI_Barrier(world);
  //if (comm->me == 0) printf("Step %d Start check splitting (fix_adaptive)\n",update->ntimestep);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    //if (element->tag[i]==681 
    //    || element->tag[i] == 683
    //    || element->tag[i] == 381
    //    || element->tag[i] == 382
    //    ) write_neighbor(ii);
    if (!(emask[i] & groupbit)) continue;
    
    if (element_shape_ids[etype[i]] != Element::HEXAHEDRON) 
      continue;

    // loop over atom neighbors

    jnum = numneighe2a[ii];
    jlist = firstneighe2a[ii];
    n = 0;
    P[0] = P[1] = P[2] = 0.0;
    A[0][0] = A[0][1] = A[0][2] = 0.0;
    A[1][0] = A[1][1] = A[1][2] = 0.0;
    A[2][0] = A[2][1] = A[2][2] = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (amask[j] & groupbit) {
        P[0] += ax[j][0];
        P[1] += ax[j][1];
        P[2] += ax[j][2];
      }
      n++;
    }

    P[0] /= n;
    P[1] /= n;
    P[2] /= n;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (amask[j] & groupbit) {
        delx = ax[j][0] - P[0];
        dely = ax[j][1] - P[1];
        delz = ax[j][2] - P[2];
        A[0][0] += delx*delx;
        A[1][1] += dely*dely;
        A[2][2] += delz*delz;
        A[0][1] += delx*dely;
        A[0][2] += delx*delz;
        A[1][2] += dely*delz;
      }
    }

    if (n >= thres) {
      A[1][0] = A[0][1];
      A[2][0] = A[0][2];
      A[2][1] = A[1][2];

      eigen_decomposition(A,eigen_values);

      // the smallest eigen_values must be much smaller than the next

      if (eigen_values[0]/eigen_values[1] < eig_frac_thres) {
        normal[0] = A[0][0];
        normal[1] = A[1][0];
        normal[2] = A[2][0];
        if (element->evec->check_split_element(i,normal,P,tol,split_width,cut,cutdim)) {
          if (nsplit == maxsplit) grow_splitlist();
          splitlist[nsplit][0] = i;
          splitlist[nsplit][1] = cut[0];
          splitlist[nsplit][2] = cutdim;
          if (debug) {
            write_split_atom(i,jnum,jlist);
            printf("Split #%d me = %d element %d from atomistic region P = %g %g %g n = %g %g %g\n",
                nsplits_total,comm->me,element->tag[i],P[0],P[1],P[2],normal[0],normal[1],normal[2]);
          }
          nsplits_total++;
          nsplit++;
        }
      }
    }

    // loop over element neighbors

    jnum = numneighe2e[ii];
    jlist = firstneighe2e[ii];

    //if (element->tag[i]>50) write_surface_intpl(i);
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jetype = etype[j];
      j &= NEIGHMASK;
      if (!(emask[j] & groupbit)) continue;

      if (element_shape_ids[jetype] == Element::HEXAHEDRON) {
        for (kk = 0; kk < nsurface[element_shape_ids[jetype]]; kk++) {
          knum = nsurface_intpl[jetype][kk];
          klist = surface_intpl[jetype][kk];
          for (k = 0; k < knum; k++) {
            jintpl = klist[k];
            jouter = is_outer[jetype][jintpl]-1;
            if (jouter < 0) error->one(FLERR,"Error is_outer array");
            if (iapattern[j][jouter] == perfect_crystal) continue;
            if (element->evec->check_split_element(i,j,kk,tol,cut,cutdim,MAXCUT)) {
              if (nsplit == maxsplit) grow_splitlist();
              splitlist[nsplit][0] = i;
              splitlist[nsplit][1] = cut[0];
              splitlist[nsplit][2] = cutdim;
              if (debug) {
                write_split_element(i,j);
                printf("Split #%d me = %d element %d from element %d surface %d\n"
                    ,nsplits_total,comm->me,element->tag[i],element->tag[j],kk);
              }
              nsplits_total++;
              nsplit++;
            }
            break;
          }
        }
      } else if (element_shape_ids[jetype] == Element::QUADRILATERAL) {
        if (domain->dimension == 3) {
          knum = nsurface_intpl[jetype][0];
          klist = surface_intpl[jetype][0];
          for (k = 0; k < knum; k++) {
            jintpl = klist[k];
            jouter = is_outer[jetype][jintpl]-1;
            if (jouter < 0) error->one(FLERR,"Error is_outer array");
            if (iapattern[j][jouter] == perfect_crystal) continue;
            int ncut = element->evec->check_split_element(i,j,0,tol,cut,cutdim,MAXCUT);
            if (ncut) { 
              if (nsplit + ncut > maxsplit) grow_splitlist();
              for (int icut = 0; icut < ncut; icut++) {
                splitlist[nsplit][0] = i;
                splitlist[nsplit][1] = cut[icut];
                splitlist[nsplit][2] = cutdim;
                nsplits_total++;
                nsplit++;
              if (debug) {
                write_split_element(i,j);
                printf("Split #%d me = %d element %d from element %d surface %d\n"
                    ,nsplits_total,comm->me,element->tag[i],element->tag[j],kk);
              }

              }
            }
            break;
          }
        } else {
          // 2D Quad element here
        }
      } else {
        // add new element types here
      }
    }
  }
  //MPI_Barrier(world);
  //if (comm->me == 0) printf("Step %d Done check splitting (fix_adaptive)\n",update->ntimestep);

  int nsplit_all;
  MPI_Allreduce(&nsplit,&nsplit_all,1,MPI_INT,MPI_MAX,world);

  // return and update next reneighbor if no plane need to split
  // this will not force reneighbor at this time step

  if (nsplit_all == 0) {
    next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    if (comm->me == 0) fprintf(screen,"No new dislocation found from fix adaptive timestep = %d\n",update->ntimestep);
    return;
  } else 
    system_changed = 1;


  pending = 1;
  int nelocal = element->nlocal;
  int nalocal = atom->nlocal;

  // split elements

  int nesplit,nesplit_all;
  nesplit = element->evec->split_element(splitlist,nsplit);
  //MPI_Barrier(world);
  //if (comm->me == 0) printf("Step %d Done splitting (fix_adaptive)\n",update->ntimestep);

  // extend tags to new atoms and elements

  atom->tag_extend();
  element->tag_extend();

  int nnew_atoms = atom->nlocal - nalocal;
  int nnew_elems = element->nlocal - nelocal;
  int nnew_atoms_all,nnew_elems_all;
  MPI_Allreduce(&nnew_atoms,&nnew_atoms_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nnew_elems,&nnew_elems_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nesplit,&nesplit_all,1,MPI_INT,MPI_SUM,world);

  // update number of elements and number of atoms

  atom->natoms += nnew_atoms_all;
  element->nelements += nnew_elems_all;

  // migrates new elements if there are any

  if (nnew_elems_all) {
    if (domain->triclinic) {
      domain->x2lamda(atom->nlocal,atom->x);
      domain->x2lamda(element->nlocal,element->x);
      domain->nodex2lamda(element->nlocal,element->nodex);
    }
    domain->reset_box();

    IrregularComm *irrcomm = new IrregularComm(cac);
    irrcomm->migrate(1);
    delete irrcomm;

    if (domain->triclinic) {
      domain->lamda2x(atom->nlocal,atom->x);
      domain->lamda2x(element->nlocal,element->x);
      domain->lamda2nodex(element->nlocal,element->nodex);
    }
  }

  if (comm->me == 0) {
    //printf("Ntimestep = %d\n",update->ntimestep);
    if (screen) fprintf(screen,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nesplit_all,nnew_atoms_all,nnew_elems_all+nesplit_all);
    if (logfile) fprintf(logfile,"Dislocation detected and %d elements need splitting.\n  %d new atoms created.\n  %d new elements created\n",nesplit_all,nnew_atoms_all,nnew_elems_all+nesplit_all);
  }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based/elem-based arrays
   ------------------------------------------------------------------------- */

void FixAdaptive::grow_splitlist()
{
  maxsplit += maxsplit/2;
  memory->grow(splitlist,maxsplit,3,"fix:splitlist");
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
  int i = list->eilist[ii];
  int jnum = list->numneighe2a[ii];
  if (jnum == 0) 
    printf("No neighbor for element %d me = %d\n",element->tag[i],comm->me); 

  int *jlist = list->firstneighe2a[ii];

  char *file = new char[100];
  sprintf(file,"dump_fix_adaptive_neighbor_%d_%d_%d.atom",element->tag[i],comm->me,update->ntimestep);
  FILE *fp = fopen(file,"w");


  fprintf(fp,"ITEM: TIMESTEP\n");

  fprintf(fp,"%d\n",update->ntimestep);

  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",jnum);

  char boundstr[9];
  domain->boundary_string(boundstr);

  if (domain->triclinic == 0) {
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
  } else {
    fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
  }

  fprintf(fp,"ITEM: ATOMS id type x y z\n");

  int j;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e\n",tag[j],type[j],x[j][0],x[j][1],x[j][2]);
  }

  fclose(fp);
  delete [] file;
}

void FixAdaptive::write_surface_intpl(int i)
{
  char *file = new char[100];
  sprintf(file,"dump_fix_adaptive_surface_intpl_%d_%d_%d.atom",element->tag[i],comm->me,update->ntimestep);
  FILE *fp = fopen(file,"w");

  int *nsurface = element->nsurface;
  int ***surface_intpl = element->surface_intpl;
  int **nsurface_intpl = element->nsurface_intpl;

  int itype = element->etype[i];
  int sum = 0;
  for (int ii = 0; ii < nsurface[itype]; ii++)
    sum += nsurface_intpl[itype][ii];

  fprintf(fp,"ITEM: TIMESTEP\n");

  fprintf(fp,"%d\n",update->ntimestep);

  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",sum);

  char boundstr[9];
  domain->boundary_string(boundstr);

  if (domain->triclinic == 0) {
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
  } else {
    fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
  }

  fprintf(fp,"ITEM: ATOMS id type x y z\n");

  int j,*jlist;
  double ***nodex = element->nodex;
  int *type = element->ctype;
  tagint *tag = element->tag;
  double coord[3];
  for (int ii = 0; ii < nsurface[itype]; ii++) {
    jlist = surface_intpl[itype][ii];
    for (int jj = 0; jj < nsurface_intpl[itype][ii]; jj++) {
      j = jlist[jj];
      element->evec->interpolate(coord,nodex,i,j,3);
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e\n",tag[j],type[j],coord[0],coord[1],coord[2]);
    }
  }
  fclose(fp);
  delete [] file;
}

void FixAdaptive::write_split_element(int i, int j)
{
  char *file = new char[100];
  sprintf(file,"dump_fix_adaptive_split_%d_me_%d.atom",nsplits_total,comm->me);
  FILE *fp = fopen(file,"w");

  int ietype = element->etype[i];
  int jetype = element->etype[j];
  int inintpl = element->nintpl[ietype];
  int jnintpl = element->nintpl[jetype];
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",inintpl+jnintpl);
  char boundstr[9];
  domain->boundary_string(boundstr);
  if (domain->triclinic == 0) {
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
  } else {
    fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
  }
  fprintf(fp,"ITEM: ATOMS id type x y z structure_id eid inode\n");

  int count = 1;
  tagint *tag = element->tag;
  double ***nodex = element->nodex;
  double **iapattern = structure->vector_intpl_atom;
  int **is_outer = element->is_outer;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  int iouter,jouter,iintg,jintg,inode,jnode;
  double coord[3];
  for (int iintpl; iintpl < inintpl; iintpl++) {
    iouter = is_outer[ietype][iintpl]-1;
    iintg = ia2i[ietype][iintpl];
    inode = -1;
    if (iintg >= 0 && iouter >= 0) 
      inode = i2n[ietype][iintg];
    element->evec->interpolate(coord,nodex,i,iintpl,3);
    if (iouter >= 0) 
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,1,coord[0],coord[1],coord[2],iapattern[i][iouter],tag[i],inode);
    else
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,1,coord[0],coord[1],coord[2],-1,tag[i],inode);
  }
  for (int jintpl; jintpl < jnintpl; jintpl++) {
    jouter = is_outer[jetype][jintpl]-1;
    jintg = ia2i[jetype][jintpl];
    jnode = -1;
    if (jintg >= 0 && jouter >= 0) 
      jnode = i2n[jetype][jintg];
    element->evec->interpolate(coord,nodex,j,jintpl,3);
    if (iouter >= 0) 
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,2,coord[0],coord[1],coord[2],iapattern[j][jouter],tag[j],jnode);
    else
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,2,coord[0],coord[1],coord[2],-1,tag[j],jnode);
  }
  fclose(fp);
  delete [] file;

}

void FixAdaptive::write_split_atom(int i, int jnum, int *jlist)
{
  char *file = new char[100];
  sprintf(file,"dump_fix_adaptive_split_%d_me_%d.atom",nsplits_total,comm->me);
  FILE *fp = fopen(file,"w");

  int ietype = element->etype[i];
  int inintpl = element->nintpl[ietype];
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",inintpl+2);
  char boundstr[9];
  domain->boundary_string(boundstr);
  if (domain->triclinic == 0) {
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
  } else {
    fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[0],domain->boxhi_bound[0],domain->xy);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[1],domain->boxhi_bound[1],domain->xz);
    fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",domain->boxlo_bound[2],domain->boxhi_bound[2],domain->yz);
  }
  fprintf(fp,"ITEM: ATOMS id type x y z structure_id eid\n");

  int count = 0;
  tagint *etag = element->tag;
  tagint *atag = atom->tag;
  double ***nodex = element->nodex;
  double **x = atom->x;
  double **iapattern = structure->vector_intpl_atom;
  double *apattern = structure->vector_atom;
  int **is_outer = element->is_outer;
  int iouter,jouter,inode,iintg;
  int **ia2i = element->ia2i;
  int **i2n = element->i2n;
  double coord[3];
  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    count = MAX(count,atag[j]);
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n"
        ,atag[j],2,x[j][0],x[j][1],x[j][2],apattern[j],-1,-1);
  }
  count++;
  for (int iintpl; iintpl < inintpl; iintpl++) {
    iouter = is_outer[ietype][iintpl]-1;
    inode = -1;
    iintg = ia2i[ietype][iintpl];
    if (iintg >= 0)
      inode = i2n[ietype][iintg];
    element->evec->interpolate(coord,nodex,i,iintpl,3);
    if (iouter >= 0) 
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,1,coord[0],coord[1],coord[2],iapattern[i][iouter],etag[i],inode);
    else
      fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %g %d %d\n",count++,1,coord[0],coord[1],coord[2],-1,etag[i],inode);
  }
  fclose(fp);
  delete [] file;
}
