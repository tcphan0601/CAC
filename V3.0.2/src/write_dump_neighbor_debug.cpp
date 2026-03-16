#include <mpi.h>
#include <string.h>
#include "write_dump_neighbor_debug.h"
#include "atom.h"
#include "atom_vec.h"
#include "element.h"
#include "element_vec.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

//enum{II, IJ};

/* ---------------------------------------------------------------------- */

WriteDumpNeighborDebug::WriteDumpNeighborDebug(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
}

/* ----------------------------------------------------------------------
   called as write_dump_neighbor_debug_neighbor_debug command in input script
   all processor write out its own file contain local+ghost atoms/elements/nodes 
   ------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Write_data_neighbor_debug command before simulation box is defined");

  if (narg < 1) error->all(FLERR, "Illegal write_dump_neighbor_debug command");

  singleflag = 0;
  if (narg > 1) {
    if (strcmp(arg[1], "gauss") == 0) {
      singleflag = 1;
      if (5 != narg) error->all(FLERR, "Illegal write_dump_neighbor_debug command");
      id = universe->inumeric(FLERR, arg[2]);
      igcell = universe->inumeric(FLERR, arg[3]);
      ibasis = universe->inumeric(FLERR, arg[4]);
    } else if (strcmp(arg[1], "atom") == 0) {
      singleflag = 2;
      if (3 != narg) error->all(FLERR, "Illegal write_dump_neighbor_debug command");
      id = universe->inumeric(FLERR, arg[2]);
    }
  } else if (narg < 1) error->all(FLERR, "Illegal write_dump_neighbor_debug command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 20;
  char *file = new char[n];

  if (singleflag != 3) {
    if ((ptr = strchr(arg[0], '*'))) {
      *ptr = '\0';
      sprintf(file, "%s" BIGINT_FORMAT "%s_proc_%d", arg[0], update->ntimestep, ptr+1, me);
    } else sprintf(file, "%s_proc_%d", arg[0], me);
  } else {
    if ((ptr = strchr(arg[0], '*'))) {
      *ptr = '\0';
      sprintf(file, "%s" BIGINT_FORMAT "%s_proc_%d.dat", arg[0], update->ntimestep, ptr+1, me);
    } else sprintf(file, "%s_proc_%d.dat", arg[0], me);
  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_dump_neighbor_debug immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  //  if (noinit == 0) {



  if (comm->me == 0 && screen)
    fprintf(screen, "System init for write_dump_neighbor_debug ...\n");
  cac->init();

  // move atoms to new processors before writing file
  // do setup_pre_exchange to force update of per-atom info if needed
  // enforce PBC in case atoms are outside box
  // call borders() to rebuild atom map since exchange() destroys map

  atom->setup();
  element->setup();
  modify->setup_pre_exchange();
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
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost, atom->x);
    domain->lamda2x(element->nlocal+element->nghost, element->x);
    domain->lamda2nodex(element->nlocal+element->nghost, element->nodex);
  }
  modify->setup_pre_neighbor();
  neighbor->build();

  if (force->pair) {
    list = force->pair->list;
  } else error->all(FLERR, "pair has not been defined");
  if (list == NULL) error->all(FLERR, "neighbor list not build yet");


  write(file);
  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
   ------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::write(char *file)
{
  // special case where reneighboring is not done in integrator
  //   on timestep data file is written (due to build_once being set)
  // if box is changing, must be reset, else data file will have
  //   wrong box size and atoms will be lost when data file is read
  // other calls to pbc and domain and comm are not made, 
  //   b/c they only make sense if reneighboring is actually performed

  //if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write data file

  bigint nalocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nalocal, &natoms, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == Thermo::ERROR)
    error->all(FLERR, "Atom count is inconsistent, cannot write data file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal, &nelements, 1, MPI_CAC_BIGINT, MPI_SUM, world);
  if (nelements != element->nelements && output->thermo->elemlostflag == Thermo::ERROR)
    error->all(FLERR, "Element count is inconsistent, cannot write data file");


  // sum up bond, angle counts
  // may be different than atom->nbonds, nangles if broken/turned-off
  /*
     if (atom->molecular == 1 && (atom->nbonds || atom->nbondtypes)) {
     nbonds_local = atom->avec->pack_bond(NULL);
     MPI_Allreduce(&nbonds_local, &nbonds, 1, MPI_CAC_BIGINT, MPI_SUM, world);
     }
     if (atom->molecular == 1 && (atom->nangles || atom->nangletypes)) {
     nangles_local = atom->avec->pack_angle(NULL);
     MPI_Allreduce(&nangles_local, &nangles, 1, MPI_CAC_BIGINT, MPI_SUM, world);
     }
     */
  // open data file

  if (singleflag == 1) {
    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }
    ielem = element->map(id);
    int *etype = element->etype;
    int *ngcell = element->ngcell;
    int *apc = element->apc;

    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }
    if (ielem < 0 || ielem >= element->nlocal) return;
    ii = igcell * apc[etype[ielem]] + ibasis + atom->nlocal;
    for (int i = 0; i < ielem; i++) {
      ii += ngcell[etype[i]] * apc[etype[i]];
    }


  } else if (singleflag == 2) {
    int mapflag = 0;
    if (atom->map_style == 0) {
      mapflag = 1;
      atom->map_init();
      atom->map_set();
    }
    iatom = atom->map(id);
    if (mapflag) {
      atom->map_delete();
      atom->map_style = 0;
    }

    if (iatom < 0 || iatom >= atom->nlocal) return;
  }


  fp = fopen(file, "w");
  if (fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open data file %s", file);
    error->one(FLERR, str);
  }

  header();
  int count = 1;
  tagint *atag = atom->tag;
  tagint *etag = element->tag;
  double **ax = atom->x;

  ElementVec *evec = element->evec;
  int i, iindex, iucell, ietype, j, jindex, jucell, jbasis;
  int jetype, japc, iapc;
  double ****nodex = element->nodex;
  double coord[3];
  int *apc = element->apc;
  int *etype = element->etype;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;
  int **g2u = element->g2u;
  int itype, jtype, itag, jtag;

  if (singleflag == 0) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      iindex = iindexlist[ii];
      if (iindex < 0) {
        itype = 1;
        coord[0] = ax[i][0];
        coord[1] = ax[i][1];
        coord[2] = ax[i][2];
        igcell = ibasis = -1;
        itag = atag[i];
      } else {
        itag = etag[i];
        itype = 2;
        ietype = etype[i];
        iapc = apc[ietype];
        igcell = iindex / iapc;
        ibasis = iindex % iapc;
        iucell = g2u[ietype][igcell];
        evec->interpolate(coord, nodex, i, ibasis, iucell, 3);
      }
      fprintf(fp, "%d %d %g %g %g %d %d %d %d\n"
          , count++, itype, coord[0], coord[1], coord[2], itag, numneigh[ii], igcell, ibasis);
    }
  } else {
    int jnum, *jlist, *jindexlist;
    if (singleflag == 1) {
      jlist = firstneigh[ii];
      jindexlist = firstneighindex[ii];
      jnum = numneigh[ii];
      evec->interpolate(coord, nodex, ielem, ibasis, iucell, 3);
    } else if (singleflag == 2) {
      jlist = firstneigh[iatom];
      jindexlist = firstneighindex[iatom];
      jnum = numneigh[iatom];
      coord[0] = ax[iatom][0];
      coord[1] = ax[iatom][1];
      coord[2] = ax[iatom][2];
    }
    fprintf(fp, "%d 1 %g %g %g %d %d %d\n"
        , count++, coord[0], coord[1], coord[2], id, -2, -2);
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      if (jindex < 0) {
        jtype = 2;
        coord[0] = ax[j][0];
        coord[1] = ax[j][1];
        coord[2] = ax[j][2];
        jucell = jbasis = -1;
        jtag = atag[j];
      } else {
        jtag = etag[j];
        jtype = 3;
        jetype = etype[j];
        japc = apc[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        evec->interpolate(coord, nodex, j, jbasis, jucell, 3);
      }
      fprintf(fp, "%d %d %g %g %g %d %d %d %d\n"
          , count++, jtype, coord[0], coord[1], coord[2], jtag, jucell, jbasis);
    }

  }
  fclose(fp);
}

/* ----------------------------------------------------------------------
   proc me writes out data file header
   ------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::header()
{
  int n = 0;
  int nnode = element->npe[element->etype[ielem]];
  if (singleflag == 0) {
    n = atom->nlocal;
    int *ngcell= element->ngcell;
    int *apc = element->apc;
    int *etype = element->etype;
    for (int i = 0; i < element->nlocal; i++)
      n += ngcell[etype[i]]*apc[etype[i]];
  } else if (singleflag == 1 && ii >= 0) {
    n = list->numneigh[ii] + 1;
  } else if (singleflag == 2 && iatom >= 0) {
    n = list->numneigh[iatom] + 1;
  }

  fprintf(fp, "ITEM: TIMESTEP\n");

  fprintf(fp, "%d\n", me);

  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, "%d\n", n);
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
  if (singleflag == 0) {
    fprintf(fp, "ITEM: ATOMS id type x y z tag numneigh igcell ibasis\n");
  } else if (singleflag == 1) {
    fprintf(fp, "ITEM: ATOMS id type x y z tag iucell ibasis\n");
  } else if (singleflag == 2) {
    fprintf(fp, "ITEM: ATOMS id type x y z tag iucell ibasis\n");
  }

}


