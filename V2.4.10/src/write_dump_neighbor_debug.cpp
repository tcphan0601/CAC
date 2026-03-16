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

//enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteDumpNeighborDebug::WriteDumpNeighborDebug(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_dump_neighbor_debug_neighbor_debug command in input script
   all processor write out its own file contain local+ghost atoms/elements/nodes 
------------------------------------------------------------------------- */

void WriteDumpNeighborDebug::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data_neighbor_debug command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_dump_neighbor_debug command");

  singleflag = 0;
  if (narg > 1) {
    if (strcmp(arg[1],"iintg") == 0) {
      singleflag = 1;
      if (4 != narg) error->all(FLERR,"Illegal write_dump_neighbor_debug command");
      id = universe->inumeric(FLERR,arg[2]);
      iintg = universe->inumeric(FLERR,arg[3]);
    } else if (strcmp(arg[1],"atom") == 0) {
      singleflag = 2;
      if (3 != narg) error->all(FLERR,"Illegal write_dump_neighbor_debug command");
      id = universe->inumeric(FLERR,arg[2]);
    } else if (strcmp(arg[1],"element") == 0) {
      singleflag = 3;
      if (3 != narg) error->all(FLERR,"Illegal write_dump_neighbor_debug command");
      id = universe->inumeric(FLERR,arg[2]);
    } else error->all(FLERR,"Illegal write_dump_neighbor_debug command");
  } else if (narg < 1) error->all(FLERR,"Illegal write_dump_neighbor_debug command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 20;
  char *file = new char[n];

  if (singleflag != 3) {
    if ((ptr = strchr(arg[0],'*'))) {
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s_proc_%d",arg[0],update->ntimestep,ptr+1,me);
    } else sprintf(file,"%s_proc_%d",arg[0],me);
  } else {
    if ((ptr = strchr(arg[0],'*'))) {
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s_proc_%d.dat",arg[0],update->ntimestep,ptr+1,me);
    } else sprintf(file,"%s_proc_%d.dat",arg[0],me);
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
    fprintf(screen,"System init for write_dump_neighbor_debug ...\n");
  cac->init();

  // move atoms to new processors before writing file
  // do setup_pre_exchange to force update of per-atom info if needed
  // enforce PBC in case atoms are outside box
  // call borders() to rebuild atom map since exchange() destroys map

  atom->setup();
  element->setup();
  modify->setup_pre_exchange();
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
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (domain->triclinic) {
    domain->lamda2x(atom->nlocal+atom->nghost,atom->x);
    domain->lamda2x(element->nlocal+element->nghost,element->x);
    domain->lamda2nodex(element->nlocal+element->nghost,element->nodex);
  }

  modify->setup_pre_neighbor();
  neighbor->build();

  if (force->pair) {
    if (force->pair->double_neighbor_layers) {
      listhalf = force->pair->listhalf;
      listfull = force->pair->listfull;
      list = listfull; 
    } else {
      list = force->pair->list;
    }
  } else error->all(FLERR,"pair has not been defined");
  if (list == NULL) error->all(FLERR,"neighbor list not build yet");


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
  MPI_Allreduce(&nalocal,&natoms,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->atomlostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  bigint nelocal = element->nlocal;
  bigint nelements;
  MPI_Allreduce(&nelocal,&nelements,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (nelements != element->nelements && output->thermo->elemlostflag == Thermo::ERROR)
    error->all(FLERR,"Element count is inconsistent, cannot write data file");


  // sum up bond,angle counts
  // may be different than atom->nbonds,nangles if broken/turned-off
  /*
     if (atom->molecular == 1 && (atom->nbonds || atom->nbondtypes)) {
     nbonds_local = atom->avec->pack_bond(NULL);
     MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_CAC_BIGINT,MPI_SUM,world);
     }
     if (atom->molecular == 1 && (atom->nangles || atom->nangletypes)) {
     nangles_local = atom->avec->pack_angle(NULL);
     MPI_Allreduce(&nangles_local,&nangles,1,MPI_CAC_BIGINT,MPI_SUM,world);
     }
     */
  // open data file

  iatom = ielem = iintg_local = -1;
  if (singleflag == 1) {
    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }
    ielem = element->map(id);
    int **i2ia = element->i2ia;
    int *e2ilist = list->e2ilist;

    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }
    if (ielem < 0 || ielem >= element->nlocal) return;

    iintg_local = e2ilist[ielem] + iintg;
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
  } else if (singleflag == 3) {
    int mapflag = 0;
    if (element->map_style == 0) {
      mapflag = 1;
      element->map_init();
      element->map_set();
    }
    ielem = element->map(id);
    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }
    if (ielem < 0 || ielem >= element->nlocal) return;

  }

  fp = fopen(file,"w");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open data file %s",file);
    error->one(FLERR,str);
  }

  header();
  int count = 1;

  int j,jintpl;
  double ***nodex = element->nodex;
  double coord[3];

  if (singleflag == 0) {
    // per atom info

    if (atom->nlocal) {
      int *tag = atom->tag;
      double **x = atom->x;
      int *numneigha2a = list->numneigha2a;
      int *numneigha2ia = list->numneigha2ia;

      if (list->ainum < atom->nlocal) error->one(FLERR,"number of neighbors for atoms not defined yet");
      for (int i = 0; i < atom->nlocal; i++) 
        fprintf(fp,"%d 1 %g %g %g %d %d %d %d -1\n",
            count++,x[i][0],x[i][1],x[i][2],tag[i],numneigha2a[i]+numneigha2ia[i],numneigha2a[i],numneigha2ia[i]);

    }

    // per element info

    if (element->nlocal) {

      int *tag = element->tag;
      int *nintg = element->nintg;
      int *etype = element->etype;
      int **i2ia = element->i2ia;
      int *numneighi2a = list->numneighi2a;
      int *numneighi2ia = list->numneighi2ia;
      int *e2ilist = list->e2ilist;

      ElementVec *evec = element->evec;

      if (list->einum < element->nlocal) error->one(FLERR,"number of neighbors for integration points not defined yet");
      for (int i = 0; i < element->nlocal; i++) {
        int itype = etype[i];
        for (int j = 0; j < nintg[itype]; j++) {
          int iintpl = i2ia[itype][j];
          iintg_local = e2ilist[i] + j;
          evec->interpolate(coord,nodex,i,iintpl,3);
          fprintf(fp,"%d 3 %g %g %g %d %d %d %d %d\n"
              ,count++,coord[0],coord[1],coord[2],tag[i],numneighi2a[iintg_local]+numneighi2ia[iintg_local],numneighi2a[iintg_local],numneighi2ia[iintg_local],j);
        }
      }
    }
  } else if (singleflag == 1) {

    int iintpl = element->i2ia[element->etype[ielem]][iintg];
    int numneigha = list->numneighi2a[iintg_local];
    int numneighia = list->numneighi2ia[iintg_local];
    int *jlist;
    int *jlist_index;
    element->evec->interpolate(coord,nodex,ielem,iintpl,3);

    fprintf(fp,"%d 1 %g %g %g %d %d\n"
        ,count++,coord[0],coord[1],coord[2],id,-2);
    if (numneigha) {
      jlist = list->firstneighi2a[iintg_local];
      for (int i = 0; i < numneigha; i++) {
        j = jlist[i];
        fprintf(fp,"%d 2 %g %g %g %d %d\n"
            ,count++,atom->x[j][0],atom->x[j][1],atom->x[j][2],atom->tag[j],-1);
      }
    }
    if (numneighia) {
      ElementVec *evec = element->evec;
      jlist = list->firstneighi2ia[iintg_local];
      jlist_index = list->firstneighi2ia_index[iintg_local];
      for (int i = 0; i < numneighia; i++) {
        j = jlist[i];
        jintpl = jlist_index[i];
        evec->interpolate(coord,nodex,j,jintpl,3);
        fprintf(fp,"%d 3 %g %g %g %d %d\n"
            ,count++,coord[0],coord[1],coord[2],element->tag[j],jintpl);
      }
    }

  } else if (singleflag == 2) {
    fprintf(fp,"%d 1 %g %g %g %d %d\n"
        ,count++,atom->x[iatom][0],atom->x[iatom][1],atom->x[iatom][2],id,-2);
    int numneigha = list->numneigha2a[iatom];
    int numneighia = list->numneigha2ia[iatom];
    int *jlist;
    int *jlist_index;
    if (numneigha) {
      jlist = list->firstneigha2a[iatom];
      for (int i = 0; i < numneigha; i++) {
        j = jlist[i];
        fprintf(fp,"%d 2 %g %g %g %d %d\n"
            ,count++,atom->x[j][0],atom->x[j][1],atom->x[j][2],atom->tag[j],-1);
      }
    }
    if (numneighia) {
      ElementVec *evec = element->evec;
      jlist = list->firstneigha2ia[iatom];
      jlist_index = list->firstneigha2ia_index[iatom];
      for (int i = 0; i < numneighia; i++) {
        j = jlist[i];
        jintpl = jlist_index[i];
        evec->interpolate(coord,nodex,j,jintpl,3);
        fprintf(fp,"%d 3 %g %g %g %d %d\n"
            ,count++,coord[0],coord[1],coord[2],element->tag[j],jintpl);
      }
    }
  } else if (singleflag == 3) {
    int ietype = element->etype[ielem];
    int *element_shape_ids = element->element_shape_ids;
    int inpe = element->npe[ietype];
    for (j = 0; j < inpe; j++)
      fprintf(fp,"%g %g %g %d %d 0\n",nodex[ielem][j][0],nodex[ielem][j][1],nodex[ielem][j][2],id,ietype);

    int numneigh = list->numneighe2e[ielem];
    int *klist = list->firstneighe2e[ielem];
    for (int kk = 0; kk < numneigh; kk++) {
      int k = klist[kk];
      int ketype = element->etype[k];
      id = element->tag[k];
      int knpe = element->npe[ketype];
      for (j = 0; j < knpe; j++)
        fprintf(fp,"%g %g %g %d %d %d\n",nodex[k][j][0],nodex[k][j][1],nodex[k][j][2],id,ketype,k>=element->nlocal);
    }
    int count;

    if (element_shape_ids[ietype] == Element::QUADRILATERAL) {
      fprintf(fp,"1 2 3 4\n");
      count = 5;
    } else if (element_shape_ids[ietype] == Element::TRIANGLE) {
      fprintf(fp,"1 2 3 3\n");
      count = 4;
    } else if (element_shape_ids[ietype] == Element::HEXAHEDRON) {
      fprintf(fp,"1 2 3 4 5 6 7 8\n");
      count = 9;
    } else if (element_shape_ids[ietype] == Element::PYRAMID) {
      fprintf(fp,"1 2 3 4 5 5 5 5\n");
      count = 6;
    } else if (element_shape_ids[ietype] == Element::TETRAHEDRON) {
      fprintf(fp,"1 2 3 3 4 4 4 4\n");
      count = 5;
    } else if (element_shape_ids[ietype] == Element::WEDGE) {
      fprintf(fp,"1 2 3 3 4 5 6 6\n");
      count = 7;
    }
    for (int kk = 0; kk < numneigh; kk++) {
      int k = klist[kk];
      int ketype = element->etype[k];
      if (element_shape_ids[ketype] == Element::QUADRILATERAL) {
        fprintf(fp,"%d %d %d %d\n",count,count+1,count+2,count+3);
        count += 4;
      } else if (element_shape_ids[ketype] == Element::TRIANGLE) {
        fprintf(fp,"%d %d %d %d\n",count,count+1,count+2,count+2);
        count += 3;
      } else if (element_shape_ids[ketype] == Element::HEXAHEDRON) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+3,count+4,count+5,count+6,count+7);
        count += 8;
      } else if (element_shape_ids[ketype] == Element::PYRAMID) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+3,count+4,count+4,count+4,count+4);
        count += 5;
      } else if (element_shape_ids[ketype] == Element::TETRAHEDRON) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+2,count+3,count+3,count+3,count+3);
        count += 4;
      } else if (element_shape_ids[ketype] == Element::WEDGE) {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n"
            ,count,count+1,count+2,count+2,count+3,count+4,count+5,count+5);
        count += 6;
      }
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
    int *nintg = element->nintg;
    int *etype = element->etype;
    for (int i = 0; i < element->nlocal; i++)
      n += nintg[etype[i]];
  } else if (singleflag == 1 && iintg_local >= 0) {
    n = list->numneighi2a[iintg_local] + list->numneighi2ia[iintg_local] + 1;
  } else if (singleflag == 2 && iatom >= 0) {
    n = list->numneigha2a[iatom] + list->numneigha2ia[iatom] + 1;
  } else if (singleflag == 3 && ielem >= 0) {
    n = list->numneighe2e[ielem]+1;
    int *jlist = list->firstneighe2e[ielem];
    for (int i = 0; i < list->numneighe2e[ielem]; i++) {
      int j = jlist[i];
      nnode += element->npe[element->etype[j]];
    }
  }

  if (singleflag != 3) {
    fprintf(fp,"ITEM: TIMESTEP\n");

    fprintf(fp,"%d\n",me);

    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",n);
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
    if (singleflag == 0) {
      fprintf(fp,"ITEM: ATOMS id type x y z tag numneigh numneigha numneighia iintg\n");
    } else if (singleflag == 1) {
      fprintf(fp,"ITEM: ATOMS id type x y z tag iintpl\n");
    } else if (singleflag == 2) {
      fprintf(fp,"ITEM: ATOMS id type x y z tag iintpl\n");
    }
  } else {
    fprintf(fp,"title = \"CAC Simulation\"");
    fprintf(fp,", variables = \"x\", \"y\", \"z\", \"tag\", \"etype\", \"ghost\"\n");
    fprintf(fp,"zone t = \"Coarse Element\"");
    fprintf(fp," n = %d e = %d",nnode,n);
    if (domain->dimension == 3) 
      fprintf(fp," datapacking = point, zonetype = febrick\n");
    else 
      fprintf(fp," datapacking = point, zonetype = fequadrilateral\n");

  }
}


