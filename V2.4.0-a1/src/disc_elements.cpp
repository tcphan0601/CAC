#include <mpi.h>
#include <string.h>
#include "disc_elements.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "irregular_comm.h"
#include "universe.h"

#include <map>

using namespace CAC_NS;
/* ---------------------------------------------------------------------- */

DiscElements::DiscElements(CAC *cac) : Pointers(cac)
{
	MPI_Comm_rank(world,&me);
}

/* ----------------------------------------------------------------------
   Called as disc_elements command in input script.
   New atoms created will have ID starting immediately
   following the largest atom ID existing before the
   disc_elements command was evoked.
   Call IrregularComm to migrate atoms since atoms might move far away
   ------------------------------------------------------------------------- */

void DiscElements::command(int narg, char **arg)
{
	disc_style = 0;
	new_etype = 0;
	if (domain->box_exist == 0)
		error->all(FLERR,"Disc_elements command before simulation box is defined");

	if (element->nelements == 0) {
		if (me == 0) error->warning(FLERR,"No element in simulation box");
		return;
	}

	if (narg < 1) error->all(FLERR,"Illegal dis_elements command");

	// flag elements for discretizing

	if (strcmp(arg[0],"region") == 0) list_disc_region(narg,arg);
	else if (strcmp(arg[0],"group") == 0) list_disc_group(narg,arg);
	else error->all(FLERR,"Illegal disc_elements command");

	element->init();


	ElementVec *evec = element->evec;

	// discretize elements in list

	if (disc_style == 0) {
		atom->init();
		int nnewatoms = 0;
		int nlocal_previous = atom->nlocal;
		int ndel = 0;
		int i = 0;
		while (i < element->nlocal)
			if (disc_flag_list[i]) {
				ndel++;
				nnewatoms += evec->element2atom(i);
				disc_flag_list[i] = disc_flag_list[element->nlocal];
			} else i++;

		// init per-atom fix/compute/variable values for created atoms

		atom->data_fix_compute_variable(nlocal_previous,atom->nlocal);

		// assign tag to new atoms

		if (atom->tag_enable) atom->tag_extend();
		atom->tag_check();

		// update new global natoms and nelements

		int totalnewatoms;
		int totaldel;
		MPI_Allreduce(&nnewatoms,&totalnewatoms,1,MPI_INT,MPI_SUM,world);
		MPI_Allreduce(&ndel,&totaldel,1,MPI_INT,MPI_SUM,world);
		atom->natoms += totalnewatoms;
		element->nelements -= totaldel;
		element->nintpls -= totalnewatoms;

		FILE *out;
		if (me == 0)
			for (int m = 0; m < 2; m++) {
				if (m == 0) out = logfile;
				else out = screen;
				if (out) {
					fprintf(out,"Discretized %d elements into %d atoms\n",totaldel,totalnewatoms);
					fprintf(out,"  New total atoms = %ld\n",atom->natoms);
					fprintf(out,"  New total elements = %ld\n",element->nelements);
				}
			}

	} else {
		int nnewelems = 0;
		int nlocal_previous = element->nlocal;
		int ndel = 0;

		//if (element->apc == 1) {

		// old elements have the same local index after discretize so just use for loop

		for (int i = 0; i < nlocal_previous; i++)
			if (disc_flag_list[i]) {
				ndel++;
				nnewelems += evec->element2element(i,new_etype,0);
			}

		// assign tags for created elements

		if (element->tag_enable) element->tag_extend();
		element->tag_check();

		//} else {
		//  int mapflag = 0;
		//  int apc = element->apc;
		//  if (element->map_style == 0) {
		//    mapflag = 1;
		//    element->map_init();
		//    element->map_set();
		//  }

		//  tagint maxtag = 0;
		//  for (int i = 0; i < element->nlocal; i++) maxtag = MAX(maxtag,element->tag[i]);
		//  tagint maxtag_all;
		//  MPI_Reduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,0,world);

		// do this in serial (not very time consuming so it's fine)
		// start from proc 0 until last proc
		// receive tag_offset from me-1 (if me > 0), discretize,
		// send new tag_offset to me+1 (if me < nprocs-1)

		// old elements have the same local index after discretize so just use for loop
		// split all elements I own in clusters

		//  tagint tag_offset;

		//  if (me) MPI_Recv(&tag_offset,1,MPI_CAC_TAGINT,me-1,0,world,MPI_STATUS_IGNORE);
		//  else tag_offset = maxtag_all + 1;
		//  for (int i = 0; i < nlocal_previous; i++)
		//    if (disc_flag_list[i]) {
		//      tagint itag = (element->tag[i] - 1)/apc*apc + 1;
		//      int n = 0;
		//      for (int k = 0; k < element->apc; k++) {
		//        int j = element->map(itag+k);
		//        if (j < 0 || j >= nlocal_previous) error->one(FLERR,"All elements in a clusters must be owned by me, try splitting elements right after creating elements");
		//        n += evec->element2element(j,new_etype,tag_offset+k);
		//        ndel++;
		//        disc_flag_list[j] = 0;
		//      }

		//      tag_offset += n;
		//      nnewelems += n;
		//    }
		//  if (me < comm->nprocs-1) MPI_Send(&tag_offset,1,MPI_CAC_TAGINT,me+1,0,world);
		//  if (mapflag) {
		//    element->map_delete();
		//    element->map_style = 0;
		//  }
		//  element->check_node_coords(0);
		//}

		// update new global nelements
		// total nintpls remains constant so no need to update

		int totalnewelems;
		int totaldel;
		MPI_Allreduce(&nnewelems,&totalnewelems,1,MPI_INT,MPI_SUM,world);
		MPI_Allreduce(&ndel,&totaldel,1,MPI_INT,MPI_SUM,world);
		element->nelements += totalnewelems - totaldel;

		FILE *out;
		if (me == 0)
			for (int m = 0; m < 2; m++) {
				if (m == 0) out = logfile;
				else out = screen;
				if (out) {
					fprintf(out,"Discretized %d elements into %d smaller elements of type %d\n",totaldel,totalnewelems,new_etype);
					fprintf(out,"  New total elements = %ld\n",element->nelements);
				}
			}
	}

	// if compress flag set
	// reset element tags to be continuous
	// set all element IDs to 0, call tag_extend()

	if (compress_flag) {
		tagint *tag = element->tag;
		int nlocal = element->nlocal;
		for (int i = 0; i < nlocal; i++) tag[i] = 0;
		element->tag_extend();
	}

	// Invoke new irregular comm, migrate atoms/elements and destroy it
	// must be done after updating new global natoms and nelements

	double **ax = atom->x;
	imageint *aimage = atom->image;
	int nalocal = atom->nlocal;
	for (int i = 0; i < nalocal; i++) domain->remap(ax[i],aimage[i]);

	double **ex = element->x;
	double ***nodex = element->nodex;
	imageint *eimage = element->image;
	int nelocal = element->nlocal;
	int *npe = element->npe;
	int *etype = element->etype;
	for (int i = 0; i < nelocal; i++) domain->remap(ex[i],nodex[i],eimage[i],npe[etype[i]]);

	if (migrate_flag) {
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
	// reset element->map and atom->map if it exists
	// set nghost to 0 so old ghosts of deleted elements and new atoms won't be mapped

	if (element->map_style) {
		element->nghost = 0;
		element->map_init();
		element->map_set();
	}

	if (atom->map_style) {
		atom->nghost = 0;
		atom->map_init();
		atom->map_set();
	}

	// clean up

	memory->destroy(disc_flag_list);
}

/* ----------------------------------------------------------------------
   flag all elements in group to be discretize, group will still exist
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_group(int narg, char **arg)
{
	if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

	int igroup = group->find(arg[1]);
	if (igroup == -1) error->all(FLERR,"Could not find disc_elements group ID");

	options(narg-2,&arg[2]);

	// allocate discretization list

	int nlocal = element->nlocal;
	memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");


	int *mask = element->mask;
	int groupbit = group->bitmask[igroup];

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) disc_flag_list[i] = 1;
		else disc_flag_list[i] = 0;
	}
}

/* ----------------------------------------------------------------------
   flag all elements in region to discretize
   ------------------------------------------------------------------------- */

void DiscElements::list_disc_region(int narg, char **arg)
{
	if (narg < 2) error->all(FLERR,"Illegal disc_elements command");

	int iregion = domain->find_region(arg[1]);
	if (iregion == -1) error->all(FLERR,"Could not find disc_elements region ID");
	domain->regions[iregion]->prematch();

	options(narg-2,&arg[2]);

	// allocate discretization list

	int nlocal = element->nlocal;
	memory->create(disc_flag_list,nlocal,"disc_elements:disc_flag_list");

	double **x = element->x;

	int *apc = element->apc;
//  if (apc == 1) {
	for (int i = 0; i < nlocal; i++) {
		if (domain->regions[iregion]->match(x[i])) disc_flag_list[i] = 1;
		else disc_flag_list[i] = 0;
	}
//  } else {
//    int mapflag = 0;
//    if (element->map_style == 0) {
//      mapflag = 1;
//      element->map_init();
//      element->map_set();
//    }
//
//    tagint itag;
//    tagint *tag = element->tag;
//    for (int i = 0; i < nlocal; i++) {
//      disc_flag_list[i] = 0;
//      itag = (tag[i]-1)/apc*apc+1;
//      for (int j = 0; j < apc; j++) {
//        int ii = element->map(itag+j);
//        if (ii < 0) error->one(FLERR,"Ghost element not correct for link");
//        if (domain->regions[iregion]->match(x[ii]))
//          disc_flag_list[i] = 1;
//      }
//    }
//
//    if (mapflag) {
//      element->map_delete();
//      element->map_style = 0;
//    }
//
//  }
}

/* ----------------------------------------------------------------------
   process command options
   ------------------------------------------------------------------------- */

void DiscElements::options(int narg, char **arg)
{
	compress_flag = 1;
	migrate_flag = 1;
	int iarg = 0;
	while (iarg < narg) {
		if (strcmp(arg[iarg],"compress") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
			if (strcmp(arg[iarg+1],"yes") == 0) compress_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) compress_flag = 0;
			else error->all(FLERR,"Illegal disc_elements command");

			iarg += 2;
		} else if (strcmp(arg[iarg],"migrate") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
			if (strcmp(arg[iarg+1],"yes") == 0) migrate_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) migrate_flag = 0;
			else error->all(FLERR,"Illegal disc_elements command");

			iarg += 2;

		} else if (strcmp(arg[iarg],"style") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal disc_elements command");
			if (strcmp(arg[iarg+1],"atom") == 0) {
				disc_style = 0;
				iarg += 2;
			} else if (strcmp(arg[iarg+1],"element") == 0) {
				if (iarg+3 > narg) error->all(FLERR,"Illegal disc_elements command");
				new_etype = universe->inumeric(FLERR,arg[iarg+2]);
				if (new_etype <= 0 || new_etype > element->netypes)
					error->all(FLERR,"Invalid element type in disc_elements command");
				disc_style = 1;
				iarg += 3;
			} else error->all(FLERR,"Illegal disc_elements command");
		} else error->all(FLERR,"Illegal disc_elements command");
	}
	// if more than one atom per node cluster, no compress element tags

	//if (element->apc != 1) compress_flag = 0;
}
