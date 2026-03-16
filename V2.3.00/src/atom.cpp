#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom_vec.h"
#include "atom.h"
#include "style_atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "update.h"
#include "atom_masks.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;
using namespace MathConst;

#define DELTA 1
#define DELTA_MEMSTR 1024
#define EPSILON 1.0e-6
#define MAXBODY 20       // max # of lines in one body, also in ReadData class

enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/*---------------------------------------------------------------------------------*/

Atom::Atom(CAC *cac) : Pointers(cac)
{
  natoms = 0;
  nlocal = nghost = nmax = 0;
  ntypes = 0;

  firstgroupname = NULL;
  sortfreq = 0;
  nextsort = 0;
  userbinsize = 0.0;
  maxbin = maxnext = 0;
  binhead = NULL;
  next = permute = NULL;

  // initialize atom arrays
  // customize by adding new array

  tag = NULL;
  type = mask = NULL;
  image = NULL;
  x = v = f = NULL;

  // custom atom arrays

  nivector = ndvector = 0;
  ivector = NULL;
  dvector = NULL;
  iname = dname = NULL;
  
   // Peridynamic scale factor

  pdscale = 1.0;

  // ntype-length arrays

  mass = NULL;
  mass_setflag = NULL;

  // callback lists & extra restart info
 
  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = NULL;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = NULL;

  //default atom ID and mapping values

  tag_enable = 1;
  map_style = map_user = 0;
  map_tag_max = -1;
  map_maxarray = map_nhash = -1;

  max_same = 0;
  sametag = NULL;
  map_array = NULL;
  map_bucket = NULL;
  map_hash = NULL;

  atom_style = NULL;
  avec = NULL;

  avec_map = new AtomVecCreatorMap();

#define ATOM_CLASS
#define AtomStyle(key,Class) \
  (*avec_map)[#key] = &avec_creator<Class>;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
}

/*----------------------------------------------------------------------------------*/

Atom::~Atom()
{

  delete [] atom_style;
  delete avec;

  delete [] firstgroupname;
  memory->destroy(binhead);
  memory->destroy(next);
  memory->destroy(permute);

  // delete atom arrays
  // customize by adding new arrays
  
  memory->destroy(tag);
  memory->destroy(type);
  memory->destroy(mask);
  memory->destroy(image);
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);

  // delete custom atom arrays

  for (int i = 0; i < nivector; i++) {
    delete [] iname[i];
    memory->destroy(ivector[i]);
  }
  for (int i = 0; i < ndvector; i++) {
    delete [] dname[i];
    memory->destroy(dvector[i]);
  }

  memory->sfree(iname);
  memory->sfree(dname);
  memory->sfree(ivector);
  memory->sfree(dvector);

  // delete per-type arrays

  delete [] mass;
  delete [] mass_setflag;

  // delete extra arrays

  memory->destroy(extra_grow);
  memory->destroy(extra_restart);
  memory->destroy(extra_border);
  memory->destroy(extra);
  
  // delete mapping data structures

  map_delete();
}

/* ----------------------------------------------------------------------
     create an AtomVec style
     called from cac.cpp, input script
  ------------------------------------------------------------------------- */

void Atom::create_avec(const char *style, int narg, char **arg)
{

  delete [] atom_style;
  if (avec) delete avec;
  atom_style = NULL;
  avec = NULL;

  // create instance of AtomVec
  // use grow() to initialize atom-based arrays to length 1
  // so that x[0][0] can always be referenced even if proc has no atoms
  
  avec = new_avec(style);
  avec->store_args(narg,arg);
  avec->process_args(narg,arg); 
  avec->grow(1); 
}

/* ----------------------------------------------------------------------
   generate an AtomVec class, first with suffix appended
  ------------------------------------------------------------------------- */

AtomVec *Atom::new_avec(const char *style)
{
  if (avec_map->find(style) != avec_map->end()) {
    AtomVecCreator avec_creator = (*avec_map)[style];
    return avec_creator(cac);
  } else error->all(FLERR,"Unknown atom style");
  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per AtomVec style in style_atom.h
------------------------------------------------------------------------- */

template <typename T>
AtomVec *Atom::avec_creator(CAC *cac)
{
  return new T(cac);
}


/* ----------------------------------------------------------------------
   decrement ptrs in callback lists to fixes beyond the deleted ifix
   happens after fix is deleted
------------------------------------------------------------------------- */

void Atom::update_callback(int ifix)
{
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix) extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix) extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix) extra_border[i]--;
}

/* ----------------------------------------------------------------------
  allocate arrays of length ntypes
  only done after ntypes is set
 ------------------------------------------------------------------------- */

void Atom::allocate_type_arrays()
{
  if (avec->mass_type) {
    mass = new double[ntypes+1];
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
}

/* ----------------------------------------------------------------------
   unpack n lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_atoms(int n, char *buf, tagint id_offset, int type_offset, 
                      int shiftflag, double *shift)
{
  int m,xptr,iptr;
  imageint imagedata;
  double xdata[3],lambda[3];
  double *coord;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_atom && nwords != avec->size_data_atom + 3)
    error->all(FLERR,"Incorrect atom format in data file");
  
  char **values = new char*[nwords];

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }
  
  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lambda[0]; subhi[0] = domain->subhi_lambda[0];
    sublo[1] = domain->sublo_lambda[1]; subhi[1] = domain->subhi_lambda[1];
    sublo[2] = domain->sublo_lambda[2]; subhi[2] = domain->subhi_lambda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }
  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // xptr = which word in line starts xyz coords
  // iptr = which word in line starts ix,iy,iz image flags

  xptr = avec->xcol_data - 1;
  int imageflag = 0;
  if (nwords > avec->size_data_atom) imageflag = 1;
  if (imageflag) iptr = nwords - 3;

  // loop over lines of atom data
  // tokenize the line into values
  // extract xyz coords and image flags
  // remap atom into simulation box

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR,"Incorrect atom format in data file");
    for (m = 1; m < nwords; m++) {
      values[m] = strtok(NULL," \t\n\r\f");
      if (values[m] == NULL)
        error->all(FLERR,"Incorrect atom format in data file");
    }

    if (imageflag)
      imagedata = ((imageint) (atoi(values[iptr]) + IMGMAX) & IMGMASK) |
        (((imageint) (atoi(values[iptr+1]) + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (atoi(values[iptr+2]) + IMGMAX) & IMGMASK) << IMG2BITS);
    else imagedata = ((imageint) IMGMAX << IMG2BITS) |
           ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    xdata[0] = atof(values[xptr]);
    xdata[1] = atof(values[xptr+1]);
    xdata[2] = atof(values[xptr+2]);
    if (shiftflag) {
      xdata[0] += shift[0];
      xdata[1] += shift[1];
      xdata[2] += shift[2];
    }

    domain->remap(xdata,imagedata);
    if (triclinic) {
      domain->x2lambda(xdata,lambda);
      coord = lambda;
    } else coord = xdata;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
      avec->data_atom(xdata,imagedata,values);
      if (id_offset) tag[nlocal-1] += id_offset;
      if (type_offset) {
        type[nlocal-1] += type_offset;
        if (type[nlocal-1] > ntypes)
          error->one(FLERR,"Invalid atom type in Atoms section of data file");
      }
    }
   
    buf = next + 1;
  }

  delete [] values;
}

/* ----------------------------------------------------------------------
   unpack N lines from Atom Velocity section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Atom::data_vels(int n, char *buf, tagint id_offset)
{
  int j,m;
  tagint tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = universe->count_words(buf);
  *next = '\n';

  if (nwords != avec->size_data_vel)
    error->all(FLERR,"Incorrect velocity format in data file");

  char **values = new char*[nwords];

  // loop over lines of atom velocities
  // tokenize the line into values
  // if I own atom tag, unpack its values

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = ATOTAGINT(values[0]) + id_offset;
    if (tagdata <= 0 || tagdata > map_tag_max)
      error->one(FLERR,"Invalid atom ID in Velocities section of data file");
    if ((m = map(tagdata)) >= 0) avec->data_vel(m,&values[1]);

    buf = next + 1;
  }

  delete [] values;
}
/* ----------------------------------------------------------------------
    set a mass and flag it as set
   called from reading of data file
   type_offset may be used when reading multiple data files
------------------------------------------------------------------------- */

void Atom::set_mass(const char *str, int type_offset)
{
  if (mass == NULL) error->all(FLERR,"Cannot set mass for this atom style");

  int itype;
  double mass_one;
  int n = sscanf(str,"%d %lg",&itype,&mass_one);
  if (n != 2) error->all(FLERR,"Invalid mass line in data file");
  itype += type_offset;

  if (itype < 1 || itype > ntypes)
    error->all(FLERR,"Invalid type for mass set");

  mass[itype] = mass_one;
  mass_setflag[itype] = 1;

  if (mass[itype] <= 0.0) error->all(FLERR,"Invalid mass value");

}

/*-----------------------------------------------------------------------
  set a mass and flag it as set
  called from EAM pair routine
-------------------------------------------------------------------------*/

void Atom::set_mass(int itype, double value)
{
  if (mass == NULL) error->all(FLERR,"Cannot set mass for this atom style");
  if (itype < 1 || itype > ntypes)
    error->all(FLERR,"Invalid type for mass set");

  mass[itype] = value;
  mass_setflag[itype] = 1;

  if (mass[itype] <=0.0) error->all(FLERR,"Invalid mass value");
}

/* ---------------------------------------------------------------------- */

void Atom::init()
{
  // delete extra array since it doesn't persist past first run

  if (nextra_store) {
    memory->destroy(extra);
    extra = NULL;
    nextra_store = 0;
  }

  // check arrays that are atom type in length

  check_mass();

  // setup of firstgroup

  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all(FLERR,"Could not find atom_modify first group ID");
  } else firstgroup = -1;

  // init AtomVec
  
  avec->init();

}

/* ----------------------------------------------------------------------
   check that all masses have been set
------------------------------------------------------------------------- */

void Atom::check_mass()
{
  if (mass == NULL) return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0) error->all(FLERR,"All masses are not set");
}

/* ---------------------------------------------------------------------- */

void Atom::setup()
{
  // setup bins for sorting
  // cannot do this in init() because uses neighbor cutoff
  if (sortfreq > 0) setup_sort_bins();
}

/* ----------------------------------------------------------------------
   setup bins for spatial sorting of atoms
------------------------------------------------------------------------- */

void Atom::setup_sort_bins()
{
  // binsize:
  // user setting if explicitly set
  // 1/2 of neighbor cutoff for non-CUDA
  // CUDA_CHUNK atoms/proc for CUDA
  // check if neighbor cutoff = 0.0

  double binsize;
  if (userbinsize > 0.0) binsize = userbinsize;
  else binsize = 0.5 * neighbor->cutneighmax;

  if (binsize == 0.0) error->all(FLERR,"Atom sorting has bin size = 0.0");

  double bininv = 1.0/binsize;

  // nbin xyz = local bins
  // bbox lo/hi = bounding box of my sub-domain

  bboxlo[0] = domain->sublo[0];
  bboxlo[1] = domain->sublo[1];
  bboxlo[2] = domain->sublo[2];
  bboxhi[0] = domain->subhi[0];
  bboxhi[1] = domain->subhi[1];
  bboxhi[2] = domain->subhi[2];

  nbinx = static_cast<int> ((bboxhi[0]-bboxlo[0]) * bininv);
  nbiny = static_cast<int> ((bboxhi[1]-bboxlo[1]) * bininv);
  nbinz = static_cast<int> ((bboxhi[2]-bboxlo[2]) * bininv);
  if (domain->dimension == 2) nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  bininvx = nbinx / (bboxhi[0]-bboxlo[0]);
  bininvy = nbiny / (bboxhi[1]-bboxlo[1]);
  bininvz = nbinz / (bboxhi[2]-bboxlo[2]);

  if (1.0*nbinx*nbiny*nbinz > INT_MAX)
    error->one(FLERR,"Too many atom sorting bins");

  nbins = nbinx*nbiny*nbinz;

  // reallocate per-bin memory if needed
  
  if (nbins > maxbin) {
    memory->destroy(binhead);
    maxbin = nbins;
    memory->create(binhead,maxbin,"atom:binhead");
  }
}

/* ----------------------------------------------------------------------
   perform spatial sort of atoms within my sub-domain
   always called between comm->exchange() and comm->borders()
   don't have to worry about clearing/setting atom->map since done in comm
------------------------------------------------------------------------- */

void Atom::sort()
{
  int i,m,n,ix,iy,iz,ibin,empty;

  // set next timestep for sorting to take place
   
  nextsort = (update->ntimestep/sortfreq)*sortfreq + sortfreq;
 
  // re-setup sort bins if needed

  if (nbins == 1) return;

  // reallocate per-atom vectors if needed

  if (nlocal > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next,maxnext,"atom:next");
    memory->create(permute,maxnext,"atom:permute");
  } 

  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  for (i = nlocal-1; i >= 0; i--) {
    ix = static_cast<int> ((x[i][0]-bboxlo[0])*bininvx);
    iy = static_cast<int> ((x[i][1]-bboxlo[1])*bininvy);
    iz = static_cast<int> ((x[i][2]-bboxlo[2])*bininvz);
    ix = MAX(ix,0);
    iy = MAX(iy,0);
    iz = MAX(iz,0);
    ix = MIN(ix,nbinx-1);
    iy = MIN(iy,nbiny-1);
    iz = MIN(iz,nbinz-1);
    ibin = iz*nbiny*nbinx + iy*nbinx + ix;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  // permute = desired permutation of atoms
  // permute[I] = J means Ith new atom will be Jth old atom

  n = 0;
  for (m = 0; m < nbins; m++) {
    i = binhead[m];
    while (i >= 0) {
      permute[n++] = i;
      i = next[i];
    }
  }

  // current = current permutation, just reuse next vector
  // current[I] = J means Ith current atom is Jth old atom

  int *current = next;
  for (i = 0; i < nlocal; i++) current[i] = i;

  // reorder local atom list, when done, current = permute
  // perform "in place" using copy() to extra atom location at end of list
  // inner while loop processes one cycle of the permutation
  // copy before inner-loop moves an atom to end of atom list
  // copy after inner-loop moves atom at end of list back into list
  // empty = location in atom list that is currently empty

  for (i = 0; i < nlocal; i++) {
    if (current[i] == permute[i]) continue;
    avec->copy(i,nlocal,0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty],empty,0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal,empty,0);
    current[empty] = permute[empty];
  }  
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   call to avec tallies per-atom vectors
   add in global to local mapping storage
------------------------------------------------------------------------- */

bigint Atom::memory_usage()
{
  memlength = DELTA_MEMSTR;
  memory->create(memstr,memlength,"atom:memstr");
  memstr[0] = '\0';
  bigint bytes = avec->memory_usage();
  memory->destroy(memstr);

  bytes += max_same*sizeof(int);
  if (map_style == 1)
    bytes += memory->usage(map_array,map_maxarray);
  else if (map_style == 2) {
    bytes += map_nbucket*sizeof(int);
    bytes += map_nhash*sizeof(HashElem);
  }
  if (maxnext) {
    bytes += memory->usage(next,maxnext);
    bytes += memory->usage(permute,maxnext);
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   accumulate per-atom vec names in memstr, padded by spaces
   return 1 if padded str is not already in memlist, else 0
------------------------------------------------------------------------- */

int Atom::memcheck(const char *str)
{
  int n = strlen(str) + 3;
  char *padded = new char[n];
  strcpy(padded," ");
  strcat(padded,str);
  strcat(padded," ");
  
  if (strstr(memstr,padded)) {
    delete [] padded;
    return 0;
  }

  if (strlen(memstr) + n >= memlength) {
    memlength += DELTA_MEMSTR;
    memory->grow(memstr,memlength,"atom:memstr");
  }

  strcat(memstr,padded);
  delete [] padded;
  return 1;
}

/* ----------------------------------------------------------------------
   add unique tags to any atoms with tag = 0
   new tags are grouped by proc and start after max current tag
   called after creating new atoms
   error if new tags will exceed MAXTAGINT
------------------------------------------------------------------------- */

void Atom::tag_extend()
{
  // maxtag_all = max tag for all atoms

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag,tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_CAC_TAGINT,MPI_MAX,world);

  // notag = # of atoms I own with no tag (tag = 0)
  // notag_sum = # of total atoms on procs <= me with no tag

  bigint notag = 0;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) notag++;

  bigint notag_total;
  MPI_Allreduce(&notag,&notag_total,1,MPI_CAC_BIGINT,MPI_SUM,world);
  if (notag_total >= MAXTAGINT)
    error->all(FLERR,"New atom IDs exceed maximum allowed ID");

  bigint notag_sum;
  MPI_Scan(&notag,&notag_sum,1,MPI_CAC_BIGINT,MPI_SUM,world);

  // itag = 1st new tag that my untagged atoms should use

  tagint itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++) if (tag[i] == 0) tag[i] = itag++;

}
