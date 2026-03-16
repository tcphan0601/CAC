#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "style_nbin.h"
#include "style_nstencil.h"
#include "style_npair.h"
#include "atom.h"
#include "element.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "universe.h"
#include "pair.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "update.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace CAC_NS;
using namespace NeighConst;

#define RQDELTA 1
#define EXDELTA 1

#define BIG 1.0e20

enum{NSQ,BIN,MULTI};     // also in NBin, NeighList, NStencil
enum{NONE,ALL,PARTIAL,TEMPLATE};
enum{ATOM,ELEMENT};

//#define NEIGH_LIST_DEBUG 1

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor(CAC *cac) : Pointers(cac),
  pairclass(NULL), pairnames(NULL), pairmasks(NULL)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  firsttime = 1;

  style = BIN;
  every = 1;
  delay = 10;
  dist_check = 1;
  pgsize = 100000;
  oneatom = 2000;
  binsizeflag = 0;
  build_once = 0;
  cluster_check = 0;
  ago = -1;

  cutneighmax = 0.0;
  cutneighsq = NULL;
  //cutneighghostsq = NULL;
  cuttype = NULL;
  cuttypesq = NULL;
  fixchecklist = NULL;

  // pairwise neighbor lists and associated data structs

  nlist = 0;
  lists = NULL;

  nbin = 0;
  neigh_bin = NULL;

  nstencil = 0;
  neigh_stencil = NULL;

  neigh_pair = NULL;

  nstencil_perpetual = 0;
  slist = NULL;

  npair_perpetual = 0;
  plist = NULL;

  nrequest = maxrequest = 0;
  requests = NULL;

  old_nrequest = 0;
  old_requests = NULL;

  old_style = style;
  old_triclinic = 0;
  old_pgsize = pgsize;
  old_oneatom = oneatom;

  zeroes = NULL;

  binclass = NULL;
  binnames = NULL;
  binmasks = NULL;
  stencilclass = NULL;
  stencilnames = NULL;
  stencilmasks = NULL;

  // coords at last neighboring

  atommaxhold = 0;
  nodemaxhold = 0;
  atomxhold = NULL;
  nodexhold = NULL;
  lastcall = -1;
  last_setup_bins = -1;
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{

  if (copymode) return;

  memory->destroy(cutneighsq);
  //memory->destroy(cutneighghostsq);
  delete [] cuttype;
  delete [] cuttypesq;
  delete [] fixchecklist;

  for (int i = 0; i < nlist; i++) delete lists[i];
  for (int i = 0; i < nbin; i++) delete neigh_bin[i];
  for (int i = 0; i < nstencil; i++) delete neigh_stencil[i];
  for (int i = 0; i < nlist; i++) delete neigh_pair[i];

  delete [] lists;
  delete [] neigh_bin;
  delete [] neigh_stencil;
  delete [] neigh_pair;

  delete [] slist;
  delete [] plist;

  for (int i = 0; i < nrequest; i++) 
    if (requests[i]) delete requests[i];
  memory->sfree(requests);
  for (int i = 0; i < old_nrequest; i++)
    if (old_requests[i]) delete old_requests[i];
  memory->sfree(old_requests);

  delete [] zeroes;

  delete [] binclass;
  delete [] binnames;
  delete [] binmasks;
  delete [] stencilclass;
  delete [] stencilnames;
  delete [] stencilmasks;
  delete [] pairclass;
  delete [] pairnames;
  delete [] pairmasks;

  memory->destroy(atomxhold);
  memory->destroy(nodexhold);
}

/* ----------------------------------------------------------------------
   update cut off due to changed element size
   called in verlet
   ------------------------------------------------------------------------- */

void Neighbor::update_cut()
{
  cut_all = cutneighmax + element->max_diag_size;
  cut_elem = cutneighmax + element->max_diag_size/2.0;
  cut_subelem = cutneighmax + element->max_diag_size/(2.0*element->esplit);
  cut_allsq = cut_all*cut_all;
  cut_elemsq = cut_elem*cut_elem;
  cut_subelemsq = cut_subelem*cut_subelem;
}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,n;

  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
  newton_pair = force->newton_pair;

  // error check

  if (delay > 0 && (delay % every) != 0)
    error->all(FLERR,"Neighbor delay must be 0 or multiple of every setting");

  if (pgsize < 10*oneatom)
    error->all(FLERR,"Neighbor page size must be >= 10x the one atom setting");

  // ------------------------------------------------------------------
  // settings

  // bbox lo/hi ptrs = bounding box of entire domain, stored by Domain

  if (triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  // set neighbor cutoffs (force cutoff + skin)
  // trigger determines when atoms migrate and neighbor lists are rebuilt
  //   needs to be non-zero for migration distance check
  //   even if pair = NULL and no neighbor lists are used
  // cutneigh = force cutoff + skin if cutforce > 0, else cutneigh = 0
  // cutneighghost = pair cutghost if it requests it, else same as cutneigh

  triggersq = 0.25*skin*skin;
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
        (dimension == 3 && domain->zperiodic)))
    boxcheck = 1;

  n = atom->ntypes;

  if (cutneighsq == NULL) {

    memory->create(cutneighsq,n+1,n+1,"neigh:cutneighsq");
    //memory->create(cutneighghostsq,n+1,n+1,"neigh:cutneighghostsq");
    cuttype = new double[n+1];
    cuttypesq = new double[n+1];
  }

  double cutoff,delta,cut;
  cutneighmin = BIG;
  cutneighmax = 0.0;

  for (i = 1; i <= n; i++) {
    cuttype[i] = cuttypesq[i] = 0.0;
    for (j = 1; j <= n; j++) {
      if (force->pair) cutoff = sqrt(force->pair->cutsq[i][j]);
      else cutoff = 0.0;
      if (cutoff > 0.0) delta = skin;
      else delta = 0.0;
      cut = cutoff + delta;

      cutneighsq[i][j] = cut*cut;
      cuttype[i] = MAX(cuttype[i],cut);
      cuttypesq[i] = MAX(cuttypesq[i],cut*cut);
      cutneighmin = MIN(cutneighmin,cut);
      cutneighmax = MAX(cutneighmax,cut);

      //if (force->pair && force->pair->ghostneigh) {
      //  cut = force->pair->cutghost[i][j] + skin;
      //  cutneighghostsq[i][j] = cut*cut;
      //} else cutneighghostsq[i][j] = cut*cut;
    }
  }
  cutneighmaxsq = cutneighmax * cutneighmax;

  // set cutoffs used for finding element neighbor

  update_cut();  

  // fixchecklist = other classes that can induce reneighboring in decide()

//  restart_check = 0;
//  if (output->restart_flag) restart_check = 1;

  delete [] fixchecklist;
  fixchecklist = NULL;
  fixchecklist = new int[modify->nfix];

  fix_check = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;

  must_check = 0;
//  if (restart_check || fix_check) must_check = 1;
  if (fix_check) must_check = 1;

  // set special_flag for 1-2, 1-3, 1-4 neighbors
  // flag[0] is not used, flag[1] = 1-2, flag[2] = 1-3, flag[3] = 1-4
  // flag = 0 if both LJ/Coulomb special values are 0.0
  // flag = 1 if both LJ/Coulomb special values are 1.0
  // flag = 2 otherwise or if KSpace solver is enabled
  // pairwise portion of KSpace solver uses all 1-2,1-3,1-4 neighbors
  // or selected Coulomb-approixmation pair styles require it

  if (force->special_lj[1] == 0.0 && force->special_coul[1] == 0.0)
    special_flag[1] = 0;
  else if (force->special_lj[1] == 1.0 && force->special_coul[1] == 1.0)
    special_flag[1] = 1;
  else special_flag[1] = 2;

  if (force->special_lj[2] == 0.0 && force->special_coul[2] == 0.0)
    special_flag[2] = 0;
  else if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0)
    special_flag[2] = 1;
  else special_flag[2] = 2;

  if (force->special_lj[3] == 0.0 && force->special_coul[3] == 0.0)
    special_flag[3] = 0;
  else if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0)
    special_flag[3] = 1;
  else special_flag[3] = 2;

  special_flag[1] = special_flag[2] = special_flag[3] = 2;

  // maxwt = max multiplicative factor on atom indices stored in neigh list

  maxwt = 0;
  if (special_flag[1] == 2) maxwt = 2;
  if (special_flag[2] == 2) maxwt = 3;
  if (special_flag[3] == 2) maxwt = 4;

  // ------------------------------------------------------------------
  // xhold arrays

  // free if not needed for this run

  if (dist_check == 0) {

    memory->destroy(atomxhold);
    atommaxhold = 0;
    atomxhold = NULL;
    memory->destroy(nodexhold);
    nodemaxhold = 0;
    nodexhold = NULL;
  }

  // first time allocation

  if (dist_check) {
    if (atommaxhold == 0) {
      atommaxhold = atom->nmax;
      memory->create(atomxhold,atommaxhold,3,"neigh:atomxhold");
    }
    
    if (nodemaxhold == 0) {
      nodemaxhold = element->nmax;
      memory->create(nodexhold,nodemaxhold,element->npe,3,"neigh:nodexhold");
    }

  }

  // create pairwise lists
  // one-time call to init_styles() to scan style files and setup
  // init_pair() creates auxiliary classes: NBin, NStencil, NPair

  if (firsttime) init_styles();
  firsttime = 0;

  int same = init_pair();

  // invoke copy_neighbor_info() in Bin,Stencil,Pair classes
  // copied once per run in case any cutoff, exclusion, special info changed

  for (i = 0; i < nbin; i++)
    neigh_bin[i]->copy_neighbor_info();

  for (i = 0; i < nstencil; i++)
    neigh_stencil[i]->copy_neighbor_info();

  for (i = 0; i < nlist; i++)
    if (neigh_pair[i]) neigh_pair[i]->copy_neighbor_info();

  if (!same && comm->me == 0) print_pairwise_info();

  // can now delete requests so next run can make new ones
  // print_pairwise_info() made use of requests
  // set of NeighLists now stores all needed info

  for (int i = 0; i < nrequest; i++) {
    delete requests[i];
    requests[i] = NULL;
  }
  nrequest = 0;

}

/* ----------------------------------------------------------------------
   create and initialize lists of Nbin, Nstencil, NPair classes
   lists have info on all classes in 3 style*.h files
   cannot do this in constructor, b/c too early to instantiate classes
   ------------------------------------------------------------------------- */

void Neighbor::init_styles()
{
  // extract info from NBin classes listed in style_nbin.h

  nbclass = 0;

#define NBIN_CLASS
#define NBinStyle(key,Class,bitmasks) nbclass++;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS

  binclass = new BinCreator[nbclass]; 
  binnames = new char*[nbclass];
  binmasks = new int[nbclass];
  nbclass = 0;

#define NBIN_CLASS
#define NBinStyle(key,Class,bitmasks) \
  binnames[nbclass] = (char *) #key; \
  binclass[nbclass] = &bin_creator<Class>; \
  binmasks[nbclass++] = bitmasks;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS

  // extract info from NStencil classes listed in style_nstencil.h

  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) nsclass++;
#include "style_nstencil.h"
#undef NStencilStyle
#undef NSTENCIL_CLASS

  stencilclass = new StencilCreator[nsclass]; 
  stencilnames = new char*[nsclass];
  stencilmasks = new int[nsclass];
  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) \
  stencilnames[nsclass] = (char *) #key; \
  stencilclass[nsclass] = &stencil_creator<Class>; \
  stencilmasks[nsclass++] = bitmasks;
#include "style_nstencil.h"
#undef NStencilStyle
#undef NSTENCIL_CLASS

  // extract info from NPair classes listed in style_npair.h

  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) npclass++;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS

  pairclass = new PairCreator[npclass]; 
  pairnames = new char*[npclass];
  pairmasks = new int[npclass];
  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) \
  pairnames[npclass] = (char *) #key; \
  pairclass[npclass] = &pair_creator<Class>; \
  pairmasks[npclass++] = bitmasks;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS
}

/* ----------------------------------------------------------------------
   create and initialize NPair classes
   ------------------------------------------------------------------------- */

int Neighbor::init_pair()
{
  int i,j,k,m;

  // test if pairwise lists need to be re-created
  // no need to re-create if:
  //   neigh style, triclinic, pgsize, oneatom have not changed
  //   current requests = old requests
  // so just return:
  //   delete requests so next run can make new ones
  //   current set of NeighLists already stores all needed info
  // requests are compared via identical() before:
  //   any requests are morphed using logic below
  //   any requests are added below, e.g. as parents of pair hybrid skip lists
  // copy them via requests_new2old() BEFORE any changes made to requests
  //   necessary b/c morphs can change requestor settings (see comment below)

  int same = 1;
  if (style != old_style) same = 0;
  if (triclinic != old_triclinic) same = 0;
  if (pgsize != old_pgsize) same = 0;
  if (oneatom != old_oneatom) same = 0;

  if (nrequest != old_nrequest) same = 0;
  else
    for (i = 0; i < nrequest; i++)
      if (requests[i]->identical(old_requests[i]) == 0) same = 0;

#ifdef NEIGH_LIST_DEBUG
  if (comm->me == 0) printf("SAME flag %d\n",same);
#endif

  if (same) return same;
  requests_new2old();

  // delete old lists since creating new ones

  for (i = 0; i < nlist; i++) delete lists[i];
  for (i = 0; i < nbin; i++) delete neigh_bin[i];
  for (i = 0; i < nstencil; i++) {
    delete neigh_stencil[i];
  }
  for (i = 0; i < nlist; i++) delete neigh_pair[i];
  delete [] lists;
  delete [] neigh_bin;
  delete [] neigh_stencil;
  delete [] neigh_pair;

  // error check on requests
  // do not allow occasional, ghost, bin list
  //   b/c it still uses variant of coord2bin() in NPair() method
  //     instead of atom2bin, this could cause error b/c stoms have
  //     moved out of proc domain by time occasional list is built
  //   solution would be to use a different NBin variant
  //     that used Npair::coord2bin(x,ix,iy,iz) (then delete it from NPair)
  //     and stored the ix,iy,iz values for all atoms (including ghosts)
  //     at time of binning when neighbor lists are rebuilt,
  //     similar to what vanilla Nbin::coord2atom() does now in atom2bin

  if (style == BIN) {
    for (i = 0; i < nrequest; i++)
      if (requests[i]->occasional && requests[i]->ghost)
        error->all(FLERR,"Cannot request an occasional binned neighbor list "
            "with ghost info");
  }

  // morph requests in various ways
  // purpose is to avoid duplicate or inefficient builds
  // may add new requests if a needed request to derive from does not exist
  // methods:
  //   (1) other = point history and rRESPA lists at their partner lists
  //   (2) skip = create any new non-skip lists needed by pair hybrid skip lists
  //   (3) granular = adjust parent and skip lists for granular onesided usage
  //   (4) h/f = pair up any matching half/full lists
  //   (5) copy = convert as many lists as possible to copy lists
  // order of morph methods matters:
  //   (1) before (2), b/c (2) needs to know history partner pairings
  //   (2) after (1), b/c (2) may also need to create new history lists
  //   (3) after (2), b/c it adjusts lists created by (2)
  //   (4) after (2) and (3), 
  //       b/c (2) may create new full lists, (3) may change them
  //   (5) last, after all lists are finalized, so all possible copies found

  int nrequest_original = nrequest;

  morph_other();
  //  morph_skip();
  //  morph_granular();     // this method can change flags set by requestor
  //  morph_halffull();
  morph_copy();

  // create new lists, one per request including added requests
  // wait to allocate initial pages until copy lists are detected
  // NOTE: can I allocate now, instead of down below?

  nlist = nrequest;

  lists = new NeighList*[nrequest];
  neigh_bin = new NBin*[nrequest];
  neigh_stencil = new NStencil*[nrequest];
  neigh_pair = new NPair*[nrequest];

  // allocate new lists
  // pass list ptr back to requestor (except for Command class)
  // only for original requests, not ones added by Neighbor class

  for (i = 0; i < nrequest; i++) {
    //    if (requests[i]->kokkos_host || requests[i]->kokkos_device) 
    //      create_kokkos_list(i);
    //    else lists[i] = new NeighList(cac);
    lists[i] = new NeighList(cac);
    lists[i]->index = i;

    if (requests[i]->pair && i < nrequest_original) {
      Pair *pair = (Pair *) requests[i]->requestor;
      pair->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->fix && i < nrequest_original) {
      Fix *fix = (Fix *) requests[i]->requestor;
      fix->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->compute && i < nrequest_original) {
      Compute *compute = (Compute *) requests[i]->requestor;
      compute->init_list(requests[i]->id,lists[i]);
    }
  }

  // invoke post_constructor() for all lists
  // copies info from requests to lists, sets ptrs to related lists

  for (i = 0; i < nrequest; i++)
    lists[i]->post_constructor(requests[i]);

  // assign Bin,Stencil,Pair style to each list

  int flag;
  for (i = 0; i < nrequest; i++) {
    flag = choose_bin(requests[i]);
    lists[i]->bin_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor bin option does not exist");

    flag = choose_stencil(requests[i]);
    lists[i]->stencil_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor stencil method does not exist");

    flag = choose_pair(requests[i]);
    lists[i]->pair_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor pair method does not exist");
  }

  // instantiate unique Bin,Stencil classes in neigh_bin & neigh_stencil vecs
  // unique = only one of its style, or request unique flag set (custom cutoff)

  nbin = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_bin = -1;
    flag = lists[i]->bin_method;
    if (flag == 0) continue;
    for (j = 0; j < nbin; j++)
      if (neigh_bin[j]->istyle == flag) break;
    if (j < nbin && !requests[i]->unique) {
      requests[i]->index_bin = j;
      continue;
    }

    BinCreator bin_creator = binclass[flag-1];
    neigh_bin[nbin] = bin_creator(cac);
    neigh_bin[nbin]->post_constructor(requests[i]);
    neigh_bin[nbin]->istyle = flag;
    requests[i]->index_bin = nbin;
    nbin++;
  }

  nstencil = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_stencil = -1;
    flag = lists[i]->stencil_method;
    if (flag == 0) continue;
    for (j = 0; j < nstencil; j++)
      // same reason above
      if (neigh_stencil[j]->istyle == flag) break;
    if (j < nstencil && !requests[i]->unique) {
      requests[i]->index_stencil = j;
      continue;
    }

    StencilCreator stencil_creator = stencilclass[flag-1];
    neigh_stencil[nstencil] = stencil_creator(cac);
    neigh_stencil[nstencil]->post_constructor(requests[i]);
    neigh_stencil[nstencil]->istyle = flag;

    // assign neigh_bin to corresponding neigh_stencil
    if (lists[i]->bin_method > 0) {
      neigh_stencil[nstencil]->nb = neigh_bin[requests[i]->index_bin];
      if (neigh_stencil[nstencil]->nb == NULL)
        error->all(FLERR,"Could not assign bin method to neighbor stencil");
    }

    requests[i]->index_stencil = nstencil;
    nstencil++;
  }

  // instantiate one Pair class per list in neigh_pair vec

  for (i = 0; i < nrequest; i++) {
    requests[i]->index_pair = -1;
    flag = lists[i]->pair_method;
    if (flag == 0) {
      neigh_pair[i] = NULL;
      continue;
    }

    PairCreator pair_creator = pairclass[flag-1];
    neigh_pair[i] = pair_creator(cac);
    neigh_pair[i]->post_constructor(requests[i]);
    neigh_pair[i]->istyle = flag;

    if (lists[i]->bin_method > 0) {
      neigh_pair[i]->nb = neigh_bin[requests[i]->index_bin];
      if (neigh_pair[i]->nb == NULL)
        error->all(FLERR,"Could not assign bin method to neighbor pair");
    }
    if (lists[i]->stencil_method > 0) {
      neigh_pair[i]->ns = neigh_stencil[requests[i]->index_stencil];
      if (neigh_pair[i]->ns == NULL)
        error->all(FLERR,"Could not assign stencil method to neighbor pair");
    }

    requests[i]->index_pair = i;
  }

  // allocate initial pages for each list, except if copy flag set
  // allocate dnum vector of zeroes if set
  
     int dnummax = 0;
     for (i = 0; i < nlist; i++) {
     if (lists[i]->copy) continue;
     lists[i]->setup_pages(pgsize,oneatom);
     dnummax = MAX(dnummax,lists[i]->dnum);
     }

     if (dnummax) {
     delete [] zeroes;
     zeroes = new double[dnummax];
     for (i = 0; i < dnummax; i++) zeroes[i] = 0.0;
     }
     
  // first-time allocation of per-atom data for lists that are built and store
  // lists that are not built: granhistory, respa inner/middle (no neigh_pair)
  // lists that do not store: copy 
  // use nmax for both grow() args
  //   i.e. grow first time to expanded size to avoid future reallocs

  int maxatom = atom->nmax;
  int maxelem = element->nmax;
  int maxintg = element->nmaxintg;
  bigint maxintpl = element->nmaxintpl;

  for (i = 0; i < nlist; i++)
    if (neigh_pair[i] && !lists[i]->copy) {
      if (lists[i]->atomlist) 
        lists[i]->grow_atom(maxatom,maxatom);
      if (lists[i]->elemlist) 
        lists[i]->grow_elem(maxelem,maxelem);
      if (lists[i]->intglist) 
        lists[i]->grow_intg(maxintg,maxintg);
      if (lists[i]->intpllist) 
        lists[i]->grow_intpl(maxintpl,maxintpl);
      if (lists[i]->threebody)
        lists[i]->grow_nia(maxintpl);
    }

  // plist = indices of perpetual NPair classes
  //         perpetual = non-occasional, re-built at every reneighboring
  // slist = indices of perpetual NStencil classes
  //         perpetual = used by any perpetual NPair class

  delete [] slist;
  delete [] plist;
  nstencil_perpetual = npair_perpetual = 0;
  slist = new int[nstencil];
  plist = new int[nlist];

  for (i = 0; i < nlist; i++) {
    if (lists[i]->occasional == 0 && lists[i]->pair_method)
      plist[npair_perpetual++] = i;
  }

  for (i = 0; i < nstencil; i++) {
    flag = 0;
    for (j = 0; j < npair_perpetual; j++)
      if (lists[plist[j]]->stencil_method == neigh_stencil[i]->istyle) 
        flag = 1;
    if (flag) slist[nstencil_perpetual++] = i;
  }

  // reorder plist vector if necessary
  // relevant for lists that are derived from a parent list:
  //   half-full,copy,skip
  // the child index must appear in plist after the parent index
  // swap two indices within plist when dependency is mis-ordered
  // start double loop check again whenever a swap is made
  // done when entire double loop test results in no swaps

  NeighList *ptr;

  int done = 0;
  while (!done) {
    done = 1;
    for (i = 0; i < npair_perpetual; i++) {
      for (k = 0; k < 3; k++) {
        ptr = NULL;
        if (k == 0) ptr = lists[plist[i]]->listcopy;
        //        if (k == 1) ptr = lists[plist[i]]->listskip;
        if (k == 2) ptr = lists[plist[i]]->listfull;
        if (ptr == NULL) continue;
        for (m = 0; m < nrequest; m++)
          if (ptr == lists[m]) break;
        for (j = 0; j < npair_perpetual; j++)
          if (m == plist[j]) break;
        if (j < i) continue;
        int tmp = plist[i];     // swap I,J indices
        plist[i] = plist[j];
        plist[j] = tmp;
        done = 0;
        break;
      }
      if (!done) break;
    }
  }

  // debug output

#ifdef NEIGH_LIST_DEBUG
  for (i = 0; i < nrequest; i++) lists[i]->print_attributes();
#endif

  return same;
}

/* ----------------------------------------------------------------------
   scan NeighRequests to set additional flags
   only for custom cutoff lists
   ------------------------------------------------------------------------- */

void Neighbor::morph_other()
{
  NeighRequest *irq;

  for (int i = 0; i < nrequest; i++) {
    irq = requests[i];

    // if cut flag set by requestor, set unique flag
    // this forces Pair,Stencil,Bin styles to be instantiated separately

    if (irq->cut) irq->unique = 1;
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests for possible copies
   if 2 requests match, turn one into a copy of the other
   ------------------------------------------------------------------------- */

void Neighbor::morph_copy()
{
  int i,j,inewton,jnewton;
  NeighRequest *irq,*jrq;

  for (i = 0; i < nrequest; i++) {
    irq = requests[i];

    // this list is already a copy list due to another morph method

    if (irq->copy) continue;

    // check all other lists

    for (j = 0; j < nrequest; j++) {
      if (i == j) continue;
      jrq = requests[j];

      // other list is already copied from this one

      if (jrq->copy && jrq->copylist == i) continue;

      // other list (jrq) to copy from must be perpetual
      // list that becomes a copy list (irq) can be perpetual or occasional
      // if both lists are perpetual, require j < i
      //   to prevent circular dependence with 3 or more copies of a list

      if (jrq->occasional) continue;
      if (!irq->occasional && j > i) continue;

      // both lists must be half, or both full

      if (irq->half != jrq->half) continue;
      if (irq->full != jrq->full) continue;

      // both lists must be newton on, or both newton off
      // IJ newton = 1 for newton on, 2 for newton off

      inewton = irq->newton;
      if (inewton == 0) inewton = force->newton_pair ? 1 : 2; 
      jnewton = jrq->newton;
      if (jnewton == 0) jnewton = force->newton_pair ? 1 : 2;
      if (inewton != jnewton) continue;

      // ok for non-ghost list to copy from ghost list, but not vice versa

      if (irq->ghost && !jrq->ghost) continue;

      // if irq requests a list, jrq must request the same
      
      if (irq->atomlist && !jrq->atomlist) continue;
      if (irq->elemlist && !jrq->elemlist) continue;
      if (irq->intglist && !jrq->intglist) continue;
      if (irq->intpllist && !jrq->intpllist) continue;

      // these flags must be same,
      //   else 2 lists do not store same pairs
      //   or their data structures are different
      // this includes custom cutoff set by requestor
      // no need to check omp b/c it stores same pairs
      // no need to check dnum b/c only set for history

      if (irq->size != jrq->size) continue;
      if (irq->cut != jrq->cut) continue;
      if (irq->cutoff != jrq->cutoff) continue;
      if (irq->sort != jrq->sort) continue;
      if (irq->sortnum != jrq->sortnum) continue;

      // 2 lists are a match

      break;
    }

    // turn list I into a copy of list J
    // do not copy a list from another copy list, but from its parent list

    if (j < nrequest) {
      irq->copy = 1;
      if (jrq->copy) irq->copylist = jrq->copylist;
      else irq->copylist = j;
    }
  }
}

/* ----------------------------------------------------------------------
   output summary of pairwise neighbor list info
   only called by proc 0
   ------------------------------------------------------------------------- */

void Neighbor::print_pairwise_info()
{
  int i,m;
  char str[128];
  NeighRequest *rq;
  FILE *out;
  double elem_cut = 0.0;
  if (element->nelements) elem_cut = element->max_diag_size;

  double cutghost[3],cut;
  cutghost[0] = cutghost[1] = cutghost[2] = MAX(cutneighmax,comm->cutghostuser);
  if (element->nelements) {
    cutghost[0] += element->max_size_array[0];
    cutghost[1] += element->max_size_array[1];
    cutghost[2] += element->max_size_array[2];
  }

  double atombinsize,elembinsize,bbox[3];
  bbox[0] =  bboxhi[0]-bboxlo[0];
  bbox[1] =  bboxhi[1]-bboxlo[1];
  bbox[2] =  bboxhi[2]-bboxlo[2];

  if (atom->natoms) {
    if (binsizeflag) atombinsize = binsize_user;
    else if (style == BIN) atombinsize = 0.5*cutneighmax;
    else atombinsize = 0.5*cutneighmin;
    if (atombinsize == 0.0) atombinsize = bbox[0];
  }

  if (element->nelements) {
    if (binsizeflag) elembinsize = binsize_user;
    else if (style == BIN) elembinsize = 0.5*(cutneighmax + elem_cut);
    else elembinsize = 0.5*(cutneighmin + elem_cut);
    if (elembinsize == 0.0) elembinsize = bbox[0];
  }

  int nperpetual = 0;
  int noccasional = 0;
  int nextra = 0;
  for (i = 0; i < nlist; i++) {
    if (lists[i]->pair_method == 0) nextra++;
    else if (lists[i]->occasional) noccasional++;
    else nperpetual++;
  }

  int esplit = element->esplit;
  for (m = 0; m < 2; m++) {
    if (m == 0) out = screen;
    else out = logfile;

    if (out) {
      fprintf(out,"Neighbor list info ...\n");
      fprintf(out,"  update every %d steps, delay %d steps, check %s\n",
          every,delay,dist_check ? "yes" : "no");

      if (element->nelements) fprintf(out,"  element split: %d x %d x %d\n",esplit,esplit,esplit);
      fprintf(out,"  max neighbors/atom: %d, page size: %d\n",
          oneatom, pgsize);
      fprintf(out,"  master list distance cutoff = %g\n",cutneighmax);
      if (element->nelements) fprintf(out,"  ghost atom/element cutoff = %g %g %g\n",cutghost[0],cutghost[1],cutghost[2]);
      else fprintf(out,"  ghost atom/element cutoff = %g\n",cut);
      if (style != NSQ) {
        if (atom->natoms) fprintf(out,"  atom bin: binsize = %g, bins = %g %g %g\n",atombinsize,
            ceil(bbox[0]/atombinsize), ceil(bbox[1]/atombinsize),
            ceil(bbox[2]/atombinsize));
        if (element->nelements) fprintf(out,"  elem bin: binsize = %g, bins = %g %g %g\n",elembinsize,
            ceil(bbox[0]/elembinsize), ceil(bbox[1]/elembinsize),
            ceil(bbox[2]/elembinsize));
      }

      fprintf(out,"  %d neighbor lists, "
          "perpetual/occasional/extra = %d %d %d\n",
          nlist,nperpetual,noccasional,nextra);

      for (i = 0; i < nlist; i++) {
        rq = requests[i];
        if (rq->pair) {
          char *pname = force->pair_match_ptr((Pair *) rq->requestor);
          sprintf(str,"  (%d) pair %s",i+1,pname);
        } else if (rq->fix) {
          sprintf(str,"  (%d) fix %s",i+1,((Fix *) rq->requestor)->style);
        } else if (rq->compute) {
          sprintf(str,"  (%d) compute %s",i+1,
              ((Compute *) rq->requestor)->style);
        } else if (rq->command) {
          sprintf(str,"  (%d) command %s",i+1,rq->command_style);
        } else if (rq->neigh) {
          sprintf(str,"  (%d) neighbor class addition",i+1);
        }
        fprintf(out,"%s",str);

        if (rq->occasional) fprintf(out,", occasional");
        else fprintf(out,", perpetual");

        // order these to get single output of most relevant

        if (rq->copy)
          fprintf(out,", copy from (%d)",rq->copylist+1);
        //        else if (rq->halffull) 
        //          fprintf(out,", half/full from (%d)",rq->halffulllist+1);
        //        else if (rq->skip)
        //          fprintf(out,", skip from (%d)",rq->skiplist+1);

        fprintf(out,"\n");

        // list of neigh list attributes

        fprintf(out,"      attributes: ");
        if (rq->half) fprintf(out,"half");
        else if (rq->full) fprintf(out,"full");

        if (rq->newton == 0) {
          if (force->newton_pair) fprintf(out,", newton on");
          else fprintf(out,", newton off");
        } else if (rq->newton == 1) fprintf(out,", newton on");
        else if (rq->newton == 2) fprintf(out,", newton off");

        if (rq->intglist) fprintf(out,", integration point");
        if (rq->intpllist) fprintf(out,", interpolated atom");
        if (rq->ghost) fprintf(out,", ghost");
        if (rq->size) fprintf(out,", size");
        if (rq->cut) fprintf(out,", cut %g",rq->cutoff);
        if (rq->sort) {
          fprintf(out,", number of sorted neighbors");
          if (rq->sortnum) fprintf(out,": %d",rq->sortnum);
          else fprintf(out,": all");
        }
        if (rq->off2on) fprintf(out,", off2on");
        fprintf(out,"\n");

        fprintf(out,"      ");
        if (lists[i]->pair_method == 0) fprintf(out,"pair build: none\n");
        else fprintf(out,"pair build: %s\n",pairnames[lists[i]->pair_method-1]);

        fprintf(out,"      ");
        if (lists[i]->stencil_method == 0) fprintf(out,"stencil: none\n");
        else fprintf(out,"stencil: %s\n",
            stencilnames[lists[i]->stencil_method-1]);

        fprintf(out,"      ");
        if (lists[i]->bin_method == 0) fprintf(out,"bin: none\n");
        else fprintf(out,"bin: %s\n",binnames[lists[i]->bin_method-1]);
      }


      fprintf(out,"  %d stencil methods\n",nstencil);
      for (i = 0; i < nstencil; i++)
        fprintf(out,"    (%d) %s\n",
            i+1,stencilnames[neigh_stencil[i]->istyle-1]);

      fprintf(out,"  %d bin methods\n",nbin);
      for (i = 0; i < nbin; i++)
        fprintf(out,"    (%d) %s\n",i+1,binnames[neigh_bin[i]->istyle-1]);

    }
  }
}

/* ----------------------------------------------------------------------
   make copy of current requests and Neighbor params
   used to compare to when next run occurs
   ------------------------------------------------------------------------- */

void Neighbor::requests_new2old()
{
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);

  old_nrequest = nrequest;
  old_requests = (NeighRequest **)
    memory->smalloc(old_nrequest*sizeof(NeighRequest *),
        "neighbor:old_requests");

  for (int i = 0; i < old_nrequest; i++) {
    old_requests[i] = new NeighRequest(cac);
    old_requests[i]->copy_request(requests[i]);
  }

  old_style = style;
  old_triclinic = triclinic;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
}

/* ----------------------------------------------------------------------
   find and return request made by classptr
   if not found or classpt = NULL, return NULL
   ------------------------------------------------------------------------- */

NeighRequest *Neighbor::find_request(void *classptr)
{
  if (classptr == NULL) return NULL;

  for (int i = 0; i < nrequest; i++)
    if (requests[i]->requestor == classptr) return requests[i];

  return NULL;
}

/* ----------------------------------------------------------------------
   assign NBin class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known Nbin classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
   ------------------------------------------------------------------------- */

int Neighbor::choose_bin(NeighRequest *rq)
{
  // no binning needed

  if (style == NSQ) return 0;

  // use request settings to match exactly one NBin class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < nbclass; i++) {
    mask = binmasks[i];

    // only one bin style
    // add requirement here

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   assign NStencil class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known NStencil classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
   ------------------------------------------------------------------------- */

int Neighbor::choose_stencil(NeighRequest *rq)
{
  // no stencil creation needed

  if (style == NSQ) return 0;
  if (rq->copy) return 0;
  //  if (rq->skip || rq->copy || rq->halffull) return 0;
  //  if (rq->history) return 0;
  //  if (rq->respainner || rq->respamiddle) return 0;

  // convert newton request to newtflag = on or off

  int newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = 1;
  else if (rq->newton == 0 && !newton_pair) newtflag = 0;
  else if (rq->newton == 1) newtflag = 1;
  else if (rq->newton == 2) newtflag = 0;

  //printf("STENCIL RQ FLAGS: hff %d %d n %d g %d s %d newtflag %d\n",
  //       rq->half,rq->full,rq->newton,rq->ghost,rq->ssa,
  //       newtflag);

  // use request and system settings to match exactly one NStencil class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < nsclass; i++) {
    mask = stencilmasks[i];

    //printf("III %d: half %d full %d newton %d newtoff %d ghost %d ssa %d\n",
    //       i,mask & NS_HALF,mask & NS_FULL,mask & NS_NEWTON,
    //       mask & NS_NEWTOFF,mask & NS_GHOST,mask & NS_SSA);

    // exactly one of half or full is set and must match

    if (rq->half) {
    //  if (!(mask & NS_HALF)) continue;
    //} else if (rq->full) {
      if (!(mask & NS_FULL)) continue;
    }

    // newtflag is on or off and must match

    if (newtflag) {
      if (!(mask & NS_NEWTON)) continue;
    } else if (!newtflag) {
      if (!(mask & NS_NEWTOFF)) continue;
    }

    // require match of these request flags and mask bits 
    // (!A != !B) is effectively a logical xor

    if (!rq->ghost != !(mask & NS_GHOST)) continue;
    //    if (!rq->ssa != !(mask & NS_SSA)) continue;

    // neighbor style is BIN or MULTI and must match

    if (style == BIN) {
      if (!(mask & NS_BIN)) continue;
    } else if (style == MULTI) {
      if (!(mask & NS_MULTI)) continue;
    }

    // dimension is 2 or 3 and must match

    if (dimension == 2) {
      if (!(mask & NS_2D)) continue;
    } else if (dimension == 3) {
      if (!(mask & NS_3D)) continue;
    }

    // domain triclinic flag is on or off and must match

    if (triclinic) {
      if (!(mask & NS_TRI)) continue;
    } else if (!triclinic) {
      if (!(mask & NS_ORTHO)) continue;
    }

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   assign NPair class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known NPair classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
   ------------------------------------------------------------------------- */

int Neighbor::choose_pair(NeighRequest *rq)
{
  // convert newton request to newtflag = on or off

  int newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = 1;
  else if (rq->newton == 0 && !newton_pair) newtflag = 0;
  else if (rq->newton == 1) newtflag = 1;
  else if (rq->newton == 2) newtflag = 0;

  // use request and system settings to match exactly one NPair class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < npclass; i++) {
    mask = pairmasks[i];

    // if copy request, no further checks needed, just return or continue

    if (rq->copy) {
      if (!(mask & NP_COPY)) continue;
      return i+1;
    }

    // exactly one of half or full is set and must match

    if (rq->half) {
      if (!(mask & NP_HALF)) continue;
    } else if (rq->full) {
      if (!(mask & NP_FULL)) continue;
    }

    // neighborlist for atoms, elements, integration points and/or interpolated atoms
    // must match

    if (rq->atomlist)
      if (!(mask & NP_ATOM)) continue;
    if (rq->elemlist)
      if (!(mask & NP_ELEM)) continue;
    if (rq->intglist) 
      if (!(mask & NP_INTG)) continue;
    if (rq->intpllist)
      if (!(mask & NP_INTPL)) continue;

    // newtflag is on or off and must match

    if (newtflag) {
      if (!(mask & NP_NEWTON)) continue;
    } else if (!newtflag) {
      if (!(mask & NP_NEWTOFF)) continue;
    }

    // threebody flag, must match
   
    if (rq->threebody)
      if (!(mask & NP_THREEBODY)) continue;

    // require match of these request flags and mask bits 
    // (!A != !B) is effectively a logical xor

    if (!rq->ghost != !(mask & NP_GHOST)) continue;

    // neighbor style is one of NSQ,BIN,MULTI and must match

    if (style == NSQ) {
      if (!(mask & NP_NSQ)) continue;
    } else if (style == BIN) {
      if (!(mask & NP_BIN)) continue;
    } else if (style == MULTI) {
      if (!(mask & NP_MULTI)) continue;
    }

    // domain triclinic flag is on or off and must match

    if (triclinic) {
      if (!(mask & NP_TRI)) continue;
    } else if (!triclinic) {
      if (!(mask & NP_ORTHO)) continue;
    }

    // neighbor sorting must match

    if (rq->sort)
      if (!(mask & NP_SORT)) continue;

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   called by other classes to request a pairwise neighbor list
   ------------------------------------------------------------------------- */

int Neighbor::request(void *requestor, int instance)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)
      memory->srealloc(requests,maxrequest*sizeof(NeighRequest *),
          "neighbor:requests");
  }

  requests[nrequest] = new NeighRequest(cac);
  requests[nrequest]->index = nrequest;
  requests[nrequest]->requestor = requestor;
  requests[nrequest]->requestor_instance = instance;
  nrequest++;
  return nrequest-1;
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_bin.h
   ------------------------------------------------------------------------- */

  template <typename T>
NBin *Neighbor::bin_creator(CAC *cac)
{
  return new T(cac);
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_stencil.h
   ------------------------------------------------------------------------- */

  template <typename T>
NStencil *Neighbor::stencil_creator(CAC *cac)
{
  return new T(cac);
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_pair.h
   ------------------------------------------------------------------------- */

  template <typename T>
NPair *Neighbor::pair_creator(CAC *cac)
{
  return new T(cac);
}

/* ----------------------------------------------------------------------
   setup neighbor binning and neighbor stencils
   called before run and every reneighbor if box size/shape changes
   only operates on perpetual lists
   build_one() operates on occasional lists
   ------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  // invoke setup_bins() for all NBin
  // actual binning is performed in build()

  for (int i = 0; i < nbin; i++) neigh_bin[i]->setup_bins(style);

  // invoke create_setup() and create() for all perpetual NStencil
  // same ops performed for occasional lists in build_one()

  for (int i = 0; i < nstencil_perpetual; i++) {
    neigh_stencil[slist[i]]->create_setup();
    neigh_stencil[slist[i]]->create();
  }

  last_setup_bins = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int Neighbor::decide()
{
  if (must_check) {
    bigint n = update->ntimestep;
    //    if (restart_check && n == output->next_restart) return 1;
    for (int i = 0; i < fix_check; i++)
      if (n == modify->fix[fixchecklist[i]]->next_reneighbor) return 1;
  }

  ago++;
  if (ago >= delay && ago % every == 0) {
    if (build_once) return 0;
    if (dist_check == 0) return 1;
    return check_distance();
  } else return 0;
}

/* ----------------------------------------------------------------------
   if any atom/node moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
   compute distance each of 8 corners of box has moved since last reneighbor
   reduce skin distance by sum of 2 largest of the 8 values
   new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
   ------------------------------------------------------------------------- */

int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;
  double delta,deltasq,delta1,delta2;

  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx*delx + dely*dely + delz*delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx*delx + dely*dely + delz*delz);
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx*delx + dely*dely + delz*delz);
        if (delta > delta1) delta1 = delta;
        else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    }
  } else deltasq = triggersq;

  double **x = atom->x;
  int nlocal = atom->nlocal;


  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - atomxhold[i][0];
    dely = x[i][1] - atomxhold[i][1];
    delz = x[i][2] - atomxhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > deltasq) flag = 1;
  }

  double ***nodex = element->nodex;
  nlocal = element->nlocal;
  int npe = element->npe;
  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < npe; k++) {
      delx = nodex[i][k][0] - nodexhold[i][k][0];
      dely = nodex[i][k][1] - nodexhold[i][k][1];
      delz = nodex[i][k][2] - nodexhold[i][k][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq > deltasq) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

/* ----------------------------------------------------------------------
   build perpetual neighbor lists
   called at setup and every few timesteps during run or minimization
   topology lists also built if topoflag = 1 (Kokkos calls with topoflag=0)
   ------------------------------------------------------------------------- */

void Neighbor::build(int topoflag)
{
  int i,j,m;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;

  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  int npe = element->npe;
  int max_nintg = element->max_nintg;
  int max_nintpl = element->max_nintpl;
  bigint nlocalintpl = nelocal*max_nintpl;
  bigint nallintpl = neall*max_nintpl;

  // check that using special bond flags will not overflow neigh lists

  if (naall > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // store current atom/node positions and box size if needed

  if (dist_check) {
    double **x = atom->x;
    double ***nodex = element->nodex;
    //    if (includegroup) nlocal = atom->nfirst;
    if (atom->nmax > atommaxhold) {
      atommaxhold = atom->nmax;
      memory->destroy(atomxhold);
      memory->create(atomxhold,atommaxhold,3,"neigh:atomxhold");
    }
    if (element->nmax > nodemaxhold) {
      nodemaxhold = element->nmax;
      memory->destroy(nodexhold);
      memory->create(nodexhold,nodemaxhold,npe,3,"neigh:nodexhold");
    }     
    for (i = 0; i < nalocal; i++) {
      atomxhold[i][0] = x[i][0];
      atomxhold[i][1] = x[i][1];
      atomxhold[i][2] = x[i][2];
    }
    for (i = 0; i < nelocal; i++) {
      for (j = 0; j < npe; j++) {
        nodexhold[i][j][0] = nodex[i][j][0];
        nodexhold[i][j][1] = nodex[i][j][1];
        nodexhold[i][j][2] = nodex[i][j][2];
      }
    }
    if (boxcheck) {
      if (triclinic == 0) {
        boxlo_hold[0] = bboxlo[0];
        boxlo_hold[1] = bboxlo[1];
        boxlo_hold[2] = bboxlo[2];
        boxhi_hold[0] = bboxhi[0];
        boxhi_hold[1] = bboxhi[1];
        boxhi_hold[2] = bboxhi[2];
      } else {
        domain->box_corners();
        corners = domain->corners;
        for (i = 0; i < 8; i++) {
          corners_hold[i][0] = corners[i][0];
          corners_hold[i][1] = corners[i][1];
          corners_hold[i][2] = corners[i][2];
        }
      }
    }
  }

  // bin atoms for all NBin instances
  // not just NBin associated with perpetual lists, also occasional lists
  // b/c cannot wait to bin occasional lists in build_one() call
  // if bin then, atoms may have moved outside of proc domain & bin extent,
  //   leading to errors or even a crash

  if (style != NSQ) {
    for (i = 0; i < nbin; i++) {
      neigh_bin[i]->bin_setup(naall,neall);
      neigh_bin[i]->bin_all();
    }
  }
  // build pairwise lists for all perpetual NPair/NeighList
  // grow() with nlocal/nall args so that only realloc if have to

  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->copy) {
      if (lists[m]->atomlist) 
        lists[m]->grow_atom(nalocal,naall);
      if (lists[m]->elemlist) 
        lists[m]->grow_elem(nelocal,neall);
      if (lists[m]->intglist) 
        lists[m]->grow_intg(nelocal*max_nintg,neall*max_nintg);
      if (lists[m]->intpllist) 
        lists[m]->grow_intpl(nlocalintpl,nallintpl);
      if (lists[m]->threebody) {
        lists[m]->grow_nia(nallintpl);
        lists[m]->grow_elem(neall,neall);
      }
    }
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }

  // build topology lists for bonds/angles/etc

}

/* ----------------------------------------------------------------------
   build a single occasional pairwise neighbor list indexed by I
   called by other classes
   ------------------------------------------------------------------------- */

void Neighbor::build_one(class NeighList *mylist, int preflag)
{
  // check if list structure is initialized

  if (mylist == NULL)
    error->all(FLERR,"Trying to build an occasional neighbor list "
        "before initialization completed");

  // build_one() should never be invoked on a perpetual list

  if (!mylist->occasional) 
    error->all(FLERR,"Neighbor build one invoked on perpetual list");

  // no need to build if already built since last re-neighbor
  // preflag is set by fix bond/create and fix bond/swap
  //   b/c they invoke build_one() on same step neigh list is re-built,
  //   but before re-build, so need to use ">" instead of ">="

  NPair *np = neigh_pair[mylist->index];

  if (preflag) {
    if (np->last_build > lastcall) return;
  } else {
    if (np->last_build >= lastcall) return;
  }

  // if this is copy list and parent is occasional list,
  // or this is halffull and parent is occasional list,
  // insure parent is current 

  if (mylist->listcopy && mylist->listcopy->occasional)
    build_one(mylist->listcopy,preflag);
  if (mylist->listfull && mylist->listfull->occasional)
    build_one(mylist->listfull,preflag);



  // create stencil if hasn't been created since last setup_bins() call

  NStencil *ns = np->ns;
  if (ns && ns->last_stencil < last_setup_bins) {
    ns->create_setup();
    ns->create();
  }

  // build the list

  if (!mylist->copy) {
    int nalocal = atom->nlocal;
    int naall = nalocal + atom->nghost;
    int nelocal = element->nlocal;
    int neall = nelocal + element->nghost;
    int max_nintg = element->max_nintg;
    int max_nintpl = element->max_nintpl;
    bigint nlocalintpl = nelocal*max_nintpl;
    bigint nallintpl = neall*max_nintpl;
    if (mylist->atomlist) 
      mylist->grow_atom(nalocal,naall);
    if (mylist->elemlist) 
      mylist->grow_elem(nelocal,neall);
    if (mylist->intglist) 
      mylist->grow_intg(nelocal*max_nintg,neall*max_nintg);
    if (mylist->intpllist) 
      mylist->grow_intpl(nlocalintpl,nallintpl);
    if (mylist->threebody) {
      mylist->grow_nia(nlocalintpl);
      mylist->grow_elem(neall,neall);
    }
  }

  np->build_setup();
  np->build(mylist);
}

/* ----------------------------------------------------------------------
   set neighbor style and skin distance
   ------------------------------------------------------------------------- */

void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command");

  skin = universe->numeric(FLERR,arg[0]);
  if (skin < 0.0) error->all(FLERR,"Illegal neighbor command");

  if (strcmp(arg[1],"nsq") == 0) {style = NSQ; error->all(FLERR,"No NSQ neighbor style yet");}
  else if (strcmp(arg[1],"bin") == 0) style = BIN;
  else if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else error->all(FLERR,"Illegal neighbor command");

}

/* ----------------------------------------------------------------------
   reset timestamps in all NeignBin, NStencil, NPair classes
   so that neighbor lists will rebuild properly with timestep change
   ditto for lastcall and last_setup_bins
   ------------------------------------------------------------------------- */

void Neighbor::reset_timestep(bigint ntimestep)
{
  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->last_bin = -1;

  for (int i = 0; i < nstencil; i++)
    neigh_stencil[i]->last_stencil = -1;

  for (int i = 0; i < nlist; i++) {
    if (!neigh_pair[i]) continue;
    neigh_pair[i]->last_build = -1;
  }

  lastcall = -1;
  last_setup_bins = -1;
}

/* ----------------------------------------------------------------------
   modify parameters of the pair-wise neighbor build
   ------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      every = universe->inumeric(FLERR,arg[iarg+1]);
      if (every <= 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      delay = universe->inumeric(FLERR,arg[iarg+1]);
      if (delay < 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) dist_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) dist_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"once") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) build_once = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) build_once = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      old_pgsize = pgsize;
      pgsize = universe->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      old_oneatom = oneatom;
      oneatom = universe->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      binsize_user = universe->numeric(FLERR,arg[iarg+1]);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cluster") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) cluster_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) cluster_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal neigh_modify command");
  } 
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint Neighbor::memory_usage()
{
  int npe = element->npe;
  bigint bytes = 0;
  bytes += memory->usage(atomxhold,atommaxhold,3);
  bytes += memory->usage(nodexhold,nodemaxhold,npe,3);

  for (int i = 0; i < nlist; i++)
    if (lists[i]) bytes += lists[i]->memory_usage();
  for (int i = 0; i < nstencil; i++) 
    bytes += neigh_stencil[i]->memory_usage();
  for (int i = 0; i < nbin; i++) 
    bytes += neigh_bin[i]->memory_usage();
  return bytes;
}
