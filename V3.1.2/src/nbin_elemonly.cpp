#include "nbin_elemonly.h"
#include "neighbor.h"
#include "atom.h"
#include "element.h"
#include "group.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "error.h"

using namespace CAC_NS;


#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100
#define BIG           1e30

/*  ----------------------------------------------------------------------  */

NBinElemonly::NBinElemonly(CAC *cac) : NBin(cac) {}

/*  ----------------------------------------------------------------------
   setup neighbor binning geometry
   bin numbering in each dimension is global:
   0 = 0.0 to binsize, 1 = binsize to 2 * binsize, etc
   nbin-1, nbin, etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
   -1, -2, etc = -binsize to 0.0, -2 * binsize to -binsize, etc
   code will work for any binsize
   since next(xyz) and stencil extend as far as necessary
   binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
   a dim must be filled exactly by integer # of bins
   in periodic, procs on both sides of PBC must see same bin boundary
   in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
   tilted simulation box cannot contain integer # of bins
   stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
   setup 2 set of bins for atom and element separately
   -------------------------------------------------------------------------  */

void NBinElemonly::setup_bins(int style)
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->max_cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->max_cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3], bsubboxlo[3], bsubboxhi[3];
  double *max_cutghost = comm->max_cutghost;
  double max_bound_box[3];
  MPI_Allreduce(element->max_element_bound_box_size, max_bound_box, 3, MPI_DOUBLE, MPI_MAX, world);

  if (triclinic == 0) {
    bsubboxlo[0] = domain->sublo[0] - max_cutghost[0] - max_bound_box[0] * 2;
    bsubboxlo[1] = domain->sublo[1] - max_cutghost[1] - max_bound_box[1] * 2;
    bsubboxlo[2] = domain->sublo[2] - max_cutghost[2] - max_bound_box[2] * 2;
    bsubboxhi[0] = domain->subhi[0] + max_cutghost[0] + max_bound_box[0] * 2;
    bsubboxhi[1] = domain->subhi[1] + max_cutghost[1] + max_bound_box[1] * 2;
    bsubboxhi[2] = domain->subhi[2] + max_cutghost[2] + max_bound_box[2] * 2;
  } else {
    double lo[3], hi[3];
    lo[0] = domain->sublo_lamda[0] - max_cutghost[0] - max_bound_box[0] * 2;
    lo[1] = domain->sublo_lamda[1] - max_cutghost[1] - max_bound_box[1] * 2;
    lo[2] = domain->sublo_lamda[2] - max_cutghost[2] - max_bound_box[2] * 2;
    hi[0] = domain->subhi_lamda[0] + max_cutghost[0] + max_bound_box[0] * 2;
    hi[1] = domain->subhi_lamda[1] + max_cutghost[1] + max_bound_box[1] * 2;
    hi[2] = domain->subhi_lamda[2] + max_cutghost[2] + max_bound_box[2] * 2;
    domain->bbox(lo, hi, bsubboxlo, bsubboxhi);
  }

  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  // setup element bins

  if (element->nelements) {
    double binsizex_optimal;
    double binsizey_optimal;
    double binsizez_optimal;

    // if triclinic box, recalculate max_box in box coords since element->max_element_bound_box_size is in lamda coords

    double max_box[3];
    if (triclinic) {
      domain->lamda2nodex(element->nlocal, element->nodex);
      double maxcoord, mincoord, max[3];
      max[2] = 0.0;
      for (int i = 0; i < element->nlocal; i++) {
        double ***inodex = element->nodex[i];
        int inpe = element->npe[element->etype[i]];
        int iapc = element->apc[element->etype[i]];
        for (int dim = 0; dim < dimension; dim++) {
          maxcoord = -BIG;
          mincoord = BIG;
          for (int j = 0; j < iapc; j++) 
            for (int k = 0; k < inpe; k++) {
              maxcoord = MAX(maxcoord, inodex[j][k][dim]);
              mincoord = MIN(mincoord, inodex[j][k][dim]);
            }
          max[dim] = MAX(maxcoord-mincoord, max[dim]);
        }
      }
      MPI_Allreduce(max, max_box, 3, MPI_DOUBLE, MPI_MAX, world);
      domain->nodex2lamda(element->nlocal, element->nodex);
    } else {
      max_box[0] = max_bound_box[0];
      max_box[1] = max_bound_box[1];
      max_box[2] = max_bound_box[2];
    }

    if (elembinsizeflag) 
      binsizex_optimal = binsizey_optimal = binsizez_optimal = elembinsize_user;
    else if (style == Neighbor::BIN) {
      binsizex_optimal = 0.5 * (cutneighmax + max_box[0]);
      binsizey_optimal = 0.5 * (cutneighmax + max_box[1]);
      binsizez_optimal = 0.5 * (cutneighmax + max_box[2]);
    } else {
      binsizex_optimal = 0.5 * (cutneighmin + max_box[0]);
      binsizey_optimal = 0.5 * (cutneighmin + max_box[1]);
      binsizez_optimal = 0.5 * (cutneighmin + max_box[2]);
    }
    if (binsizex_optimal == 0.0) binsizex_optimal = bbox[0];
    if (binsizey_optimal == 0.0) binsizey_optimal = bbox[1];
    if (binsizez_optimal == 0.0) binsizez_optimal = bbox[2];
    double binsizexinv = 1.0/binsizex_optimal;
    double binsizeyinv = 1.0/binsizey_optimal;
    double binsizezinv = 1.0/binsizez_optimal;

    // test for too many global bins in any dimension due to huge global domain

    if (bbox[0] * binsizexinv > MAXSMALLINT || 
        bbox[1] * binsizeyinv > MAXSMALLINT ||
        bbox[2] * binsizeyinv > MAXSMALLINT)
      error->all(FLERR, "Domain too large for neighbor element bins");

    // create actual bins
    // always have one bin even if cutoff > bbox
    // for 2d, nbinz = 1

    nelembinx = static_cast<int> (bbox[0] * binsizexinv);
    nelembiny = static_cast<int> (bbox[1] * binsizeyinv);
    if (dimension == 3) nelembinz = static_cast<int> (bbox[2] * binsizezinv);
    else nelembinz = 1;

    if (nelembinx == 0) nelembinx = 1;
    if (nelembiny == 0) nelembiny = 1;
    if (nelembinz == 0) nelembinz = 1;

    // compute actual bin size for nbins to fit into box exactly
    // error if actual bin size << cutoff, since will create a zillion bins
    // this happens when nbin = 1 and box size << cutoff
    // typically due to non-periodic, flat system in a particular dim
    // in that extreme case, should use NSQ not BIN neighbor style

    elembinsizex = bbox[0]/nelembinx;
    elembinsizey = bbox[1]/nelembiny;
    elembinsizez = bbox[2]/nelembinz;

    elembininvx = 1.0 / elembinsizex;
    elembininvy = 1.0 / elembinsizey;
    elembininvz = 1.0 / elembinsizez;

    if (binsizex_optimal * elembininvx > CUT2BIN_RATIO ||
        binsizey_optimal * elembininvy > CUT2BIN_RATIO ||
        binsizez_optimal * elembininvz > CUT2BIN_RATIO)
      error->all(FLERR, "Cannot use neighbor element bins - box size << cutoff");

    // mbinlo/hi = lowest and highest global bins my ghost elements could be in
    // coord = lowest and highest values of coords for my ghost elements
    // static_cast(-1.5) = -1, so subract additional -1
    // add in SMALL for round-off safety

    int melembinxhi, melembinyhi, melembinzhi;
    double coord;

    coord = bsubboxlo[0] - SMALL * bbox[0];
    melembinxlo = static_cast<int> ((coord-bboxlo[0]) * elembininvx);
    if (coord < bboxlo[0]) melembinxlo--;
    coord = bsubboxhi[0] + SMALL * bbox[0];
    melembinxhi = static_cast<int> ((coord-bboxlo[0]) * elembininvx);

    coord = bsubboxlo[1] - SMALL * bbox[1];
    melembinylo = static_cast<int> ((coord-bboxlo[1]) * elembininvy);
    if (coord < bboxlo[1]) melembinylo--;
    coord = bsubboxhi[1] + SMALL * bbox[1];
    melembinyhi = static_cast<int> ((coord-bboxlo[1]) * elembininvy);

    if (dimension == 3) {
      coord = bsubboxlo[2] - SMALL * bbox[2];
      melembinzlo = static_cast<int> ((coord-bboxlo[2]) * elembininvz);
      if (coord < bboxlo[2]) melembinzlo--;
      coord = bsubboxhi[2] + SMALL * bbox[2];
      melembinzhi = static_cast<int> ((coord-bboxlo[2]) * elembininvz);
    }

    // extend bins by 2 to insure stencil extent is included
    // for 2d, only 1 bin in z

    melembinxlo -= 2;
    melembinxhi += 2;
    melembinx = melembinxhi - melembinxlo + 1;

    melembinylo -= 2;
    melembinyhi += 2;
    melembiny = melembinyhi - melembinylo + 1;

    if (dimension == 3) {
      melembinzlo -= 2;
      melembinzhi += 2;
    } else melembinzlo = melembinzhi = 0;
    melembinz = melembinzhi - melembinzlo + 1;

    bigint bbin = ((bigint) melembinx) * ((bigint) melembiny) * ((bigint) melembinz) + 1;
    if (bbin > MAXSMALLINT) error->one(FLERR, "Too many neighbor element bins");
    melembins = bbin;
  } else {
    melembins = 0;
  }

  //if (comm->me == 166)
  //  printf("Setup_bin step = %d me = %d istyle = %d melembin = %d %d %d\n", update->ntimestep, comm->me, istyle, melembinx, melembiny, melembinz);
  // setup atom bins

  if (atom->natoms) {

    // create actual bins
    // always have one bin even if cutoff > bbox
    // for 2d, nbinz = 1

    natombinx = nelembinx;
    natombiny = nelembiny;
    natombinz = nelembinz;

    atombinsizex = elembinsizex;
    atombinsizey = elembinsizey;
    atombinsizez = elembinsizez;

    atombininvx = elembininvx;
    atombininvy = elembininvy;
    atombininvz = elembininvz;

    matombinx = melembinx;
    matombiny = melembiny;
    matombinz = melembinz;

    matombinxlo = melembinxlo;
    matombinylo = melembinylo;
    matombinzlo = melembinzlo;
    matombins = melembins;
  } else {
    matombins = 0; 
  }

}

/*  ----------------------------------------------------------------------
   bin owned and ghost atoms & elements
   -------------------------------------------------------------------------  */

void NBinElemonly::bin_all()
{

  int i, ibin, nall;
  double **x;
  last_bin = update->ntimestep;
  for (i = 0; i < matombins; i++) 
    atombinhead[i] = -1;
  for (i = 0; i < melembins; i++) 
    elembinhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms/elements at end of list, which is necessary

  x = atom->x;
  nall = atom->nlocal + atom->nghost;
  for (i = nall-1; i >= 0; i--) {
    ibin = coord2atombin(x[i][0], x[i][1], x[i][2]);
    atom2bin[i] = ibin;
    atombins[i] = atombinhead[ibin];
    atombinhead[ibin] = i;
  }

  x = element->x;
  nall = element->nlocal + element->nghost;
  for (i = nall-1; i >= 0; i--) {
    ibin = coord2elembin(x[i][0], x[i][1], x[i][2]);
    elem2bin[i] = ibin;
    elembins[i] = elembinhead[ibin];
    elembinhead[ibin] = i;
  }
}

