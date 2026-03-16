#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_eam.h"
#include "atom.h"
#include "force.h"
#include "element.h"
#include "comm.h"
#include "neighbor.h"
#include "universe.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace CAC_NS;

#define MAXLINE 1024

/* -------------------------------------------------------------------------------------------- */

PairEAM::PairEAM(CAC *cac) : Pair(cac)
{
  restartinfo = 0;
  manybody_flag = 1;
  no_virial_fdotr_compute = 1;

  nemax = 0;
  namax = 0;
  atomrho = NULL;
  atomfp = NULL;
  noderho = NULL;

  map = NULL;
  type2frho = NULL;

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;
  scale = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;

  // set comm size needed by this pair

  comm_atom_forward = 1;
  comm_elem_forward = element->maxnpe * element->maxapc;
  comm_atom_reverse = 1;
  comm_elem_reverse = element->maxnpe * element->maxapc;
}

/* ----------------------------------------------------------------------------------------
  check if allocated, since class can be destructed when incomplete
  ------------------------------------------------------------------------------------------ */

PairEAM::~PairEAM()
{
  if (copymode) return;

  memory->destroy(atomrho);
  memory->destroy(atomfp);
  memory->destroy(noderho);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    map = NULL;
    type2frho = NULL;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
    memory->destroy(scale);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete [] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
    }
    memory->sfree(funcfl);
    funcfl = NULL;
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
    setfl = NULL;
  }

  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
    fs = NULL;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}


/*  ----------------------------------------------------------------------
   global settings
-------------------------------------------------------------------------  */

void PairEAM::settings(int narg, char **arg)
{
  //int iarg = 0;
  //edensity_style = 0;
  //while (iarg < narg) {
  //  if (strcmp(arg[iarg], "edensity") == 0) {
  //    if (iarg+2 > narg) error->all(FLERR, "Illegal pair_style command");
  //    if (strcmp(arg[iarg+1], "enhanced")) 
  //      edensity_style = 1;
  //    else if (strcmp(arg[iarg+1], "normal")) 
  //      edensity_style = 0;
  //    iarg += 2;
  //  } else error->all(FLERR, "Illegal pair_style command");
  //}
  if (narg > 0) error->all(FLERR, "Illegal pair_style command");
}


/*  ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
-------------------------------------------------------------------------  */

void PairEAM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo, ihi, jlo, jhi;
  universe->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  universe->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2], funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *)
      memory->srealloc(funcfl, nfuncfl * sizeof(Funcfl), "pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file, arg[2]);
  }

  // set setflag and map only for i, i type pairs
  // set mass of atom type if i = j

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      if (i == j) {
        setflag[i][i] = 1;
        map[i] = ifuncfl;
        atom->set_mass(FLERR, i, funcfl[ifuncfl].mass);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*  ----------------------------------------------------------------------
   allocate all arrays
-------------------------------------------------------------------------  */

void PairEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor, n+1, n+1, "pair:type2rhor");
  memory->create(type2z2r, n+1, n+1, "pair:type2z2r");
  memory->create(scale, n+1, n+1, "pair:scale");
}

/* ---------------------------------------------------------------------------------
  read potential values from a DYNAMO single element funcfl file
  ----------------------------------------------------------------------------------- */

void PairEAM::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl-1];

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str, "Cannot open EAM potential file %s", filename);
      error->one(FLERR, str);
    }
  }

  int tmp, nwords;
  if (me == 0) {
    fgets(line, MAXLINE, fptr);
    fgets(line, MAXLINE, fptr);
    sscanf(line, "%d %lg", &tmp, &file->mass);
    fgets(line, MAXLINE, fptr);
    nwords = sscanf(line, "%d %lg %d %lg %lg", 
           &file->nrho, &file->drho, &file->nr, &file->dr, &file->cut);
  }

  MPI_Bcast(&nwords, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->mass, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  if ((nwords != 5) || (file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
    error->all(FLERR, "Invalid EAM potential file");

  memory->create(file->frho, (file->nrho+1), "pair:frho");
  memory->create(file->rhor, (file->nr+1), "pair:rhor");
  memory->create(file->zr, (file->nr+1), "pair:zr");

  if (me == 0) grab(fptr, file->nrho, &file->frho[1]);
  MPI_Bcast(&file->frho[1], file->nrho, MPI_DOUBLE, 0, world);

  if (me == 0) grab(fptr, file->nr, &file->zr[1]);
  MPI_Bcast(&file->zr[1], file->nr, MPI_DOUBLE, 0, world);

  if (me == 0) grab(fptr, file->nr, &file->rhor[1]);
  MPI_Bcast(&file->rhor[1], file->nr, MPI_DOUBLE, 0, world);

  if (me == 0) fclose(fptr);
}

/*  ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
-------------------------------------------------------------------------  */

void PairEAM::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line, MAXLINE, fptr);
    ptr = strtok(line, " \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL, " \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/* ----------------------------------------------------------------------------------
  init specific to this pair style
  ------------------------------------------------------------------------------------ */

void PairEAM::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this, instance_me);
}

/*  ----------------------------------------------------------------------
   init for one type pair i, j and corresponding j, i
-------------------------------------------------------------------------  */

double PairEAM::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax, funcfl[m].cut);
  } else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;

  cutforcesq = cutmax * cutmax;

  return cutmax;
}

/*  ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
-------------------------------------------------------------------------  */

void PairEAM::file2array()
{
  int i, j, k, m, n;
  int ntypes = atom->ntypes;
  double sixth = 1.0/6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr, file->dr);
    drho = MAX(drho, file->drho);
    rmax = MAX(rmax, (file->nr-1) * file->dr);
    rhomax = MAX(rhomax, (file->nrho-1) * file->drho);
  }

  // set nr, nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho+1, "pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r, p, cof1, cof2, cof3, cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1) * drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k, file->nrho-2);
      k = MAX(k, 2);
      p -= k;
      p = MIN(p, 2.0);
      cof1 = -sixth * p * (p-1.0) * (p-2.0);
      cof2 = 0.5 * (p * p-1.0) * (p-2.0);
      cof3 = -0.5 * p * (p+1.0) * (p-2.0);
      cof4 = sixth * p * (p * p-1.0);
      frho[n][m] = cof1 * file->frho[k-1] + cof2 * file->frho[k] +
        cof3 * file->frho[k+1] + cof4 * file->frho[k+2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nr+1, "pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1) * dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k, file->nr-2);
      k = MAX(k, 2);
      p -= k;
      p = MIN(p, 2.0);
      cof1 = -sixth * p * (p-1.0) * (p-2.0);
      cof2 = 0.5 * (p * p-1.0) * (p-2.0);
      cof3 = -0.5 * p * (p+1.0) * (p-2.0);
      cof4 = sixth * p * (p * p-1.0);
      rhor[n][m] = cof1 * file->rhor[k-1] + cof2 * file->rhor[k] +
        cof3 * file->rhor[k+1] + cof4 * file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I, J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N * (N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl * (nfuncfl+1)/2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nr+1, "pair:z2r");

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  double zri, zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
        r = (m-1) * dr;

        p = r/ifile->dr + 1.0;
        k = static_cast<int> (p);
        k = MIN(k, ifile->nr-2);
        k = MAX(k, 2);
        p -= k;
        p = MIN(p, 2.0);
        cof1 = -sixth * p * (p-1.0) * (p-2.0);
        cof2 = 0.5 * (p * p-1.0) * (p-2.0);
        cof3 = -0.5 * p * (p+1.0) * (p-2.0);
        cof4 = sixth * p * (p * p-1.0);
        zri = cof1 * ifile->zr[k-1] + cof2 * ifile->zr[k] +
          cof3 * ifile->zr[k+1] + cof4 * ifile->zr[k+2];

        p = r/jfile->dr + 1.0;
        k = static_cast<int> (p);
        k = MIN(k, jfile->nr-2);
        k = MAX(k, 2);
        p -= k;
        p = MIN(p, 2.0);
        cof1 = -sixth * p * (p-1.0) * (p-2.0);
        cof2 = 0.5 * (p * p-1.0) * (p-2.0);
        cof3 = -0.5 * p * (p+1.0) * (p-2.0);
        cof4 = sixth * p * (p * p-1.0);
        zrj = cof1 * jfile->zr[k-1] + cof2 * jfile->zr[k] +
          cof3 * jfile->zr[k+1] + cof4 * jfile->zr[k+2];

        z2r[n][m] = 27.2 * 0.529 * zri * zrj;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow, icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}


/*  ----------------------------------------------------------------------  */

void PairEAM::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline, nfrho, nrho+1, 7, "pair:frho");
  memory->create(rhor_spline, nrhor, nr+1, 7, "pair:rhor");
  memory->create(z2r_spline, nz2r, nr+1, 7, "pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho, drho, frho[i], frho_spline[i]);

  for (int i = 0; i < nrhor; i++) 
    interpolate(nr, dr, rhor[i], rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr, dr, z2r[i], z2r_spline[i]);
}

/*  ----------------------------------------------------------------------  */

void PairEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0 * (spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0 * (spline[m+1][6]-spline[m][6]) -
      2.0 * spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0 * (spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0 * spline[m][4]/delta;
    spline[m][0] = 3.0 * spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------------------------- */

void PairEAM::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, m, n, ictype, jctype, ietype, jetype, japc, jnpe, jbasis;
  int inode, node, jucell, iucell, igcell, inpe, iapc, ibasis, jgcell, jnode;
  double r, rhoip, rhojp, z2, z2p, recip, phi, phip, psip, fpair;
  int jnum;
  double evdwl, iescale, ivscale, jescale, jvscale;
  double p, *coeff;
  int *jlist, *jindexlist;
  int iindex, jindex;
  double *ieptr, *jeptr, *ivptr, *jvptr, *irhoptr, *jrhoptr;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, rtmp;
  double *fi, *fj;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;
  int newton_pair = force->newton_pair;
  int *npe = element->npe;
  int *apc = element->apc;
  tagint *atag = atom->tag;
  tagint *etag = element->tag;
  double **ex = element->x;
  double ****nodex = element->nodex;
  double ****gaussf = element->gaussf;
  int *etype = element->etype;
  int **ctype = element->ctype;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  int **g2u = element->g2u;
  int **u2g = element->u2g;
  int **g2n = element->g2n;
  int **is_outer = element->is_outer;
  double ***shape_array = element->shape_array;
  double ***weighted_shape_array = element->weighted_shape_array;
  double *nodal_weight = element->nodal_weight;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *iindexlist = list->iindexlist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **firstneighindex = list->firstneighindex;

  // rhoi: the electron density at gaussian cell
  // rhoj: the electron density at the neighbor site of the gaussian cell

  double rhoi, rhoj;

  //fpi: derivative of embedding energy at gaussian cell
  //fpj: derivative of embedding energy at the neighbor site of the gaussian cell

  double fpi, fpj;

  // grow energy and fp arrays if necessary
  // need to be at least atom->nmax or element->nmax in length

  if (atom->nmax > namax) {
    memory->destroy(atomrho);
    memory->destroy(atomfp);
    namax = atom->nmax;
    memory->create(atomrho, namax, "pair:atomrho");
    memory->create(atomfp, namax, "pair:fp");
  }

  if (element->nmax > nemax) {
    memory->destroy(noderho);
    nemax = element->nmax;
    memory->create(noderho, nemax, element->maxapc, element->maxnpe, "pair:noderho");
  }

  // zero out density

  if (newton_pair) {
    for (i = 0; i < naall; i++) atomrho[i] = 0.0;
    for (i = 0; i < neall; i++) 
      for (j = 0; j < element->maxapc; j++) 
        for (k = 0; k < element->maxnpe; k++) 
          noderho[i][j][k] = 0.0;
  } else {
    for (i = 0; i < nalocal; i++) atomrho[i] = 0.0;
    for (i = 0; i < nelocal; i++) 
      for (j = 0; j < element->maxapc; j++) 
        for (k = 0; k < element->maxnpe; k++) 
          noderho[i][j][k] = 0.0;
  }

  /* ------------------------------------------------------------------- 
   *CALCULATE ELECTRON DENSITY
   *----------------------------------------------------------------- */ 

  // atomrho = electron density at each atom
  // noderho = electron density at each element node

  // calculate ATOM's electron density
  // same as LAMMPS
  // loop over neighbors (atoms and virtual atoms) of my atoms to calculate electron density at atom I

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    iindex = iindexlist[ii];

    // i is atom

    if (iindex < 0) {
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
      ictype = atype[i];
      irhoptr = &atomrho[i];
    } 

    // i is gauss point
    // skip if i is not node

    else {
      ietype = etype[i];
      iapc = apc[ietype];
      igcell = iindex / iapc;
      inode = g2n[ietype][igcell];
      if (inode < 0) continue;
      ibasis = iindex % iapc;
      ictype = ctype[i][ibasis];
      irhoptr = &noderho[i][ibasis][inode];
      xtmp = nodex[i][ibasis][inode][0];
      ytmp = nodex[i][ibasis][inode][1];
      ztmp = nodex[i][ibasis][inode][2];
    }

    // loop over neighbors

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];
    for ( jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      jrhoptr = NULL;
      delx = xtmp; dely = ytmp; delz = ztmp;

      // j is atom

      if (jindex < 0) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
        jctype = atype[j];
        if (newton_pair || j < nalocal) jrhoptr = &atomrho[j];
      } 

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jnpe = npe[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = ctype[j][jbasis];
        for (node = 0; node < jnpe; node++) {
          delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
          dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
          delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
        }
        if (newton_pair || j < nelocal) {
          jgcell = u2g[jetype][jucell];
          if (jgcell >= 0) {
            jnode = g2n[jetype][jgcell];
            if (jnode >= 0) 
              jrhoptr = &noderho[j][jbasis][jnode];
          }
        }
      }
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutforcesq) {
        p = sqrt(rsq) * rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m, nr-1);
        p -= m;
        p = MIN(p, 1.0);
        coeff = rhor_spline[type2rhor[jctype][ictype]][m];
        *irhoptr += ((coeff[3] * p+coeff[4]) * p+coeff[5]) * p+coeff[6];

        if (jrhoptr) {
          coeff = rhor_spline[type2rhor[ictype][jctype]][m];
          *jrhoptr += ((coeff[3] * p+coeff[4]) * p+coeff[5]) * p+coeff[6];
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms), 
  //  will exceed table, so add linear term to conserve energy

  for (i = 0; i < nalocal; i++) {
    ictype = atype[i];
    atomfp[i] = compute_fp(atomrho[i], coeff, p, ictype);
    if (eflag) {
      phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      if (atomrho[i] > rhomax) phi += atomfp[i] * (atomrho[i]-rhomax);
      phi *= scale[ictype][ictype];
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function at atom sites
  // communicate electron density at element node sites

  comm->forward_comm_pair(this);

  /* ----------------------------------------------------
   *CALCULATE FORCE
   *------------------------------------------------- */

  // rhoip = derivative of (density at atom j due to atom i)
  // rhojp = derivative of (density at atom j due to atom j)
  // phi = pair potential energy
  // phip = phi'
  // z2 = phi * r
  // z2p = (phi * r)' = (phi' r) + phi
  // psip needs both fp[i] and fp[j] terms since r_ij appears in two 
  //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
  //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    iindex = iindexlist[ii];
    ieptr = ivptr = NULL;
    iescale = ivscale = 0.0;

    // i is atom

    if (iindex < 0) {
      ictype = atype[i];
      xtmp = ax[i][0];
      ytmp = ax[i][1];
      ztmp = ax[i][2];
      fpi = atomfp[i];
      if (eflag_global) iescale = 0.5;
      if (vflag_global) ivscale = 0.5;
      if (eflag_atom) ieptr = &eatom[i];
      if (vflag_atom) ivptr = vatom[i];
      fi = af[i];
    }

    // i is gauss point

    else {
      ietype = etype[i];
      iapc = apc[ietype];
      inpe = npe[ietype];
      igcell = iindex / iapc;
      ibasis = iindex % iapc;
      ictype = ctype[i][ibasis];
      inode = g2n[ietype][igcell];
      fi = gaussf[i][ibasis][igcell];
      if (inode < 0) {
        iucell = g2u[ietype][igcell];
        rhoi = xtmp = ytmp = ztmp = 0.0;
        for (node = 0; node < inpe; node++) {
          xtmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][0];
          ytmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][1];
          ztmp += shape_array[ietype][iucell][node] * nodex[i][ibasis][node][2];
          rhoi += shape_array[ietype][iucell][node] * noderho[i][ibasis][node];
        }
      } else {
        xtmp = nodex[i][ibasis][inode][0];
        ytmp = nodex[i][ibasis][inode][1];
        ztmp = nodex[i][ibasis][inode][2];
        rhoi = noderho[i][ibasis][inode];
        if (eflag_global) {
          if (nodal_energy_weight < 0) iescale = nodal_weight[ietype] / 2.0;
          else iescale = nodal_energy_weight / 2.0;
        }
        if (vflag_global) ivscale = nodal_weight[ietype] / 2.0;
        if (eflag_atom) ieptr = &enode[i][ibasis][inode];
        if (vflag_atom) ivptr = vnode[i][ibasis][inode];
      }

      // calculate fp at gaussian point i

      fpi = compute_fp(rhoi, coeff, p, ictype);

      if (eflag && inode >= 0) {
        phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
        if (rhoi > rhomax) phi += fpi * (rhoi-rhomax);
        phi *= scale[ictype][ictype];
        if (eflag_global) eng_vdwl += phi * iescale * 2.0;
        if (eflag_atom) enode[i][ibasis][inode] += phi;
      }
    }

    // loop over neighbors

    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    jindexlist = firstneighindex[ii];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jindex = jindexlist[jj];
      delx = xtmp; dely = ytmp; delz = ztmp;
      jescale = jvscale = 0.0;
      fj = jeptr = jvptr = NULL;

      // j is atom

      if (jindex < 0) {
        delx -= ax[j][0];
        dely -= ax[j][1];
        delz -= ax[j][2];
        fpj = atomfp[j];
        jctype = atype[j];
        if (newton_pair || j < nalocal) {
          fj = af[j];
          if (eflag_atom) jeptr = &eatom[j];
          if (vflag_atom) jvptr = vatom[j];
          if (eflag_global) jescale = 0.5;
          if (vflag_global) jvscale = 0.5;
        }
      } 

      // j is virtual atom

      else {
        jetype = etype[j];
        japc = apc[jetype];
        jnpe = npe[jetype];
        jucell = jindex / japc;
        jbasis = jindex % japc;
        jctype = ctype[j][jbasis];
        rhoj = 0.0;
        for (node = 0; node < jnpe; node++) {
          delx -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][0];
          dely -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][1];
          delz -= shape_array[jetype][jucell][node] * nodex[j][jbasis][node][2];
          rhoj += shape_array[jetype][jucell][node] * noderho[j][jbasis][node];
        }

        // calculate fp at virtual atom j

        fpj = compute_fp(rhoj, coeff, p, jctype);

        if (newton_pair || j < nelocal) {
          jgcell = u2g[jetype][jucell];
          if (jgcell >= 0) {
            jnode = g2n[jetype][jgcell];
            fj = gaussf[j][jbasis][jgcell];
            if (jnode >= 0) {
              if (eflag_atom) jeptr = &enode[j][jbasis][jnode];
              if (vflag_atom) jvptr = vnode[j][jbasis][jnode];
              if (eflag_global) {
                if (nodal_energy_weight < 0) jescale = nodal_weight[jetype] / 2.0;
                else jescale = nodal_energy_weight / 2.0;
              }
              if (vflag_global) jvscale = nodal_weight[jetype] / 2.0;
            }
          }
        }
      }

      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutforcesq) {

        // add force contribution from site J

        r = sqrt(rsq);
        p = r * rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m, nr-1);
        p -= m;
        p = MIN(p, 1.0);
        coeff = rhor_spline[type2rhor[ictype][jctype]][m];
        rhoip = (coeff[0] * p + coeff[1]) * p + coeff[2];
        coeff = rhor_spline[type2rhor[jctype][ictype]][m];
        rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];
        coeff = z2r_spline[type2z2r[ictype][jctype]][m];
        z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
        z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
        recip = 1.0/r;
        phi = z2 * recip;
        phip = z2p * recip - phi * recip;
        psip = fpi * rhojp + fpj * rhoip + phip;
        fpair = -scale[ictype][jctype] * psip * recip;

        fi[0] += delx * fpair;
        fi[1] += dely * fpair;
        fi[2] += delz * fpair;
        if (fj != NULL) {
          fj[0] -= delx * fpair;
          fj[1] -= dely * fpair;
          fj[2] -= delz * fpair;
        }

        if (eflag) evdwl = scale[ictype][jctype] * phi;
        if (evflag) ev_tally(iescale + jescale, ivscale + jvscale, 
            ieptr, jeptr, ivptr, jvptr, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------------------------------- */

int PairEAM::pack_atom_forward_comm(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = atomfp[j];
  }
  return m;
}

/* ------------------------------------------------------------------------------------------ */

void PairEAM::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) atomfp[i] = buf[m++];
}


/* --------------------------------------------------------------------------------------------- */

int PairEAM::pack_elem_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, k, l, m;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < apc[etype[j]]; k++) 
      for (l = 0; l < npe[etype[j]]; l++) 
        buf[m++] = noderho[j][k][l];
  }
  return m;
}

/* ---------------------------------------------------------------------------------------------- */

void PairEAM::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i, j, k, m, last;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) 
        noderho[i][j][k] = buf[m++];
}

/*  ----------------------------------------------------------------------  */

int PairEAM::pack_atom_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = atomrho[i];
  return m;
}

/*  ----------------------------------------------------------------------  */

void PairEAM::unpack_atom_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    atomrho[j] += buf[m++];
  }
}

/*  ----------------------------------------------------------------------  */

int PairEAM::pack_elem_reverse_comm(int n, int first, double *buf)
{
  int i, j, k, m, last;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
    for (j = 0; j < apc[etype[i]]; j++) 
      for (k = 0; k < npe[etype[i]]; k++) 
        buf[m++] = noderho[i][j][k];
  return m;
}

/*  ----------------------------------------------------------------------  */

void PairEAM::unpack_elem_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, l, m;
  int *npe = element->npe;
  int *apc = element->apc;
  int *etype = element->etype;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < apc[etype[j]]; k++) 
      for (l = 0; l < npe[etype[j]]; l++) 
        noderho[j][k][l] += buf[m++];
  }
}

/*  ----------------------------------------------------------------------
    compute fp
    -------------------------------------------------------------------------  */

double PairEAM::compute_fp(double rho, double *coeff, double &p, int type)
{
  p = rho * rdrho + 1.0;
  int m = static_cast<int> (p);
  m = MAX(1, MIN(m, nrho-1));
  p -= m;
  p = MIN(p, 1.0);
  coeff = frho_spline[type2frho[type]][m];
  return (coeff[0] * p + coeff[1]) * p + coeff[2];
}

/*  ----------------------------------------------------------------------
    memory usage of local atom/elem-based arrays
    -------------------------------------------------------------------------  */

double PairEAM::memory_usage()
{
  int maxnpe = element->maxnpe;
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom * 6 * sizeof(double);
  bytes += maxeelem * maxnpe * sizeof(double);
  bytes += maxvelem * maxnpe * 6 * sizeof(double);
  bytes += 2 * namax * sizeof(double);
  bytes += nemax * maxnpe * sizeof(double);
  return bytes;
}
