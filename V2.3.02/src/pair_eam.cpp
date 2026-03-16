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

/*--------------------------------------------------------------------------------------------*/

PairEAM::PairEAM(CAC *cac) : Pair(cac)
{
  restartinfo = 0;
  manybody_flag = 1;

  nemax = 0;
  namax = 0;
  atomrho = NULL;
  atomfp = NULL;
  noderho = NULL;

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
  comm_elem_forward = npe;
}

/*----------------------------------------------------------------------------------------
  check if allocated, since class can be destructed when incomplete
  ------------------------------------------------------------------------------------------*/

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
    for (int i = 0; i < nfuncfl;i++) {
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


/*----------------------------------------------------------------------------------------------*/

void PairEAM::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR, "Illegal pair_style command");
}


/*------------------------------------------------------------------------------------------------
  set coeff for one or more type pairs
  read DYNAMO funcfl file
  -------------------------------------------------------------------------------------------------*/

void PairEAM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
  // for now the default results for ilo ihi jlo jhi are 1 1 1 1

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2],funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *)
      memory->srealloc(funcfl,nfuncfl*sizeof(Funcfl),"pair:funcfl");
    read_file(arg[2]);
    int n = strlen(arg[2]) + 1;
    funcfl[ifuncfl].file = new char[n];
    strcpy(funcfl[ifuncfl].file,arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i=j

  int count = 0;
  for (int i = ilo; i <= ihi;i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (i==j) {
        setflag[i][i] = 1;
        map[i] = ifuncfl;
        atom->set_mass(i,funcfl[ifuncfl].mass);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/*------------------------------------------------------------------------------
  allocate all arrays
  -------------------------------------------------------------------------------*/

void PairEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j<= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
  memory->create(scale,n+1,n+1,"pair:scale");
}

/*---------------------------------------------------------------------------------
  read potential values from a DYNAMO single element funcfl file
  -----------------------------------------------------------------------------------*/

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
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  int tmp;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg", &tmp, &file->mass);
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %lg",
        &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->mass,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  memory->create(file->frho,(file->nrho+1),"pair:frho");
  memory->create(file->rhor,(file->nr+1),"pair:rhor");
  memory->create(file->zr,(file->nr+1),"pair:zr");

  if (me == 0) grab(fptr,file->nrho,&file->frho[1]);
  MPI_Bcast(&file->frho[1],file->nrho,MPI_DOUBLE,0,world);

  // if (me=0) grab(fptr,file->nr,&file->zr[1]);
  // MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) grab(fptr,file->nr,&file->rhor[1]);
  MPI_Bcast(&file->rhor[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) grab(fptr,file->nr,&file->zr[1]);
  MPI_Bcast(&file->zr[1],file->nr,MPI_DOUBLE,0,world);

  if (me == 0) fclose(fptr); 
}

/*------------------------------------------------------------------------------
  grab n values from file fp and put them in list
  values can be several to a line
  only called by proc 0
  --------------------------------------------------------------------------------*/

void PairEAM::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/*----------------------------------------------------------------------------------
  init specific to this pair style
  ------------------------------------------------------------------------------------*/

void PairEAM::init_style()
{

  // convert read-in file(s) to arrays and spline them
  
  file2array();
  array2spline();
  neighbor->request(this,instance_me);
}

/*----------------------------------------------------------------------------------
  init for one type pair i,j and corresponding j,i
  ------------------------------------------------------------------------------------*/

double PairEAM::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setflag[i][j] == 0) scale [i][j] = 1.0;
  scale[j][i] = scale[i][j];

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax,funcfl[m].cut);
  } else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;

  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/*-----------------------------------------------------------------------------------
  convert read_in funcfl potential(s) to standard array format
  interpolate all file values to a single grid and cutoff
  -------------------------------------------------------------------------------------*/

void PairEAM::file2array()
{
  int i,j,k,m,n;
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
    dr = MAX(dr,file->dr);
    drho = MAX(drho,file->drho);
    rmax = MAX(rmax,(file->nr-1)*file->dr);
    rhomax = MAX(rhomax,(file->nrho-1)*file->drho);
  }

  // set nr, nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);


  //--------------------------------------------------------------------
  // setup frho arrays
  //--------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r, p, cof1,cof2,cof3,cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/file->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      frho[n][m] = cof1*file->frho[k-1] + cof2*file->frho[k] +
        cof3*file->frho[k+1] + cof4*file->frho[k+2];
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

  //--------------------------------------------------------------------------
  // setup rhor arrays
  //--------------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/file->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,file->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -sixth*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = sixth*p*(p*p-1.0);
      rhor[n][m] = cof1*file->rhor[k-1] + cof2*file->rhor[k] +
        cof3*file->rhor[k+1] + cof4*file->rhor[k+2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  //--------------------------------------------------------------------
  // setup z2r arrays
  //--------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl*(nfuncfl+1)/2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r,nr+1,"pair:z2r");

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff

  double zri,zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
        r = (m-1)*dr;
        p = r/jfile->dr + 1.0;
        k = static_cast<int> (p);
        k = MIN(k,jfile->nr-2);
        k = MAX(k,2);
        p -= k;
        p = MIN(p,2.0);
        cof1 = -sixth*p*(p-1.0)*(p-2.0);
        cof2 = 0.5*(p*p-1.0)*(p-2.0);
        cof3 = -0.5*p*(p+1.0)*(p-2.0);
        cof4 = sixth*p*(p*p-1.0);
        zrj = cof1*jfile->zr[k-1] + cof2*jfile->zr[k] +
          cof3*jfile->zr[k+1] + cof4*jfile->zr[k+2];
        z2r[n][m] = zrj;

      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid)
  //  type2z2r is not used by non-opt
  //  but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol ==-1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n+=m+1;
      n+=icol;
      type2z2r[i][j] = n;
    }
  }

  //error->all(FLERR, "test from file2array");
}

/*----------------------------------------------------------------------------------------*/

void PairEAM::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr, dr,rhor[i],rhor_spline[i]);

  for (int i=0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/*-----------------------------------------------------------------------------------------*/

void PairEAM::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5*(spline[3][6] - spline[1][6]);
  spline[n-1][5] = 0.5*(spline[n][6] - spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6] - spline[m+2][6]) +
        8.0*(spline[m+1][6]-spline[m-1][6]))/12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6] - spline[m][6])-
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] - 
      2.0*(spline[m+1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m<=n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/*-----------------------------------------------------------------------------------------*/

void PairEAM::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,m,n,ictype,jctype,ietype,jetype;
  int iinode,inode,jnode,jintpl,iintpl,iintg,iintpl_local,iintg_local;
  double r,rhoip,rhojp,z2,z2p,recip,phi,phip,psip,fpair;
  bigint jnum;
  double evdwl;
  double p,*coeff;
  int *jlist,*jindexlist;
  double ix,iy,iz,jx,jy,jz,delx,dely,delz,rsq;
  double fx,fy,fz;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  double **ax = atom->x;
  double **af = atom->f;
  int *atype = atom->type;
  int nalocal = atom->nlocal;
  int naall = nalocal + atom->nghost;
  int newton_pair = force->newton_pair;
  tagint *atag = atom->tag;
  tagint *etag = element->tag;
  double **ex = element->x;
  double ***nodex = element->nodex;
  double ***nodef = element->nodef;
  double **weight = element->weight;
  int *nintg = element->nintg;
  int *nintpl = element->nintpl;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int nelocal = element->nlocal;
  int neall = nelocal + element->nghost;
  int **i2ia = element->i2ia;
  int **i2n = element->i2n;
  double ***shape_array = element->shape_array;
  double *fscale = new double[element->netypes+1];
  for (i = 1; i <= element->netypes; i++)
    fscale[i] = (static_cast<double>(npe))/nintpl[i];

  int ainum = list->ainum;
  int einum = list->einum;
  int iinum = list->iinum;
  int *ailist = list->ailist;
  int *eilist = list->eilist;
  int **n2ilist = list->n2ilist;
  int *e2ilist = list->e2ilist;
  int *numneigha2a = list->numneigha2a;
  int *numneigha2e = list->numneigha2e;
  int *numneigha2ia = list->numneigha2ia;
  int *numneighi2a = list->numneighi2a;
  int *numneighe2e = list->numneighe2e;
  int *numneighi2ia = list->numneighi2ia;
  int **firstneigha2a = list->firstneigha2a;
  int **firstneigha2e = list->firstneigha2e;
  int **firstneigha2ia = list->firstneigha2ia;
  int **firstneigha2ia_index = list->firstneigha2ia_index;
  int **firstneighi2a = list->firstneighi2a;
  int **firstneighe2e = list->firstneighe2e;
  int **firstneighi2ia = list->firstneighi2ia;
  int **firstneighi2ia_index = list->firstneighi2ia_index;


  // rhoi: the electron density at integration point
  // rhoj: the electron density at the neighbor site of the integration point
  
  double rhoi, rhoj;

  //fpi: derivative of embedding energy at integration point
  //fpj: derivative of embedding energy at the neighbor site of the integration point
  
  double fpi, fpj;

  // grow energy and fp arrays if necessary
  // need to be at least atom->nmax or element->nmax in length

  if (atom->nmax > namax) {
    memory->destroy(atomrho);
    memory->destroy(atomfp);
    namax = atom->nmax;
    memory->create(atomrho,namax,"pair:atomrho");
    memory->create(atomfp,namax,"pair:fp");
  }

  if (element->nmax > nemax){
    memory->destroy(noderho);
    nemax = element->nmax;
    memory->create(noderho,nemax,npe,"pair:noderho");
  }

  // zero out density

  if (newton_pair) {
    for (i = 0; i < naall; i++) atomrho[i] = 0.0;
    for (i = 0; i < neall; i++) {
      for (j = 0; j < npe; j++) {
        noderho[i][j] = 0.0;
      }
    }
  } else {
    for (i = 0; i < nalocal; i++) atomrho[i] = 0.0;
    for (i = 0; i < nelocal; i++) {
      for (j = 0; j < npe; j++) {
        noderho[i][j] = 0.0;
      }
    }
  }


  /*------------------------------------------------------------------- 
   * CALCULATE ELECTRON DENSITY
   *-----------------------------------------------------------------*/ 

  // atomrho = electron density at each atom
  // noderho = electron density at each element node

  // calculate ATOM's electron density
  // same as LAMMPS
  // loop over neighbors (atoms and interpolated atoms) of my atoms to calculate electron density at atom i


  for (ii = 0; ii < ainum; ii++) {

    i = ailist[ii];

    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];
    ictype = atype[i];

    // loop over neighboring atoms
   
    jnum = numneigha2a[i];
    if (jnum) {
      jlist = firstneigha2a[i];
      for ( jj = 0; jj < jnum;jj++) {
        j = jlist[jj];
        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutforcesq) {
          jctype = atype[j];
          p = sqrt(rsq)*rdr + 1.0;
          m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,1.0);
          coeff = rhor_spline[type2rhor[jctype][ictype]][m];
          atomrho[i] += ((coeff[3]*p+coeff[4])*p+coeff[5])*p+coeff[6];
        }
      }
    }

    // loop over neighboring interpolated atoms 

    jnum = numneigha2ia[i];
    if (jnum) {
      jlist = firstneigha2ia[i];
      jindexlist = firstneigha2ia_index[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jetype = etype[j];
        jintpl = jindexlist[jj];

        // interpolate positions at site J

        jx = jy = jz = 0.0;
        for (jnode = 0; jnode < npe; jnode++) {
          jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
          jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
          jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
        }

        delx = ix - jx; 
        dely = iy - jy; 
        delz = iz - jz; 
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutforcesq) {
          jctype = ctype[j];
          p = sqrt(rsq)*rdr + 1.0;
          m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,1.0);
          coeff = rhor_spline[type2rhor[jctype][ictype]][m];
          atomrho[i] += ((coeff[3]*p+coeff[4])*p+coeff[5])*p+coeff[6];
        }
      }
    }
  }

  // calculate ELEMENT NODE's electron density
  // 1. Loop over neighbor of element nodes, 
  //    use iintglist to get the corresponding integration 
  //    points and full the corresponding neighbor list
  // 2. Calculate electron density at nodes (noderho)
  // loop over neighbors (atoms and interpolated atoms) of my element nodes to calculate electron density at node k in element i

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii]; 
    ietype = etype[i];
    ictype = ctype[i];

    for (k = 0; k < npe;k++) {

      ix = nodex[i][k][0];
      iy = nodex[i][k][1];
      iz = nodex[i][k][2];

      iintg = n2ilist[i][k]; 

      // loop over neighboring atoms j to add electron density contribution to node k in element i

      jnum = numneighi2a[iintg];    
      if (jnum) {
        jlist = firstneighi2a[iintg];
        for (jj = 0; jj < jnum;jj++){
          j = jlist[jj];
          delx = ix - ax[j][0];
          dely = iy - ax[j][1];
          delz = iz - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutforcesq) {
            jctype = atype[j];
            p = sqrt(rsq)*rdr + 1.0;
            m = static_cast<int> (p);
            m = MIN(m,nr-1);
            p -= m;
            p = MIN(p,1.0);
            coeff = rhor_spline[type2rhor[jctype][ictype]][m];
            noderho[i][k] += ((coeff[3]*p+coeff[4])*p+coeff[5])*p+coeff[6];
          }
        }
      }

      // loop over neighboring interpolated atoms jintpl in element j to add electron density contribution to node k in element i

      jnum = numneighi2ia[iintg];
      if (jnum) {
        jlist = firstneighi2ia[iintg];
        jindexlist = firstneighi2ia_index[iintg];
        for (jj = 0; jj < jnum;jj++){
          j = jlist[jj];
          jintpl = jindexlist[jj];
          jetype = etype[j];

          // interpolate positions at site J

          jx = jy = jz = 0.0;
          for (jnode = 0; jnode < npe; jnode++) {
            jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
            jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
            jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
          }

          delx = ix - jx; 
          dely = iy - jy; 
          delz = iz - jz; 
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutforcesq) {
            jctype = ctype[j];
            p = sqrt(rsq)*rdr + 1.0;
            m = static_cast<int> (p);
            m = MIN(m,nr-1);
            p -= m;
            p = MIN(p,1.0);
            coeff = rhor_spline[type2rhor[jctype][ictype]][m];
            noderho[i][k] += ((coeff[3]*p+coeff[4])*p+coeff[5])*p+coeff[6];
          }
        }
      }
    }
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //  will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];
    ictype = atype[i];
    p = atomrho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[ictype]][m];
    atomfp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
  }


  // communicate derivative of embedding function at atom sites
  // communicate electron density at element node sites

  comm->forward_comm_pair(this);


  /*----------------------------------------------------
   * CALCULATE FORCE
   * -------------------------------------------------*/

  // rhoip = derivative of (density at atom j due to atom i)
  // rhojp = derivative of (density at atom j due to atom j)
  // phi = pair potential energy
  // phip = phi'
  // z2 = phi*r
  // z2p = (phi * r)' = (phi' r) + phi
  // psip needs both fp[i] and fp[j] terms since r_ij appears in two 
  //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
  //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip


  // compute forces on each ATOM
  // loop over neighbors of my atoms
  // same as LAMMPS

  for (ii = 0; ii < ainum; ii++) {
    i = ailist[ii];
    ictype = atype[i];
    ix = ax[i][0];
    iy = ax[i][1];
    iz = ax[i][2];
    fpi = atomfp[i];

    // loop over neighboring atoms

    jnum = numneigha2a[i];
    if (jnum) {
      jlist = firstneigha2a[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        delx = ix - ax[j][0];
        dely = iy - ax[j][1];
        delz = iz - ax[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutforcesq) {
          fpj = atomfp[j];
          jctype = atype[j];

          // add force contribution from site J

          r = sqrt(rsq);
          p = r*rdr + 1.0;
          m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,1.0);
          coeff = rhor_spline[type2rhor[ictype][jctype]][m];
          rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
          coeff = rhor_spline[type2rhor[jctype][ictype]][m];
          rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
          coeff = z2r_spline[type2z2r[ictype][jctype]][m];
          z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
          z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          recip = 1.0/r;
          phi = z2*recip;
          phip = z2p*recip - phi*recip;
          psip = fpi*rhojp + fpj*rhoip + phip;
          fpair = -scale[ictype][jctype]*psip*recip;

          af[i][0] += delx*fpair;
          af[i][1] += dely*fpair;
          af[i][2] += delz*fpair;

          if (vflag) atom_v_tally(i,fpair,delx,dely,delz);
        }
      }
    }

    // loop over neighboring interpolated atoms

    jnum = numneigha2ia[i];
    if (jnum){
      jlist = firstneigha2ia[i];
      jindexlist = firstneigha2ia_index[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jintpl = jindexlist[jj];
        jetype = etype[j];

        // interpolate positions at site J from nodex

        jx = jy = jz = 0.0;
        for (jnode = 0; jnode < npe; jnode++) {
          jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
          jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
          jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
        }

        delx = ix - jx; 
        dely = iy - jy; 
        delz = iz - jz; 
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutforcesq) {
          jctype = ctype[j];

          // interpolate rho at site J from noderho

          rhoj = 0.0;
          for (jnode = 0; jnode < npe; jnode++) 
            rhoj += shape_array[jetype][jintpl][jnode]*noderho[j][jnode];

          // calculate fp at site J

          p = rhoj*rdrho + 1.0;
          m = static_cast<int> (p);
          m = MAX(1,MIN(m,nrho-1));
          p -= m;
          p = MIN(p,1.0);
          coeff = frho_spline[type2frho[jctype]][m];
          fpj = (coeff[0]*p+coeff[1])*p + coeff[2];

          // add force contribution from site J

          r = sqrt(rsq);
          p = r*rdr + 1.0;
          m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,1.0);
          coeff = rhor_spline[type2rhor[ictype][jctype]][m];
          rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
          coeff = rhor_spline[type2rhor[jctype][ictype]][m];
          rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
          coeff = z2r_spline[type2z2r[ictype][jctype]][m];
          z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
          z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          recip = 1.0/r;
          phi = z2*recip;
          phip = z2p*recip - phi*recip;
          psip = fpi*rhojp + fpj*rhoip + phip;
          fpair = -scale[ictype][jctype]*psip*recip;

          af[i][0] += delx*fpair;
          af[i][1] += dely*fpair;
          af[i][2] += delz*fpair;

          if (vflag) atom_v_tally(i,fpair,delx,dely,delz);
        }
      } 
    }
  }

  // compute forces on each NODE
  // compute forces on each integration points and distribute to nodes
  // loop over my integration points iing in element i
  // fx fy fz = force at current integration point.
  // use e2ilist to get the local index of the first integration point
  //   of element i
  // 1. Interpolate electron density from nodes to current integration point
  // 2. Loop over neighboring atoms. Add force contribution to fx fy fz.
  // 3. Loop over neighboring interpolated atoms. Add force contribution to fx fy fz
  // 4. After getting total forces at all integration points in element i, 
  //    distribute forces at integration points to nodef in element i.
  // 5. Scale down forces on nodes 

  for (ii = 0; ii < einum; ii++) {
    i = eilist[ii];
    ietype = etype[i];
    ictype = ctype[i];

    for (iintg = 0; iintg < nintg[ietype]; iintg++){

      // initialize force values

      fx = fy = fz = 0.0;

      // interpolate postions and rho of integration point from node values

      iintpl = i2ia[ietype][iintg];
      iinode = i2n[ietype][iintg];
      ix = iy = iz = 0.0;
      rhoi = 0.0;
      for (inode = 0; inode < npe; inode++) {
        ix += shape_array[ietype][iintpl][inode]*nodex[i][inode][0];
        iy += shape_array[ietype][iintpl][inode]*nodex[i][inode][1];
        iz += shape_array[ietype][iintpl][inode]*nodex[i][inode][2];
        rhoi += shape_array[ietype][iintpl][inode]*noderho[i][inode];
      }

      // calculate fp at integration point

      p = rhoi*rdrho + 1.0;
      m = static_cast<int> (p);
      m = MAX(1,MIN(m,nrho-1));
      p -= m;
      p = MIN(p,1.0);
      coeff = frho_spline[type2frho[ictype]][m];
      fpi = (coeff[0]*p+coeff[1])*p + coeff[2];

      // loop over neighboring atoms and add force contribution

      iintg_local = e2ilist[ii] + iintg;     
      jnum = numneighi2a[iintg_local];
      if (jnum) {
        jlist = firstneighi2a[iintg_local];
        for (jj = 0; jj < jnum; jj++){
          j = jlist[jj];
          delx = ix - ax[j][0];
          dely = iy - ax[j][1];
          delz = iz - ax[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutforcesq) {
            jctype = atype[j];
            fpj = atomfp[j];

            // add force contribution from site J

            r = sqrt(rsq);
            p = r*rdr + 1.0;
            m = static_cast<int> (p);
            m = MIN(m,nr-1);
            p -= m;
            p = MIN(p,1.0);
            coeff = rhor_spline[type2rhor[ictype][jctype]][m];
            rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
            coeff = rhor_spline[type2rhor[jctype][ictype]][m];
            rhojp = (coeff[0]*p + coeff[1])*p+coeff[2];
            coeff = z2r_spline[type2z2r[ictype][jctype]][m];
            z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
            z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
            recip = 1.0/r;
            phi = z2*recip;
            phip = z2p*recip - phi*recip;
            psip = fpi*rhojp + fpj*rhoip + phip;
            fpair = -scale[ictype][jctype]*psip*recip;

            fx += delx*fpair;
            fy += dely*fpair;
            fz += delz*fpair;

            // tally virial if integration point is a node
            
            if (vflag && iinode >= 0) node_v_tally(i,iinode,fpair,delx,dely,delz,fscale[ietype]);
          }
        }
      }

      // loop over neighboring interpolated atoms and add force contribution

      jnum = numneighi2ia[iintg_local];
      if (jnum){
        jlist = firstneighi2ia[iintg_local]; 
        jindexlist = firstneighi2ia_index[iintg_local];
        for (jj = 0; jj < jnum; jj++){
          j = jlist[jj];
          jintpl = jindexlist[jj];
          jetype = etype[j];

          // interpolating positions at site J from nodex

          jx = jy = jz = 0.0;
          for (jnode = 0; jnode < npe; jnode++) {
            jx += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][0];
            jy += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][1];
            jz += shape_array[jetype][jintpl][jnode]*nodex[j][jnode][2];
          }
          delx = ix - jx; 
          dely = iy - jy; 
          delz = iz - jz; 
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutforcesq) {
            jctype = ctype[j];

            // interpolate rho at site J from noderho

            rhoj = 0.0;
            for (jnode = 0; jnode < npe; jnode++) 
              rhoj += shape_array[jetype][jintpl][jnode]*noderho[j][jnode];

            // calculate fp at site J

            p = rhoj*rdrho + 1.0;
            m = static_cast<int> (p);
            m = MAX(1,MIN(m,nrho-1));
            p -= m;
            p = MIN(p,1.0);
            coeff = frho_spline[type2frho[jctype]][m];
            fpj = (coeff[0]*p+coeff[1])*p + coeff[2];

            // add force contribution from site J

            r = sqrt(rsq);
            p = r*rdr + 1.0;
            m = static_cast<int> (p);
            m = MIN(m,nr-1);
            p -= m;
            p = MIN(p,1.0);
            coeff = rhor_spline[type2rhor[ictype][jctype]][m];
            rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
            coeff = rhor_spline[type2rhor[jctype][ictype]][m];
            rhojp = (coeff[0]*p + coeff[1])*p+coeff[2];
            coeff = z2r_spline[type2z2r[ictype][jctype]][m];
            z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
            z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
            recip = 1.0/r;
            phi = z2*recip;
            phip = z2p*recip - phi*recip;
            psip = fpi*rhojp + fpj*rhoip + phip;
            fpair = -scale[ictype][jctype]*psip*recip;

            fx += delx*fpair;
            fy += dely*fpair;
            fz += delz*fpair;

            // tally virial if integration point is a node
            
            if (vflag && iinode >= 0) node_v_tally(i,iinode,fpair,delx,dely,delz,fscale[ietype]);
          }
        }
      }

      // distribute forces at integration point to forces at nodes 

      for (inode = 0; inode < npe; inode++) {
        nodef[i][inode][0] += weight[ietype][iintg]*shape_array[ietype][iintpl][inode]*fx;
        nodef[i][inode][1] += weight[ietype][iintg]*shape_array[ietype][iintpl][inode]*fy;
        nodef[i][inode][2] += weight[ietype][iintg]*shape_array[ietype][iintpl][inode]*fz;
      }
    }

    // scale back forces on node

    for (inode = 0; inode < npe; inode++) {
      nodef[i][inode][0] *= fscale[ietype];
      nodef[i][inode][1] *= fscale[ietype];
      nodef[i][inode][2] *= fscale[ietype];
    }
  }

  // clean up
  
  delete [] fscale;
}

/*-----------------------------------------------------------------------------------------------*/

int PairEAM::pack_atom_forward_comm(int n, int *list, double *buf,
    int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = atomfp[j];
  }
  return m;
}
/*------------------------------------------------------------------------------------------*/

void PairEAM::unpack_atom_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) atomfp[i] = buf[m++];
}


/*---------------------------------------------------------------------------------------------*/

int PairEAM::pack_elem_forward_comm(int n, int *list, double *buf,
    int pbc_flag, int *pbc)
{
  int i,j,m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (int k = 0; k < npe; k++) buf[m++] = noderho[j][k];
  }
  return m;
}

/*----------------------------------------------------------------------------------------------*/

void PairEAM::unpack_elem_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (int k = 0; k < npe; k++) noderho[i][k] = buf[m++];
}

/* ----------------------------------------------------------------------
   compute force from atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairEAM::pair_a2a(int i, int j, double *f)
{
  int m;
  double r,rsq,rhoip,rhojp,z2,z2p,recip,phi,phip,psip,fpair;
  double p,*coeff;
  int ictype = atom->type[i];
  int jctype = atom->type[j];
  double **x = atom->x;
  double delx,dely,delz;

  // check if distance < force cutoff
  // return 0 force otherwise
  
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutforcesq) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate force
  
  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);
  coeff = rhor_spline[type2rhor[ictype][jctype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jctype][ictype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[ictype][jctype]][m];
  z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = atomfp[i]*rhojp + atomfp[j]*rhoip + phip;
  fpair = -scale[ictype][jctype]*psip*recip;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;
}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on atom I
   ------------------------------------------------------------------------- */

void PairEAM::pair_a2ia(int i, int j, int jintpl, double *jx, double *f)
{
  int m,node;
  double r,rsq,rhoip,rhojp,z2,z2p,recip,phi,phip,psip,fpair;
  double p,*coeff,fpj;
  double rhoj;
  double delx,dely,delz;
  int ictype = atom->type[i];
  int jctype = element->ctype[j];
  int jetype = element->etype[j];
  double ***shape_array = element->shape_array;
  int npe = element->npe;

  // calculate rho at site J from noderho

  for (node = 0; node < npe; node++) 
    rhoj += shape_array[jetype][jintpl][node]*noderho[j][node];

  // check if distance < force cutoff
  // return 0 force otherwise
  
  delx = atom->x[i][0] - jx[0];
  dely = atom->x[i][1] - jx[1];
  delz = atom->x[i][2] - jx[2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutforcesq) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate fpj
  
  p = rhoj*rdrho + 1.0;
  m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  p = MIN(p,1.0);
  coeff = frho_spline[type2frho[jctype]][m];
  fpj = (coeff[0]*p+coeff[1])*p + coeff[2];

  // calculate force 
  
  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);
  coeff = rhor_spline[type2rhor[ictype][jctype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jctype][ictype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[ictype][jctype]][m];
  z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = atomfp[i]*rhojp + fpj*rhoip + phip;
  fpair = -scale[ictype][jctype]*psip*recip;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;

}

/* ----------------------------------------------------------------------
   compute force from interpolated atom J acting on interpolated atom I
   ------------------------------------------------------------------------- */

void PairEAM::pair_ia2ia(int i, int iintpl, double *ix, int j, int jintpl, double *jx, double *f)
{
  int m,node;
  double r,rsq,rhoip,rhojp,z2,z2p,recip,phi,phip,psip,fpair;
  double p,*coeff,fpi,fpj;
  double rhoi,rhoj;
  double delx,dely,delz;
  int ictype = element->ctype[i];
  int ietype = element->etype[i];
  int jctype = element->ctype[j];
  int jetype = element->etype[j];
  double ***shape_array = element->shape_array;
  int npe = element->npe;

  // calculate rho at site I & J from noderh
  
  rhoi = rhoj = 0.0;
  for (node = 0; node < npe; node++) {
    rhoi += shape_array[ietype][iintpl][node]*noderho[i][node];
    rhoj += shape_array[jetype][jintpl][node]*noderho[j][node];
  }

  // check if distance < force cutoff
  // return 0 force otherwise

  delx = ix[0] - jx[0];
  dely = ix[1] - jx[1];
  delz = ix[2] - jx[2];
  rsq = delx*delx+dely*dely+delz*delz;
  if (rsq >= cutforcesq) {
    f[0] = f[1] = f[2] = 0.0;
    return;
  }

  // calculate fpi & fpj

  p = rhoi*rdrho + 1.0;
  m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  p = MIN(p,1.0);
  coeff = frho_spline[type2frho[ictype]][m];
  fpi = (coeff[0]*p+coeff[1])*p + coeff[2];

  p = rhoj*rdrho + 1.0;
  m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  p = MIN(p,1.0);
  coeff = frho_spline[type2frho[jctype]][m];
  fpj = (coeff[0]*p+coeff[1])*p + coeff[2];

  // calculate force 

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);
  coeff = rhor_spline[type2rhor[ictype][jctype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jctype][ictype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[ictype][jctype]][m];
  z2p = (coeff[0]*p + coeff[1])*p+coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fpi*rhojp + fpj*rhoip + phip;
  fpair = -scale[ictype][jctype]*psip*recip;

  f[0] = delx*fpair;
  f[1] = dely*fpair;
  f[2] = delz*fpair;
}

/* ----------------------------------------------------------------------
   memory usage of local atom/elem-based arrays
   ------------------------------------------------------------------------- */

double PairEAM::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom * 6 * sizeof(double);
  bytes += 2 * namax * sizeof(double);
  bytes += nemax * npe *  sizeof(double);
  return bytes;
}
