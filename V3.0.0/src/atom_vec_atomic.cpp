#include <stdlib.h>
#include "atom_vec_atomic.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace CAC_NS;

/* ---------------------------------------------------------------- */

AtomVecAtomic::AtomVecAtomic(CAC *cac) : AtomVec(cac)
{
  mass_type = 1;
  comm_x_only = comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 7;
  size_velocity = 3;
  size_data_atom = 5;
  size_tecplot_atom = 6;
  size_data_vel = 4;
  xcol_data = 3;
}

/*  ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
   -------------------------------------------------------------------------  */

void AtomVecAtomic::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR, "Per-processor system is too big");

  tag = memory->grow(atom->tag, nmax, "atom:tag");
  type = memory->grow(atom->type, nmax, "atom:type");
  mask = memory->grow(atom->mask, nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");
  x = memory->grow(atom->x, nmax, 3, "atom:x");
  if (atom->atom_strain_flag)
    x_current = memory->grow(atom->x_current, nmax, 3, "atom:x_current");
  v = memory->grow(atom->v, nmax, 3, "atom:v");
  f = memory->grow(atom->f, nmax * comm->nthreads, 3, "atom:f");
  grain_tag = memory->grow(atom->grain_tag, nmax, "atom:grain_tag");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_atom_arrays(nmax);
}

/*  ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
   -------------------------------------------------------------------------  */

void AtomVecAtomic::data_atom(double *coord, imageint imagetmp, int grain_id, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (tag[nlocal] <= 0)
    error->one(FLERR, "Invalid atom ID in Atoms section of data file");
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR, "Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;
  grain_tag[nlocal] = grain_id;

  // 1 equal to the all group, 3 equal to atom group

  mask[nlocal] = 1|4;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/*  ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
   -------------------------------------------------------------------------  */

void AtomVecAtomic::data_atom_strain(double *coord, double *current, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (tag[nlocal] <= 0)
    error->one(FLERR, "Invalid atom ID in Atoms section of data file");
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR, "Invalid atom type in Atoms section of data file");
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  x_current[nlocal][0] = current[0];
  x_current[nlocal][1] = current[1];
  x_current[nlocal][2] = current[2];
  image[nlocal] = imagetmp;

  // 1 equal to the all group, 3 equal to atom group

  mask[nlocal] = 1|4;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/*  ----------------------------------------------------------------------
   copy atom I info to atom J
   -------------------------------------------------------------------------  */

void AtomVecAtomic::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  if (atom->atom_strain_flag) {
    x_current[j][0] = x_current[i][0];
    x_current[j][1] = x_current[i][1];
    x_current[j][2] = x_current[i][2];
  }
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2]; 

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_atom_arrays(i, j, delflag);

}

/*  ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
   -------------------------------------------------------------------------  */

int AtomVecAtomic::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  if (atom->atom_strain_flag) {
    buf[m++] = x_current[i][0];
    buf[m++] = x_current[i][1];
    buf[m++] = x_current[i][2];
  }
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = ubuf(grain_tag[i]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_atom_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);


  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  if (atom->atom_strain_flag) {
    x_current[nlocal][0] = buf[m++];
    x_current[nlocal][1] = buf[m++];
    x_current[nlocal][2] = buf[m++];
  }
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  grain_tag[nlocal] = (int) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_atom_exchange(nlocal, &buf[m]);

  atom->nlocal++;
  return m;
}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::pack_border_vel(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;
  double dx, dy, dz, dvx, dvy, dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    }
  }
  return m;

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_atom_border(n, list, &buf[m]);

}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::pack_border(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;
  double dx, dy, dz;
  double dx_current, dy_current, dz_current;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(grain_tag[j]).d;
      if (atom->atom_strain_flag) {
        buf[m++] = x_current[j][0];
        buf[m++] = x_current[j][1];
        buf[m++] = x_current[j][2];
        buf[m++] = ubuf(image[j]).d;
      }

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
      dx_current = pbc[0] * domain->xprd_current;
      dy_current = pbc[1] * domain->yprd_current;
      dz_current = pbc[2] * domain->zprd_current;
    } else {
      dx = dx_current = pbc[0];
      dy = dy_current = pbc[1];
      dz = dz_current = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = ubuf(grain_tag[j]).d;
      if (atom->atom_strain_flag) {
        buf[m++] = x_current[j][0] + dx_current;
        buf[m++] = x_current[j][1] + dy_current;
        buf[m++] = x_current[j][2] + dz_current;
        buf[m++] = ubuf(image[j]).d;
      }

    }
  }

  return m;

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_atom_border(n, list, &buf[m]);
}

/*  ----------------------------------------------------------------------  */

void AtomVecAtomic::unpack_border(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    while (i >= nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    grain_tag[i] = (int) ubuf(buf[m++]).i;
    if (atom->atom_strain_flag) {
      x_current[i][0] = buf[m++];
      x_current[i][1] = buf[m++];
      x_current[i][2] = buf[m++];
      image[i] = (imageint) ubuf(buf[m++]).i;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_atom_border(n, first, &buf[m]);
}

/*  ----------------------------------------------------------------------  */

void AtomVecAtomic::unpack_border_vel(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    while (i >= nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_atom_border(n, first, &buf[m]);
}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::pack_comm(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
    }
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::pack_comm_vel(int n, int *list, double *buf, 
    int pbc_flag, int *pbc)
{
  int i, j, m;
  double dx, dy, dz, dvx, dvy, dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

void AtomVecAtomic::unpack_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}

/*  ----------------------------------------------------------------------  */

void AtomVecAtomic::unpack_comm_vel(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/*  ----------------------------------------------------------------------  */

int AtomVecAtomic::pack_reverse(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/*  ----------------------------------------------------------------------  */

void AtomVecAtomic::unpack_reverse(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/*  ----------------------------------------------------------------------
   return # of bytes of allocated memory
   -------------------------------------------------------------------------  */

bigint AtomVecAtomic::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag, nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type, nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask, nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image, nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x, nmax, 3);
  if (atom->memcheck("v")) bytes += memory->usage(v, nmax, 3);
  if (atom->memcheck("f")) bytes += memory->usage(f, nmax * comm->nthreads, 3);
  if (atom->memcheck("grain_tag")) bytes += memory->usage(grain_tag, nmax);

  return bytes;
}

/*  ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
   -------------------------------------------------------------------------  */

void AtomVecAtomic::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e\n", 
        (tagint) ubuf(buf[i][0]).i, (int) ubuf(buf[i][1]).i, 
        buf[i][2], buf[i][3], buf[i][4]);
}

/*  ----------------------------------------------------------------------
   pack atom info for data file including
   -------------------------------------------------------------------------  */

void AtomVecAtomic::pack_data(double **buf)
{
  for (int i = 0; i < atom->nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = x[i][0];
    buf[i][3] = x[i][1];
    buf[i][4] = x[i][2];
  }
}

/*  ----------------------------------------------------------------------
   pack atom info for data file from write_data_elem command
   -------------------------------------------------------------------------  */

void AtomVecAtomic::pack_data_elem(double **buf, int ghostflag)
{
  if (!ghostflag) {
    for (int i = 0; i < atom->nlocal; i++) {
      buf[i][0] = ubuf(tag[i]).d;
      buf[i][1] = ubuf(type[i]).d;
      buf[i][2] = x[i][0];
      buf[i][3] = x[i][1];
      buf[i][4] = x[i][2];
    }
  } else {
    for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
      buf[i][0] = ubuf(i + 1).d;
      buf[i][1] = ubuf(type[i]).d;
      buf[i][2] = x[i][0];
      buf[i][3] = x[i][1];
      buf[i][4] = x[i][2];
    }
  }
}


/*  ----------------------------------------------------------------------
   write atom info to dat file
   -------------------------------------------------------------------------  */

void AtomVecAtomic::write_tecplot(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp, "%-1.16e %-1.16e %-1.16e " TAGINT_FORMAT " %d %d\n", 
        buf[i][0], buf[i][1], buf[i][2], (tagint) ubuf(buf[i][3]).i
        , (int) ubuf(buf[i][4]).i, (int) ubuf(buf[i][5]).i);
}


/*  ----------------------------------------------------------------------
   pack atom info to dat file 
   -------------------------------------------------------------------------  */

void AtomVecAtomic::pack_tecplot(double **buf)
{
  for (int i = 0; i < atom->nlocal; i++) {
    buf[i][0] = x[i][0];
    buf[i][1] = x[i][1];
    buf[i][2] = x[i][2];
    buf[i][3] = ubuf(tag[i]).d;
    buf[i][4] = ubuf(-1).d;
    buf[i][5] = ubuf(type[i]).d;
  }
}

/*  ----------------------------------------------------------------------
   pack atom info to plt or szplt file 
   -------------------------------------------------------------------------  */

void AtomVecAtomic::pack_tecplot_binary(double *buf, int ival)
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (ival < 3 && ival >= 0) buf[i] = x[i][ival];
    else if (ival == 3) buf[i] = tag[i];
    else if (ival == 4) buf[i] = -1;
    else if (ival == 5) buf[i] = type[i];
    else error->one(FLERR, "Invalid ival in pack_tecplot_binary");
  }
}


/*  ----------------------------------------------------------------------
   add a new atom
   called by disc_elements function to discretize element to atoms
   called by create_atoms command
   -------------------------------------------------------------------------  */

void AtomVecAtomic::create_atom(double *coord, int itype, tagint itag, int gtag)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  type[nlocal] = itype;
  tag[nlocal] = itag;
  mask[nlocal] = 1|4; 
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  grain_tag[nlocal] = gtag;
  v[nlocal][0] = 0.0; 
  v[nlocal][1] = 0.0; 
  v[nlocal][2] = 0.0; 
  atom->nlocal++;
}

