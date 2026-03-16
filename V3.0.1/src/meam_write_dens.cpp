#include "meam.h"
#include <cstddef>
#include "memory.h"
#include "error.h"
#include "element.h"
#include "element_vec.h"
#include "atom.h"
#include "domain.h"
#include "cac.h"


using namespace CAC_NS;

/*  ------------------------------------------------------------------  */
void MEAM::write_dens()
{
  FILE *fp = fopen("dump_densities.atom", "w");
  if (fp == NULL) {
    char str[128];
    sprintf(str, "Cannot open data file dump_densities.atom");
    error->one(FLERR, str);
  }
  class Domain *domain = cac->domain;
  fprintf(fp, "ITEM: TIMESTEP\n");

  fprintf(fp, "%d\n", 0);

  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, BIGINT_FORMAT "\n", atom->natoms + element->nelements * element->ngcell[1] * element->apc[1]);
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
  fprintf(fp, "ITEM: ATOMS id type x y z tag igcell inode rho0 rho1 rho2 rho3 frhop dgamma1 dgamma2 dgamma3 arho1[1] arho1[2] arho1[3] arho2[1] arho2[2] arho2[3] arho2b arho2[4] arho2[5] arho2[6] arho3[1] arho3[2] arho3[3] arho3[4] arho3[5] arho3[6] arho3[7] arho3[8] arho3[9] arho3[10] arho3b[1] arho3b[2] arho3b[3] t_ave[1] t_ave[2] t_ave[3] tsq_ave[1] tsq_ave[2] tsq_ave[3]\n");
  double densities[37];
  for (int i = 0; i < atom->nlocal; i++) {
    fprintf(fp, "%d 3 %g %g %g %d %d %d "
        , atom->tag[i]
        , atom->x[i][0]
        , atom->x[i][1]
        , atom->x[i][2]
        , -1, -1, -1);
    copy_densities(densities, i, -1, -1, -1);
    for (int k = 0; k < 37; k++)
      fprintf(fp, "%g ", densities[k]);
    fprintf(fp, "\n");
  }
  int count = atom->natoms+1;
  double xtmp, ytmp, ztmp;
  double ****nodex = element->nodex;
  for (int i = 0; i < element->nlocal; i++) 
    for (int ibasis = 0; ibasis < element->apc[element->etype[i]]; ibasis++) 
      for (int igcell = 0; igcell < element->ngcell[element->etype[i]]; igcell++) {
        int inode = element->g2n[element->etype[i]][igcell];
        int iucell = element->g2u[element->etype[i]][igcell];

        xtmp = element->evec->interpolate(nodex,i,ibasis,iucell,0);
        ytmp = element->evec->interpolate(nodex,i,ibasis,iucell,1);
        ztmp = element->evec->interpolate(nodex,i,ibasis,iucell,2);

        int tmp;
        if (inode>=0) tmp = 1;
        else tmp = 2;

        fprintf(fp, "%d %d %g %g %g %d %d %d ",
            count++,
            tmp,
            xtmp,
            ytmp,
            ztmp,
            element->tag[i], igcell, inode);
        copy_densities(densities, i, inode, iucell, ibasis);
        for (int k = 0; k < 37; k++)
          fprintf(fp, "%g ", densities[k]);
        fprintf(fp, "\n");
      }

  fclose(fp);
}
