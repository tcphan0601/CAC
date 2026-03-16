#include "meam.h"
#include "atom.h"
#include "element.h"

using namespace CAC_NS;

void MEAM::meam_dens_final(int eflag_either, int eflag_global, int eflag_atom, double *eng_vdwl,
                      double *eatom, double **enode, int *fmap, double **scale, int &errorflag)
{
  int i, elti;
  int m;
  double rhob, G, dG, Gbar, dGbar, gam, shp[3], Z;
  double denom, rho_bkgd, Fl;
  double scaleii;
  int inode;
  int *atype = atom->type;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *npe = element->npe;

  //     Complete the calculation of density for atoms

  for (i = 0; i < atom->nlocal; i++) {
    elti = fmap[atype[i]];
    if (elti >= 0) {
      scaleii = scale[atype[i]][atype[i]];
      atomrho1[i] = 0.0;
      atomrho2[i] = -1.0 / 3.0 * atomarho2b[i] * atomarho2b[i];
      atomrho3[i] = 0.0;
      for (m = 0; m < 3; m++) {
        atomrho1[i] += atomarho1[i][m] * atomarho1[i][m];
        atomrho3[i] -= 3.0 / 5.0 * atomarho3b[i][m] * atomarho3b[i][m];
      }
      for (m = 0; m < 6; m++) {
        atomrho2[i] += this->v2D[m] * atomarho2[i][m] * atomarho2[i][m];
      }
      for (m = 0; m < 10; m++) {
        atomrho3[i] += this->v3D[m] * atomarho3[i][m] * atomarho3[i][m];
      }

      if (atomrho0[i] > 0.0) {
        if (this->ialloy == 1) {
          atomt_ave[i][0] = fdiv_zero(atomt_ave[i][0], atomtsq_ave[i][0]);
          atomt_ave[i][1] = fdiv_zero(atomt_ave[i][1], atomtsq_ave[i][1]);
          atomt_ave[i][2] = fdiv_zero(atomt_ave[i][2], atomtsq_ave[i][2]);
        } else if (this->ialloy == 2) {
          atomt_ave[i][0] = this->t1_meam[elti];
          atomt_ave[i][1] = this->t2_meam[elti];
          atomt_ave[i][2] = this->t3_meam[elti];
        } else {
          atomt_ave[i][0] /= atomrho0[i];
          atomt_ave[i][1] /= atomrho0[i];
          atomt_ave[i][2] /= atomrho0[i];
        }
      }

      atomgamma[i] = atomt_ave[i][0] * atomrho1[i] + atomt_ave[i][1] * atomrho2[i] + atomt_ave[i][2] * atomrho3[i];

      if (atomrho0[i] > 0.0) {
        atomgamma[i] /= (atomrho0[i] * atomrho0[i]);
      }

      Z = get_Zij(this->lattce_meam[elti][elti]);

      G = G_gam(atomgamma[i], this->ibar_meam[elti], errorflag);
      if (errorflag != 0)
        return;

      get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);

      if (this->ibar_meam[elti] <= 0) {
        Gbar = 1.0;
        dGbar = 0.0;
      } else {
        if (this->mix_ref_t == 1) {
          gam = (atomt_ave[i][0] * shp[0] + atomt_ave[i][1] * shp[1] + atomt_ave[i][2] * shp[2]) / (Z * Z);
        } else {
          gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
                (Z * Z);
        }
        Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
      }
      atomrho[i] = atomrho0[i] * G;

      if (this->mix_ref_t == 1) {
        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          gam = (atomt_ave[i][0] * shp[0] + atomt_ave[i][1] * shp[1] + atomt_ave[i][2] * shp[2]) / (Z * Z);
          Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
        }
        rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
      } else {
        if (this->bkgd_dyn == 1) {
          rho_bkgd = this->rho0_meam[elti] * Z;
        } else {
          rho_bkgd = this->rho_ref_meam[elti];
        }
      }
      rhob = atomrho[i] / rho_bkgd;
      denom = 1.0 / rho_bkgd;

      G = dG_gam(atomgamma[i], this->ibar_meam[elti], dG);

      atomdgamma1[i] = (G - 2 * dG * atomgamma[i]) * denom;

      if (!iszero(atomrho0[i])) {
        atomdgamma2[i] = (dG / atomrho0[i]) * denom;
      } else {
        atomdgamma2[i] = 0.0;
      }

      //     dgamma3 is nonzero only if we are using the "mixed" rule for
      //     computing t in the reference system (which is not correct, but
      //     included for backward compatibility
      if (this->mix_ref_t == 1) {
        atomdgamma3[i] = atomrho0[i] * G * dGbar / (Gbar * Z * Z) * denom;
      } else {
        atomdgamma3[i] = 0.0;
      }

      Fl = embedding(this->A_meam[elti], this->Ec_meam[elti][elti], rhob, atomfrhop[i]);

      if (eflag_either != 0) {
        Fl *= scaleii;
        if (eflag_global != 0) {
          *eng_vdwl = *eng_vdwl + Fl;
        }
        if (eflag_atom != 0) {
          eatom[i] = eatom[i] + Fl;
        }
      }
    }
  }

  //     Complete the calculation of density for nodes
  
  for (i = 0; i < element->nlocal; i++) {
    elti = fmap[ctype[i]];
    if (elti >= 0) {
      for (inode = 0; inode < npe[etype[i]]; inode++) {
        scaleii = scale[ctype[i]][ctype[i]];
        noderho1[i][inode] = 0.0;
        noderho2[i][inode] = -1.0 / 3.0 * nodearho2b[i][inode] * nodearho2b[i][inode];
        noderho3[i][inode] = 0.0;
        for (m = 0; m < 3; m++) {
          noderho1[i][inode] += nodearho1[i][inode][m] * nodearho1[i][inode][m];
          noderho3[i][inode] -= 3.0 / 5.0 * nodearho3b[i][inode][m] * nodearho3b[i][inode][m];
        }
        for (m = 0; m < 6; m++) {
          noderho2[i][inode] += this->v2D[m] * nodearho2[i][inode][m] * nodearho2[i][inode][m];
        }
        for (m = 0; m < 10; m++) {
          noderho3[i][inode] += this->v3D[m] * nodearho3[i][inode][m] * nodearho3[i][inode][m];
        }

        if (noderho0[i][inode] > 0.0) {
          if (this->ialloy == 1) {
            nodet_ave[i][inode][0] = fdiv_zero(nodet_ave[i][inode][0], nodetsq_ave[i][inode][0]);
            nodet_ave[i][inode][1] = fdiv_zero(nodet_ave[i][inode][1], nodetsq_ave[i][inode][1]);
            nodet_ave[i][inode][2] = fdiv_zero(nodet_ave[i][inode][2], nodetsq_ave[i][inode][2]);
          } else if (this->ialloy == 2) {
            nodet_ave[i][inode][0] = this->t1_meam[elti];
            nodet_ave[i][inode][1] = this->t2_meam[elti];
            nodet_ave[i][inode][2] = this->t3_meam[elti];
          } else {
            nodet_ave[i][inode][0] /= noderho0[i][inode];
            nodet_ave[i][inode][1] /= noderho0[i][inode];
            nodet_ave[i][inode][2] /= noderho0[i][inode];
          }
        }

        nodegamma[i][inode] = nodet_ave[i][inode][0] * noderho1[i][inode] 
          + nodet_ave[i][inode][1] * noderho2[i][inode] 
          + nodet_ave[i][inode][2] * noderho3[i][inode];

        if (noderho0[i][inode] > 0.0) {
          nodegamma[i][inode] /= (noderho0[i][inode] * noderho0[i][inode]);
        }

        Z = get_Zij(this->lattce_meam[elti][elti]);

        G = G_gam(nodegamma[i][inode], this->ibar_meam[elti], errorflag);
        if (errorflag != 0)
          return;

        get_shpfcn(this->lattce_meam[elti][elti], this->stheta_meam[elti][elti], this->ctheta_meam[elti][elti], shp);

        if (this->ibar_meam[elti] <= 0) {
          Gbar = 1.0;
          dGbar = 0.0;
        } else {
          if (this->mix_ref_t == 1) {
            gam = (nodet_ave[i][inode][0] * shp[0] + nodet_ave[i][inode][1] * shp[1] + nodet_ave[i][inode][2] * shp[2]) / (Z * Z);
          } else {
            gam = (this->t1_meam[elti] * shp[0] + this->t2_meam[elti] * shp[1] + this->t3_meam[elti] * shp[2]) /
              (Z * Z);
          }
          Gbar = G_gam(gam, this->ibar_meam[elti], errorflag);
        }
        noderho[i][inode] = noderho0[i][inode] * G;

        if (this->mix_ref_t == 1) {
          if (this->ibar_meam[elti] <= 0) {
            Gbar = 1.0;
            dGbar = 0.0;
          } else {
            gam = (nodet_ave[i][inode][0] * shp[0] 
                + nodet_ave[i][inode][1] * shp[1] 
                + nodet_ave[i][inode][2] * shp[2]) / (Z * Z);
            Gbar = dG_gam(gam, this->ibar_meam[elti], dGbar);
          }
          rho_bkgd = this->rho0_meam[elti] * Z * Gbar;
        } else {
          if (this->bkgd_dyn == 1) {
            rho_bkgd = this->rho0_meam[elti] * Z;
          } else {
            rho_bkgd = this->rho_ref_meam[elti];
          }
        }
        rhob = noderho[i][inode] / rho_bkgd;
        denom = 1.0 / rho_bkgd;

        G = dG_gam(nodegamma[i][inode], this->ibar_meam[elti], dG);

        nodedgamma1[i][inode] = (G - 2 * dG * nodegamma[i][inode]) * denom;

        if (!iszero(noderho0[i][inode])) {
          nodedgamma2[i][inode] = (dG / noderho0[i][inode]) * denom;
        } else {
          nodedgamma2[i][inode] = 0.0;
        }

        //     dgamma3 is nonzero only if we are using the "mixed" rule for
        //     computing t in the reference system (which is not correct, but
        //     included for backward compatibility
        if (this->mix_ref_t == 1) {
          nodedgamma3[i][inode] = noderho0[i][inode] * G * dGbar / (Gbar * Z * Z) * denom;
        } else {
          nodedgamma3[i][inode] = 0.0;
        }

        Fl = embedding(this->A_meam[elti], this->Ec_meam[elti][elti], rhob, nodefrhop[i][inode]);

        if (eflag_either != 0) {
          Fl *= scaleii;
          if (eflag_global != 0) {
            *eng_vdwl = *eng_vdwl + Fl;
          }
          if (eflag_atom != 0) {
            enode[i][inode] = enode[i][inode] + Fl;
          }
        }
      }
    }
  }
}

