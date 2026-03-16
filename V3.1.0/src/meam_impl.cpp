#include "meam.h"
#include <cstddef>
#include "memory.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

MEAM::MEAM(Memory * mem, Atom * atm, Element * elem, Error * err)
  : memory(mem), atom(atm), element(elem), error(err)
{
  phir = phirar = phirar1 = phirar2 = phirar3 = phirar4 = phirar5 = phirar6 = nullptr;

  debug = 0;
  namax = 0;
  nemax = 0;
  atomrho = atomrho0 = atomrho1 = atomrho2 = atomrho3 = atomfrhop = nullptr;
  atomgamma = atomdgamma1 = atomdgamma2 = atomdgamma3 = atomarho2b = nullptr;
  atomarho1 = atomarho2 = atomarho3 = atomarho3b = atomt_ave = atomtsq_ave = nullptr;
  noderho = noderho0 = noderho1 = noderho2 = noderho3 = nodefrhop = nullptr;
  nodegamma = nodedgamma1 = nodedgamma2 = nodedgamma3 = nodearho2b = nullptr;
  nodearho1 = nodearho2 = nodearho3 = nodearho3b = nodet_ave = nodetsq_ave = nullptr;


  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = nullptr;

  neltypes = 0;
  for (int i = 0; i < maxelt; i++) {
    A_meam[i] = rho0_meam[i] = beta0_meam[i] =
      beta1_meam[i]= beta2_meam[i] = beta3_meam[i] =
      t0_meam[i] = t1_meam[i] = t2_meam[i] = t3_meam[i] =
      rho_ref_meam[i] = ibar_meam[i] = ielt_meam[i] = 0.0;
    for (int j = 0; j < maxelt; j++) {
      lattce_meam[i][j] = FCC;
      Ec_meam[i][j] = re_meam[i][j] = alpha_meam[i][j] = delta_meam[i][j] = ebound_meam[i][j] = attrac_meam[i][j] = repuls_meam[i][j] = 0.0;
      nn2_meam[i][j] = zbl_meam[i][j] = eltind[i][j] = 0;
    }
  }
  // DEBUG
  /* 
  fp_pair = fopen("pair_cac.txt", "w");
  fp_three_body = fopen("threebody_cac.txt", "w");
  fp_three_body_k = fopen("threebody_k_cac.txt", "w");
  */ 
  //----------------------------------

}

MEAM::~MEAM()
{
  // DEBUG
  /* 
  fclose(fp_pair);
  fclose(fp_three_body);
  fclose(fp_three_body_k);
   */
  //----------------------------------
  memory->destroy(this->phirar6);
  memory->destroy(this->phirar5);
  memory->destroy(this->phirar4);
  memory->destroy(this->phirar3);
  memory->destroy(this->phirar2);
  memory->destroy(this->phirar1);
  memory->destroy(this->phirar);
  memory->destroy(this->phir);

  memory->destroy(this->atomrho);
  memory->destroy(this->atomrho0);
  memory->destroy(this->atomrho1);
  memory->destroy(this->atomrho2);
  memory->destroy(this->atomrho3);
  memory->destroy(this->atomfrhop);
  memory->destroy(this->atomgamma);
  memory->destroy(this->atomdgamma1);
  memory->destroy(this->atomdgamma2);
  memory->destroy(this->atomdgamma3);
  memory->destroy(this->atomarho2b);

  memory->destroy(this->atomarho1);
  memory->destroy(this->atomarho2);
  memory->destroy(this->atomarho3);
  memory->destroy(this->atomarho3b);
  memory->destroy(this->atomt_ave);
  memory->destroy(this->atomtsq_ave);

  memory->destroy(this->noderho);
  memory->destroy(this->noderho0);
  memory->destroy(this->noderho1);
  memory->destroy(this->noderho2);
  memory->destroy(this->noderho3);
  memory->destroy(this->nodefrhop);
  memory->destroy(this->nodegamma);
  memory->destroy(this->nodedgamma1);
  memory->destroy(this->nodedgamma2);
  memory->destroy(this->nodedgamma3);
  memory->destroy(this->nodearho2b);

  memory->destroy(this->nodearho1);
  memory->destroy(this->nodearho2);
  memory->destroy(this->nodearho3);
  memory->destroy(this->nodearho3b);
  memory->destroy(this->nodet_ave);
  memory->destroy(this->nodetsq_ave);

  memory->destroy(this->scrfcn);
  memory->destroy(this->dscrfcn);
  memory->destroy(this->fcpair);
}
