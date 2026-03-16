#include <mpi.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair.h"
#include "atom.h"
#include "element.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "math_const.h"

using namespace CAC_NS;
using namespace MathConst;

//enum{NONE,RLINEAR,RSQ,BMP};

// allocate space for static class instance variable and initialize it
int Pair::instance_total = 0;

/* ---------------------------------------------------------------------- */

Pair::Pair(CAC *cac) : Pointers(cac)
{
        instance_me = instance_total++;

        eng_vdwl = eng_coul = 0.0;

        comm_atom_forward = comm_elem_forward = comm_atom_reverse = comm_elem_reverse = comm_atom_reverse_off = comm_elem_reverse_off = 0;

        single_enable = 1;
        restartinfo = 1;
        respa_enable = 0;
        one_coeff = 0;
        no_virial_fdotr_compute = 0;
        writedata = 0;
        ghostneigh = 0;

        nextra = 0;
        pvector = NULL;
        single_extra = 0;
        svector = NULL;

        ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = dipoleflag = 0;
        reinitflag = 1;

        // pair_modify settingsx

        compute_flag = 1;
        manybody_flag = 0;
        threebody_flag = 0;
        offset_flag = 0;
        mix_flag = GEOMETRIC;
        tail_flag = 0;
        etail = ptail = etail_ij = ptail_ij = 0.0;
        ncoultablebits = 12;
        ndisptablebits = 12;
        tabinner = sqrt(2.0);
        tabinner_disp = sqrt(2.0);

        allocated = 0;
//  suffix_flag = Suffix::NONE;

        maxeatom = maxvatom = 0;
        maxeelem = maxvelem = 0;
        eatom = NULL;
        vatom = NULL;
        enode = NULL;
        vnode = NULL;

        datamask = ALL_MASK;
        datamask_ext = ALL_MASK;

//  execution_space = Host;
//  datamask_read = ALL_MASK;
//  datamask_modify = ALL_MASK;

        copymode = 0;
        dvalue = 0;
        elemscale = NULL;
}

Pair::~Pair()
{
        if (copymode) return;
        memory->destroy(eatom);
        memory->destroy(vatom);
        memory->destroy(enode);
        memory->destroy(vnode);
        delete [] elemscale;

}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
        int i,j;

        if (offset_flag && tail_flag)
                error->all(FLERR,"Cannot have both pair_modify shift and tail set to yes");
        if (tail_flag && domain->dimension == 2)
                error->all(FLERR,"Cannot use pair tail corrections with 2d simulations");
        if (tail_flag && domain->nonperiodic && comm->me == 0)
                error->warning(FLERR,"Using pair tail corrections with nonperiodic system");
        if (!compute_flag && tail_flag)
                error->warning(FLERR,"Using pair tail corrections with compute set to no");
        if (!compute_flag && offset_flag)
                error->warning(FLERR,"Using pair potential shift with compute set to no");

        // for manybody potentials
        // check if bonded exclusions could invalidate the neighbor list

        // I,I coeffs must be set\
        // init_one() will check if I,J is set explicitly or inferred by mixing

        if (!allocated) error->all(FLERR,"All pair coeffs are not set");

        for (i = 1; i <= atom->ntypes; i++)
                if (setflag[i][i] == 0) error->all(FLERR,"All pair coeffs are not set");

        // style-specific initialization

        init_style();

        // call init_one() for each I,J
        // set cutsq for each I,J, used to neighbor
        // cutforce = max of all I,J cutoffs

        cutforce = 0.0;
        etail = ptail = 0.0;
        double cut;

        for (i = 1; i <= atom->ntypes; i++)
                for (j = i; j <= atom->ntypes; j++) {
                        cut = init_one(i,j);
                        cutsq[i][j] = cutsq[j][i] = cut*cut;
                        cutforce = MAX(cutforce,cut);
                        if (tail_flag) {
                                etail += etail_ij;
                                ptail += ptail_ij;
                                if (i != j) {
                                        etail += etail_ij;
                                        ptail += ptail_ij;
                                }
                        }
                }

        delete [] elemscale;
        int netypes = element->netypes;
        int npe = element->npe;
        int *nintpl = element->nintpl;
        elemscale = new double[netypes+1];
        for (i = 1; i <= netypes; i++)
                elemscale[i] = (static_cast<double>(npe))/nintpl[i];

}

/* ----------------------------------------------------------------------
   init specific to a pair style
   specific pair style can override this function
     if needs its own error checks
     if needs another kind of neighbor list
   request default neighbor list = half list
   ------------------------------------------------------------------------- */

void Pair::init_style()
{
        neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
    mixing of pair potential prefactors (epsilon)
   ------------------------------------------------------------------------- */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
        if (mix_flag == GEOMETRIC)
                return sqrt(eps1*eps2);
        else if (mix_flag == ARITHMETIC)
                return sqrt(eps1*eps2);
        else if (mix_flag == SIXTHPOWER)
                return (2.0 * sqrt(eps1*eps2) *
                        pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) + pow(sig2,6.0)));
        else return 0.0;
}

/* ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
   ------------------------------------------------------------------------- */

double Pair::mix_distance(double sig1, double sig2)
{
        if (mix_flag == GEOMETRIC)
                return sqrt(sig1*sig2);
        else if (mix_flag == ARITHMETIC)
                return (0.5 * (sig1+sig2));
        else if (mix_flag == SIXTHPOWER)
                return pow((0.5 * (pow(sig1,6.0) + pow(sig2,6.0))),1.0/6.0);
        else return 0.0;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   specific pair style can override this function
   ------------------------------------------------------------------------- */

void Pair::init_list(int which, NeighList *ptr)
{
        list = ptr;
}

/* ---------------------------------------------------------------------- */

void Pair::compute_dummy(int eflag, int vflag)
{
        evflag = 0;
}

/* ---------------------------------------------------------------------- */

double Pair::memory_usage()
{
        int npe = element->npe;
        int n = comm->nthreads;
        double bytes = n * maxeatom * sizeof(double);
        bytes += n * maxvatom * 6 * sizeof(double);
        bytes += n * maxeelem * npe * sizeof(double);
        bytes += n * maxvelem * npe * 6 * sizeof(double);
        return bytes;
}
/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
   ------------------------------------------------------------------------- */

void Pair::ev_setup(int eflag, int vflag, int alloc)
{
        int npe = element->npe;
        int i,j,n;

        evflag = 1;

        eflag_either = eflag;
        eflag_global = eflag % 2;
        eflag_atom = eflag / 2;

        vflag_either = vflag;
        vflag_global = vflag % 4;
        vflag_atom = vflag / 4;

        // reallocate per-atom/per-node arrays if necessary

        if (eflag_atom) {
                if (atom->nmax > maxeatom) {
                        maxeatom = atom->nmax;
                        if (alloc) {
                                memory->destroy(eatom);
                                memory->create(eatom,comm->nthreads*maxeatom,"pair:eatom");
                        }
                }
                if (element->nmax > maxeelem) {
                        maxeelem = element->nmax;
                        if (alloc) {
                                memory->destroy(enode);
                                memory->create(enode,comm->nthreads*maxeelem,npe,"pair:enode");
                        }
                }
        }

        if (vflag_atom) {
                if (atom->nmax > maxvatom) {
                        maxvatom = atom->nmax;
                        if (alloc) {
                                memory->destroy(vatom);
                                memory->create(vatom,comm->nthreads*maxvatom,6,"pair:vatom");
                        }
                }
                if (element->nmax > maxvelem) {
                        maxvelem = element->nmax;
                        if (alloc) {
                                memory->destroy(vnode);
                                memory->create(vnode,comm->nthreads*maxvelem,npe,6,"pair:vnode");
                        }
                }
        }

        // zero accumulators
        // use force->newton instead of newton_pair
        //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info

        if (eflag_global) eng_vdwl = eng_coul = 0.0;
        if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
        if (eflag_atom && alloc) {
                n = atom->nlocal;
                if (force->newton) n += atom->nghost;
                for (i = 0; i < n; i++) eatom[i] = 0.0;
                n = element->nlocal;
                if (force->newton) n += element->nghost;
                for (i = 0; i < n; i++)
                        for (j = 0; j < npe; j++) enode[i][j] = 0.0;
        }
        if (vflag_atom && alloc) {
                n = atom->nlocal;
                if (force->newton) n += atom->nghost;
                for (i = 0; i < n; i++) {
                        vatom[i][0] = 0.0;
                        vatom[i][1] = 0.0;
                        vatom[i][2] = 0.0;
                        vatom[i][3] = 0.0;
                        vatom[i][4] = 0.0;
                        vatom[i][5] = 0.0;
                }
                n = element->nlocal;
                if (force->newton) n += element->nghost;
                for (i = 0; i < n; i++)
                        for (j = 0; j < npe; j++) {
                                vnode[i][j][0] = 0.0;
                                vnode[i][j][1] = 0.0;
                                vnode[i][j][2] = 0.0;
                                vnode[i][j][3] = 0.0;
                                vnode[i][j][4] = 0.0;
                                vnode[i][j][5] = 0.0;
                        }
        }

        // if vflag_global = 2 and pair::compute() calls virial_fdotr_compute() and no element in simulation
        // compute global virial via (F dot r) instead of via pairwise summation
        // unset other flags as appropriate

        if (vflag_global == 2 && no_virial_fdotr_compute == 0 && element->nelements == 0) {
                vflag_fdotr = 1;
                vflag_global = 0;
                if (vflag_atom == 0) vflag_either = 0;
                if (vflag_either == 0 && eflag_either == 0) evflag = 0;
        } else vflag_fdotr = 0;
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators for atom-intpl pairs
   need i < nlocal test since called by bond_quartic and dihedral_charmm
   ------------------------------------------------------------------------- */

void Pair::atom_ev_tally(int i, int j, int nlocal, int newton_pair,
                         double evdwl, double ecoul, double fpair,
                         double delx, double dely, double delz)
{
        double evdwlhalf,ecoulhalf,epairhalf,v[6];

        if (eflag_either) {
                if (eflag_global) {
                        if (newton_pair) {
                                eng_vdwl += evdwl;
                                eng_coul += ecoul;
                        } else {
                                evdwlhalf = 0.5*evdwl;
                                ecoulhalf = 0.5*ecoul;
                                if (i < nlocal) {
                                        eng_vdwl += evdwlhalf;
                                        eng_coul += ecoulhalf;
                                }
                                if (j < nlocal) {
                                        eng_vdwl += evdwlhalf;
                                        eng_coul += ecoulhalf;
                                }
                        }
                }
                if (eflag_atom) {
                        epairhalf = 0.5 * (evdwl + ecoul);
                        if (newton_pair || i < nlocal) eatom[i] += epairhalf;
                        if (newton_pair || j < nlocal) eatom[j] += epairhalf;
                }
        }

        if (vflag_either) {
                v[0] = delx*delx*fpair;
                v[1] = dely*dely*fpair;
                v[2] = delz*delz*fpair;
                v[3] = delx*dely*fpair;
                v[4] = delx*delz*fpair;
                v[5] = dely*delz*fpair;

                if (vflag_global) {
                        if (newton_pair) {
                                virial[0] += v[0]; virial[1] += v[1];
                                virial[2] += v[2]; virial[3] += v[3];
                                virial[4] += v[4]; virial[5] += v[5];
                        } else {
                                if (i < nlocal) {
                                        virial[0] += 0.5*v[0]; virial[1] += 0.5*v[1];
                                        virial[2] += 0.5*v[2]; virial[3] += 0.5*v[3];
                                        virial[4] += 0.5*v[4]; virial[5] += 0.5*v[5];
                                }
                                if (j < nlocal) {
                                        virial[0] += 0.5*v[0]; virial[1] += 0.5*v[1];
                                        virial[2] += 0.5*v[2]; virial[3] += 0.5*v[3];
                                        virial[4] += 0.5*v[4]; virial[5] += 0.5*v[5];
                                }
                        }
                }

                if (vflag_atom) {
                        if (newton_pair || i < nlocal) {
                                vatom[i][0] += 0.5*v[0]; vatom[i][1] += 0.5*v[1];
                                vatom[i][2] += 0.5*v[2]; vatom[i][3] += 0.5*v[3];
                                vatom[i][4] += 0.5*v[4]; vatom[i][5] += 0.5*v[5];
                        }
                        if (newton_pair || j < nlocal) {
                                vatom[j][0] += 0.5*v[0]; vatom[j][1] += 0.5*v[1];
                                vatom[j][2] += 0.5*v[2]; vatom[j][3] += 0.5*v[3];
                                vatom[j][4] += 0.5*v[4]; vatom[j][5] += 0.5*v[5];
                        }
                }
        }

//  if (num_tally_compute > 0) {
//    for (int k = 0; k < num_tally_compute; ++k) {
//      Compute *c = list_tally_compute[k];
//      c->pair_tally_callback(i, j, nlocal, newton_pair,
//                             evdwl, ecoul, fpair, delx, dely, delz);
//    }
//  }
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by EAM potential
   ------------------------------------------------------------------------- */

void Pair::atom_ev_tally_full(int i, double evdwl, double ecoul, double fpair, double delx, double dely, double delz)
{
        double v[6];
        if (eflag_either) {
                if (eflag_global) {
                        eng_vdwl += 0.5*evdwl;
                        eng_coul += 0.5*ecoul;
                }
                if (eflag_atom)
                        eatom[i] += 0.5*(evdwl + ecoul);
        }

        if (vflag_either) {
                v[0] = 0.5*delx*delx*fpair;
                v[1] = 0.5*dely*dely*fpair;
                v[2] = 0.5*delz*delz*fpair;
                v[3] = 0.5*delx*dely*fpair;
                v[4] = 0.5*delx*delz*fpair;
                v[5] = 0.5*dely*delz*fpair;

                if (vflag_global) {
                        virial[0] += v[0]; virial[1] += v[1];
                        virial[2] += v[2]; virial[3] += v[3];
                        virial[4] += v[4]; virial[5] += v[5];
                }
                if (vflag_atom) {
                        vatom[i][0] += v[0]; vatom[i][1] += v[1];
                        vatom[i][2] += v[2]; vatom[i][3] += v[3];
                        vatom[i][4] += v[4]; vatom[i][5] += v[5];
                }
        }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-node accumulators
   skip if inode < 0
   ------------------------------------------------------------------------- */

void Pair::node_ev_tally(int i, int inode, double evdwl, double ecoul, double fpair, double delx, double dely, double delz)
{
        //if (inode < 0) return;
        double v[6];
        double *tempscale = force->pair->elemscale;
        double iscale = tempscale[element->etype[i]];
        if (eflag_either) {
                if (eflag_global) {
                        eng_vdwl += 0.5*evdwl/iscale;
                        eng_coul += 0.5*ecoul/iscale;
                }
                if (eflag_atom)
                        enode[i][inode] += 0.5*(evdwl + ecoul);
        }

        if (vflag_either) {
                v[0] = 0.5*delx*delx*fpair;
                v[1] = 0.5*dely*dely*fpair;
                v[2] = 0.5*delz*delz*fpair;
                v[3] = 0.5*delx*dely*fpair;
                v[4] = 0.5*delx*delz*fpair;
                v[5] = 0.5*dely*delz*fpair;
                if (vflag_global) {
                        virial[0] += v[0]/iscale; virial[1] += v[1]/iscale;
                        virial[2] += v[2]/iscale; virial[3] += v[3]/iscale;
                        virial[4] += v[4]/iscale; virial[5] += v[5]/iscale;
                }
                if (vflag_atom) {
                        vnode[i][inode][0] += v[0]; vnode[i][inode][1] += v[1];
                        vnode[i][inode][2] += v[2]; vnode[i][inode][3] += v[3];
                        vnode[i][inode][4] += v[4]; vnode[i][inode][5] += v[5];
                }
        }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW potentials when triple involves 3 atoms, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
   ------------------------------------------------------------------------- */

void Pair::atom_ev_tally3(int i, int j, int k, double evdwl, double ecoul,
                          double *fj, double *fk, double delxji, double delyji, double delzji, double delxki, double delyki, double delzki)
{
        double evdwlthird,ecoulthird,epairthird,v[6];

        if (eflag_either) {
                if (eflag_global) {
                        eng_vdwl += evdwl;
                        eng_coul += ecoul;
                }
                if (eflag_atom) {
                        epairthird = THIRD * (evdwl + ecoul);
                        eatom[i] += epairthird;
                        eatom[j] += epairthird;
                        eatom[k] += epairthird;
                }
        }

        if (vflag_either) {
                v[0] = delxji*fj[0] + delxki*fk[0];
                v[1] = delyji*fj[1] + delyki*fk[1];
                v[2] = delzji*fj[2] + delzki*fk[2];
                v[3] = delxji*fj[1] + delxki*fk[1];
                v[4] = delxji*fj[2] + delxki*fk[2];
                v[5] = delyji*fj[2] + delyki*fk[2];

                if (vflag_global) {
                        virial[0] += v[0]; virial[1] += v[1];
                        virial[2] += v[2]; virial[3] += v[3];
                        virial[4] += v[4]; virial[5] += v[5];
                }

                if (vflag_atom) {
                        vatom[i][0] += THIRD*v[0]; vatom[i][1] += THIRD*v[1];
                        vatom[i][2] += THIRD*v[2]; vatom[i][3] += THIRD*v[3];
                        vatom[i][4] += THIRD*v[4]; vatom[i][5] += THIRD*v[5];
                        vatom[j][0] += THIRD*v[0]; vatom[j][1] += THIRD*v[1];
                        vatom[j][2] += THIRD*v[2]; vatom[j][3] += THIRD*v[3];
                        vatom[j][4] += THIRD*v[4]; vatom[j][5] += THIRD*v[5];
                        vatom[k][0] += THIRD*v[0]; vatom[k][1] += THIRD*v[1];
                        vatom[k][2] += THIRD*v[2]; vatom[k][3] += THIRD*v[3];
                        vatom[k][4] += THIRD*v[4]; vatom[k][5] += THIRD*v[5];
                }
        }
}


/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW potentials when triple involves 2 atoms, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
   ------------------------------------------------------------------------- */

void Pair::atom_ev_tally3(int i, int j, double evdwl, double ecoul,
                          double *fj, double *fk, double delxji, double delyji, double delzji, double delxki, double delyki, double delzki)
{
        double evdwlthird,ecoulthird,epairthird,v[6];

        if (eflag_either) {
                evdwlthird = THIRD*evdwl;
                ecoulthird = THIRD*ecoul;
                if (eflag_global) {
                        eng_vdwl += 2.0*evdwlthird;
                        eng_coul += 2.0*ecoulthird;
                }
                if (eflag_atom) {
                        epairthird = evdwlthird + ecoulthird;
                        eatom[i] += epairthird;
                        eatom[j] += epairthird;
                }
        }

        if (vflag_either) {
                v[0] = THIRD*(delxji*fj[0] + delxki*fk[0]);
                v[1] = THIRD*(delyji*fj[1] + delyki*fk[1]);
                v[2] = THIRD*(delzji*fj[2] + delzki*fk[2]);
                v[3] = THIRD*(delxji*fj[1] + delxki*fk[1]);
                v[4] = THIRD*(delxji*fj[2] + delxki*fk[2]);
                v[5] = THIRD*(delyji*fj[2] + delyki*fk[2]);

                if (vflag_global) {
                        virial[0] += 2.0*v[0]; virial[1] += 2.0*v[1];
                        virial[2] += 2.0*v[2]; virial[3] += 2.0*v[3];
                        virial[4] += 2.0*v[4]; virial[5] += 2.0*v[5];
                }

                if (vflag_atom) {
                        vatom[i][0] += v[0]; vatom[i][1] += v[1];
                        vatom[i][2] += v[2]; vatom[i][3] += v[3];
                        vatom[i][4] += v[4]; vatom[i][5] += v[5];
                        vatom[j][0] += v[0]; vatom[j][1] += v[1];
                        vatom[j][2] += v[2]; vatom[j][3] += v[3];
                        vatom[j][4] += v[4]; vatom[j][5] += v[5];

                }
        }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW potentials when triple involves 1 atom, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
   ------------------------------------------------------------------------- */

void Pair::atom_ev_tally3(int i, double evdwl, double ecoul,
                          double *fj, double *fk, double delxji, double delyji, double delzji, double delxki, double delyki, double delzki)
{
        double evdwlthird,ecoulthird,epairthird,v[6];

        if (eflag_either) {
                evdwlthird = THIRD*evdwl;
                ecoulthird = THIRD*ecoul;
                if (eflag_global) {
                        eng_vdwl += evdwlthird;
                        eng_coul += ecoulthird;
                }
                if (eflag_atom)
                        eatom[i] += evdwlthird + ecoulthird;
        }

        if (vflag_either) {
                v[0] = THIRD*(delxji*fj[0] + delxki*fk[0]);
                v[1] = THIRD*(delyji*fj[1] + delyki*fk[1]);
                v[2] = THIRD*(delzji*fj[2] + delzki*fk[2]);
                v[3] = THIRD*(delxji*fj[1] + delxki*fk[1]);
                v[4] = THIRD*(delxji*fj[2] + delxki*fk[2]);
                v[5] = THIRD*(delyji*fj[2] + delyki*fk[2]);

                if (vflag_global) {
                        virial[0] += v[0]; virial[1] += v[1];
                        virial[2] += v[2]; virial[3] += v[3];
                        virial[4] += v[4]; virial[5] += v[5];
                }

                if (vflag_atom) {
                        vatom[i][0] += v[0]; vatom[i][1] += v[1];
                        vatom[i][2] += v[2]; vatom[i][3] += v[3];
                        vatom[i][4] += v[4]; vatom[i][5] += v[5];
                }
        }
}


/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-node accumulators
   check if inode,jnode, or knode >= 0 (= -1 if not a node)
   called by SW potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
   ------------------------------------------------------------------------- */

void Pair::node_ev_tally3(int i, int inode, int j, int jnode,
                          int k, int knode, double evdwl, double ecoul, double *fj,
                          double *fk, double delxji, double delyji, double delzji, double delxki, double delyki, double delzki)
{
        double evdwlthird,ecoulthird,epairthird,v[6];

        double *tempscale = force->pair->elemscale;
        double iscale = tempscale[element->etype[i]];
        double jscale = tempscale[element->etype[j]];
        double kscale = tempscale[element->etype[k]];
        if (eflag_either) {
                evdwlthird = THIRD * evdwl;
                ecoulthird = THIRD * ecoul;
                if (eflag_global) {
                        if (inode >= 0) {
                                eng_vdwl += evdwlthird/iscale;
                                eng_coul += ecoulthird/iscale;
                        }
                        if (jnode >= 0) {
                                eng_vdwl += evdwlthird/jscale;
                                eng_coul += ecoulthird/jscale;
                        }
                        if (knode >= 0) {
                                eng_vdwl += evdwlthird/kscale;
                                eng_coul += ecoulthird/kscale;
                        }
                }
                if (eflag_atom) {
                        epairthird = evdwlthird + ecoulthird;
                        if (inode >= 0) enode[i][inode] += epairthird;
                        if (jnode >= 0) enode[j][jnode] += epairthird;
                        if (knode >= 0) enode[k][knode] += epairthird;
                }
        }

        if (vflag_either) {
                v[0] = THIRD*(delxji*fj[0] + delxki*fk[0]);
                v[1] = THIRD*(delyji*fj[1] + delyki*fk[1]);
                v[2] = THIRD*(delzji*fj[2] + delzki*fk[2]);
                v[3] = THIRD*(delxji*fj[1] + delxki*fk[1]);
                v[4] = THIRD*(delxji*fj[2] + delxki*fk[2]);
                v[5] = THIRD*(delyji*fj[2] + delyki*fk[2]);

                if (vflag_global) {
                        if (inode >= 0) {
                                virial[0] += v[0]/iscale; virial[1] += v[1]/iscale;
                                virial[2] += v[2]/iscale; virial[3] += v[3]/iscale;
                                virial[4] += v[4]/iscale; virial[5] += v[5]/iscale;
                        }
                        if (jnode >= 0) {
                                virial[0] += v[0]/jscale; virial[1] += v[1]/jscale;
                                virial[2] += v[2]/jscale; virial[3] += v[3]/jscale;
                                virial[4] += v[4]/jscale; virial[5] += v[5]/jscale;
                        }
                        if (knode >= 0) {
                                virial[0] += v[0]/kscale; virial[1] += v[1]/kscale;
                                virial[2] += v[2]/kscale; virial[3] += v[3]/kscale;
                                virial[4] += v[4]/kscale; virial[5] += v[5]/kscale;
                        }
                }

                if (vflag_atom) {
                        if (inode >= 0) {
                                vnode[i][inode][0] += v[0]; vnode[i][inode][1] += v[1];
                                vnode[i][inode][2] += v[2]; vnode[i][inode][3] += v[3];
                                vnode[i][inode][4] += v[4]; vnode[i][inode][5] += v[5];
                        }
                        if (jnode >= 0) {
                                vnode[j][jnode][0] += v[0]; vnode[j][jnode][1] += v[1];
                                vnode[j][jnode][2] += v[2]; vnode[j][jnode][3] += v[3];
                                vnode[j][jnode][4] += v[4]; vnode[j][jnode][5] += v[5];
                        }
                        if (knode >= 0) {
                                vnode[k][knode][0] += v[0]; vnode[k][knode][1] += v[1];
                                vnode[k][knode][2] += v[2]; vnode[k][knode][3] += v[3];
                                vnode[k][knode][4] += v[4]; vnode[k][knode][5] += v[5];
                        }
                }
        }
}

/* ----------------------------------------------------------------------
   compute global pair virial via summing F dot r over own & ghost atoms
   at this point, only pairwise forces have been accumulated in atom->f
   only called when no element in simulation
   ------------------------------------------------------------------------- */

void Pair::virial_fdotr_compute()
{
        double **x = atom->x;
        double **f = atom->f;

        // sum over force on all particles including ghosts

        int nall = atom->nlocal + atom->nghost;
        for (int i = 0; i < nall; i++) {
                virial[0] += f[i][0]*x[i][0];
                virial[1] += f[i][1]*x[i][1];
                virial[2] += f[i][2]*x[i][2];
                virial[3] += f[i][1]*x[i][0];
                virial[4] += f[i][2]*x[i][0];
                virial[5] += f[i][2]*x[i][1];
        }

        // prevent multiple calls to update the virial
        // when a hybrid pair style uses both a gpu and non-gpu pair style
        // or when respa is used with gpu pair styles

        vflag_fdotr = 0;
}
