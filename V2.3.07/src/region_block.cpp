#include "stdlib.h"
#include "string.h"
#include "region_block.h"
#include "domain.h"
#include "comm.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;

#define BIG 1.0e20

/*-----------------------------------------------------------------------------------*/

RegBlock::RegBlock(CAC *cac, int narg, char **arg) : Region(cac, narg, arg)
{
   options(narg-8,&arg[8]);

   if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {  
     if (domain->box_exist == 0)
       error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
     else if (domain->triclinic == 0) xlo = domain->boxlo[0];
     else xlo = domain->boxlo_bound[0];
   } else xlo = xscale*universe->numeric(FLERR,arg[2]);

   if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
     if (domain->box_exist == 0)
       error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[3],"INF") == 0) xhi = BIG;
     else if (domain->triclinic == 0) xhi = domain->boxhi[0];
     else xhi = domain->boxhi_bound[0];
   } else xhi = xscale*universe->numeric(FLERR,arg[3]);

   if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
     if (domain->box_exist == 0)
       error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
     else if (domain->triclinic == 0) ylo = domain->boxlo[1];
     else ylo = domain->boxlo_bound[1];
   } else ylo = yscale*universe->numeric(FLERR,arg[4]);

   if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
     if (domain->box_exist == 0)
       error->all(FLERR, "Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[5],"INF") == 0) yhi = BIG;
     else if (domain->triclinic == 0) yhi = domain->boxhi[1];
     else yhi = domain->boxhi_bound[1];
   } else yhi = yscale*universe->numeric(FLERR,arg[5]);

   if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
     if (domain->box_exist == 0)
       error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
     else if (domain->triclinic== 0) zlo = domain->boxlo[2];
     else zlo = domain->boxlo_bound[2];
   } else zlo = zscale*universe->numeric(FLERR,arg[6]);

   if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
     if (domain->box_exist == 0)
       error->all(FLERR, "Cannot use region INF or EDGE when box does not exist");
     if (strcmp(arg[7],"INF") == 0) zhi = BIG;
     else if (domain->triclinic == 0) zhi = domain->boxhi[2];
     else zhi = domain->boxhi_bound[2];
   } else zhi = zscale*universe->numeric(FLERR,arg[7]);

   // error check

   if (xlo > xhi || ylo > yhi || zlo > zhi)
     error->all(FLERR, "Illegal region block command");

   // extent of block

   if (interior) {
     bboxflag = 1;
     extent_xlo = xlo;
     extent_xhi = xhi;
     extent_ylo = ylo;
     extent_yhi = yhi;
     extent_zlo = zlo;
     extent_zhi = zhi;
   } else bboxflag = 0;

   // particle could be close to all 6 planes

   cmax = 6;
   contact = new Contact[cmax];

}

/*-------------------------------------------------------------*/

RegBlock::~RegBlock()
{
  delete [] contact;
}

/*-----------------------------------------------------------------------*/

int RegBlock::inside(double x, double y, double z)
{
  if (x >= xlo && x <= xhi && 
      y >= ylo && y <= yhi && 
      z >= zlo && z <= zhi)
    return 1;
  return 0;
}
