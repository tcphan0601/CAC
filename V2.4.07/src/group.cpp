#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "group.h"
#include "domain.h"
#include "region.h"
#include "atom.h"
#include "element.h"
#include "element_vec.h"
#include "comm.h"
#include "force.h"
#include "input.h"
#include "dump.h"
#include "modify.h"
#include "output.h"
#include "compute.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "math_extra.h"

#include <map>

using namespace CAC_NS;
using namespace MathExtra;

#define MAX_GROUP 32
#define BIG 1e30
#define EPSILON 1.0e-6

enum{ETYPE,CTYPE,ID};
enum{CENTER,ALLNODE,ONENODE,ONEATOM};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN,NOTBETWEEN};

Group::Group(CAC *cac) : Pointers(cac)
{
  MPI_Comm_rank(world,&me);

  names = new char*[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++) {
    names[i] = NULL;
    bitmask[i] = 1 << i;
    inversemask[i] = bitmask[i] ^ ~0;
    dynamic[i] = 0;
  }

  int n;

  // create "all" group
  
  char *str1 = (char *) "all";
  n = strlen(str1) + 1;
  names[0] = new char[n];
  strcpy(names[0],str1);

  // create "element" group containing all elements and their nodes
  
  char *str2 = (char *) "element";
  n = strlen(str2) + 1;
  names[1] = new char[n];
  strcpy(names[1],str2);

  // create "atom" group containing all atoms
  
  char *str3 = (char *) "atom";
  n = strlen(str3) + 1;
  names[2] = new char[n];
  strcpy(names[2],str3);

  ngroup = 3;
}

/* ----------------------------------------------------------------------
   free all memory
   ------------------------------------------------------------------------- */

Group::~Group()
{
  for (int i = 0; i < MAX_GROUP; i++) delete [] names[i];
  delete [] names;
  delete [] bitmask;
  delete [] inversemask;
  delete [] dynamic;
}

/* ----------------------------------------------------------------------
   return group index if name matches existing group, -1 if no such group
   ------------------------------------------------------------------------- */

int Group::find(const char *name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && strcmp(name,names[igroup]) == 0) 
      return igroup;
  return -1;
}
/* ----------------------------------------------------------------------
   find group with name or create group if it doesn't exist
   return group index
------------------------------------------------------------------------- */

int Group::find_or_create(const char *name)
{
  int igroup = find(name);
  if (igroup >= 0) return igroup;

  if (ngroup == MAX_GROUP) error->all(FLERR,"Too many groups");
  igroup = find_unused();
  int n = strlen(name) + 1;
  names[igroup] = new char[n];
  strcpy(names[igroup],name);
  ngroup++;

  return igroup;
}


/*------------------------------------------------------------------------
  assign atoms & element nodes to a new or existing group
  --------------------------------------------------------------------------*/

void Group::assign(int narg, char **arg)
{

  if (domain->box_exist == 0)
    error->all(FLERR,"Group command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal group command");

  int igroup = find(arg[0]);

  // delete the group if not being used elsewhere
  // clear mask of each atom/element/node assigned to this group

  if (strcmp(arg[1],"delete") == 0) {
    if (igroup == -1) error->all(FLERR,"Could not find group delete group ID");
    if (igroup == 0) error->all(FLERR,"Cannot delete group all");
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->igroup == igroup)
        error->all(FLERR,"Cannot delete group currently used by a fix");
    for (int i = 0; i < modify->ncompute; i++)
      if (modify->compute[i]->igroup == igroup)
        error->all(FLERR,"Cannot delete group currently used by a compute");
    for (int i = 0; i < output->ndump; i++)
      if (output->dump[i]->igroup == igroup)
        error->all(FLERR,"Cannot delete group currently used by a dump");
    //if (atom->firstgroupname && strcmp(arg[0],atom->firstgroupname) == 0)
    //  error->all(FLERR,
    //             "Cannot delete group currently used by atom_modify first");

    clear(igroup);
    
    if (dynamic[igroup]) {
      int n = strlen("GROUP_") + strlen(names[igroup]) + 1;
      char *fixID = new char[n];
      sprintf(fixID,"GROUP_%s",names[igroup]);
      modify->delete_fix(fixID);
      delete [] fixID;
    }

    delete [] names[igroup];
    names[igroup] = NULL;
    dynamic[igroup] = 0;
    ngroup--;

    return;
  }

  // clear the group

  if (strcmp(arg[1],"clear") == 0) {

    if (igroup == -1) error->all (FLERR,"Could not find group clear group ID");
    if (igroup == 0) error->all (FLERR,"Cannot clear group all");

    clear(igroup);

    return;
  }

  if (igroup == -1) {
    if (ngroup == MAX_GROUP) error->all(FLERR,"Too many groups");
    igroup = find_unused(); 
    int n = strlen(arg[0]) + 1;
    names[igroup] = new char[n];
    strcpy(names[igroup],arg[0]);
    ngroup++;
  }

  double **ax = atom->x;
  double **ex = element->x;
  double ***nodex = element->nodex;
  int *npe = element->npe;
  int *amask = atom->mask;
  int *emask = element->mask;
  int **nodemask = element->nodemask;
  int nalocal = atom->nlocal;
  int nelocal = element->nlocal;
  int bit = bitmask[igroup];
  double ***shape_array = element->shape_array; 
  int *nintpl = element->nintpl;
  int *etype = element->etype;

  // style = region
  // add to group if atom, element center, or element node is in region

  if (strcmp(arg[1],"region") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal group command");

    // option args

    // group element style determine how element is 
    // considered to be in region
    // 0 if element center is in region (default)
    // 1 if all element nodes are in region
    // 2 if at least 1 element node is in region
    // 3 if at least one interpolated atom 
    // in element is in region (i.e. region cut through element)

    int group_element_style = CENTER;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"style") == 0) {
        if (iarg + 2 > narg) error->all(FLERR,"Illegal group command");
        if (strcmp(arg[iarg+1],"center") == 0) 
          group_element_style = CENTER;
        else if (strcmp(arg[iarg+1],"allnode") == 0) 
          group_element_style = ALLNODE;
        else if (strcmp(arg[iarg+1],"onenode") == 0) 
          group_element_style = ONENODE;
        else if (strcmp(arg[iarg+1],"oneatom") == 0) 
          group_element_style = ONEATOM;
        else error->all(FLERR,"Illegal group command");
        iarg += 2;
      } else error->all(FLERR,"Illegal group command");
    }

    int iregion = domain->find_region(arg[2]);
    Region *region = domain->regions[iregion];

    if (iregion == -1) error->all(FLERR,"Group region ID does not exist");
    region->init();
    region->prematch();

    for(int i = 0; i < nalocal; i++) 
      if (region->match(ax[i]))
        amask[i] |= bit;

    int ok;
    double x[3];
    int ietype;

    // if multiple atoms per node cluster
    // assign all elements in clusters to group if at least one match 
    
    int mapflag = 0;
//    int apc = element->apc;
//    if (apc > 1) {
//      if (element->map_style == 0) {
//        mapflag = 1;
//        element->map_init();
//        element->map_set();
//      }
//    }

    for (int i = 0; i < nelocal; i++) {
      if (group_element_style == CENTER) {
        if (region->match(ex[i]))
          ok = 1;
        else ok = 0;
        for (int j = 0; j < npe[etype[i]]; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
        }
      } else if (group_element_style == ALLNODE) {
        ok = 1;
        for (int j = 0; j < npe[etype[i]]; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
          else ok = 0;
        }
      } else if (group_element_style == ONENODE) {
        ok = 0;
        for (int j = 0; j < npe[etype[i]]; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x)) {
            nodemask[i][j] |= bit;
            ok = 1;
          }
        }
      } else if (group_element_style == ONEATOM) {
        ok = 0;
        ietype = etype[i];
        for (int j = 0; j < nintpl[ietype]; j++) {
          element->evec->interpolate(x,nodex,i,j,3);
          domain->remap(x);
          if (region->match(x)) {
            ok = 1;
            break;
          }
        }

        for (int j = 0; j < npe[etype[i]]; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
        }
      }
      if (ok) {
          emask[i] |= bit;
      }
    }

    if (mapflag) {
      element->map_delete();
      element->map_style = 0;
    }

    // style = subtract

  } else if (strcmp(arg[1],"subtract") == 0) {

    if (narg < 4) error->all(FLERR,"Illegal group command");

    int length = narg-2;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 2; iarg < narg; iarg++) {
      jgroup = find(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      list[iarg-2] = jgroup;
    }

    // add to group if in 1st group in list

    int otherbit = bitmask[list[0]];

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & otherbit) amask[i] |= bit;

    for (int i = 0; i < nelocal; i++) {
      if (emask[i] & otherbit) emask[i] |= bit;
      for (int j = 0; j < npe[etype[i]]; j++)
        if (nodemask[i][j] & otherbit) nodemask[i][j] |= bit;
    }

    //remove atoms/elements/nodes if they are in any of the other groups
    //AND with inverse mask removes the atom/element/node from group

    int inverse = inversemask[igroup];

    for (int ilist = 1; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      for (int i = 0; i < nalocal; i++)
        if (amask[i] & otherbit) amask[i] &= inverse;
      for (int i = 0; i < nelocal; i++) {
        if (emask[i] & otherbit) emask[i] &= inverse;
        for (int j = 0; j < npe[etype[i]]; j++)
          if (nodemask[i][j] & otherbit) nodemask[i][j] &= inverse;
      }
    }

    delete [] list;

    // style = union

  } else if (strcmp(arg[1],"union") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal group command");

    int length = narg-2;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 2; iarg < narg; iarg++) {
      jgroup = find(arg[iarg]);
      if (jgroup==-1) error->all(FLERR,"Group ID does not exist");
      list[iarg-2] = jgroup;
    }

    // add to group if in any other group in list

    int otherbit;

    for (int ilist = 0; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      for (int i = 0; i < nalocal; i++)
        if (amask[i] & otherbit) amask[i] |= bit;

      for (int i = 0; i < nelocal; i++){
        if (emask[i] & otherbit) emask[i] |= bit;
        for (int j = 0; j < npe[etype[i]]; j++)
          if (nodemask[i][j] & otherbit) nodemask[i][j] |= bit;
      }
    }

    delete [] list;

    // style = intersect

  } else if (strcmp(arg[1],"intersect") == 0) {

    if (narg < 4) error->all(FLERR,"Illegal group command");

    int length = narg-2;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 2; iarg < narg; iarg++) {
      jgroup = find(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      if (dynamic[jgroup])
        error->all(FLERR,"Cannot intersect groups using a dynamic group");
      list[iarg-2] = jgroup;
    }

    // add to group if in all groups in list

    int otherbit,ok,ilist;

    for (int i = 0; i < nalocal; i++) {
      ok = 1;
      for (int ilist = 0; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        if ((amask[i] & otherbit) == 0) ok = 0;
      }
      if (ok) amask[i] |= bit;
    }

    for (int i = 0; i < nelocal; i++) {
      ok = 1;
      for (int ilist = 0; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        if ((emask[i] & otherbit) == 0) ok = 0;
      }
      if (ok) emask[i] |= bit;

      for (int j = 0; j < npe[etype[i]]; j++){
        ok = 1;
        for (int ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((nodemask[i][j] & otherbit) == 0) ok = 0;
        }
        if (ok) nodemask[i][j] |= bit;
      }

    }

    delete [] list;

    // style = ctype, etype, id
    // no node added, only add to group if either atom or element matches type/id or condition

  } else if (strcmp(arg[2],"ctype") == 0 ||
      strcmp(arg[2],"etype") == 0 ||
      strcmp(arg[2],"id") == 0) {

    if (narg < 4) error->all(FLERR,"Illegal group command");
    int element_flag;
    int *mask;
    int nlocal;
    if (strcmp(arg[1],"element") == 0) {
      element_flag = 1;
      mask = emask; 
      nlocal = nelocal;
    } else if (strcmp(arg[1],"atom") == 0) {
      element_flag = 0;
      mask = amask;
      nlocal = nalocal;
    } else error->all(FLERR,"Illegar group command");

    int category;

    if (strcmp(arg[2],"ctype") == 0) category = CTYPE;
    else if (strcmp(arg[2],"id") == 0) category = ID;
    else if (strcmp(arg[2],"etype") == 0 &&
        element_flag) category = ETYPE;
    else error->all(FLERR,"Illegal group command");

    int *attribute = NULL;
    tagint *tattribute = NULL;

    if (category == ETYPE) attribute = element->etype;
    else if (category == CTYPE) {
      if (element_flag) attribute = element->ctype;
      else attribute = atom->type;
    } else if (category == ID) {
      if (element_flag) tattribute = element->tag;
      else tattribute = atom->tag;
    }

    // args = logical condition

    if (narg > 4 &&
        (strcmp(arg[3],"<") == 0 || 
         strcmp(arg[3],">") == 0 ||
         strcmp(arg[3],"<=") == 0 || 
         strcmp(arg[3],">=") == 0 ||
         strcmp(arg[3],"==") == 0 || 
         strcmp(arg[3],"!=") == 0 ||
         strcmp(arg[3],"<>") == 0)) {

      int condition = -1;
      if (strcmp(arg[3],"<") == 0) condition = LT;
      else if (strcmp(arg[3],"<=") == 0) condition = LE;
      else if (strcmp(arg[3],">") == 0) condition = GT;
      else if (strcmp(arg[3],">=") == 0) condition = GE;
      else if (strcmp(arg[3],"==") == 0) condition = EQ;
      else if (strcmp(arg[3],"!=") == 0) condition = NEQ;
      else if (strcmp(arg[3],"<>") == 0) condition = BETWEEN;
      else if (strcmp(arg[3],"><") == 0) condition = NOTBETWEEN;
      else error->all(FLERR,"Illegal group command");

      tagint bound1,bound2;
      bound1 = universe->tnumeric(FLERR,arg[4]);
      bound2 = -1;

      if (condition == BETWEEN || condition == NOTBETWEEN) {
        if (narg != 6) error->all(FLERR,"Illegal group command");
        bound2 = universe->tnumeric(FLERR,arg[5]);
        if (bound1 > bound2) error->all(FLERR,"Illegal group command");
      } else if (narg != 5) error->all(FLERR,"Illegal group command");

      // add to group if meets condition

      if (attribute) {
        if (condition == LT) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] < bound1) mask[i] |= bit;
        } else if (condition == LE) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] <= bound1) mask[i] |= bit;
        } else if (condition == GT) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] > bound1) mask[i] |= bit;
        } else if (condition == GE) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] >= bound1) mask[i] |= bit;
        } else if (condition == EQ) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] == bound1) mask[i] |= bit;
        } else if (condition == NEQ) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] != bound1) mask[i] |= bit;
        } else if (condition == BETWEEN) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] >= bound1 && attribute[i] <= bound2)
              mask[i] |= bit;
        } else if (condition == NOTBETWEEN) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] <= bound1 || attribute[i] >= bound2)
              mask[i] |= bit;
        }
      } else {
        if (condition == LT) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] < bound1) mask[i] |= bit;
        } else if (condition == LE) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] <= bound1) mask[i] |= bit;
        } else if (condition == GT) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] > bound1) mask[i] |= bit;
        } else if (condition == GE) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] >= bound1) mask[i] |= bit;
        } else if (condition == EQ) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] == bound1) mask[i] |= bit;
        } else if (condition == NEQ) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] != bound1) mask[i] |= bit;
        } else if (condition == BETWEEN) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] >= bound1 && tattribute[i] <= bound2)
              mask[i] |= bit;
        } else if (condition == NOTBETWEEN) {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] <= bound1 || tattribute[i] >= bound2)
              mask[i] |= bit;
        }
      }

      // args = list of values

    } else {
      char *ptr;
      tagint start,stop,delta;

      for (int iarg = 3; iarg < narg; iarg++) {
        delta = 1;
        if (strchr(arg[iarg],':')) {
          ptr = strtok(arg[iarg],":");
          start = universe->tnumeric(FLERR,ptr);
          ptr = strtok(NULL,":");
          stop = universe->tnumeric(FLERR,ptr);
          ptr = strtok(NULL,":");
          if (ptr) delta = universe->tnumeric(FLERR,ptr);
        } else {
          start = stop = universe->tnumeric(FLERR,arg[iarg]);
        }
        if (delta < 1)
          error->all(FLERR,"Illegal range increment value");

        // add to group if attribute matches value or sequence

        if (attribute) {
          for (int i = 0; i < nlocal; i++)
            if (attribute[i] >= start && attribute[i] <= stop &&
                (attribute[i]-start) % delta == 0) mask[i] |= bit;
        } else {
          for (int i = 0; i < nlocal; i++)
            if (tattribute[i] >= start && tattribute[i] <= stop &&
                (tattribute[i]-start) % delta == 0) mask[i] |= bit;
        }
      }
    } 
  } else error->all(FLERR,"Illegal group command");

  // print stats for changed group

  int n,m,k;
  n = m = k = 0;

  for (int i = 0; i < nalocal; i++) if (amask[i] & bit) n++;
  for (int i = 0; i < nelocal; i++) {
    if (emask[i] & bit) k++;
    for (int j = 0; j < npe[etype[i]]; j++)
      if (nodemask[i][j] & bit) m++;
  }

  double atomlocal = n;
  double atomall;

  MPI_Allreduce(&atomlocal,&atomall,1,MPI_DOUBLE,MPI_SUM,world);

  double elemlocal = k;
  double elemall;

  MPI_Allreduce(&elemlocal,&elemall,1,MPI_DOUBLE,MPI_SUM,world);

  double nodelocal = m;
  double nodeall;

  MPI_Allreduce(&nodelocal,&nodeall,1,MPI_DOUBLE,MPI_SUM,world);

  FILE *out;
  if (me == 0) 
    for (int i = 0; i < 2; i++) {
      if (i == 0) out = screen;
      else out = logfile;
      if (out) {
        if (atomall) fprintf(out,"%.15g atoms in group %s\n",atomall,names[igroup]);
        if (elemall) fprintf(out,"%.15g elements in group %s\n",elemall,names[igroup]);
        if (nodeall) fprintf(out,"%.15g nodes in group %s\n",nodeall,names[igroup]);
      }
    }
}

/*-------------------------------------------------------------------------
  return index of first available group
  should never be called when group limit has been reached
  ---------------------------------------------------------------------------*/

int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == NULL) return igroup;
  return -1;
}

/*-------------------------------------------------------------------------
  clear mask of each atom/element/node assigned to igroup
  ---------------------------------------------------------------------------*/

void Group::clear(int igroup) 
{
  int *amask = atom->mask;
  int *emask = element->mask;
  int **nmask = element->nodemask;
  int bits = inversemask[igroup];
  int *npe = element->npe;
  int *etype = element->etype;
  for (int i = 0; i < atom->nlocal; i++) amask[i] &= bits;
  for (int i = 0; i < element->nlocal; i++) {
    emask[i] &= bits;
    for (int j = 0; j < npe[etype[i]]; j++)
      nmask[i][j] &= bits;
  }
}

// ----------------------------------------------------------------------
// computations on a group of atoms
// ----------------------------------------------------------------------


/* ----------------------------------------------------------------------
   count atoms in group and region
   ------------------------------------------------------------------------- */

bigint Group::count_atom(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i])) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count nodes in group and region using element mask
   ------------------------------------------------------------------------- */

bigint Group::count_element_node(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = element->x;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i])) n += npe[etype[i]];

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count nodes in group and region using nodemask
   ------------------------------------------------------------------------- */

bigint Group::count_node(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double ***nodex = element->nodex;
  int **nodemask = element->nodemask;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < npe[etype[i]]; j++)
      if (nodemask[i][j] & groupbit && region->match(nodex[i][j])) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count interpolated atoms in group and region
   ------------------------------------------------------------------------- */

bigint Group::count_intpl(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = element->x;
  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *nintpl = element->nintpl;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i])) n += nintpl[etype[i]];

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count elements in group and region
   ------------------------------------------------------------------------- */

bigint Group::count_elem(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = element->x;
  int *mask = element->mask;
  int nlocal = element->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i])) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count atoms in group
   ------------------------------------------------------------------------- */

bigint Group::count_atom(int igroup)
{
  int groupbit = bitmask[igroup];

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count elements in group
   ------------------------------------------------------------------------- */

bigint Group::count_elem(int igroup)
{
  int groupbit = bitmask[igroup];

  int *mask = element->mask;
  int nlocal = element->nlocal;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count nodes in group using element mask
   ------------------------------------------------------------------------- */

bigint Group::count_element_node(int igroup)
{
  int groupbit = bitmask[igroup];

  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) n += npe[etype[i]];

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count nodes in group using nodemask
   ------------------------------------------------------------------------- */

bigint Group::count_node(int igroup)
{
  int groupbit = bitmask[igroup];

  int **mask = element->nodemask;
  int nlocal = element->nlocal;
  int *npe = element->npe;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < npe[etype[i]]; j++)
      if (mask[i][j] & groupbit) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   count interpolated atoms in group
   ------------------------------------------------------------------------- */

bigint Group::count_intpl(int igroup)
{
  int groupbit = bitmask[igroup];

  int *mask = element->mask;
  int nlocal = element->nlocal;
  int *nintpl = element->nintpl;
  int *etype = element->etype;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n += nintpl[etype[i]];

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   compute the coordinate bounds of the group of atoms/elements
   periodic images are not considered, so atoms/elements are NOT unwrapped
------------------------------------------------------------------------- */

void Group::bounds(int igroup, double *minmax)
{
  int groupbit = bitmask[igroup];

  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;

  double **x,***nodex;
  int *mask,*etype,*npe;
  int nlocal,i,j;

  x = atom->x;
  mask = atom->mask;
  nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }

  nodex = element->nodex;
  etype = element->etype;
  npe = element->npe;
  mask = element->mask;
  nlocal = element->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (j = 0; j < npe[etype[i]]; j++) {
        extent[0] = MIN(extent[0],nodex[i][j][0]);
        extent[1] = MAX(extent[1],nodex[i][j][0]);
        extent[2] = MIN(extent[2],nodex[i][j][1]);
        extent[3] = MAX(extent[3],nodex[i][j][1]);
        extent[4] = MIN(extent[4],nodex[i][j][2]);
        extent[5] = MAX(extent[5],nodex[i][j][2]);
      }
    }
  }

  // compute extent across all procs
  // flip sign of MIN to do it in one Allreduce MAX
  // set box by extent in shrink-wrapped dims

  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];

  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);

  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}

/* ----------------------------------------------------------------------
   compute the coordinate bounds of the group of atoms/elements in region
   periodic images are not considered, so atoms/elements are NOT unwrapped
   ------------------------------------------------------------------------- */

void Group::bounds(int igroup, double *minmax, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;

  double **x,***nodex;
  int *mask,*etype,*npe;
  int nlocal,i,j;

  x = atom->x;
  mask = atom->mask;
  nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }

  x = element->x;
  nodex = element->nodex;
  etype = element->etype;
  npe = element->npe;
  mask = element->mask;
  nlocal = element->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      for (j = 0; j < npe[etype[i]]; j++) {
        extent[0] = MIN(extent[0],nodex[i][j][0]);
        extent[1] = MAX(extent[1],nodex[i][j][0]);
        extent[2] = MIN(extent[2],nodex[i][j][1]);
        extent[3] = MAX(extent[3],nodex[i][j][1]);
        extent[4] = MIN(extent[4],nodex[i][j][2]);
        extent[5] = MAX(extent[5],nodex[i][j][2]);
      }
    }
  }

  // compute extent across all procs
  // flip sign of MIN to do it in one Allreduce MAX
  // set box by extent in shrink-wrapped dims

  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];

  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);

  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}

/* ----------------------------------------------------------------------
   compute the center-of-mass coords of group of atoms
   masstotal = total mass
   return center-of-mass coords in cm[]
   must unwrap atoms to compute center-of-mass correctly
   ------------------------------------------------------------------------- */

void Group::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype  = element->etype;
  int *ectype = element->ctype;
  int *nintpl = element->nintpl;
  double **ex = element->x;
  imageint *eimage = element->image;

  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;

  double massone;
  double unwrap[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      cmone[0] += unwrap[0] * massone;
      cmone[1] += unwrap[1] * massone;
      cmone[2] += unwrap[2] * massone;
    }

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit) {
      massone = mass[ectype[i]] * nintpl[etype[i]];
      domain->unmap(ex[i],eimage[i],unwrap);
      cmone[0] += unwrap[0] * massone;
      cmone[1] += unwrap[1] * massone;
      cmone[2] += unwrap[2] * massone;
    }

  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass coords of group of atoms in region
   mastotal = total mass
   return center-of-mass coords in cm[]
   must unwrap atoms to compute center-of-mass correctly
   ------------------------------------------------------------------------- */

void Group::xcm(int igroup, double masstotal, double *cm, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype  = element->etype;
  int *ectype = element->ctype;
  int *nintpl = element->nintpl;
  double **ex = element->x;
  imageint *eimage = element->image;

  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;

  double massone;
  double unwrap[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      massone = mass[type[i]];
      domain->unmap(x[i],image[i],unwrap);
      cmone[0] += unwrap[0] * massone;
      cmone[1] += unwrap[1] * massone;
      cmone[2] += unwrap[2] * massone;
    }

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit && region->match(ex[i][0],ex[i][1],ex[i][2])) {
      massone = mass[ectype[i]] * nintpl[etype[i]];
      domain->unmap(ex[i],eimage[i],unwrap);
      cmone[0] += unwrap[0] * massone;
      cmone[1] += unwrap[1] * massone;
      cmone[2] += unwrap[2] * massone;
    }

  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms
   use either per-type mass 
   ------------------------------------------------------------------------- */

double Group::mass(int igroup)
{
  int groupbit = bitmask[igroup];

  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double one = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) one += mass[type[i]];

  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit) 
      one += mass[ctype[i]] * nintpl[etype[i]];

  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms in region
   use either per-type mass 
   ------------------------------------------------------------------------- */

double Group::mass(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double one = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
      one += mass[type[i]];

  double **ex = element->x;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
      one += mass[ctype[i]] * nintpl[etype[i]];

  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   compute the total kinetic energy of group of atoms and return it
   ------------------------------------------------------------------------- */

double Group::ke(int igroup)
{
  int groupbit = bitmask[igroup];

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double one = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
        mass[type[i]];

  double ***nodev = element->nodev;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  int *npe = element->npe;
  double *nodal_weight = element->nodal_weight;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit)
      for (int j = 0; j < npe[etype[i]]; j++)
        one += (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]) *
          mass[ctype[i]] * nodal_weight[etype[i]];

  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}

/* ----------------------------------------------------------------------
   compute the total kinetic energy of group of atoms in region and return it
   ------------------------------------------------------------------------- */

double Group::ke(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double one = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
      one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
        mass[type[i]];

  double ***nodev = element->nodev;
  double **ex = element->x;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  int *npe = element->npe;
  double *nodal_weight = element->nodal_weight;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit && region->match(ex[i][0],ex[i][1],ex[i][2]))
      for (int j = 0; j < npe[etype[i]]; j++)
        one += (nodev[i][j][0]*nodev[i][j][0] + nodev[i][j][1]*nodev[i][j][1] + nodev[i][j][2]*nodev[i][j][2]) *
          mass[ctype[i]] * nodal_weight[etype[i]];

  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms
   masstotal = total mass
   return center-of-mass velocity in cm[]
   ------------------------------------------------------------------------- */

void Group::vcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      p[0] += v[i][0]*massone;
      p[1] += v[i][1]*massone;
      p[2] += v[i][2]*massone;
    }

  double ***nodev = element->nodev;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  int *npe = element->npe;
  double *nodal_weight = element->nodal_weight;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit) 
      for (int j = 0; j < npe[etype[i]]; j++) {
        massone = mass[ctype[i]] * nodal_weight[etype[i]];
        p[0] += nodev[i][j][0]*massone;
        p[1] += nodev[i][j][1]*massone;
        p[2] += nodev[i][j][2]*massone;
      }

  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms in region
   masstotal = total mass
   return center-of-mass velocity in cm[]
   ------------------------------------------------------------------------- */

void Group::vcm(int igroup, double masstotal, double *cm, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      massone = mass[type[i]];
      p[0] += v[i][0]*massone;
      p[1] += v[i][1]*massone;
      p[2] += v[i][2]*massone;
    }

  double **ex = element->x;
  double ***nodev = element->nodev;
  int *emask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  int *npe = element->npe;
  double *nodal_weight = element->nodal_weight;

  for (int i = 0; i < nelocal; i++)
    if (emask[i] & groupbit && region->match(ex[i][0],ex[i][1],ex[i][2])) {
      for (int j = 0; j < npe[etype[i]]; j++) {
        massone = mass[ctype[i]] * nodal_weight[etype[i]];
        p[0] += nodev[i][j][0]*massone;
        p[1] += nodev[i][j][1]*massone;
        p[2] += nodev[i][j][2]*massone;
      }
    }
  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the angular momentum L (lmom) of group
   around center-of-mass cm
   must unwrap atoms and nodes to compute L correctly
------------------------------------------------------------------------- */

void Group::angmom(int igroup, double *cm, double *lmom) 
{
  int groupbit = bitmask[igroup];
  int i,j;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double dx,dy,dz,massone;
  double unwrap[3];

  double p[3];
  p[0] = p[1] = p[2] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      massone = mass[type[i]];
      p[0] += massone * (dy * v[i][2] - dz * v[i][1]);
      p[1] += massone * (dz * v[i][0] - dx * v[i][2]);
      p[2] += massone * (dx * v[i][1] - dy * v[i][0]);
    }

  mask = element->mask;
  int nelocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  image = element->image;
  double coord[3],vel[3];

  for (i = 0; i < nelocal; i++) {
    massone = mass[ctype[i]];
    if (mask[i] & groupbit) {
      for (j = 0; j < nintpl[etype[i]]; j++) {
        element->evec->interpolate(coord,nodex,i,j,3);
        element->evec->interpolate(vel,nodev,i,j,3);
        domain->unmap(coord,image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        p[0] += massone * (dy*vel[2] - dz*vel[1]);
        p[1] += massone * (dz*vel[0] - dx*vel[2]);
        p[2] += massone * (dx*vel[1] - dy*vel[0]);
      }
    }
  }

  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   compute the angular momentum L (lmom) of group
   around center-of-mass cm
   must unwrap atoms and nodes to compute L correctly
------------------------------------------------------------------------- */

void Group::angmom(int igroup, double *cm, double *lmom, int iregion) 
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  int i,j;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double dx,dy,dz,massone;
  double unwrap[3];

  double p[3];
  p[0] = p[1] = p[2] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      massone = mass[type[i]];
      p[0] += massone * (dy * v[i][2] - dz * v[i][1]);
      p[1] += massone * (dz * v[i][0] - dx * v[i][2]);
      p[2] += massone * (dx * v[i][1] - dy * v[i][0]);
    }

  mask = element->mask;
  nlocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  double ***nodex = element->nodex;
  double ***nodev = element->nodev;
  image = element->image;
  double coord[3],vel[3];

  for (i = 0; i < nlocal; i++) {
    massone = mass[ctype[i]];
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      for (j = 0; j < nintpl[etype[i]]; j++) {
        element->evec->interpolate(coord,nodex,i,j,3);
        element->evec->interpolate(vel,nodev,i,j,3);
        domain->unmap(coord,image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        p[0] += massone * (dy*vel[2] - dz*vel[1]);
        p[1] += massone * (dz*vel[0] - dx*vel[2]);
        p[2] += massone * (dx*vel[1] - dy*vel[0]);
      }
    }
  }

  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}


/* ----------------------------------------------------------------------
   compute moment of inertia tensor around center-of-mass cm of group
   must unwrap atoms and elments to compute itensor correctly
   ------------------------------------------------------------------------- */

void Group::inertia(int igroup, double *cm, double itensor[3][3]) 
{
  int i,j;

  int groupbit = bitmask[igroup];

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double dx,dy,dz,massone;
  double unwrap[3];

  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }

  mask = element->mask;
  nlocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  double ***nodex = element->nodex;
  image = element->image;
  double coord[3];

  for (i = 0; i < nlocal; i++) {
    massone = mass[ctype[i]];
    if (mask[i] & groupbit) {
      for (j = 0; j < nintpl[etype[i]]; j++) {
        element->evec->interpolate(coord,nodex,i,j,3);
        domain->unmap(coord,image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        ione[0][0] += massone * (dy*dy + dz*dz);
        ione[1][1] += massone * (dx*dx + dz*dz);
        ione[2][2] += massone * (dx*dx + dy*dy);
        ione[0][1] -= massone * dx*dy;
        ione[1][2] -= massone * dy*dz;
        ione[0][2] -= massone * dx*dz;
      }
    }
  }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];

  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   compute moment of inertia tensor around center-of-mass cm of group
   must unwrap atoms and elments to compute itensor correctly
   ------------------------------------------------------------------------- */

void Group::inertia(int igroup, double *cm, double itensor[3][3], int iregion) 
{
  int i,j;

  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double dx,dy,dz,massone;
  double unwrap[3];

  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }

  x = element->x;
  mask = element->mask;
  nlocal = element->nlocal;
  int *etype = element->etype;
  int *ctype = element->ctype;
  int *nintpl = element->nintpl;
  double ***nodex = element->nodex;
  image = element->image;
  double coord[3];

  for (i = 0; i < nlocal; i++) {
    massone = mass[ctype[i]];
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      for (j = 0; j < nintpl[etype[i]]; j++) {
        element->evec->interpolate(coord,nodex,i,j,3);
        domain->unmap(coord,image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        ione[0][0] += massone * (dy*dy + dz*dz);
        ione[1][1] += massone * (dx*dx + dz*dz);
        ione[2][2] += massone * (dx*dx + dy*dy);
        ione[0][1] -= massone * dx*dy;
        ione[1][2] -= massone * dy*dz;
        ione[0][2] -= massone * dx*dz;
      }
    }
  }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];

  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   compute angular velocity omega from L and I
   ------------------------------------------------------------------------- */

void Group::omega(double *angmom, double inertia[3][3], double *w) 
{
  double idiag[3], ex[3], ey[3], ez[3], cross[3];
  double evectors[3][3], inverse[3][3];

  // determinant = triple product of rows of inertia matrix

  double determinant = inertia[0][0] * inertia[1][1] * inertia[2][2] +
    inertia[0][1] * inertia[1][2] * inertia[2][0] +
    inertia[0][2] * inertia[1][0] * inertia[2][1] -
    inertia[0][0] * inertia[1][2] * inertia[2][1] -
    inertia[0][1] * inertia[1][0] * inertia[2][2] -
    inertia[2][0] * inertia[1][1] * inertia[0][2];

  // non-singular I matrix
  // use L = Iw, inverting I to solve for w
  // this should give exact zeroing of angular momentum by velocity command

  if (determinant > EPSILON) {

    inverse[0][0] =
      inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1];
    inverse[0][1] =
      -(inertia[0][1] * inertia[2][2] - inertia[0][2] * inertia[2][1]);
    inverse[0][2] =
      inertia[0][1] * inertia[1][2] - inertia[0][2] * inertia[1][1];

    inverse[1][0] =
      -(inertia[1][0] * inertia[2][2] - inertia[1][2] * inertia[2][0]);
    inverse[1][1] =
      inertia[0][0] * inertia[2][2] - inertia[0][2] * inertia[2][0];
    inverse[1][2] =
      -(inertia[0][0] * inertia[1][2] - inertia[0][2] * inertia[1][0]);

    inverse[2][0] =
      inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0];
    inverse[2][1] =
      -(inertia[0][0] * inertia[2][1] - inertia[0][1] * inertia[2][0]);
    inverse[2][2] =
      inertia[0][0] * inertia[1][1] - inertia[0][1] * inertia[1][0];

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        inverse[i][j] /= determinant;

    w[0] = inverse[0][0] * angmom[0] + inverse[0][1] * angmom[1] +
      inverse[0][2] * angmom[2];
    w[1] = inverse[1][0] * angmom[0] + inverse[1][1] * angmom[1] +
      inverse[1][2] * angmom[2];
    w[2] = inverse[2][0] * angmom[0] + inverse[2][1] * angmom[1] +
      inverse[2][2] * angmom[2];

    // handle (nearly) singular I matrix
    // typically due to 2-atom group or linear molecule
    // use jacobi() and angmom_to_omega() to calculate valid omega
    // less exact answer than matrix inversion, due to iterative Jacobi method

  } else {
    int ierror = MathExtra::jacobi(inertia, idiag, evectors);
    if (ierror)
      error->all(FLERR, "Insufficient Jacobi rotations for group::omega");

    ex[0] = evectors[0][0];
    ex[1] = evectors[1][0];
    ex[2] = evectors[2][0];
    ey[0] = evectors[0][1];
    ey[1] = evectors[1][1];
    ey[2] = evectors[2][1];
    ez[0] = evectors[0][2];
    ez[1] = evectors[1][2];
    ez[2] = evectors[2][2];

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed

    MathExtra::cross3(ex, ey, cross);
    if (MathExtra::dot3(cross, ez) < 0.0)
      MathExtra::negate3(ez);

    // if any principal moment < scaled EPSILON, set to 0.0

    double max;
    max = MAX(idiag[0], idiag[1]);
    max = MAX(max, idiag[2]);

    if (idiag[0] < EPSILON * max)
      idiag[0] = 0.0;
    if (idiag[1] < EPSILON * max)
      idiag[1] = 0.0;
    if (idiag[2] < EPSILON * max)
      idiag[2] = 0.0;

    // calculate omega using diagonalized inertia matrix

    MathExtra::angmom_to_omega(angmom, ex, ey, ez, idiag, w);
  }
}
