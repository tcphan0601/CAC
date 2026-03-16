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

enum{ETYPE,CTYPE,ID};
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

  // create "element" group containing nodes of all elements
  
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
  int npe = element->npe;
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

    // group element_style determine how element is 
    // considered to be in region
    // 0 if element center is in region (default)
    // 1 if all element nodes are in region
    // 2 if at least 1 element node is in region
    // 3 if at least one interpolated atom 
    // in element is in region (i.e. region cut through element)

    int group_element_style = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"style") == 0) {
        if (iarg + 2 > narg) error->all(FLERR,"Illegal group command");
        if (strcmp(arg[iarg+1],"center") == 0) 
          group_element_style = 0;
        else if (strcmp(arg[iarg+1],"allnode") == 0) 
          group_element_style = 1;
        else if (strcmp(arg[iarg+1],"onenode") == 0) 
          group_element_style = 2;
        else if (strcmp(arg[iarg+1],"oneatom") == 0) 
          group_element_style = 3;
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

    for(int i = 0; i < nelocal; i++) {
      if (group_element_style == 0) {
        if (region->match(ex[i]))
          ok = 1;
        else ok = 0;
        for (int j = 0; j < npe; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
        }
      } else if (group_element_style == 1) {
        ok = 1;
        for (int j = 0; j < npe; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
          else ok = 0;
        }
      } else if (group_element_style == 2) {
        ok = 0;
        for (int j = 0; j < npe; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x)) {
            nodemask[i][j] |= bit;
            ok = 1;
          }
        }
      } else if (group_element_style == 3) {
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

        for (int j = 0; j < npe; j++) {
          copy3(nodex[i][j],x);
          domain->remap(x);
          if (region->match(x))
            nodemask[i][j] |= bit;
        }
      }
      if (ok) emask[i] |= bit;
    }

    // style = subtract

  } else if (strcmp(arg[1],"subtract") == 0) {

    if (narg < 4) error->all(FLERR,"Illegal group command");

    int length = narg-2;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 2; iarg < narg ;iarg++) {
      jgroup = find(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      list[iarg-2] = jgroup;
    }

    // add to group if in 1st group in list

    int otherbit = bitmask[list[0]];

    for (int i = 0; i < nalocal; i++)
      if (amask[i] & otherbit) amask[i] |= bit;

    for (int i = 0; i < nelocal;i++) {
      if (emask[i] & otherbit) emask[i] |= bit;
      for (int j = 0; j < npe;j++)
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
        for (int j = 0; j < npe; j++)
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

    for (int ilist = 0; ilist < length;ilist++) {
      otherbit = bitmask[list[ilist]];
      for (int i = 0; i < nalocal; i++)
        if (amask[i] & otherbit) amask[i] |= bit;

      for (int i = 0; i < nelocal; i++){
        if (emask[i] & otherbit) emask[i] |= bit;
        for (int j = 0; j < npe; j++)
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

      for (int j = 0; j < npe; j++){
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
    for (int j = 0; j < npe; j++)
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
  for (int i = 0; i < atom->nlocal; i++) amask[i] &= bits;
  for (int i = 0; i < element->nlocal; i++) {
    emask[i] &= bits;
    for (int j = 0; j < element->npe; j++)
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
   count nodes in group and region
------------------------------------------------------------------------- */

bigint Group::count_node(int igroup, int iregion)
{
  int groupbit = bitmask[igroup];
  Region *region = domain->regions[iregion];
  region->prematch();

  double ***nodex = element->nodex;
  int **nodemask = element->nodemask;
  int nlocal = element->nlocal;
  int npe = element->npe;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < npe; j++)
      if (nodemask[i][j] & groupbit && region->match(nodex[i][j])) n++;

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
   count nodes in group
------------------------------------------------------------------------- */

bigint Group::count_node(int igroup)
{
  int groupbit = bitmask[igroup];

  int **mask = element->nodemask;
  int nlocal = element->nlocal;
  int npe = element->npe;

  int n = 0;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < npe; j++)
      if (mask[i][j] & groupbit) n++;

  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_CAC_BIGINT,MPI_SUM,world);
  return nall;
}


