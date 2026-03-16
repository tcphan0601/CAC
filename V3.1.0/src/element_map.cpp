#include <math.h>
#include "element.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace CAC_NS;

#define EXTRA 1000

/*  ----------------------------------------------------------------------
   allocate and initialize array or hash table for global -> local map
   for array option:
     array length = 1 to map_tag_max
     set entire array to -1 as initial values
   for hash option:
     map_nhash = length of hash table
     map_nbucket = # of hash buckets, prime larger than map_nhash * 2
       so buckets will only be filled with 0 or 1 elements on average
-------------------------------------------------------------------------  */

void Element::map_init(int check)
{
  // check for new map style if max elementID changed (check = 1 = default)
  // recreate = 1 if must delete old map and create new map
  // recreate = 0 if can re-use old map w/out realloc and just adjust settings
  // map_maxarray/map_nhash initially -1, to force recreate even when no elements

  int recreate = 0;
  if (check) recreate = map_style_set();

  if (map_style == 1 && map_tag_max > map_maxarray) recreate = 1;
  else if (map_style == 2 && nlocal+nghost > map_nhash) recreate = 1;

  // if not recreating:
  // for array, initialize current map_tag_max values
  // for hash, set all buckets to empty, put all entries in free list

  if (!recreate) {
    if (map_style == 1) {
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
    } else {
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
      map_hash[map_nhash-1].next = -1;
    }

  // recreating: delete old map and create new one for array or hash

  } else {
    map_delete();

    if (map_style == 1) {
      map_maxarray = map_tag_max;
      memory->create(map_array, map_maxarray+1, "Element:map_array");
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;

    } else {

      // map_nhash = max # of elements that can be hashed on this proc
      // set to max of ave elements/proc or elements I can store
      // multiply by 2, require at least 1000
      // doubling means hash table will need to be re-init only rarely

      int nper = static_cast<int> (nelements/comm->nprocs);
      map_nhash = MAX(nper, nmax);
      map_nhash *= 2;
      map_nhash = MAX(map_nhash, 1000);

      // map_nbucket = prime just larger than map_nhash
      // next_prime() should be fast enough, 
      //   about 10% of odd integers are prime above 1M

      map_nbucket = universe->next_prime(map_nhash);

      // set all buckets to empty
      // set hash to map_nhash in length
      // put all hash entries in free list and point them to each other

      map_bucket = new int[map_nbucket];
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

      map_hash = new HashElem[map_nhash];
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
      map_hash[map_nhash-1].next = -1;
    }
  }
}

/*  ----------------------------------------------------------------------
   clear global -> local map for all of my own and ghost elements
   for hash table option:
     global ID may not be in table if image element was already cleared
-------------------------------------------------------------------------  */

void Element::map_clear()
{
  if (map_style == 1) {
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      sametag[i] = -1;
      map_array[tag[i]] = -1;
    }

  } else {
    int previous, ibucket, index;
    tagint global;
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      sametag[i] = -1;

      // search for key
      // if don't find it, done

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index == -1) continue;

      // delete the hash entry and add it to free list
      // special logic if entry is 1st in the bucket

      if (previous == -1) map_bucket[ibucket] = map_hash[index].next;
      else map_hash[previous].next = map_hash[index].next;

      map_hash[index].next = map_free;
      map_free = index;
      map_nused--;
    }
  }
}

/*  ----------------------------------------------------------------------
   set global -> local map for all of my own and ghost elements
   loop in reverse order so that nearby images take precedence over far ones
     and owned elements take precedence over images
   for hash table option:
     if hash table too small, re-init
     global ID may already be in table if image element was set
-------------------------------------------------------------------------  */

void Element::map_set()
{
  int nall = nlocal + nghost;

  if (map_style == 1) {

    // possible reallocation of sametag must come before loop over elements
    // since loop sets sametag

    if (nall > max_same) {
      max_same = nall + EXTRA;
      memory->destroy(sametag);
      memory->create(sametag, max_same, "Element:sametag");
    }

    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
    }

  } else {

    // if this proc has more elements than hash table size, call map_init()
    //   call with 0 since max elementID in system has not changed
    // possible reallocation of sametag must come after map_init(), 
    //   b/c map_init() may invoke map_delete(), whacking sametag

    if (nall > map_nhash) map_init(0);
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memory->destroy(sametag);
      memory->create(sametag, max_same, "Element:sametag");
    }

    int previous, ibucket, index;
    tagint global;

    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_find_hash(tag[i]);

      // search for key
      // if found it, just overwrite local value with index

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index > -1) {
        map_hash[index].local = i;
        continue;
      }

      // take one entry from free list
      // add the new global/local pair as entry at end of bucket list
      // special logic if this entry is 1st in bucket

      index = map_free;
      map_free = map_hash[map_free].next;
      if (previous == -1) map_bucket[ibucket] = index;
      else map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }
  }
}

/*  ----------------------------------------------------------------------
   set global to local map for one element
   for hash table option:
     global ID may already be in table if element was already set
   called by Special class
-------------------------------------------------------------------------  */

void Element::map_one(tagint global, int local)
{
  if (map_style == 1) map_array[global] = local;
  else {
    // search for key
    // if found it, just overwrite local value with index

    int previous = -1;
    int ibucket = global % map_nbucket;
    int index = map_bucket[ibucket];
    while (index > -1) {
      if (map_hash[index].global == global) break;
      previous = index;
      index = map_hash[index].next;
    }
    if (index > -1) {
      map_hash[index].local = local;
      return;
    }

    // take one entry from free list
    // add the new global/local pair as entry at end of bucket list
    // special logic if this entry is 1st in bucket

    index = map_free;
    map_free = map_hash[map_free].next;
    if (previous == -1) map_bucket[ibucket] = index;
    else map_hash[previous].next = index;
    map_hash[index].global = global;
    map_hash[index].local = local;
    map_hash[index].next = -1;
    map_nused++;
  }
}

/*  ----------------------------------------------------------------------
   set map style to array or hash based on user request or max elementID
   set map_tag_max = max element ID (may be larger than nelements)
   called whenever map_init() called with new total # of elements
   return 1 if map_style changed, else 0
-------------------------------------------------------------------------  */

int Element::map_style_set()
{
  if (tag_enable == 0)
    error->all(FLERR, "Cannot create an element map unless elements have IDs");

  // map_tag_max = max ID of any element that will be in new map
  // map_tag_max = -1 if no elements

  tagint max = -1;
  for (int i = 0; i < nlocal; i++) max = MAX(max, tag[i]);
 
  MPI_Allreduce(&max, &map_tag_max, 1, MPI_CAC_TAGINT, MPI_MAX, world);

  // set map_style for new map
  // if user-selected, use that setting
  // else if map_tag_max > 1M, use hash
  // else use array

  int map_style_old = map_style;
  if (map_user) map_style = map_user;
  else if (map_tag_max > 1000000) map_style = 2;
  else map_style = 1;

  // recreate = 1 if must create new map b/c map_style changed

  int recreate = 0;
  if (map_style != map_style_old) recreate = 1;
  return recreate;
}

/*  ----------------------------------------------------------------------
   free the array or hash table for global to local mapping
-------------------------------------------------------------------------  */

void Element::map_delete()
{
  memory->destroy(sametag);
  sametag = nullptr;
  max_same = 0;

  if (map_style == 1) {
    memory->destroy(map_array);
    map_array = nullptr;
  } else {
    if (map_nhash) {
      delete [] map_bucket;
      delete [] map_hash;
      map_bucket = nullptr;
      map_hash = nullptr;
    }
    map_nhash = 0;
  }
}

/*  ----------------------------------------------------------------------
   lookup global ID in hash table, return local index
   called by map() in element.h
-------------------------------------------------------------------------  */

int Element::map_find_hash(tagint global)
{
  int local = -1;
  int index = map_bucket[global % map_nbucket];
  while (index > -1) {
    if (map_hash[index].global == global) {
      local = map_hash[index].local;
      break;
    }
    index = map_hash[index].next;
  }
  return local;
}

