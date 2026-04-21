#ifndef CAC_MEMORY_H
#define CAC_MEMORY_H

#include "pointers.h"
#include <unordered_map>
#include <vector>
#include "universe.h"

namespace CAC_NS {

class Memory : protected Pointers {
 public:
  Memory(class CAC *);

  void *smalloc(bigint n, const char *);
  void *srealloc(void *, bigint n, const char *);
  void sfree(void *);
  void fail(const char *);

/*  ----------------------------------------------------------------------
   create/grow/destroy vecs and multidim arrays with contiguous memory blocks
   only use with primitive data types, e.g. 1d vec of ints, 2d array of doubles
   fail() prevents use with pointers, 
     e.g. 1d vec of int * , due to mismatched destroy
   avoid use with non-primitive data types to avoid code bloat
   for these other cases, use smalloc/srealloc/sfree directly
-------------------------------------------------------------------------  */

/*  ----------------------------------------------------------------------
   create a 1d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE *create(TYPE *&array, int n, const char *name)
  {

    bigint nbytes = ((bigint) sizeof(TYPE)) * n;
    array = (TYPE *) smalloc(nbytes, name);
    return array;
  }

  template <typename TYPE>
  TYPE **create(TYPE **& /* array */, int /* n */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   grow or shrink 1d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE *grow(TYPE *&array, int n, const char *name)
  {
    if (array == nullptr) return create(array, n, name);

    bigint nbytes = ((bigint) sizeof(TYPE)) * n;
    array = (TYPE *) srealloc(array, nbytes, name);
    return array;
  }

  template <typename TYPE>
  TYPE **grow(TYPE **& /* array */, int /* n */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 1d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy(TYPE *&array) {
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 1d array with index from nlo to nhi inclusive
   cannot grow it
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE *create1d_offset(TYPE *&array, int nlo, int nhi, const char *name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * (nhi-nlo+1);
    array = (TYPE *) smalloc(nbytes, name);
    array -= nlo;
    return array;
  }

  template <typename TYPE>
  TYPE **create1d_offset(TYPE **& /* array */, int /* nlo */, int /* nhi */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 1d array with index offset
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy1d_offset(TYPE *&array, int offset)
  {
    if (array) sfree(&array[offset]);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 2d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE **create(TYPE **&array, int n1, int n2, const char *name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2;
    TYPE *data = (TYPE *) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes, name);

    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  template <typename TYPE>
  TYPE ***create(TYPE ***& /* array */, int /* n1 */, int /* n2 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE **grow(TYPE **&array, int n1, int n2, const char *name)
  {
    if (array == nullptr) return create(array, n1, n2, name);

    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2;
    TYPE *data = (TYPE *) srealloc(array[0], nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) srealloc(array, nbytes, name);

    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  template <typename TYPE>
  TYPE ***grow(TYPE ***& /* array */, int /* n1 */, int /* n2 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 2d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy(TYPE **&array)
  {
    if (array == nullptr) return;
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 2d array with a ragged 2nd dimension
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE **create_ragged(TYPE **&array, int n1, int *n2, const char *name)
  {
    bigint n2sum = 0;
    for (int i = 0; i < n1; i++) n2sum += n2[i];

    bigint nbytes = ((bigint) sizeof(TYPE)) * n2sum;
    TYPE *data = (TYPE *) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes, name);

    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2[i];
    }
    return array;
  }

  template <typename TYPE>
  TYPE ***create_ragged(TYPE ***& /* array */, int /* n1 */, int */* n2 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   create a 2d array with 2nd index from n2lo to n2hi inclusive
   cannot grow it
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE **create2d_offset(TYPE **&array, int n1, int n2lo, int n2hi, 
                         const char *name)
  {
    int n2 = n2hi - n2lo + 1;
    create(array, n1, n2, name);
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    return array;
  }

  template <typename TYPE>
  TYPE ***create2d_offset(TYPE ***& /* array */, int /* n1 */, int /* n2lo */, int /* n2hi */, 
                          const char *name) {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 2d array with 2nd index offset
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy2d_offset(TYPE **&array, int offset)
  {
    if (array == nullptr) return;
    sfree(&array[0][offset]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 3d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3;
    TYPE *data = (TYPE *) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1 * n2;
    TYPE **plane = (TYPE **) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1;
    array = (TYPE ***) smalloc(nbytes, name);

    int i, j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }

  template <typename TYPE>
  TYPE ****create(TYPE ****& /* array */, int /* n1 */, int /* n2 */, int /* n3 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name)
  {
    if (array == nullptr) return create(array, n1, n2, n3, name);

    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3;
    TYPE *data = (TYPE *) srealloc(array[0][0], nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1 * n2;
    TYPE **plane = (TYPE **) srealloc(array[0], nbytes, name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1;
    array = (TYPE ***) srealloc(array, nbytes, name);

    int i, j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }

  template <typename TYPE>
  TYPE ****grow(TYPE ****& /* array */, int /* n1 */, int /* n2 */, int /* n3 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 3d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy(TYPE ***&array)
  {
    if (array == nullptr) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 3d array with 1st index from n1lo to n1hi inclusive
   cannot grow it
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, 
                          int n2, int n3, const char *name)
  {
    int n1 = n1hi - n1lo + 1;
    create(array, n1, n2, n3, name);
    array -= n1lo;
    return array;
  }

  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****& /* array */, int /* n1lo */, int /* n1hi */, 
                           int /* n2 */, int /* n3 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   free a 3d array with 1st index offset
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy3d_offset(TYPE ***&array, int offset)
  {
    if (array == nullptr) return;
    sfree(&array[offset][0][0]);
    sfree(&array[offset][0]);
    sfree(&array[offset]);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 3d array with
   1st index from n1lo to n1hi inclusive, 
   2nd index from n2lo to n2hi inclusive, 
   3rd index from n3lo to n3hi inclusive
   cannot grow it
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, 
                          int n2lo, int n2hi, int n3lo, int n3hi, 
                          const char *name)
  {
    int n1 = n1hi - n1lo + 1;
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    create(array, n1, n2, n3, name);

    bigint m = ((bigint) n1) * n2;
    for (bigint i = 0; i < m; i++) array[0][i] -= n3lo;
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    array -= n1lo;
    return array;
  }

  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****& /* array */, int /* n1lo */, int /* n1hi */, 
                           int /* n2lo */, int /* n2hi */, int /* n3lo */, int /* n3hi */, 
                           const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   free a 3d array with all 3 indices offset
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy3d_offset(TYPE ***&array, 
                        int n1_offset, int n2_offset, int n3_offset)
  {
    if (array == nullptr) return;
    sfree(&array[n1_offset][n2_offset][n3_offset]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(&array[n1_offset]);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 4d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4, 
                  const char *name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3 * n4;
    TYPE *data = (TYPE *) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1 * n2 * n3;
    TYPE **cube = (TYPE **) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1 * n2;
    TYPE ***plane = (TYPE ***) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE ***)) * n1;
    array = (TYPE ****) smalloc(nbytes, name);

    int i, j, k;
    bigint m1, m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2 * n3 + j * n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2 * n3 + j * n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }

  template <typename TYPE>
  TYPE *****create(TYPE *****& /* array */, int /* n1 */, int /* n2 */, int /* n3 */, int /* n4 */, 
                   const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
  grow or shrink 1st dim of a 4d array
  last 3 dims must stay the same
  -------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, int n4, 
          const char *name)
  {
          if (array == nullptr) return create(array, n1, n2, n3, n4, name);

          bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3 * n4;
          TYPE *data = (TYPE *)srealloc(array[0][0][0], nbytes, name);
          nbytes = ((bigint) sizeof(TYPE *)) * n1 * n2 * n3;
          TYPE **cube = (TYPE **)srealloc(array[0][0], nbytes, name);
          nbytes = ((bigint) sizeof(TYPE **)) * n1 * n2;
          TYPE ***plane = (TYPE ***)srealloc(array[0], nbytes, name);
          nbytes = ((bigint) sizeof(TYPE ***)) * n1;
          array = (TYPE ****)srealloc(array, nbytes, name);

          int i, j, k;
          bigint m1, m2;
          bigint n = 0;
          for (i = 0; i < n1; i++) {
                  m2 = ((bigint)i) * n2;
                  array[i] = &plane[m2];
                  for (j = 0; j < n2; j++) {
                          m1 = ((bigint)i) * n2 + j;
                          m2 = ((bigint)i) * n2 * n3 + j * n3;
                          plane[m1] = &cube[m2];
                          for (k = 0; k < n3; k++) {
                                  m1 = ((bigint)i) * n2 * n3 + j * n3 + k;
                                  cube[m1] = &data[n];
                                  n += n4;
                          }
                  }
          }
          return array;
  }

  template <typename TYPE>
  TYPE *****grow(TYPE *****& /* array */, int /* n1 */, int /* n2 */, int /* n3 */, int /* n4 */, 
          const char *name)
  {
          fail(name); return nullptr;
  }

/*  ----------------------------------------------------------------------
   destroy a 4d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy(TYPE ****&array)
  {
    if (array == nullptr) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 4d array with indices
   2nd index from n2lo to n2hi inclusive
   3rd index from n3lo to n3hi inclusive
   4th index from n4lo to n4hi inclusive
   cannot grow it
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE ****create4d_offset(TYPE ****&array, int n1, int n2lo, int n2hi, 
                           int n3lo, int n3hi, int n4lo, int n4hi, 
                           const char *name)
  {
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    int n4 = n4hi - n4lo + 1;
    create(array, n1, n2, n3, n4, name);

    bigint m = ((bigint) n1) * n2 * n3;
    for (bigint i = 0; i < m; i++) array[0][0][i] -= n4lo;
    m = ((bigint) n1) * n2;
    for (bigint i = 0; i < m; i++) array [0][i] -= n3lo;
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    return array;
  }

  template <typename TYPE>
  TYPE ****create4d_offset(TYPE *****& /* array */, int /* n1 */, int /* n2lo */, int /* n2hi */, 
                           int /* n3lo */, int /* n3hi */, int /* n4lo */, int /* n4hi */, 
                           const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   free a 4d array with indices 2, 3, and 4 offset
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy4d_offset(TYPE ****&array, 
                        int n2_offset, int n3_offset, int n4_offset)
  {
    if (array == nullptr) return;
    sfree(&array[0][n2_offset][n3_offset][n4_offset]);
    sfree(&array[0][n2_offset][n3_offset]);
    sfree(&array[0][n2_offset]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   create a 5d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4, 
                   int n5, const char *name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3 * n4 * n5;
    TYPE *data = (TYPE *) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1 * n2 * n3 * n4;
    TYPE **level4 = (TYPE **) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1 * n2 * n3;
    TYPE ***level3 = (TYPE ***) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE ***)) * n1 * n2;
    TYPE ****level2 = (TYPE ****) smalloc(nbytes, name);
    nbytes = ((bigint) sizeof(TYPE ****)) * n1;
    array = (TYPE *****) smalloc(nbytes, name);

    int i, j, k, l;
    bigint m1, m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &level2[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2 * n3 +  ((bigint) j) * n3;
        level2[m1] = &level3[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2 * n3 +  ((bigint) j) * n3 + k;
          m2 = ((bigint) i) * n2 * n3 * n4 +
            ((bigint) j) * n3 * n4 + ((bigint) k) * n4;
          level3[m1] = &level4[m2];
          for (l = 0; l < n4; l++) {
            m1 = ((bigint) i) * n2 * n3 * n4 +
              ((bigint) j) * n3 * n4 + ((bigint) k) * n4 + l;
            level4[m1] = &data[n];
            n += n5;
          }
        }
      }
    }
    return array;
  }

  template <typename TYPE>
  TYPE ******create(TYPE ******& /* array */, int /* n1 */, int /* n2 */, int /* n3 */, int /* n4 */, 
                    int /* n5 */, const char *name)
  {fail(name); return nullptr;}

/*  ----------------------------------------------------------------------
   destroy a 5d array
-------------------------------------------------------------------------  */

  template <typename TYPE>
  void destroy(TYPE *****&array)
  {
    if (array == nullptr) return;
    sfree(array[0][0][0][0]);
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }

/*  ----------------------------------------------------------------------
   memory usage of arrays, including pointers
-------------------------------------------------------------------------  */

  template <typename TYPE>
  bigint usage(TYPE *array, int n)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n;
    return bytes;
  }

  template <typename TYPE>
  bigint usage(TYPE **array, int n1, int n2)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1 * n2;
    bytes += ((bigint) sizeof(TYPE *)) * n1;
    return bytes;
  }

  template <typename TYPE>
  bigint usage(TYPE ***array, int n1, int n2, int n3)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3;
    bytes += ((bigint) sizeof(TYPE *)) * n1 * n2;
    bytes += ((bigint) sizeof(TYPE **)) * n1;
    return bytes;
  }

  template <typename TYPE>
  bigint usage(TYPE ****array, int n1, int n2, int n3, int n4)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1 * n2 * n3 * n4;
    bytes += ((bigint) sizeof(TYPE *)) * n1 * n2 * n3;
    bytes += ((bigint) sizeof(TYPE **)) * n1 * n2;
    bytes += ((bigint) sizeof(TYPE ***)) * n1;
    return bytes;
  }

  size_t estimate_unordered_vector2d_map_memory(const std::unordered_map<Vector2D, int, Vector2DHasher> &map) {
    size_t total_memory_bytes = 0;

    // Estimate the size of a single map node (Key + Value + Node Overhead)
    const size_t SIZE_OF_POINTER = sizeof(void*); // 8 bytes on 64-bit
    const size_t KEY_SIZE = sizeof(Vector2D);    // 8 bytes (2 ints stored inline)
    const size_t VALUE_SIZE = sizeof(int);      // 4 bytes
    
    // Node Overhead: Includes the 'next' pointer for the linked list and padding.
    // A common estimate is Pointer (8 bytes) + ~8 bytes for padding/bookkeeping.
    const size_t NODE_OVERHEAD = SIZE_OF_POINTER + SIZE_OF_POINTER; // Conservative 16 bytes overhead

    // Total memory required for one element node in the map.
    const size_t TOTAL_NODE_SIZE = KEY_SIZE + VALUE_SIZE + NODE_OVERHEAD;
    
    // 1. Estimate memory for the bucket array (hash table)
    // Each bucket entry is a pointer (or node handle).
    total_memory_bytes += map.bucket_count() * SIZE_OF_POINTER;

    // 2. Estimate memory for all stored key-value pairs (nodes)
    // Since Vector2D stores its data inline, we calculate the total memory 
    // by multiplying the fixed node size by the element count.
    total_memory_bytes += map.size() * TOTAL_NODE_SIZE;

    // 3. Add small, fixed overhead for the map object itself (e.g., size, load factor info).
    // Typically small (e.g., 24-32 bytes) and often neglected in large estimations.
    // total_memory_bytes += sizeof(std::unordered_map<...>); 

    return total_memory_bytes;

    /*
    // 1. Estimate memory for the bucket array (the hash table itself).
    // Each bucket is roughly the size of a pointer.
    total_memory_bytes += map.bucket_count() * sizeof(void*);

    // 2. Estimate memory for each stored key-value pair (node).
    for (const auto &pair : map) {
      // --- Memory for the key (std::vector<int>) ---
      // a) The dynamically allocated data array inside the vector.
      total_memory_bytes += pair.first.capacity() * sizeof(int); 

      // b) The vector object itself (control block: pointer, size, capacity).
      total_memory_bytes += sizeof(pair.first); 

      // --- Memory for the value (int) ---
      total_memory_bytes += sizeof(pair.second);

      // --- Overhead for the map node ---
      // (e.g., the 'next' pointer in the bucket's linked list for collisions)
      total_memory_bytes += sizeof(void*); 
    }

    return total_memory_bytes;
    */
  }
};

}

#endif


