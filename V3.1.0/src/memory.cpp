#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"
#include "error.h"

using namespace CAC_NS;

/*  ----------------------------------------------------------------------  */

Memory::Memory(CAC *cac) : Pointers(cac) {}

/*  ----------------------------------------------------------------------
   safe malloc
-------------------------------------------------------------------------  */

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return nullptr;

#if defined(CAC_MEMALIGN)
  void *ptr;
  int retval = posix_memalign(&ptr, CAC_MEMALIGN, nbytes);
  if (retval) ptr = nullptr;
#else
  void *ptr = malloc(nbytes);
#endif
  if (ptr == nullptr) {
    char str[128];
    sprintf(str, "Failed to allocate " BIGINT_FORMAT " bytes for array %s", 
            nbytes, name);
    error->one(FLERR, str);
  }
  return ptr;
}

/*  ----------------------------------------------------------------------
   safe realloc
-------------------------------------------------------------------------  */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return nullptr;
  }

  ptr = realloc(ptr, nbytes);
  if (ptr == nullptr) {
    char str[128];
    sprintf(str, "Failed to reallocate " BIGINT_FORMAT " bytes for array %s", 
            nbytes, name);
   error->one(FLERR, str);
  }
  return ptr;
}

/*  ----------------------------------------------------------------------
   safe free
-------------------------------------------------------------------------  */

void Memory::sfree(void *ptr)
{
  if (ptr == nullptr) return;
  free(ptr);
}

/*  ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
-------------------------------------------------------------------------  */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str, "Cannot create/grow a vector/array of pointers for %s", name);
  error->one(FLERR, str);
}
