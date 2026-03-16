#ifndef CAC_MY_PAGE_H
#define CAC_MY_PAGE_H

#include <stdlib.h>
namespace CAC_NS {

template<class T>
class MyPage {
 public:
  int ndatum;      // total # of stored datums
  int nchunk;      // total # of stored chunks

  MyPage() {
    ndatum = nchunk = 0;
    pages = NULL;
    npage = 0;
    errorflag = 0;
  }

  // (re)initialize allocation params
  // also allocate first page(s)

  int init(int user_maxchunk = 1, int user_pagesize = 1024, 
           int user_pagedelta = 1) {
    maxchunk = user_maxchunk;
    pagesize = user_pagesize;
    pagedelta = user_pagedelta;

    if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0) return 1;
    if (maxchunk > pagesize) return 1;

    // free any previously allocated pages

    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);

    // initial page allocation

    ndatum = nchunk = 0;
    pages = NULL;
    npage = 0;
    allocate();
    if (errorflag) return 2;
    ipage = index = 0;
    page = pages[ipage];
    return 0;
  }

  // free all allocated pages

  ~MyPage() {
    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);
  }

  // get ptr to one datum
  // return NULL if run out of memory

  T *get() {
    ndatum++;
    nchunk++;
    if (index < pagesize) return &page[index++];
    ipage++;
    if (ipage == npage) {
      allocate();
      if (errorflag) return NULL;
    }
    page = pages[ipage];
    index = 0;
    return &page[index++];
  }

  // get ptr to location that can store N datums
  // error if N > maxchunk
  // return NULL if run out of memory

  T *get(int n) {
    if (n > maxchunk) {
      errorflag = 1;
      return NULL;
    }
    ndatum += n;
    nchunk++;
    if (index+n <= pagesize) {
      int start = index;
      index += n;
      return &page[start];
    }
    ipage++;
    if (ipage == npage) {
      allocate();
      if (errorflag) return NULL;
    }
    page = pages[ipage];
    index = n;
    return &page[0];
  }

  // get ptr to location that can store maxchunk datums
  // will return same ptr as previous call if vgot() not called
  // return NULL if run out of memory
  
  T *vget() {
    if (index+maxchunk <= pagesize) return &page[index];
    ipage++;
    if (ipage == npage) {
      allocate();
      if (errorflag) return NULL;
    }
    page = pages[ipage];
    index = 0;
    return &page[index];
  }

  // increment by N = # of values stored in loc returned by vget()
  // OK to not call if vget() ptr was not used
  // error if N > maxchunk

  void vgot(int n) {
    if (n > maxchunk) errorflag = 1;
    ndatum += n;
    nchunk++;
    index += n;
  }

  // clear all pages, without freeing any memory

  void reset() {
    ndatum = nchunk = 0;
    index = ipage = 0;
    page = pages[ipage];
  }

  // return total size of allocated pages

  int size() const {
    return npage*pagesize*sizeof(T);
  }

  // return error status

  int status() const {
    return errorflag;
  }

 private:
  T **pages;      // list of allocated pages
  T *page;        // ptr to current page
  int npage;      // # of allocated pages
  int ipage;      // index of current page
  int index;      // current index on current page
  
  int maxchunk;   // max # of datums in one requested chunk
  int pagesize;   // # of datums in one page, default = 1024
  int pagedelta;  // # of pages to allocate at once, default = 1

  int errorflag;  // flag > 0 if error has occurred
                  // 1 = chunk size exceeded maxchunk
                  // 2 = memory allocation error

  void allocate() {
    npage += pagedelta;
    pages = (T **) realloc(pages,npage*sizeof(T *));
    if (!pages) {
      errorflag = 2;
      return;
    }

    for (int i = npage-pagedelta; i < npage; i++) {
#if defined(CAC_MEMALIGN)
      void *ptr;
      if (posix_memalign(&ptr, CAC_MEMALIGN, pagesize*sizeof(T)))
        errorflag = 2;
      pages[i] = (T *) ptr;
#else
      pages[i] = (T *) malloc(pagesize*sizeof(T));
      if (!pages[i]) errorflag = 2;
#endif
    }
  }
};

}

#endif
