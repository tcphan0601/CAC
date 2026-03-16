#ifdef PAIR_CLASS

PairStyle(table,PairTable)

#else

#ifndef CAC_PAIR_TABLE_H
#define CAC_PAIR_TABLE_H

#include "pair.h"

namespace CAC_NS {

class PairTable : public Pair {
 public:
  PairTable(class CAC *);
  virtual ~PairTable();

  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  //void write_restart(FILE *);
  //void read_restart(FILE *);
  //void write_restart_settings(FILE *);
  //void read_restart_settings(FILE *);
  //virtual double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  enum{LOOKUP,LINEAR,SPLINE,BITMAP};

 protected:
  int tabstyle,tablength;
  struct Table {
    int ninput,rflag,fpflag,match,ntablebits;
    int nshiftbits,nmask;
    double rlo,rhi,fplo,fphi,cut;
    double *rfile,*efile,*ffile;
    double *e2file,*f2file;
    double innersq,delta,invdelta,deltasq6;
    double *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
  };
  int ntables;
  Table *tables;

  int **tabindex;

  virtual void allocate();
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  virtual void compute_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  static void spline(double *, double *, int, double, double, double *);
  static double splint(double *, double *, double *, int, double);
};

}

#endif
#endif

