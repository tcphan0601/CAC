#include <mpi.h>
#include "cac.h"
#include "input.h"
#include "error.h"
#include <stdio.h>
#include <stdlib.h>

using namespace CAC_NS;

/*  ----------------------------------------------------------------------
   main program to drive CAC
-------------------------------------------------------------------------  */

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  CAC *cac = new CAC(argc, argv, MPI_COMM_WORLD);
  cac->input->file();
  
  delete cac;

  int me;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize(); 
}
