/*****************************************************************************
 *
 *  main.c
 *
 *  Main driver code. See ludwig.c for details of timestepping etc.
 *
 *  $Id$
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#include <stdio.h>
#include "gaspi_utils.h"

#include "pe.h"
#include "ludwig.h"
#ifdef PETSC
  #include "petscksp.h"
#endif

/*****************************************************************************
 *
 *  main
 *
 *****************************************************************************/

int main(int argc, char ** argv) {

  char inputfile[FILENAME_MAX] = "input";

  gaspi_rank_t grank, gsize;

  /* MPI Init */
  MPI_Init(&argc, &argv);
  /* GASPI Init*/
  GASPIERROR(gaspi_proc_init(GASPI_BLOCK));

#ifdef PETSC
  PetscInitialize(&argc, &argv, (char*) 0, NULL); 
#endif 
  if (argc > 1) sprintf(inputfile, "%s", argv[1]);

  ludwig_run(inputfile);

#ifdef PETSC
  PetscFinalize();
#endif
  /* MPI Finalize*/
  MPI_Finalize();

  /* GASPI Finalize*/
  GASPIERROR(gaspi_proc_term(GASPI_BLOCK));


  return 0;
}
