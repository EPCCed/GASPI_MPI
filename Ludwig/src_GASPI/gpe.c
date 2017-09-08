/*****************************************************************************
 *
 *  gpe.c
 *
 *  The GASPI parallel environment.
 *
 *  This is responsible for initialisation and finalisation of
 *  the GASPI parallel environment.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Luis Cebamanos (l.cebamanos@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "gpe.h"


struct gpe_s{
  gaspi_rank_t gaspi_rank;           /* Rank */
  gaspi_rank_t gaspi_size;           /* Size */
};

static gpe_t * gpe = NULL;

static gaspi_segment_id_t gaspi_segment_id = 0;
static gaspi_pointer_t _vptr;
/*****************************************************************************
 *
 *  gpe_init
 *
 *  Initialise the GASPI model.
 *****************************************************************************/

void gpe_init(void) {

  /* GASPI create comm */
  gpe_create();
  
  /* GASPI get rank values */
  GASPIERROR(gaspi_proc_rank(&gpe->gaspi_rank));
  /* GASPI get size values */
  GASPIERROR(gaspi_proc_num(&gpe->gaspi_size));
 
  info("Welcome to Ludwig v%d.%d.%d (%s version running on %d process%s)\n\n",
       LUDWIG_MAJOR_VERSION, LUDWIG_MINOR_VERSION, LUDWIG_PATCH_VERSION,
       (gpe->gaspi_size > 1) ? "GASPI" : "Serial", gpe->gaspi_size,
       (gpe->gaspi_size == 1) ? "" : "es");

}
/*****************************************************************************
 *
 *  gpe_alloc_segment
 *
 *  Creates a GASPI segment of a given size and zero it. At the moment it creates
 *  only a double size array. If we require other types, it may be better of
 *  just pass the bytes size as parameters.
 *  Return: id segment that has just been created
 *****************************************************************************/

void gpe_alloc_segment(int nx, int ny, int nz, int data) {
  int i;
  // Size is multiply by 2 becase we want to receive in the same buffer
  // therefore we need twice as much space
  int XZsize = 2*nx*nz*data;
  int YZsize = 2*ny*nz*data;
  int XYsize = 2*nx*ny*data;

  // We create 2 segements per direction; one to the left and one the right
  for (i=0; i<2;i++){
    
    gaspi_segment_id=i;
    // create segments
    GASPIERROR(gaspi_segment_create(gaspi_segment_id, XZsize*sizeof(double),
				    GASPI_GROUP_ALL, GASPI_BLOCK,
				    GASPI_ALLOC_DEFAULT));
    //GASPI_MEM_INITIALIZED));
    GASPIERROR(gaspi_segment_create(gaspi_segment_id+2, YZsize*sizeof(double),
				    GASPI_GROUP_ALL, GASPI_BLOCK,
				    GASPI_ALLOC_DEFAULT));
    //GASPI_MEM_INITIALIZED));
    GASPIERROR(gaspi_segment_create(gaspi_segment_id+4, XYsize*sizeof(double),
				    GASPI_GROUP_ALL, GASPI_BLOCK,
				    GASPI_ALLOC_DEFAULT));
    //    GASPI_MEM_INITIALIZED));
  }

}
/*****************************************************************************
 *
 *  gpe_get_gaspi_pointer
 *
 *  Creates a GASPI segment of a given size and zero it. At the moment it creates
 *  only a double size array. If we require other types, it may be better of
 *  just pass the bytes size as parameters.
 *  Return: pointer to id segment
 *****************************************************************************/
double* gpe_get_gaspi_pointer(){
  double* mem = (double*)_vptr;
  return mem;
}

/*****************************************************************************
 *
 *  gpe_num_segment
 *
 *  Return the number of GASPI segments that have been created until now
 *****************************************************************************/
int gpe_num_segment() {
  
  return (int) gaspi_segment_id;	
  
}

/*****************************************************************************
 *
 *  gpe_create
 *
 *  Allocate memory for the GASPI communicator
 *****************************************************************************/
void gpe_create(){
  gpe = (gpe_t *) calloc(1, sizeof(gpe_t));
  if (gpe == NULL){
    printf("calloc(gpe_t) failed\n");
    exit(0);
  }
}

/*****************************************************************************
 *
 *  gpe_rank
 *
 *  Returns a GASPI rank
 *****************************************************************************/

int gpe_rank(){
  return gpe->gaspi_rank;
}
/*****************************************************************************
 *
 *  gpe_size
 *
 *  Returns GASPI proc size
 *****************************************************************************/

int gpe_size(){
  return gpe->gaspi_size;
}
