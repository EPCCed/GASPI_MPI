/*****************************************************************************
 *
 *  gpe.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *  Luis Cebamanos (l.cebamanos@epcc.ed.ac.uk)
 *
 *****************************************************************************/
#ifndef GPE_H
#define GPE_H

#include "../version.h"
#include "gaspi_utils.h"
#include "targetDP.h"
#include "pe.h"
typedef struct gpe_s gpe_t;

__targetHost__ void gpe_init(void);
__targetHost__ void gpe_create(void);
__targetHost__ int  gpe_rank(void);
__targetHost__ int  gpe_size(void);
__targetHost__ void  gpe_alloc_segment(int,int,int,int);
__targetHost__ int  gpe_num_segment(void);
__targetHost__ double* gpe_get_gaspi_pointer(void);

#endif
