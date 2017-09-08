#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "gpe.h"
#include "coords.h"
#include "model.h"
#include "lb_model_s.h"
#include "targetDP.h"
#include "io_harness.h"
#include "control.h"
#include "gaspi_utils.h"

/*****************************************************************************
 *
 *  lb_halo_via_copy_nonblocking
 *
 *  A version of the halo swap which uses a flat buffer to copy the
 *  relevant data rathar than MPI data types. This version is NON-
 *  BLOCKING and sends 26 messages instead of 6: 6 planes, 12 edges
 *  and 8 corners.
 *
 *  It works for both MODEL and MODEL_R, but the loop order favours
 *  MODEL_R.
 *
 *****************************************************************************/

int lb_halo_via_copy_nonblocking_start_GASPI(lb_t * lb) {

  /* Send messages to 6 neighbouring planes
     12 neighbouring edges
     8 neighbouring corners */
  halo_planes_GASPI(lb);
  halo_edges_GASPI(lb);
  halo_corners_GASPI(lb);

  return 0;
}

int lb_halo_via_copy_nonblocking_end_GASPI(lb_t * lb) {

  int n;
  int recvcount;
  int sendcount;

  /* for (n=0; n<6; n++){ */
  /*   MPI_Waitany(6, lb->hl.recvreqP, &recvcount, lb->hl.statusP); */
  /* } */

  /* for (n=0; n<2; n++){ */
  /*   MPI_Waitany(2, lb->hl.recvreqP2, &recvcount, lb->hl.statusP2); */
  /* } */

  for (n=0; n<12; n++){
    MPI_Waitany(12, lb->hl.recvreqE, &recvcount, lb->hl.statusE);
  }

  for (n=0; n<8; n++){
    MPI_Waitany(8, lb->hl.recvreqC, &recvcount, lb->hl.statusC);
  }

  /* Copy data from MPI buffers */
  unpack_halo_buffers_GASPI(lb);

  /* for (n=0; n<6; n++){ */
  /*   MPI_Waitany(6, lb->hl.sendreqP, &sendcount, lb->hl.statusP); */
  /* } */
  /* for (n=0; n<2; n++){ */
  /*   MPI_Waitany(2, lb->hl.sendreqP2, &sendcount, lb->hl.statusP2); */
  /* } */

  for (n=0; n<12; n++){
    MPI_Waitany(12, lb->hl.sendreqE, &sendcount, lb->hl.statusE);
  }

  for (n=0; n<8; n++){
    MPI_Waitany(8, lb->hl.sendreqC, &sendcount, lb->hl.statusC);
  }


  return 0;
}



/*****************************************************************************
 *
 *  unpack_halo_buffers
 *
 *  Unpacks buffers from MPI messages into Halo values for next step of simulation
 *
 *****************************************************************************/

void unpack_halo_buffers_GASPI(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo;
  int count;
  int nlocal[3];

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  int nsendXZ, nsendYZ, nsendXY;
  nsendYZ = NVEL*lb->ndist*nlocal[Y]*nlocal[Z];
  nsendXZ = NVEL*lb->ndist*nlocal[X]*nlocal[Z];
  nsendXY = NVEL*lb->ndist*nlocal[X]*nlocal[Y];

  int nsendX, nsendY, nsendZ;
  nsendX = NVEL*lb->ndist*nlocal[X];
  nsendY = NVEL*lb->ndist*nlocal[Y];
  nsendZ = NVEL*lb->ndist*nlocal[Z];

  int nsend;
  nsend = NVEL*lb->ndist*1;

  /**/
  int YZ_size = nsendYZ;
  const gaspi_segment_id_t seg_id_YZ_L = 0;
  const gaspi_segment_id_t seg_id_YZ_R = 1;

  int offset = YZ_size*sizeof(double);
  gaspi_pointer_t gptr_YZ_L, gptr_YZ_R;

  /* pointer to the right */
  GASPIERROR(gaspi_segment_ptr(seg_id_YZ_L, &gptr_YZ_L));
  double* ptr_YZ_L =  (double*)gptr_YZ_L;
  
  GASPIERROR(gaspi_segment_ptr(seg_id_YZ_R, &gptr_YZ_R));
  double* ptr_YZ_R =  (double*)gptr_YZ_R;
  

  double* recv  = ptr_YZ_L+YZ_size;
  double* recv2 = ptr_YZ_R+YZ_size;  
  /**/

  /* Unpack Planes */

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(0, jc, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvbackYZ[count];
	  fptr[indexhalo] = recv2[count];

          index = coords_index(nlocal[X] + 1, jc, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvforwYZ[count];
	  fptr[indexhalo] = recv[count];
	  
          ++count;
        }
      }
    }
  }
  assert(count == nsendYZ);
  
  int XZ_size =  nsendXZ;
  offset = XZ_size*sizeof(double);    
  
  /* segment ids */
  const gaspi_segment_id_t seg_id_XZ_L=2;
  /* segment ids */
  const gaspi_segment_id_t seg_id_XZ_R=3;
  
  /* GASPI pointers to the GASPI segment */
  gaspi_pointer_t gptr_XZ_L, gptr_XZ_R;
  
  /* pointer to the right */
  GASPIERROR(gaspi_segment_ptr(seg_id_XZ_L, &gptr_XZ_L));
  double* ptr_XZ_L =  (double*)gptr_XZ_L;
  
  GASPIERROR(gaspi_segment_ptr(seg_id_XZ_R, &gptr_XZ_R));
  double* ptr_XZ_R =  (double*)gptr_XZ_R;
  
  recv   = ptr_XZ_L + XZ_size;
  recv2  = ptr_XZ_R + XZ_size;
  

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(ic, 0, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvbackXZ[count];
	  fptr[indexhalo] = recv2[count];
	  
          index = coords_index(ic, nlocal[Y] + 1, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvforwXZ[count];
	  fptr[indexhalo] = recv[count];
          ++count;
        }
      }
    }
  }
  assert(count == nsendXZ);

  int XY_size =  nsendXY;
  offset = XY_size*sizeof(double);

  /* segment ids */
  const gaspi_segment_id_t seg_id_XY_L=4;
  /* segment ids */
  const gaspi_segment_id_t seg_id_XY_R=5;

  /* GASPI pointers to the GASPI segment */
  gaspi_pointer_t gptr_XY_L, gptr_XY_R;

  /* pointer to the right */
  GASPIERROR(gaspi_segment_ptr(seg_id_XY_L, &gptr_XY_L));
  double* ptr_XY_L =  (double*)gptr_XY_L;

  GASPIERROR(gaspi_segment_ptr(seg_id_XY_R, &gptr_XY_R));
  double* ptr_XY_R =  (double*)gptr_XY_R;
 
  recv   = ptr_XY_L + XY_size;
  recv2  = ptr_XY_R + XY_size;

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {

          index = coords_index(ic, jc, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvbackXY[count];
	  fptr[indexhalo] = recv2[count];

          index = coords_index(ic, jc, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //fptr[indexhalo] = lb->hl.recvforwXY[count];
	  fptr[indexhalo] = recv[count];

          ++count;
        }
      }
    }
  }
  assert(count == nsendXY);

  /* Free memory for planes buffers */
  free(lb->hl .sendforwYZ);
  free(lb->hl .recvforwYZ);
  free(lb->hl .sendbackYZ);
  free(lb->hl .recvbackYZ);
  free(lb->hl .sendforwXZ);
  free(lb->hl .recvforwXZ);
  free(lb->hl .sendbackXZ);
  free(lb->hl .recvbackXZ);
  free(lb->hl .sendforwXY);
  free(lb->hl .recvforwXY);
  free(lb->hl .sendbackXY);
  free(lb->hl .recvbackXY);

  /* Unpack Edges */

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {

        index = coords_index(ic, nlocal[Y] + 1, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXpp[count];

        index = coords_index(ic, nlocal[Y] + 1, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXpn[count];

        index = coords_index(ic, 0, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXnp[count];

        index = coords_index(ic, 0, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXnn[count];
        ++count;
      }
    }
  }
  assert(count == nsendX);

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {

        index = coords_index(nlocal[X] + 1, jc, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYpp[count];

        index = coords_index(nlocal[X] + 1, jc, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYpn[count];

        index = coords_index(0, jc, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYnp[count];

        index = coords_index(0, jc, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYnn[count];
        ++count;
      }
    }
  }
  assert(count == nsendY);

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZpp[count];

        index = coords_index(nlocal[X] + 1, 0, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZpn[count];

        index = coords_index(0, nlocal[Y] + 1, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZnp[count];

        index = coords_index(0, 0, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZnn[count];
        ++count;
      }
    }
  }
  assert(count == nsendZ);

  /* Free memory for planes buffers */
  free(lb->hl .sendXnn);
  free(lb->hl .recvXnn);
  free(lb->hl .sendXnp);
  free(lb->hl .recvXnp);
  free(lb->hl .sendXpn);
  free(lb->hl .recvXpn);
  free(lb->hl .sendXpp);
  free(lb->hl .recvXpp);
  free(lb->hl .sendYnn);
  free(lb->hl .recvYnn);
  free(lb->hl .sendYnp);
  free(lb->hl .recvYnp);

  free(lb->hl .sendYpn);
  free(lb->hl .recvYpn);
  free(lb->hl .sendYpp);
  free(lb->hl .recvYpp);
  free(lb->hl .sendZnn);
  free(lb->hl .recvZnn);
  free(lb->hl .sendZnp);
  free(lb->hl .recvZnp);
  free(lb->hl .sendZpn);
  free(lb->hl .recvZpn);
  free(lb->hl .sendZpp);
  free(lb->hl .recvZpp);

  /* Unpack corners buffers */

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {

      index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvppp[count];

      index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvppn[count];

      index = coords_index(nlocal[X] + 1, 0, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvpnp[count];

      index = coords_index(nlocal[X] + 1, 0, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvpnn[count];

      index = coords_index(0, nlocal[Y] + 1, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnpp[count];

      index = coords_index(0, nlocal[Y] + 1, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnpn[count];

      index = coords_index(0, 0, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnnp[count];

      index = coords_index(0, 0, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnnn[count];
      count++;
    }
  }

  /* Free memory for corner buffers */
  free(lb->hl .sendnnn);
  free(lb->hl .sendnnp);
  free(lb->hl .sendnpn);
  free(lb->hl .sendnpp);
  free(lb->hl .sendpnn);
  free(lb->hl .sendpnp);
  free(lb->hl .sendppn);
  free(lb->hl .sendppp);

  free(lb->hl .recvnnn);
  free(lb->hl .recvnnp);
  free(lb->hl .recvnpn);
  free(lb->hl .recvnpp);
  free(lb->hl .recvpnn);
  free(lb->hl .recvpnp);
  free(lb->hl .recvppn);
  free(lb->hl .recvppp);

  return;
}
