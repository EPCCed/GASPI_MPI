#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gpe.h"
#include "gaspi_utils.h"
#include "pe.h"
#include "coords.h"
#include "model.h"
#include "lb_model_s.h"
#include "targetDP.h"
#include "io_harness.h"
#include "control.h"

/*****************************************************************************
 *
 *  halo_planes
 *
 *  Sends 6 GASPI messages (the planes) for the Non-Blocking version.
 *
 *****************************************************************************/

void halo_planes_GASPI(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagf = 900;
  const int tagb = 901;

  /* The ranks of neighbouring planes */
  int pforwX, pbackX, pforwY, pbackY, pforwZ, pbackZ;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  /* Allocate size of sendplanes and number of elements send with each plane */
  int nsendXZ, nsendYZ, nsendXY;
  nsendYZ = NVEL*lb->ndist*nlocal[Y]*nlocal[Z];
  nsendXZ = NVEL*lb->ndist*nlocal[X]*nlocal[Z];
  nsendXY = NVEL*lb->ndist*nlocal[X]*nlocal[Y];

  /* Allocate message sizes for plane send/receives */
  lb->hl .sendforwYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .sendbackYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .recvforwYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .recvbackYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .sendforwXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .sendbackXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .recvforwXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .recvbackXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .sendforwXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .sendbackXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .recvforwXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .recvbackXY = (double *) malloc(nsendXY*sizeof(double));

  /* Receive planes in the X-direction */
  /* PPM, NMM, P=Positive, M=Middle, N=Negative for the XYZ directions respectively */
  pforwX = nonblocking_cart_neighb(PMM);
  pbackX = nonblocking_cart_neighb(NMM);

  /* MPI_Irecv(&lb->hl.recvforwYZ[0], nsendYZ, MPI_DOUBLE, pforwX, tagb, comm,  */
  /*        &lb->hl.recvreqP[0]); */
  /* MPI_Irecv(&lb->hl .recvbackYZ[0], nsendYZ, MPI_DOUBLE, pbackX, tagf, comm,  */
  /*        &lb->hl.recvreqP[1]); */

  /* Receive planes in the Y-direction */
  pforwY = nonblocking_cart_neighb(MPM);
  pbackY = nonblocking_cart_neighb(MNM);

  /* MPI_Irecv(&lb->hl .recvforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagb, comm,  */
  /*        &lb->hl.recvreqP[2]); */
  /* MPI_Irecv(&lb->hl .recvbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagf, comm,  */
  /*        &lb->hl.recvreqP[3]); */
  /* MPI_Irecv(&lb->hl .recvforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagb, comm, */
  /*        &lb->hl.recvreqP2[0]); */
  /* MPI_Irecv(&lb->hl .recvbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagf, comm, */
  /*        &lb->hl.recvreqP2[1]); */

  /* Receive planes in the Z-direction */
  pforwZ = nonblocking_cart_neighb(MMP);
  pbackZ = nonblocking_cart_neighb(MMN);

  /* MPI_Irecv(&lb->hl .recvforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagb, comm,  */
  /*        &lb->hl.recvreqP[4]); */
  /* MPI_Irecv(&lb->hl .recvbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagf, comm,  */
  /*        &lb->hl.recvreqP[5]); */
  /* MPI_Irecv(&lb->hl .recvforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagb, comm, */
  /*           &lb->hl.recvreqP2[0]); */
  /* MPI_Irecv(&lb->hl .recvbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagf, comm, */
  /*           &lb->hl.recvreqP2[1]); */

  /* Send in the X-direction (YZ plane) */
  /**************************************/
  int YZ_size = nsendYZ;
  const gaspi_size_t seg_size =  2 * YZ_size * sizeof(double);

  /* segment ids */
  const gaspi_segment_id_t seg_id_YZ_L = 0;
  const gaspi_segment_id_t seg_id_YZ_R = 1;

  int offset = YZ_size*sizeof(double);

  gaspi_pointer_t gptr_YZ_L, gptr_YZ_R;

  /* pointer to the right */
  GASPIERROR(gaspi_segment_ptr(seg_id_YZ_L, &gptr_YZ_L));
  double* ptr_YZ_L =  (double*)gptr_YZ_L;
  /* pointer to the left */
  GASPIERROR(gaspi_segment_ptr(seg_id_YZ_R, &gptr_YZ_R));
  double* ptr_YZ_R =  (double*)gptr_YZ_R;

  /**************************************/


  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(nlocal[X], jc, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl.sendforwYZ[count] = fptr[indexreal];
          ptr_YZ_R[count] = fptr[indexreal];

          index = coords_index(1, jc, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl.sendbackYZ[count] = fptr[indexreal];
          ptr_YZ_L[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendYZ);

  gaspi_notification_id_t notif_id_L = (gaspi_notification_id_t)0;
  gaspi_notification_t notif_val_L = (gaspi_notification_t)1;
  gaspi_notification_id_t notif_id_R = (gaspi_notification_id_t)2;
  gaspi_notification_t notif_val_R = (gaspi_notification_t)3;


  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));

  GASPIERROR(gaspi_write_notify( seg_id_YZ_L, 0, pbackX,
                                 seg_id_YZ_L, offset, offset, notif_id_L,
                                 notif_val_L, 0, GASPI_BLOCK));
  GASPIERROR(gaspi_write_notify( seg_id_YZ_R, 0, cart_neighb(FORWARD,X),
                                 seg_id_YZ_R, offset, offset, notif_id_R,
                                 notif_val_R, 1, GASPI_BLOCK));


  wait_or_die(seg_id_YZ_L, notif_id_L, notif_val_L);
  wait_or_die(seg_id_YZ_R, notif_id_R, notif_val_R);


  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));

  /* MPI_Issend(&lb->hl .sendbackYZ [0], nsendYZ, MPI_DOUBLE, pbackX, tagb, comm, */
  /*         &lb->hl.sendreqP[0]); */
  /* MPI_Issend(&lb->hl .sendforwYZ [0], nsendYZ, MPI_DOUBLE, pforwX, tagf, comm, */
  /*         &lb->hl.sendreqP2[0]); */

  /* Send in the Y-direction (XZ plane) */

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



  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(ic, nlocal[Y], kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl .sendforwXZ[count] = fptr[indexreal];
	  ptr_XZ_R[count] = fptr[indexreal];
          
	  index = coords_index(ic, 1, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl .sendbackXZ[count] = fptr[indexreal];
	  ptr_XZ_L[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendXZ);
  
  notif_id_L = (gaspi_notification_id_t)4;
  notif_val_L = (gaspi_notification_t)5;
  notif_id_R = (gaspi_notification_id_t)6;
  notif_val_R = (gaspi_notification_t)7;
  
  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));
  
  
  GASPIERROR( gaspi_write_notify( seg_id_XZ_L, 0, cart_neighb(BACKWARD,Y),
				  seg_id_XZ_L, offset, offset,notif_id_L,
				  notif_val_L, 0, GASPI_BLOCK) );
  GASPIERROR( gaspi_write_notify( seg_id_XZ_R, 0, cart_neighb(FORWARD,Y),
				  seg_id_XZ_R, offset, offset, notif_id_R,
				  notif_val_R, 1, GASPI_BLOCK) );
  
  wait_or_die(seg_id_XZ_L, notif_id_L, notif_val_L);
  wait_or_die(seg_id_XZ_R, notif_id_R, notif_val_R);
  
  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));
  

  /*  MPI_Issend(&lb->hl .sendbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagb, comm,  */
  /*         &lb->hl.sendreqP[2]); */
  /* MPI_Issend(&lb->hl .sendforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagf, comm,  */
  /*         &lb->hl.sendreqP[3]); */
  /* MPI_Issend(&lb->hl .sendbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagb, comm,  */
  /*         &lb->hl.sendreqP2[0]); */
  /* MPI_Issend(&lb->hl .sendforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagf, comm,  */
  /*         &lb->hl.sendreqP2[1]); */


  /* Finally, Send in the Z-direction (XY plane) */

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
  
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {

          index = coords_index(ic, jc, nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl .sendforwXY[count] = fptr[indexreal];
	  ptr_XY_R[count] = fptr[indexreal];

          index = coords_index(ic, jc, 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          //lb->hl .sendbackXY[count] = fptr[indexreal];
	  ptr_XY_L[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendXY);

  notif_id_L = (gaspi_notification_id_t)8;
  notif_val_L = (gaspi_notification_t)9;
  notif_id_R = (gaspi_notification_id_t)10;
  notif_val_R = (gaspi_notification_t)11;
  
  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));
  
  
  GASPIERROR( gaspi_write_notify( seg_id_XY_L, 0, cart_neighb(BACKWARD,Z),
				  seg_id_XY_L, offset, offset,notif_id_L,
				  notif_val_L, 0, GASPI_BLOCK) );
  GASPIERROR( gaspi_write_notify( seg_id_XY_R, 0, cart_neighb(FORWARD,Z),
				  seg_id_XY_R, offset, offset, notif_id_R,
				  notif_val_R, 1, GASPI_BLOCK) );
  
  wait_or_die(seg_id_XY_L, notif_id_L, notif_val_L);
  wait_or_die(seg_id_XY_R, notif_id_R, notif_val_R);
  
  GASPIERROR(gaspi_wait (0, GASPI_BLOCK));
  GASPIERROR(gaspi_wait (1, GASPI_BLOCK));

  /* MPI_Issend(&lb->hl .sendbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagb, comm,  */
  /*         &lb->hl.sendreqP[4]); */
  /* MPI_Issend(&lb->hl .sendforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagf, comm,  */
  /*         &lb->hl.sendreqP[5]); */
  /* MPI_Issend(&lb->hl .sendbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagb, comm, */
  /*            &lb->hl.sendreqP2[0]); */
  /* MPI_Issend(&lb->hl .sendforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagf, comm, */
  /*            &lb->hl.sendreqP2[1]); */

  return;

}
