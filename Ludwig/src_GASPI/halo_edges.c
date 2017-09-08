
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
 *  halo_edges
 *
 *  Sends 12 MPI messages (the edges) for the Non-Blocking version.
 *
 *****************************************************************************/


void halo_edges_GASPI(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagnn = 903;
  const int tagnp = 904;
  const int tagpn = 905;
  const int tagpp = 906;

  /* Ranks of neighbouring edges.
     Xpp, X is the direction, parallel to which edges are
     sent pp refers to the YZ directions respectively*/
  int Xpp, Xpn, Xnp, Xnn, Ypp, Ypn, Ynp, Ynn, Zpp, Zpn, Znp, Znn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  int nsendX, nsendY, nsendZ;
  nsendX = NVEL*lb->ndist*nlocal[X];
  nsendY = NVEL*lb->ndist*nlocal[Y];
  nsendZ = NVEL*lb->ndist*nlocal[Z];

  /* Allocate message sizes for edges send/receives */
  lb->hl .sendXnn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXnn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXnp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXnp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXpn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXpn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXpp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXpp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendYnn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYnn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendYnp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYnp = (double *) malloc(nsendY*sizeof(double));

  lb->hl .sendYpn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYpn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendYpp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYpp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendZnn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZnn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZnp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZnp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZpn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZpn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZpp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZpp = (double *) malloc(nsendZ*sizeof(double));

  /* Receive edges parallel to x-direction*/
  Xnn = nonblocking_cart_neighb(MNN);
  Xnp = nonblocking_cart_neighb(MNP);
  Xpn = nonblocking_cart_neighb(MPN);
  Xpp = nonblocking_cart_neighb(MPP);

  MPI_Irecv(&lb->hl.recvXnn[0], nsendX, MPI_DOUBLE, Xnn, tagpp, comm, &lb->hl.recvreqE[0]);
  MPI_Irecv(&lb->hl.recvXnp[0], nsendX, MPI_DOUBLE, Xnp, tagpn, comm, &lb->hl.recvreqE[1]);
  MPI_Irecv(&lb->hl.recvXpn[0], nsendX, MPI_DOUBLE, Xpn, tagnp, comm, &lb->hl.recvreqE[2]);
  MPI_Irecv(&lb->hl.recvXpp[0], nsendX, MPI_DOUBLE, Xpp, tagnn, comm, &lb->hl.recvreqE[3]);

  /* Receive edges parallel to y-direction*/
  Ynn = nonblocking_cart_neighb(NMN);
  Ynp = nonblocking_cart_neighb(NMP);
  Ypn = nonblocking_cart_neighb(PMN);
  Ypp = nonblocking_cart_neighb(PMP);

  MPI_Irecv(&lb->hl.recvYnn[0], nsendY, MPI_DOUBLE, Ynn, tagpp, comm, &lb->hl.recvreqE[4]);
  MPI_Irecv(&lb->hl.recvYnp[0], nsendY, MPI_DOUBLE, Ynp, tagpn, comm, &lb->hl.recvreqE[5]);
  MPI_Irecv(&lb->hl.recvYpn[0], nsendY, MPI_DOUBLE, Ypn, tagnp, comm, &lb->hl.recvreqE[6]);
  MPI_Irecv(&lb->hl.recvYpp[0], nsendY, MPI_DOUBLE, Ypp, tagnn, comm, &lb->hl.recvreqE[7]);

  /* Receive edges parallel to z-direction*/
  Znn = nonblocking_cart_neighb(NNM);
  Znp = nonblocking_cart_neighb(NPM);
  Zpn = nonblocking_cart_neighb(PNM);
  Zpp = nonblocking_cart_neighb(PPM);

  MPI_Irecv(&lb->hl.recvZnn[0], nsendZ, MPI_DOUBLE, Znn, tagpp, comm, &lb->hl.recvreqE[8]);
  MPI_Irecv(&lb->hl.recvZnp[0], nsendZ, MPI_DOUBLE, Znp, tagpn, comm, &lb->hl.recvreqE[9]);
  MPI_Irecv(&lb->hl.recvZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagnp, comm, &lb->hl.recvreqE[10]);
  MPI_Irecv(&lb->hl.recvZpp[0], nsendZ, MPI_DOUBLE, Zpp, tagnn, comm, &lb->hl.recvreqE[11]);

  /* Send edges parallel to x-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {

        index = coords_index(ic, 1, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXnn[count] = fptr[indexreal];

        index = coords_index(ic, 1, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXnp[count] = fptr[indexreal];

        index = coords_index(ic, nlocal[Y], 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXpn[count] = fptr[indexreal];

        index = coords_index(ic, nlocal[Y], nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendX);

  MPI_Issend(&lb->hl .sendXpp[0], nsendX, MPI_DOUBLE, Xpp, tagpp, comm, &lb->hl.sendreqE[0]);
  MPI_Issend(&lb->hl .sendXpn[0], nsendX, MPI_DOUBLE, Xpn, tagpn, comm, &lb->hl.sendreqE[1]);
  MPI_Issend(&lb->hl .sendXnp[0], nsendX, MPI_DOUBLE, Xnp, tagnp, comm, &lb->hl.sendreqE[2]);
  MPI_Issend(&lb->hl .sendXnn[0], nsendX, MPI_DOUBLE, Xnn, tagnn, comm, &lb->hl.sendreqE[3]);

  /* Send edges parallel to y-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {

        index = coords_index(1, jc, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYnn[count] = fptr[indexreal];

        index = coords_index(1, jc, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYnp[count] = fptr[indexreal];

        index = coords_index(nlocal[X], jc, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYpn[count] = fptr[indexreal];

        index = coords_index(nlocal[X], jc, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendY);

  MPI_Issend(&lb->hl .sendYpp[0], nsendY, MPI_DOUBLE, Ypp, tagpp, comm, &lb->hl.sendreqE[4]);
  MPI_Issend(&lb->hl .sendYpn[0], nsendY, MPI_DOUBLE, Ypn, tagpn, comm, &lb->hl.sendreqE[5]);
  MPI_Issend(&lb->hl .sendYnp[0], nsendY, MPI_DOUBLE, Ynp, tagnp, comm, &lb->hl.sendreqE[6]);
  MPI_Issend(&lb->hl .sendYnn[0], nsendY, MPI_DOUBLE, Ynn, tagnn, comm, &lb->hl.sendreqE[7]);

  /* Send edges parallel to z-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = coords_index(1, 1, kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZnn[count] = fptr[indexreal];

        index = coords_index(1, nlocal[Y], kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZnp[count] = fptr[indexreal];

        index = coords_index(nlocal[X], 1, kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZpn[count] = fptr[indexreal];

        index = coords_index(nlocal[X], nlocal[Y], kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendZ);

  MPI_Issend(&lb->hl .sendZpp[0] , nsendZ, MPI_DOUBLE, Zpp, tagpp, comm, &lb->hl.sendreqE[8]);
  MPI_Issend(&lb->hl .sendZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagpn, comm, &lb->hl.sendreqE[9]);
  MPI_Issend(&lb->hl .sendZnp[0], nsendZ, MPI_DOUBLE, Znp, tagnp, comm, &lb->hl.sendreqE[10]);
  MPI_Issend(&lb->hl .sendZnn[0], nsendZ, MPI_DOUBLE, Znn, tagnn, comm, &lb->hl.sendreqE[11]);

  return;
}
