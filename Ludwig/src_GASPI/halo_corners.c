
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
 *  halo_corners
 *
 *  Sends 8 MPI messages (the corners) for the Non-Blocking version.
 *
 *****************************************************************************/

void halo_corners_GASPI(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int count;
  int index, indexhalo, indexreal;
  int nlocal[3];

  const int tagnnn = 907;
  const int tagnnp = 908;
  const int tagnpn = 909;
  const int tagnpp = 910;
  const int tagpnn = 911;
  const int tagpnp = 912;
  const int tagppn = 913;
  const int tagppp = 914;

  /* Ranks of neighbouring corners XYZ direction*/
  int ppp, ppn, pnp, pnn, npp, npn, nnp, nnn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr;
  if (get_step())
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;

  int nsend;
  nsend = NVEL*lb->ndist*1;

  /* Allocate message sizes for plane send/receives */
  lb->hl .sendnnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnpn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnpp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendpnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendpnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendppn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendppp = (double *) malloc(nsend*sizeof(double));

  lb->hl .recvnnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnpn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnpp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvpnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvpnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvppn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvppp = (double *) malloc(nsend*sizeof(double));

  nnn = nonblocking_cart_neighb(NNN);
  nnp = nonblocking_cart_neighb(NNP);
  npn = nonblocking_cart_neighb(NPN);
  npp = nonblocking_cart_neighb(NPP);
  pnn = nonblocking_cart_neighb(PNN);
  pnp = nonblocking_cart_neighb(PNP);
  ppn = nonblocking_cart_neighb(PPN);
  ppp = nonblocking_cart_neighb(PPP);

  MPI_Irecv(&lb->hl.recvnnn[0], nsend, MPI_DOUBLE, nnn, tagppp, comm, &lb->hl.recvreqC[0]);
  MPI_Irecv(&lb->hl.recvnnp[0], nsend, MPI_DOUBLE, nnp, tagppn, comm, &lb->hl.recvreqC[1]);
  MPI_Irecv(&lb->hl.recvnpn[0], nsend, MPI_DOUBLE, npn, tagpnp, comm, &lb->hl.recvreqC[2]);
  MPI_Irecv(&lb->hl.recvnpp[0], nsend, MPI_DOUBLE, npp, tagpnn, comm, &lb->hl.recvreqC[3]);
  MPI_Irecv(&lb->hl.recvpnn[0], nsend, MPI_DOUBLE, pnn, tagnpp, comm, &lb->hl.recvreqC[4]);
  MPI_Irecv(&lb->hl.recvpnp[0], nsend, MPI_DOUBLE, pnp, tagnpn, comm, &lb->hl.recvreqC[5]);
  MPI_Irecv(&lb->hl.recvppn[0], nsend, MPI_DOUBLE, ppn, tagnnp, comm, &lb->hl.recvreqC[6]);
  MPI_Irecv(&lb->hl.recvppp[0], nsend, MPI_DOUBLE, ppp, tagnnn, comm, &lb->hl.recvreqC[7]);

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {

      index = coords_index(1, 1, 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnnn[count] = fptr[indexreal];

      index = coords_index(1, 1, nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnnp[count] = fptr[indexreal];

      index = coords_index(1, nlocal[Y], 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnpn[count] = fptr[indexreal];

      index = coords_index(1, nlocal[Y], nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnpp[count] = fptr[indexreal];

      index = coords_index(nlocal[X], 1, 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendpnn[count] = fptr[indexreal];

      index = coords_index(nlocal[X], 1, nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendpnp[count] = fptr[indexreal];

      index = coords_index(nlocal[X], nlocal[Y], 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendppn[count] = fptr[indexreal];

      index = coords_index(nlocal[X], nlocal[Y], nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendppp[count] = fptr[indexreal];
      count++;
    }
  }

  MPI_Issend(&lb->hl .sendppp[0], nsend, MPI_DOUBLE, ppp, tagppp, comm, &lb->hl.sendreqC[0]);
  MPI_Issend(&lb->hl .sendppn[0], nsend, MPI_DOUBLE, ppn, tagppn, comm, &lb->hl.sendreqC[1]);
  MPI_Issend(&lb->hl .sendpnp[0], nsend, MPI_DOUBLE, pnp, tagpnp, comm, &lb->hl.sendreqC[2]);
  MPI_Issend(&lb->hl .sendpnn[0], nsend, MPI_DOUBLE, pnn, tagpnn, comm, &lb->hl.sendreqC[3]);
  MPI_Issend(&lb->hl .sendnpp[0], nsend, MPI_DOUBLE, npp, tagnpp, comm, &lb->hl.sendreqC[4]);
  MPI_Issend(&lb->hl .sendnpn[0], nsend, MPI_DOUBLE, npn, tagnpn, comm, &lb->hl.sendreqC[5]);
  MPI_Issend(&lb->hl .sendnnp[0], nsend, MPI_DOUBLE, nnp, tagnnp, comm, &lb->hl.sendreqC[6]);
  MPI_Issend(&lb->hl .sendnnn[0], nsend, MPI_DOUBLE, nnn, tagnnn, comm, &lb->hl.sendreqC[7]);

  return;
}

