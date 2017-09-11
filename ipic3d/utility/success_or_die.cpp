#include "success_or_die.h"
#include <GASPI.h>
#include "assert.h"
#include <stdlib.h>
#include "debug.h"

gaspi_rank_t GPIUtil::gpi_rank = -1;
gaspi_rank_t GPIUtil::gpi_nproc;
gaspi_segment_id_t GPIUtil::segment_id_send;
gaspi_segment_id_t GPIUtil::segment_id_recv;
gaspi_pointer_t GPIUtil::segment_ptr_send;
gaspi_pointer_t GPIUtil:: segment_ptr_recv;
gaspi_queue_id_t GPIUtil::queue_id_send;
gaspi_queue_id_t GPIUtil:: queue_id_recv;
int GPIUtil::isGPIInitialized = 0;

GPIUtil& GPIUtil::instance()
{
  assert(isGPIInitialized);
  static GPIUtil* instance = new GPIUtil;
  return *instance;
}


void GPIUtil::init(int gridXC, int gridYC)
{
  assert(!isGPIInitialized);
  SUCCESS_OR_DIE(gaspi_proc_init(GASPI_BLOCK));
  isGPIInitialized = 1;
  SUCCESS_OR_DIE(gaspi_proc_rank(&gpi_rank));
  SUCCESS_OR_DIE(gaspi_proc_num(&gpi_nproc));
  const gaspi_size_t tot_exchange_size = (2 * (gridXC + gridYC) + 4) * sizeof(double);
  segment_id_send = 0;
  queue_id_send = 0;
  segment_id_recv = 1;
  queue_id_recv = 1;
  SUCCESS_OR_DIE(gaspi_segment_create(segment_id_send, tot_exchange_size, GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT));
  SUCCESS_OR_DIE(gaspi_segment_create(segment_id_recv, tot_exchange_size, GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT));
  SUCCESS_OR_DIE(gaspi_segment_ptr(segment_id_send, &segment_ptr_send));
  SUCCESS_OR_DIE(gaspi_segment_ptr(segment_id_recv, &segment_ptr_recv));
}

void GPIUtil::finalize()
{
  SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
  SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
}

gaspi_queue_id_t GPIUtil::getQueueIdSend()
{
  return instance().queue_id_send;
}

gaspi_queue_id_t GPIUtil::getQueueIdRecv()
{
  return instance().queue_id_recv;
}

gaspi_segment_id_t GPIUtil::getSegIdSend()
{
  return instance().segment_id_send;
}

gaspi_segment_id_t GPIUtil::getSegIdRecv()
{
  return instance().segment_id_recv;
}

gaspi_pointer_t GPIUtil::getSegPtrSend()
{
  return instance().segment_ptr_send;
}

gaspi_pointer_t GPIUtil::getSegPtrRecv()
{
  return instance().segment_ptr_recv;
}

gaspi_rank_t GPIUtil::getGPIRank()
{
  return instance().gpi_rank;
}

gaspi_rank_t GPIUtil::getGPINprocs()
{
  return instance().gpi_nproc;
}
