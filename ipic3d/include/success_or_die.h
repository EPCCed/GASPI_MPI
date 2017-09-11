#ifndef SUCCESS_OR_DIE_H
#define SUCCESS_OR_DIE_H

#include <GASPI.h>
#include "assert.h"
#include <stdlib.h>

#define SUCCESS_OR_DIE(f...)                                            \
  do                                                                    \
  {									\
    const gaspi_return_t r = f;						\
    if (r != GASPI_SUCCESS)						\
    {									\
      printf("Error: '%s' [%s:%i]: %i\n",#f,__FILE__,__LINE__,r);       \
      exit (EXIT_FAILURE);                                              \
    }									\
  } while (0)

class GPIUtil
{
 public:
  static GPIUtil& instance();

 private:
  ~GPIUtil(){}
  GPIUtil(){}

 public:
  static void init(int gridXC,int gridYC);
  static void finalize();
  static gaspi_queue_id_t getQueueIdSend();
  static gaspi_queue_id_t getQueueIdRecv();
  static gaspi_segment_id_t getSegIdSend();
  static gaspi_segment_id_t getSegIdRecv();
  static gaspi_pointer_t getSegPtrSend();
  static gaspi_pointer_t getSegPtrRecv();
  static gaspi_rank_t getGPIRank();
  static gaspi_rank_t getGPINprocs();

 private:
  static gaspi_rank_t gpi_rank, gpi_nproc;
  static gaspi_segment_id_t segment_id_send, segment_id_recv;
  static gaspi_pointer_t segment_ptr_send, segment_ptr_recv;
  static gaspi_queue_id_t queue_id_send, queue_id_recv;
  static int isGPIInitialized;
};

#endif
