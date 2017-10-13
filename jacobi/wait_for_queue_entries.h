#ifndef WAIT_FOR_QUEUE_ENTRIES_H_
#define WAIT_FOR_QUEUE_ENTRIES_H_

#include <stdlib.h>
#include "success_or_die.h"

static void
  wait_for_queue_entries ( gaspi_queue_id_t* queue
                         , gaspi_number_t n_requested_entries
                         )
{
  gaspi_number_t queue_size_max;
  gaspi_number_t queue_size;
  gaspi_number_t queue_num;

  SUCCESS_OR_DIE (gaspi_queue_size_max, &queue_size_max);
  SUCCESS_OR_DIE (gaspi_queue_size, *queue, &queue_size);
  SUCCESS_OR_DIE (gaspi_queue_num, &queue_num);

  if (! (queue_size + n_requested_entries <= queue_size_max))
  {
    *queue = (*queue + 1) % queue_num;
    SUCCESS_OR_DIE (gaspi_wait, *queue, GASPI_BLOCK);
  }
}

#endif /* WAIT_FOR_QUEUE_ENTRIES_H_ */
