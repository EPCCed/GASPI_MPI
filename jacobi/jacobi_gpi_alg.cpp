#include <GASPI.h>
#include "jacobi_gpi_alg.h"
#include "success_or_die.h"
#include "wait_for_queue_entries.h"
#include <assert.h>
#include <algorithm>
#include "dist.h"
#include <iostream>

void jacobi ( int n
            , int n_local_rows
            , double* local_a
            , double* local_b
            , double** p
            , double* x_new
            , int max_iter
            , double tol
            )
{
  double* x = *p;

  gaspi_rank_t my_gaspi_rank, n_gaspi_procs;
  SUCCESS_OR_DIE (gaspi_proc_rank, &my_gaspi_rank);
  SUCCESS_OR_DIE (gaspi_proc_num, &n_gaspi_procs);

  gaspi_segment_id_t segment_id_from = 0;
  gaspi_segment_id_t segment_id_to = 1;

  SUCCESS_OR_DIE ( gaspi_segment_use
                 , segment_id_from
                 , x
                 , n*sizeof (double)
                 , GASPI_GROUP_ALL
                 , GASPI_BLOCK
                 , 0
                 );

  SUCCESS_OR_DIE ( gaspi_segment_use
                 , segment_id_to
                 , x_new
                 , n*sizeof (double)
                 , GASPI_GROUP_ALL
                 , GASPI_BLOCK
                 , 0
                 );

  int offset = my_gaspi_rank * (n / n_gaspi_procs)
             + std::min ((int)my_gaspi_rank, n % n_gaspi_procs);

  gaspi_queue_id_t queue = 0;

  int n_iterations = 0;
  do
  {
    if (n_iterations > 0)
    {
      std::swap (segment_id_from, segment_id_to);
      std::swap (x, x_new);
    }

    #pragma omp parallel for
    for (int i = 0; i < n_local_rows; i++)
    {
      double d = local_b[i];
      for (int j = 0; j < i + offset; j++)
        d -= local_a[i*n+j] * x[j];
      for (int j = i + offset +1; j < n; j++)
        d -= local_a[i*n+j] * x[j];

      x_new[i+offset] = d/local_a[i*n + i + offset];
    }

    wait_for_queue_entries (&queue, n_gaspi_procs);

    for (int dest = 0; dest < n_gaspi_procs ; dest++)
    {
      if (dest == my_gaspi_rank)
        continue;

      SUCCESS_OR_DIE
        ( gaspi_write_notify
        , segment_id_to
        , offset * sizeof(double)
        , dest
        , segment_id_to
        , offset * sizeof(double)
        , n_local_rows * sizeof (double)
        , (gaspi_notification_id_t) (my_gaspi_rank)
        , (gaspi_notification_t) (my_gaspi_rank + 1)
        , queue
        , GASPI_BLOCK
        );
    }

    int completed = 0;
    while (completed < n_gaspi_procs - 1)
    {
      gaspi_notification_id_t received_notification;
      gaspi_return_t rv = gaspi_notify_waitsome
        ( segment_id_to
        , 0
        , n_gaspi_procs
        , &received_notification
        , GASPI_BLOCK
        );

      gaspi_notification_t value;
      SUCCESS_OR_DIE ( gaspi_notify_reset
                     , segment_id_to
                     , received_notification
                     , &value
                     );
      if (value)
      {
        completed++;
      }
    }
  } while (n_iterations++ < max_iter && error (x_new, x, n) >= tol);
}
