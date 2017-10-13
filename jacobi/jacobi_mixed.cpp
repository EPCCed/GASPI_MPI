#include "mpi.h"
#include <GASPI.h>
#include "success_or_die.h"
#include "wait_for_queue_entries.h"
#include "init_data.h"
#include <assert.h>
#include <cstring>
#include "dist.h"

int main(int argc, char* argv[])
{
  double *a, *x, *x_new, *b, *local_b;
  int  i, j, n;
  int my_mpi_rank, n_mpi_procs;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &n_mpi_procs);

  gaspi_rank_t my_gaspi_rank, n_gaspi_procs;
  SUCCESS_OR_DIE (gaspi_proc_init, GASPI_BLOCK);
  SUCCESS_OR_DIE (gaspi_proc_rank, &my_gaspi_rank);
  SUCCESS_OR_DIE (gaspi_proc_num, &n_gaspi_procs);

  assert(my_mpi_rank == my_gaspi_rank);
  assert(n_mpi_procs == n_gaspi_procs);

  if (my_mpi_rank == 0)
    n = init_input_data (&a, &b);

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  x = new double[n];
  if (my_mpi_rank == 0)
    init_solution (x, n);

  int* sendcounts = new int[n_mpi_procs];
  int* displs = new int[n_mpi_procs];

  comp_counts_and_displs
    (n, n, n_mpi_procs, &sendcounts, &displs);

  int n_local_rows = get_num_rows (my_mpi_rank, n, n_mpi_procs);

  double* local_a = new double[n_local_rows * n];
  MPI_Scatterv ( a, sendcounts, displs, MPI_DOUBLE, local_a
               , n_local_rows*n, MPI_DOUBLE, 0, MPI_COMM_WORLD
               );

  comp_counts_and_displs
    (n, 1, n_mpi_procs, &sendcounts, &displs);

  local_b = new double[n_local_rows];
  MPI_Scatterv ( b, sendcounts, displs, MPI_DOUBLE, local_b
               , n_local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD
               );

  delete[] sendcounts;
  delete[] displs;

  MPI_Bcast (x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  gaspi_segment_id_t segment_id_from = 0;
  gaspi_segment_id_t segment_id_to = 1;

  SUCCESS_OR_DIE ( gaspi_segment_create
                 , segment_id_from
                 , n*sizeof (double)
                 , GASPI_GROUP_ALL
                 , GASPI_BLOCK
                 , GASPI_MEM_UNINITIALIZED
                 );

  SUCCESS_OR_DIE ( gaspi_segment_create
                 , segment_id_to
                 , n*sizeof (double)
                 , GASPI_GROUP_ALL
                 , GASPI_BLOCK
                 , GASPI_MEM_INITIALIZED
                 );

  int offset = get_offset (my_gaspi_rank, n, n_gaspi_procs);
  gaspi_queue_id_t queue = 0;

  gaspi_pointer_t gaspi_ptr_from, gaspi_ptr_to;
  SUCCESS_OR_DIE (gaspi_segment_ptr, segment_id_from, &gaspi_ptr_from);
  SUCCESS_OR_DIE (gaspi_segment_ptr, segment_id_to, &gaspi_ptr_to);

  memmove (gaspi_ptr_from, x, n*sizeof(double));
  delete[] x;

  x = (double*)gaspi_ptr_from;
  x_new = (double*)gaspi_ptr_to;

  check (x, n);

  int n_iterations = 0;
  do
  {
    if (n_iterations > 0)
    {
      std::swap (segment_id_from, segment_id_to);
      std::swap (gaspi_ptr_from, gaspi_ptr_to);
      std::swap (x, x_new);
    }

    #pragma omp parallel for
    for (i = 0; i < n_local_rows; i++)
    {
      double d = local_b[i];
      for (j = 0; j < i + offset; j++)
        d -= local_a[i*n+j] * x[j];
      for (j = i + offset +1; j < n; j++)
        d -= local_a[i*n+j] * x[j];

      x_new[i + offset] = d/local_a[i*n + i + offset];
    }

    wait_for_queue_entries (&queue, n_gaspi_procs);

    #pragma omp parallel for
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
  } while (n_iterations++ < MAX_ITER && error (x_new, x, n) >= TOL);

  if (my_gaspi_rank == 0)
  {
    for (int k=0; k<n; k++)
      gaspi_printf("rank: %d, x[%d] = %f\n", my_mpi_rank, k, x[k]);

    delete[] a;
    delete[] b;
  }

  delete[] local_a;
  delete[] local_b;

  SUCCESS_OR_DIE (gaspi_proc_term, GASPI_BLOCK);

  MPI_Finalize();

  return 0;
}
