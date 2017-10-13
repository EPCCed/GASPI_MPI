#include "mpi.h"
#include "init_data.h"
#include <assert.h>
#include "dist.h"
#include "omp.h"

int main(int argc, char* argv[])

{
  double *a, *b;
  int n, my_mpi_rank, n_mpi_procs;

  int mpisupport;

  MPI_Init_thread
    (&argc, &argv, MPI_THREAD_MULTIPLE, &mpisupport);

  MPI_Comm_rank (MPI_COMM_WORLD, &my_mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &n_mpi_procs);

  if (my_mpi_rank == 0)
    n = init_input_data (&a, &b);

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double* x = new double [n];

  if (my_mpi_rank == 0)
    init_solution (x, n);

  int* counts = new int[n_mpi_procs];
  int* displs = new int[n_mpi_procs];

  comp_counts_and_displs
    (n, n, n_mpi_procs, &counts, &displs);

  int n_local_rows = get_num_rows (my_mpi_rank, n, n_mpi_procs);

  double* local_a = new double[n_local_rows * n];

  MPI_Scatterv ( a, counts, displs
               , MPI_DOUBLE, local_a
               , n_local_rows*n, MPI_DOUBLE
               , 0, MPI_COMM_WORLD
               );

  comp_counts_and_displs
    (n, 1, n_mpi_procs, &counts, &displs);

  double* local_b = new double[n_local_rows];

  MPI_Scatterv ( b, counts, displs
               , MPI_DOUBLE, local_b
               , n_local_rows, MPI_DOUBLE
               , 0, MPI_COMM_WORLD
               );

  MPI_Bcast (x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int offset = get_offset (my_mpi_rank, n, n_mpi_procs);
  double* x_new = new double[n];

  MPI_Request* send_reqs = new MPI_Request[n_mpi_procs-1];
  MPI_Request* recv_reqs = new MPI_Request[n_mpi_procs-1];
  MPI_Status* statuses = new MPI_Status[n_mpi_procs-1];

  int n_iterations = 0;
  do
  {
    if (n_iterations > 0)
      std::swap (x, x_new);

    #pragma omp parallel for
    for (int i = 0; i < n_local_rows; i++)
    {
      double d = local_b[i];

      for (int j = 0; j<i + offset; j++)
        d -= local_a[i*n+j] * x[j];
      for (int j = i + offset + 1; j < n; j++)
        d -= local_a[i*n+j] * x[j];
      x_new[i + offset] = d/local_a[i*n+i+offset];
    }

    #pragma omp parallel for
    for (int rank = 0; rank < n_mpi_procs; rank++)
    {
      if (rank == my_mpi_rank)
        continue;

      MPI_Isend
        ( x_new + get_offset(my_mpi_rank,n, n_mpi_procs)
        , get_num_rows (my_mpi_rank, n, n_mpi_procs)
        , MPI_DOUBLE,rank, 100, MPI_COMM_WORLD
        , &send_reqs[rank<my_mpi_rank?rank:rank-1]
        );
    }

    #pragma omp parallel for
    for (int rank = 0; rank < n_mpi_procs; rank++)
    {
      if (rank == my_mpi_rank)
        continue;

      MPI_Irecv
        ( x_new + get_offset (rank, n, n_mpi_procs)
        , get_num_rows(rank,n,n_mpi_procs)
        , MPI_DOUBLE
        , rank, 100, MPI_COMM_WORLD
        , &recv_reqs[rank<my_mpi_rank?rank:rank-1]
        );
    }

    MPI_Waitall(n_mpi_procs-1,recv_reqs,statuses);
    MPI_Waitall(n_mpi_procs-1,send_reqs,statuses);

  } while (n_iterations++<MAX_ITER && error(x, x_new, n) >= TOL);

  if (my_mpi_rank == 0)
  {
    for (int k=0; k<n; k++)
      printf("rank: %d, x[%d] = %f\n", my_mpi_rank, k, x[k]);

    delete[] a;
    delete[] b;
  }

  delete[] counts;
  delete[] displs;
  delete[] send_reqs;
  delete[] recv_reqs;
  delete[] local_a;
  delete[] local_b;
  delete[] x;
  delete[] x_new;

  MPI_Finalize();

  return 0;
}
