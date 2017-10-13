#include "mpi.h"
#include <GASPI.h>
#include "success_or_die.h"
#include "jacobi_gpi_alg.h"
#include "init_data.h"
#include <assert.h>
#include "dist.h"

int main(int argc, char* argv[])
{
  double *a, *x, *b, *local_b;
  int n, my_mpi_rank, n_mpi_procs;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &n_mpi_procs);

  SUCCESS_OR_DIE (gaspi_proc_init, GASPI_BLOCK);

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

  double* x_new = new double[n];
  jacobi (n, n_local_rows, local_a, local_b, &x, x_new, MAX_ITER, TOL);

  if (my_mpi_rank == 0)
  {
    for (int k=0; k<n; k++)
      gaspi_printf("rank: %d, x[%d] = %f\n", my_mpi_rank, k, x[k]);

    delete[] a;
    delete[] b;
  }

  delete[] local_a;
  delete[] local_b;

  SUCCESS_OR_DIE (gaspi_proc_term, GASPI_BLOCK);
  delete[] x_new;
  delete[] x;

  MPI_Finalize();

  return 0;
}
