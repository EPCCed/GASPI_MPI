#ifndef INIT_DATA_H
#define INIT_DATA_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <assert.h>

int init_input_data (double** a, double** b)
{
  int n;
  printf("Enter the dimension of the square matrix: ");
  scanf("%d", &n);

  *a = (double*) calloc (n*n, sizeof (double));
  *b = (double*) calloc (n, sizeof (double));

  int i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      if (i == j)
        (*a)[i*n+j] = n;
      else
        (*a)[i*n+j] = 1.0;
    }

    (*b)[i] = 1.0;
  }

  return n;
}

void init_solution (double* x, int n)
{
  for (int i =0; i<n; i++)
    x[i] = i;
}

int get_num_rows (int k, int n, int n_procs)
{
  return n/n_procs + (k<n%n_procs ? 1 : 0);
}

int get_offset (int k, int n, int n_procs)
{
  return k * (n/n_procs) + std::min (k, n % n_procs);
}

void comp_counts_and_displs (int n, int m, int n_procs, int** counts, int** displs)
{
  for (int k = 0; k<n_procs; k++)
  {
    (*counts)[k] = m * (n/n_procs + (k<n%n_procs ? 1 : 0));
    (*displs)[k] = m * (k * (n/n_procs) + std::min (k, n % n_procs));
  }
}

void check (double* x, int n)
{
  for (int i=0; i<n; i++)
    assert (x[i] == i);
}
#endif
