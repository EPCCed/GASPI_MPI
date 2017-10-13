#ifndef DIST_H_
#define DIST_H_

#include <cmath>

const int MAX_ITER = 1000;
const double TOL = 0.01;

static double error (double* x, double* y, int n)
{
  double sum_1 = 0.0, sum_2 = 0.0;

  for (int i = 0; i < n; i++)
  {
    sum_1 += (x[i] - y[i])*(x[i] - y[i]);
    sum_2 += x[i]*x[i];
  }

  return sqrt (sum_1/sum_2);
}

#endif
