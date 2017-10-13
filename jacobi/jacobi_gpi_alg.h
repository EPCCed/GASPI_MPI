#ifndef JACOBI_H_
#define JACOBI_H_

#define MIN(a,b) (((a)<(b))?(a):(b))

void jacobi ( int n
            , int n_local_rows
            , double* local_a
            , double* local_b
            , double** x
            , double* x_new
            , int max_iter
            , double tol
            );

#endif /* JACOBI_H_ */
