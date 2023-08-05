# Jacobi_solver

The code defines a struct named matrix_t to represent matrices and contains their dimensions and elements.

The compute_gold.c file contains the reference implementation of the Jacobi method (not shown here).

The main function takes two command-line arguments, matrix-size and num-threads, specifying the size of the square matrix and the number of threads to use for parallel execution.

The code generates a diagonally dominant matrix (A) and a constant vector (B) of the specified size.

It then computes the Jacobi solution using the reference code and displays the results (the solution matrix and statistics).

Two versions of the Jacobi method are implemented using pthreads:
a. compute_using_pthreads_v1: This version divides the rows of the matrix into equal-sized chunks and assigns each chunk to a different thread for computation.
b. compute_using_pthreads_v2: This version uses striding, where each thread computes the elements of the solution vector with a specific stride.

The two pthread versions are called, and the solutions are displayed along with statistics.

The program then frees the allocated memory and exits.

Note: The code contains a few issues that need to be addressed:

There is a type mismatch in the helper_v2 function. The thread function receives myThreadStruct2* as an argument, but it is incorrectly cast as myThreadStruct1*. This should be fixed by changing the casting to myThreadStruct2*.
The thread function helper_v2 contains a for-loop that increments tid_count, but this should be modified to be the stride of the loop.
The variable THRESHOLD is not defined in the code, which may lead to compilation errors.
The compute_gold.c file, which contains the reference implementation, is not included in the provided code. It should be included in the compilation.
