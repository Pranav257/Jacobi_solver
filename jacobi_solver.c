/* Code for the Jacobi method of solving a system of linear equations 
 * by iteration.

 * Author: Naga Kandasamy
 * Date modified: APril 26, 2023
 *
 * Student name(s): Baran Arig, Pranav Reddy
 * Date modified: 5/10/23
 *
 * Compile as follows:
 * gcc -o jacobi_solver jacobi_solver.c compute_gold.c -Wall -O3 -lpthread -lm 
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "jacobi_solver.h"
#define _POSIX_C_SOURCE 200112L /* Or higher */


/* Uncomment the line below to spit out debug information */ 
//#define DEBUG

pthread_mutex_t mutex;
pthread_barrier_t barrier, barrier2;


typedef struct myThreadStruct1{
	int num_threads;
	int tid;
	int max_iter;
	matrix_t A;
	matrix_t B;
	matrix_t mt_sol_x_v1;
	float chunkSize;
} myThreadStruct1;

typedef struct myThreadStruct2{
	int num_threads;
	int tid;
	int max_iter;
	matrix_t A;
	matrix_t B;
	matrix_t mt_sol_x_v2;
	int num_elements;
	int stride;
} myThreadStruct2;

//signatures for fxns
void helper_v1(void* ts1_void);
void helper_v2(void* ts1_void);


int main(int argc, char **argv) 
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s matrix-size num-threads\n", argv[0]);
        fprintf(stderr, "matrix-size: width of the square matrix\n");
        fprintf(stderr, "num-threads: number of worker threads to create\n");
		exit(EXIT_FAILURE);
	}

    int matrix_size = atoi(argv[1]);
    int num_threads = atoi(argv[2]);

    matrix_t  A;                    /* N x N constant matrix */
	matrix_t  B;                    /* N x 1 b matrix */
	matrix_t reference_x;           /* Reference solution */ 
    matrix_t mt_solution_x_v1;      /* Solution computed by pthread code using chunking */
    matrix_t mt_solution_x_v2;      /* Solution computed by pthread code using striding */

	/* Generate diagonally dominant matrix */
    fprintf(stderr, "\nCreating input matrices\n");
	srand(time(NULL));
	A = create_diagonally_dominant_matrix(matrix_size, matrix_size);
	if (A.elements == NULL) {
        fprintf(stderr, "Error creating matrix\n");
        exit(EXIT_FAILURE);
	}
	
    /* Create other matrices */
    B = allocate_matrix(matrix_size, 1, 1);
	reference_x = allocate_matrix(matrix_size, 1, 0);
	mt_solution_x_v1 = allocate_matrix(matrix_size, 1, 0);
    mt_solution_x_v2 = allocate_matrix(matrix_size, 1, 0);

#ifdef DEBUG
	print_matrix(A);
	print_matrix(B);
	print_matrix(reference_x);
#endif

    /* Compute Jacobi solution using reference code */
	fprintf(stderr, "Generating solution using reference code\n");
    int max_iter = 100000; /* Maximum number of iterations to run */
    compute_gold(A, reference_x, B, max_iter);
    display_jacobi_solution(A, reference_x, B); /* Display statistics */
	
	/* Compute the Jacobi solution using pthreads. 
     * Solutions are returned in mt_solution_x_v1 and mt_solution_x_v2.
     * */
    fprintf(stderr, "\nPerforming Jacobi iteration using pthreads using chunking\n");
	compute_using_pthreads_v1(A, mt_solution_x_v1, B, max_iter, num_threads);
    display_jacobi_solution(A, mt_solution_x_v1, B); /* Display statistics */
    
    fprintf(stderr, "\nPerforming Jacobi iteration using pthreads using striding\n");
	compute_using_pthreads_v2(A, mt_solution_x_v2, B, max_iter, num_threads);
    display_jacobi_solution(A, mt_solution_x_v2, B); /* Display statistics */

    free(A.elements); 
	free(B.elements); 
	free(reference_x.elements); 
	free(mt_solution_x_v1.elements);
    free(mt_solution_x_v2.elements);
	
    exit(EXIT_SUCCESS);
}

/* FIXME: Complete this function to perform the Jacobi calculation using pthreads using chunking. 
 * Result must be placed in mt_sol_x_v1. */
void compute_using_pthreads_v1(const matrix_t A, matrix_t mt_sol_x_v1, const matrix_t B, int max_iter, int num_threads)
{
pthread_t* threads = (pthread_t *)malloc(num_threads*sizeof(pthread_t));
myThreadStruct1* ts1;
int i;
pthread_barrier_init(&barrier, NULL, num_threads);
pthread_barrier_init(&barrier2, NULL, num_threads);
pthread_mutex_init(&mutex,NULL);


for( i =0; i<num_threads; i++){
	ts1 = (myThreadStruct1*) malloc(sizeof(myThreadStruct1));
	ts1 -> tid = i;
	ts1 -> A = A;
	ts1->max_iter = max_iter;
	ts1->B = B;
	ts1->mt_sol_x_v1 = mt_sol_x_v1;
	ts1->num_threads = num_threads;
	ts1->chunkSize = floor(B.num_rows/num_threads);

	pthread_create(&threads[i], NULL, (void*)helper_v1, ts1);

}
for(i=0; i<num_threads; i++){
	pthread_join(threads[i], NULL);
}
free(threads);
pthread_mutex_destroy(&mutex);
pthread_barrier_destroy(&barrier);
pthread_barrier_destroy(&barrier2);


}


void helper_v1(void* ts1_void){
	myThreadStruct1*ts1 = ts1_void;

	//define variables for ease
	//float chunksize = ts1->chunkSize;
	matrix_t A = ts1->A;
	matrix_t B = ts1->B;
	matrix_t mt_sol_x_v1 = ts1-> mt_sol_x_v1;
	//float new_x[A.num_rows];

	int max_iter = ts1->max_iter;

	matrix_t new_x = allocate_matrix(A.num_rows, 1, 0);
	int ret;
	int num_iter = 0;

	int begin = ts1->tid * ts1->chunkSize;
	int end = ts1->tid * ts1->chunkSize + ts1->chunkSize;
	
	
	int i, j;
	for (i = begin; i < end; i++){ /* Initialize unknowns */ 
		mt_sol_x_v1.elements[i] = B.elements[i];
		}
	//ret = pthread_barrier_wait(&barrier);
	int done = 0;
	
	//printf("%d %d %d\n", ts1->tid, begin, end);
	
	
	//ret = pthread_barrier_wait(&barrier);

	//exit(1);
	while (!done) {
		//ret = pthread_barrier_wait(&barrier);

		pthread_mutex_lock(&mutex);	

		for (i = begin; i < end; i++){ /* Loop iterates over each unknown */ 
			double sum = 0;
			sum = -A.elements[i * A.num_columns + i] * mt_sol_x_v1.elements[i];
			for ( j = begin; j < end; j++) { /* Implement Equation 1 */ 
					sum += A.elements[i*A.num_columns+j] * mt_sol_x_v1.elements[j];
				
			}

			new_x.elements[i] = (B.elements[i] - sum) / A.elements[i*A.num_columns+i];

		}

		pthread_mutex_unlock(&mutex);	

		ret = pthread_barrier_wait(&barrier);

			/* Update unknown values and test for convergence */

		double ssd = 0;
		pthread_mutex_lock(&mutex);

		for ( i = begin; i < end; i++) {
				ssd += (mt_sol_x_v1.elements[i] - new_x.elements[i]) * (mt_sol_x_v1.elements[i] - new_x.elements[i]);
				mt_sol_x_v1.elements[i] = new_x.elements[i];
		}
		
		pthread_mutex_unlock(&mutex);	

		ret = pthread_barrier_wait(&barrier2);


		pthread_mutex_lock(&mutex);

		double mse = sqrt(ssd); /* Mean squared error. */
		fprintf(stderr, "MSE = %f\n", mse); 

		if (mse < THRESHOLD)
			done = 1;
		
		pthread_mutex_unlock(&mutex);	
		ret = pthread_barrier_wait(&barrier);


	} /* End while. */
}
/* FIXME: Complete this function to perform the Jacobi calculation using pthreads using striding. 
 * Result must be placed in mt_sol_x_v2. */
void compute_using_pthreads_v2(const matrix_t A, matrix_t mt_sol_x_v2, const matrix_t B, int max_iter, int num_threads)
{
pthread_t* threads = (pthread_t *)malloc(num_threads*sizeof(pthread_t));
myThreadStruct2* ts2;
int i;
pthread_barrier_init(&barrier, NULL, num_threads);
pthread_barrier_init(&barrier2, NULL, num_threads);
pthread_mutex_init(&mutex,NULL);


for( i =0; i<num_threads; i++){
	ts2 = (myThreadStruct1*) malloc(sizeof(myThreadStruct1));
	ts2 -> tid = i;
	ts2 -> A = A;
	ts2->max_iter = max_iter;
	ts2->B = B;
	ts2->mt_sol_x_v2 = mt_sol_x_v2;
	ts2->num_threads = num_threads;
	ts2->num_elements = B.num_rows;
	ts2->stride = num_threads;

	pthread_create(&threads[i], NULL, (void*)helper_v2, ts2);

}
for(i=0; i<num_threads; i++){
	pthread_join(threads[i], NULL);
}
free(threads);
pthread_mutex_destroy(&mutex);
pthread_barrier_destroy(&barrier);
pthread_barrier_destroy(&barrier2);

}

void helper_v2(void* ts2_void){
	myThreadStruct2*ts2 = ts2_void;

	int num_threads = ts2->num_threads;
	int tid = ts2->tid;
	matrix_t A = ts2->A;
	matrix_t B = ts2->B;
	matrix_t mt_sol_x_v2 = ts2-> mt_sol_x_v2;
	int num_elements = ts2->num_elements;
	int stride = ts2->stride;
	int ret;

	matrix_t new_x = allocate_matrix(A.num_rows, 1, 0);
/*
    while(ts2->tid < ts2-> num_elements){
        ts2->y[ts2->tid] = ts2->a * ts2->x[ts2->tid] + ts2->y[ts2->tid];
        ts2->tid = ts2->tid + ts2->stride;
    }
*/
int tid_count = ts2->tid;
while(tid_count < num_elements){
//for (i = 0; i < num_threads; i++) /* Initialize unknowns */ 
	mt_sol_x_v2.elements[tid_count] = B.elements[tid_count];
	tid_count = tid_count + stride;
}
int done = 0;

ret = pthread_barrier_wait(&barrier);


while (!done) {
int tid_count2 = ts2->tid;
while(tid_count2 < num_elements){
//for (i = 0; i < n; i++){ /* Loop iterates over each unknown */ sum = 0;
	double sum = 0;
	int tid_count3 = ts2->tid;
	while(tid_count3 < num_elements){
	//for (j = 0; j < n; j++)  /* Implement Equation 1 */ 
		if (tid_count2 != tid_count3)
        sum += A.elements[tid_count2*A.num_columns + tid_count3] * mt_sol_x_v2.elements[tid_count3];
		tid_count3 += stride;
	}

	new_x.elements[tid_count2] = (B.elements[tid_count2] - sum) / A.elements[tid_count2*A.num_columns+tid_count2];
	tid_count2 += stride;
}
ret = pthread_barrier_wait(&barrier2);

    /* Update unknown values and test for convergence */
double ssd = 0;
int tid_count4 = ts2->tid;

pthread_mutex_lock(&mutex);	

while(tid_count4 < num_elements){
//for (i = 0; i < n; i++) {
		ssd += (mt_sol_x_v2.elements[tid_count4] - new_x.elements[tid_count4]) * (mt_sol_x_v2.elements[tid_count4] - new_x.elements[tid_count4]);
        mt_sol_x_v2.elements[tid_count4] = new_x.elements[tid_count4];
		tid_count4 += stride;
    }

pthread_mutex_unlock(&mutex);	


if (sqrt(ssd) < TOLERANCE) done = 1;
} /* End while. */
}




/* Allocate a matrix of dimensions height * width.
   If init == 0, initialize to all zeroes.  
   If init == 1, perform random initialization.
*/
matrix_t allocate_matrix(int num_rows, int num_columns, int init)
{
    int i;    
    matrix_t M;
    M.num_columns = num_columns;
    M.num_rows = num_rows;
    int size = M.num_rows * M.num_columns;
		
	M.elements = (float *)malloc(size * sizeof(float));
	for (i = 0; i < size; i++) {
		if (init == 0) 
            M.elements[i] = 0; 
		else
            M.elements[i] = get_random_number(MIN_NUMBER, MAX_NUMBER);
	}
    
    return M;
}	

/* Print matrix to screen */
void print_matrix(const matrix_t M)
{
    int i, j;
	for (i = 0; i < M.num_rows; i++) {
        for (j = 0; j < M.num_columns; j++) {
			fprintf(stderr, "%f ", M.elements[i * M.num_columns + j]);
        }
		
        fprintf(stderr, "\n");
	} 
	
    fprintf(stderr, "\n");
    return;
}

/* Return a floating-point value between [min, max] */
float get_random_number(int min, int max)
{
    float r = rand ()/(float)RAND_MAX;
	return (float)floor((double)(min + (max - min + 1) * r));
}

/* Check if matrix is diagonally dominant */
int check_if_diagonal_dominant(const matrix_t M)
{
    int i, j;
	float diag_element;
	float sum;
	for (i = 0; i < M.num_rows; i++) {
		sum = 0.0; 
		diag_element = M.elements[i * M.num_rows + i];
		for (j = 0; j < M.num_columns; j++) {
			if (i != j)
				sum += abs(M.elements[i * M.num_rows + j]);
		}
		
        if (diag_element <= sum)
			return -1;
	}

	return 0;
}

/* Create diagonally dominant matrix */
matrix_t create_diagonally_dominant_matrix(int num_rows, int num_columns)
{
	matrix_t M;
	M.num_columns = num_columns;
	M.num_rows = num_rows; 
	int size = M.num_rows * M.num_columns;
	M.elements = (float *)malloc(size * sizeof(float));

    int i, j;
	fprintf(stderr, "Generating %d x %d matrix with numbers between [-.5, .5]\n", num_rows, num_columns);
	for (i = 0; i < size; i++)
        M.elements[i] = get_random_number(MIN_NUMBER, MAX_NUMBER);
	
	/* Make diagonal entries large with respect to the entries on each row. */
    float row_sum;
	for (i = 0; i < num_rows; i++) {
		row_sum = 0.0;		
		for (j = 0; j < num_columns; j++) {
			row_sum += fabs(M.elements[i * M.num_rows + j]);
		}
		
        M.elements[i * M.num_rows + i] = 0.5 + row_sum;
	}

    /* Check if matrix is diagonal dominant */
	if (check_if_diagonal_dominant(M) < 0) {
		free(M.elements);
		M.elements = NULL;
	}
	
    return M;
}



