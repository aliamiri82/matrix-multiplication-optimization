#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void fill(double* x, int n) {
    int i;
    for (i = 0, n = n * n; i < n; i++, x++)
        *x = ((double)(1 + rand() % 12345)) / ((double)(1 + rand() % 6789));
}

void matrix_mult_index(int n, double* a, double* b, double* c) {
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            c[i * n + j] = 0;
            for (k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
        }
}

void matrix_mult_ptr_reg(int n, double* a, double* b, double* c) {
    register double cij;
    register double *at, *bt;
    register int i, j, k;
    for (i = 0; i < n; i++, a += n)
        for (j = 0; j < n; j++, c++) {
            cij = 0;
            for (k = 0, at = a, bt = &b[j]; k < n; k++, at++, bt += n)
                cij += *at * *bt;
            *c = cij;
        }
}

void matrix_mult_ptr_no_reg(int n, double* a, double* b, double* c) {
    double cij;
    double *at, *bt;
    int i, j, k;
    for (i = 0; i < n; i++, a += n)
        for (j = 0; j < n; j++, c++) {
            cij = 0;
            for (k = 0, at = a, bt = &b[j]; k < n; k++, at++, bt += n)
                cij += *at * *bt;
            *c = cij;
        }
}
void matrix_mult_blocked(int n, int BS, double* a, double* b, double* c) {
    register double *at, *bt;
    register int i, j, k, ii, kk, jj, counter_a , counter_b;
    for (ii = 0; ii < n; ii += BS, a += n * BS)
        for (jj = 0; jj < n; jj += BS, b -= (BS * n) * counter_b){
                counter_b = 0;
            for (kk = 0; kk < n; kk += BS, a -= counter_a * n, b += BS * n){
                counter_b++;
                counter_a = 0;
                for (i = ii; i < ii + BS && i < n; i++ , a += n){
                    counter_a++;
                    for (j = jj; j < jj + BS && j < n; j++)
                        for (k = kk , at = a ,bt = b + j; k < kk + BS && k < n; k++ ,at++ ,bt += n)
                            c[i*n + j] += *at * *bt;
                }
            }
        }
}

void matrix_mult_ptr_reg_transpose(int n, double* a, double* b, double* c) {
    register double cij;
    register double *at, *bt;
    register int i, j, k;
    for (i = 0; i < n; i++, a += n){
        for(j  =0; j<n;j++,c++){
            cij = 0;
            for (k = 0,at = a, bt = &b[j]; k < n ; k++, at++, bt++){
                cij += *at * *bt;
            }
            *c = cij;
    }
}
}


void transpose(int n, double* b, double* transposeB) {
    register int i, j;
    register double cij;
    register double *bt;
    for (i = 0; i < n; i++){
        for (j = 0, bt = &b[j*n + i]; j < n; j++ , bt+= n , transposeB++){
            cij = *bt;
            *transposeB = cij;
        }
    }
}




int main() {
    for(int l = 0;l < 20; l++){
    clock_t t0, t1;
    int n;

    printf("Input size of matrix, n = ");
    scanf("%d", &n);

    double *A = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *B = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *C = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *transposeB = (double*)_aligned_malloc(n * n * sizeof(double), 64);

    if (A == NULL || B == NULL || C == NULL) {
        printf("Memory Allocation Error\n");
        return -1;
    }

    srand(time(NULL));
    fill(A, n);
    fill(B, n);


    t0 = clock();
    transpose(n, B, transposeB);
    matrix_mult_blocked(n, 2, A, B,C);
    t1 = clock();
    printf("\nExecution time of matrix_mult_ptr_reg_transpose = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
    for ( int i = 0; i < n*n; i++)
    {
        printf("%f ",C[i]);
    }
    
    
    
    t0 = clock();
    matrix_mult_ptr_reg_transpose(n, A, transposeB, C);
    t1 = clock();
    printf("\nExecution time of matrix_mult_ptr_reg_transpose = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
    for ( int i = 0; i < n*n; i++)
    {
        printf("%f ",C[i]);
    }
//    t0 = clock();
//    transpose(n, B, transposeB);
//    matrix_mult_ptr_reg_transpose(n, A, transposeB, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_ptr_reg_block = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//
//    
//    // Execute matrix_mult_index
//    t0 = clock();
//    matrix_mult_index(n, A, B, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_index = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//
//
//    // Execute matrix_mult_ptr_reg
//    t0 = clock();
//    matrix_mult_ptr_reg(n, A, B, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_ptr_reg = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//
//
//    // Execute matrix_mult_ptr_no_reg
//    t0 = clock();
//    matrix_mult_ptr_no_reg(n, A, B, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_ptr_no_reg = %0.2f s\n", (float)(t1 - t0) / CLOCKS_PER_SEC);



    _aligned_free(A);
    _aligned_free(B);
    _aligned_free(C);
    _aligned_free(transposeB);

    }
    return 0;
}
