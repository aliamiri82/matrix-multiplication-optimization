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

void matrix_mult_blocked(int n, int BS, double* A, double* B, double* C) {
    register int ii, jj, kk, i, j, k;
    register double *c_block;
    register double *a_row, *b_col;

    for (ii = 0; ii < n; ii += BS)
        for (jj = 0; jj < n; jj += BS)
            for (kk = 0; kk < n; kk += BS) {

                c_block = C + ii * n + jj;

                for (i = 0; i < BS && (ii + i) < n; i++) {
                    a_row = A + (ii + i) * n + kk;

                    for (j = 0; j < BS && (jj + j) < n; j++) {
                        b_col = B + (kk * n) + (jj + j);
                        register double cij = *(c_block + i * n + j);

                        for (k = 0; k < BS && (kk + k) < n; k++) {
                            cij += (*(a_row + k)) * (*(b_col + k * n));
                        }
                        *(c_block + i * n + j) = cij;
                    }
                }
            }
}



void matrix_mult_ptr_reg_transpose(int n, double* a, double* b, double* c) {
    register double cij;
    register double *at, *bt;
    register int i, j, k;
    for (i = 0; i < n; i++, a += n, b -= n*n){
        for(j  =0; j<n;j++,c++, b += n){
            cij = 0;
            for (k = 0,at = a, bt = b; k < n ; k++, at++, bt++){
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
    int n;

    printf("Input size of matrix, n = ");
    scanf("%d", &n);
    for(int l = 0;l < 20; l++){
    clock_t t0, t1;


    double *A = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *B = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *C = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *transposeB = (double*)_aligned_malloc(n * n * sizeof(double), 64);
for (int i = 0; i < n*n; i++)
    C[i] = 0.0;

    if (A == NULL || B == NULL || C == NULL) {
        printf("Memory Allocation Error\n");
        return -1;
    }

    srand(time(NULL));
    fill(A, n);
    fill(B, n);
    int v;


    scanf("%d", &v);
    t0 = clock();
    matrix_mult_blocked(n, v, A, B,C);
    t1 = clock();
    printf("Execution time of matrix_mult_ptr_reg_blocked = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);



    // t0 = clock();
    // transpose(n, B, transposeB);
    // matrix_mult_ptr_reg_transpose(n, A, transposeB, C);
    // t1 = clock();
    // printf("Execution time of matrix_mult_ptr_reg_transpose = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
    // for (int i =0 ; i < n*n; i++) {
    //     printf("%lf ", C[i]);
    // }



    printf("\n");
    printf("\n");
    //Execute matrix_mult_index
//    t0 = clock();
//    matrix_mult_index(n, A, B, C);
//    t1 = clock();
    //printf("\nExecution time of matrix_mult_index = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//        for (int i =0 ; i < n*n; i++) {
//            printf("%lf ", C[i]);
//        }

    // Execute matrix_mult_ptr_reg
//    t0 = clock();
//    matrix_mult_ptr_reg(n, A, B, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_ptr_reg = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");
//    printf("\n");
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
