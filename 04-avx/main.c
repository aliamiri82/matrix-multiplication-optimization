#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h> // For AVX
#define BLOCK_SIZE 16
#include <stdio.h>

void matrix_transpose(int n, double* matrix) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double temp = matrix[i * n + j];
            matrix[i * n + j] = matrix[j * n + i];
            matrix[j * n + i] = temp;
        }
    }
}

void matrix_mult_transpose_vector1(int n, double* a, double* b, double* c) {
    typedef double v4df __attribute__((vector_size(32)));
    v4df citj;
    register v4df cij;
    register v4df *at, *bt;
    register int i, j, k;

    matrix_transpose(n, b);

    for (i = 0; i < n; i++, a += n) {
        for (j = 0, bt = (v4df *)b; j < n; j++, c++) {
            cij = (v4df){0, 0, 0, 0};
            for (k = 0, at = (v4df *)a; k < n; k += 4, at++, bt++) {
                cij += *at * *bt;
            }
            citj = cij;
            *c = citj[0] + citj[1] + citj[2] + citj[3];
        }
    }
    matrix_transpose(n, b);
}

void matrix_mult_transpose_vector_unroll_2_b(int n, double* a, double* b, double* c) {
    typedef double v4df __attribute__((vector_size(32)));
    v4df citj;
    register v4df cij;
    register v4df *at, *bt;
    register int i , j , k;
    /////////////////////
    v4df citj1;
    register v4df cij1;
    register v4df *bt1;

    matrix_transpose(n, b);

    for (i = 0; i < n; i++, a += n, c += n/2) {
        for (j = 0, bt = (v4df *)b,bt1 = (v4df *)(b + n/2 * n) ; j < n; j+=2, c++) {
            cij = (v4df){0, 0, 0, 0};
            cij1 =(v4df){0, 0, 0, 0};
            for (k = 0, at = (v4df *)a; k < n; k += 4, at++, bt++, bt1++) {
                cij += *at * *bt;
                cij1+= *at * *bt1;
            }
            citj = cij;
            citj1 = cij1;
            *c = citj[0] + citj[1] + citj[2] + citj[3];
            *(c + n/2) = citj1[0] + citj1[1] + citj1[2] + citj1[3];
        }
    }
    matrix_transpose(n, b);
}

void matrix_mult_transpose_vector_unroll_2_a(int n, double* a, double* b, double* c) {
    typedef double v4df __attribute__((vector_size(32)));
    v4df citj;
    register v4df cij;
    register v4df *at, *bt;
    int i , j;
    register int k;
    /////////////////////
    v4df citj1;
    register v4df cij1;
    register v4df *at1;

    matrix_transpose(n, b);

    for (i = 0; i < n; i+=2, a += 2*n, c += n) {
        for (j = 0, bt = (v4df *)b ; j < n; j++, c++) {
            cij = (v4df){0, 0, 0, 0};
            cij1 =(v4df){0, 0, 0, 0};
            for (k = 0, at = (v4df *)a, at1 = (v4df *)(a + n); k < n; k += 4, at++, bt++, at1++) {
                cij += *at * *bt;
                cij1+= *at1 * *bt;
            }
            citj = cij;
            citj1 = cij1;
            *c = citj[0] + citj[1] + citj[2] + citj[3];
            *(c + n) = citj1[0] + citj1[1] + citj1[2] + citj1[3];
        }
    }
    matrix_transpose(n, b);
}

void my_matrix_mult_vector(int n, double* a, double* b, double* c) {
    typedef double v4df __attribute__((vector_size(32)));
    v4df cij, cij1, cij2, cij3;
    register v4df *at, *at1, *at2, *at3, *bt;
    register int k;
    int i, j;

    matrix_transpose(n, b);

    for (i = 0; i < n; i += 4 ,a += 4 * n) {
        for (j = 0; j < n; j++) {
            cij  = (v4df){0, 0, 0, 0};
            cij1 = (v4df){0, 0, 0, 0};
            cij2 = (v4df){0, 0, 0, 0};
            cij3 = (v4df){0, 0, 0, 0};
            at  = (v4df *)(a);
            at1 = (v4df *)(a + n);
            at2 = (v4df *)(a + 2*n);
            at3 = (v4df *)(a + 3*n);


            bt  = (v4df *)(b + j * n);

            for (k = 0 ; k < n; k += 4, at++, at1++, at2++, at3++, bt++) {
                cij  += *at  * *bt;
                cij1 += *at1 * *bt;
                cij2 += *at2 * *bt;
                cij3 += *at3 * *bt;
            }
            c[(i + 0) * n + j] = cij[0]  + cij[1]  + cij[2]  + cij[3];
            c[(i + 1) * n + j] = cij1[0] + cij1[1] + cij1[2] + cij1[3];
            c[(i + 2) * n + j] = cij2[0] + cij2[1] + cij2[2] + cij2[3];
            c[(i + 3) * n + j] = cij3[0] + cij3[1] + cij3[2] + cij3[3];
        }
    }
    matrix_transpose(n, b);
}

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



void matrix_mult_ptr_reg_transpose(int n, double* a, double* b, double* c) {
    register double cij;
    register double *at, *bt;
    register int i, j, k;
    for (i = 0; i < n; i++, a += n, b -= n*n){
        for(j  =0; j<n; j++, c++, b += n){
            cij = 0;
            for (k = 0,at = a, bt = b; k < n ; k++, at++, bt++){
                cij += *at * *bt;
            }
            *c = cij;
        }
    }
}
void matrix_mult_ptr_reg_transpose_loop_unrool_2(int n, double* a, double* b, double* c) {
    register double cij ,cij_next;
    register double *at, *bt;
    register int i, j, k;
    for (i = 0; i < n; i++, a += n, b -= n*n){
        for(j  =0; j< n/2 ;j++,c += 2, b += 2*n){
            cij = 0;
            cij_next = 0;
            for (k = 0, at = a, bt = b; k < n ; k++, at++, bt++){
                cij += *at * *bt;
                cij_next += *at * *(bt + n);
            }
            *c = cij;
            *(c + 1) = cij_next;
        }
        if (n & 1) {
            cij = 0.0;
            for (k = 0, at = a, bt = b; k < n; ++k, ++at, ++bt)
                cij += *at * *bt;
            *c = cij;

            ++c;
            b += n;
        }
    }
}
void matrix_mult_ptr_reg_transpose_loop_unrool_4(int n, double *a, double *b, double *c) {
    register double cij, cij_next, next, next2;
    register double *at, *bt;
    register int i, j, k;
    const int quadCols = n / 4;
    const int rem = n & 3;

    for (i = 0; i < n; ++i, a += n, b -= n * n) {
        for (j = 0; j < quadCols; ++j, c += 4, b += 4 * n) {
            cij = cij_next = next = next2 = 0.0;
            for (k = 0, at = a, bt = b; k < n; ++k, ++at, ++bt) {
                double a_val = *at;
                cij      += a_val * bt[0];
                cij_next += a_val * bt[n];
                next     += a_val * bt[2 * n];
                next2    += a_val * bt[3 * n];
            }
            c[0] = cij;
            c[1] = cij_next;
            c[2] = next;
            c[3] = next2;
        }

        if (rem) {
            for (int r = 0; r < rem; ++r, ++c, b += n) {
                double sum = 0.0;
                for (k = 0, at = a, bt = b; k < n; ++k, ++at, ++bt)
                    sum += *at * *bt;
                *c = sum;
            }
        }
    }
}

void matrix_mult_ptr_reg_transpose_loop_unrool_8(int n, double *a, double *b, double *c) {
    register double c0, c1, c2, c3, c4, c5, c6, c7;
    register double *at, *bt;
    register int i, j, k;
    const int octCols = n / 8;
    const int rem = n & 7;

    for (i = 0; i < n; ++i, a += n, b -= n * n) {
        for (j = 0; j < octCols; ++j, c += 8, b += 8 * n) {
            c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;
            for (k = 0, at = a, bt = b; k < n; ++k, ++at, ++bt) {
                double a_val = *at;
                c0 += a_val * bt[0];
                c1 += a_val * bt[n];
                c2 += a_val * bt[2 * n];
                c3 += a_val * bt[3 * n];
                c4 += a_val * bt[4 * n];
                c5 += a_val * bt[5 * n];
                c6 += a_val * bt[6 * n];
                c7 += a_val * bt[7 * n];
            }
            c[0] = c0;
            c[1] = c1;
            c[2] = c2;
            c[3] = c3;
            c[4] = c4;
            c[5] = c5;
            c[6] = c6;
            c[7] = c7;
        }

        if (rem) {
            for (int r = 0; r < rem; ++r, ++c, b += n) {
                double sum = 0.0;
                for (k = 0, at = a, bt = b; k < n; ++k, ++at, ++bt)
                    sum += *at * *bt;
                *c = sum;
            }
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
    for(int l = 0;l < 20; l++){
        scanf("%d", &n);
        clock_t t0, t1;


    double *A = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *B = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *C = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *D = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *E = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *F = (double*)_aligned_malloc(n * n * sizeof(double), 64);
    double *transposeB = (double*)_aligned_malloc(n * n * sizeof(double), 64);
for (int i = 0; i < n*n; i++){
    C[i] = 0.0;
    E[i] = 0.0;
    F[i] = 0.0;
}
    if (A == NULL || B == NULL || C == NULL) {
        printf("Memory Allocation Error\n");
        return -1;
    }

    srand(time(NULL));
    fill(A, n);
    fill(B, n);
    int v;


//    scanf("%d", &v);
//    t0 = clock();
//    matrix_mult_blocked(n, v, A, B,C);
//    t1 = clock();
//    printf("Execution time of matrix_mult_ptr_reg_blocked = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");


//    t0 = clock();
//    matrix_mult_transpose_vector_unroll_2_a(n, A, B, E);
//    t1 = clock();
//    printf("matrix_mult_transpose_vector_unroll_2_a = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");
//    for (int i =0 ; i < n*n; i++) {
//        printf("%lf \n", E[i] );
//    }


//    t0 = clock();
//    matrix_mult_transpose_vector_unroll_4_a(n, A, B, E);
//    t1 = clock();
//    printf("matrix_mult_transpose_vector_unroll_4_a = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");
//    for (int i =0 ; i < n*n; i++) {
//        printf("%lf \n", E[i] );
//    }
//    t0 = clock();
//    matrix_mult_transpose_vector_unroll_2_b(n, A, B, E);
//    t1 = clock();
//    printf("matrix_mult_transpose_vector_unroll_2_b = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");

//    t0 = clock();
//    matrix_mult_transpose_vector1(n, A, B, D);
//    t1 = clock();
//    printf("correct_answer = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    printf("\n");
//        for (int i =0 ; i < n*n; i++) {
//        printf("%lf \n", D[i] );
//    }

//    t0 = clock();
//    transpose(n, B, transposeB);
//    matrix_mult_ptr_reg_transpose_loop_unrool_2(n, A, transposeB, D);
//    t1 = clock();
//    printf("Execution time of matrix_mult_ptr_reg_transpose_2 = %0.2f s", (float)(t1 - t0) / CLOCKS_PER_SEC);
    printf("\n");
//            for (int i =0 ; i < n*n; i++) {
//            printf("%lf ", C[i]);
//        }




//    t0 = clock();
//    transpose(n, B, transposeB);
//    matrix_mult_ptr_reg_transpose_loop_unrool_8(n, A, transposeB, F);
//    t1 = clock();
//    printf("Execution time of matrix_mult_ptr_reg_transpose_8 = %0.2f s \n", (float)(t1 - t0) / CLOCKS_PER_SEC);
//    for (int i =0 ; i < n*n; i++) {
//        printf("%lf ", C[i]);
//    }

//    //Execute matrix_mult_index

//    t0 = clock();
//    matrix_mult_index(n, A, B, C);
//    t1 = clock();
//    printf("\nExecution time of matrix_mult_index = %0.2f s \n", (float)(t1 - t0) / CLOCKS_PER_SEC);

//            for (int i =0 ; i < n*n; i++) {
//        printf("%lf \n ", C[i]);
//    }
        int counter = 0;
        for (int i = 0 ; i < n*n; i++) {

if (fabs(E[i] - D[i]) > pow(10, -6))
    counter++;

        }
        printf("\n");
        printf("Errors:");
        printf("%d\n", counter);
          printf("\n");





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
    _aligned_free(D);
    _aligned_free(E);
    _aligned_free(F);
    _aligned_free(transposeB);

    }
    return 0;
}



