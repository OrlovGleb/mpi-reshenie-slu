#include "function.h"
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// double const const9 = 1e+9;
int const dvadz = 20;
int const c5 = 5;
int const c4 = 4;
int const big = 5000;
const double e92 = 1e9;
const double e93 = 1e+9;
void vvod_b(double *b, const double *A, int n, int first_row, int last_row)
{ //ввод матрицы б
    double a = 0;
    for (int i = first_row; i <= last_row; i++) {
        for (int j = 0; j < (n - 1) / 2 + 1; j++) {
            a = a + A[n * (i - first_row) + 2 * j];
        }
        b[i - first_row] = a;
        a = 0;
    }
}

double currentTimeNano1() //функция для подсчета времени
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)(t.tv_sec * e92 + t.tv_nsec);
}

int main(int argc, char *argv[])
{
    int e = 0;
    int res1 = 0;
    int res2 = 0;
    char filename[dvadz];
    int n = 0;
    int m = 0;
    int k = 0;
    int i;
    double start = 0.0;
    double end = 0.0;
    int first_row = 0;
    int last_row = 0;
    int p; // номер потока
    int num; // число потоков
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    if (argc < c4) {

        return -1;
    }
    if (argc > c5) {

        return -1;
    }
    for (i = 1; i < c4; i++) {
        for (k = 0; k < 1; k++) {
            if (isdigit(argv[i][k]) == 0) {

                return -1;
            }
        }
    }
    if (argv[1] == NULL) {
        return -1;
    }
    if (argv[2] == NULL) {
        return -1;
    }
    if (argv[3] == NULL) {
        return -1;
    }
    n = (atoi(argv[1]));
    m = (atoi(argv[2]));
    k = (atoi(argv[3]));
    if (k == 0) {
        if (argv[4] == NULL) {
            MPI_Finalize();
            return -1;
        }
        strcpy(filename, argv[4]);
    }

    if (n < 1) {
        MPI_Finalize();
        return -1;
    }
    if (n > big) {
        MPI_Finalize();
        return -2;
    }

    if (m > n) {
        MPI_Finalize();
        return -1;
    }
    if (m < 1) {
        MPI_Finalize();
        return -1;
    }
    if (k > 4) {
        MPI_Finalize();
        return -1;
    }

    if (argc == c5) {
    } else if (argc == 4) {
        if (k == 0) {
            MPI_Finalize();
            return -1;
        }
    }

    first_row = n * p / num;
    last_row = n * (p + 1) / num - 1;
    double *A = (double *)malloc(n * (last_row - first_row + 1) * sizeof(double));
    double *b = (double *)malloc((last_row - first_row + 1) * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    double *per1 = (double *)malloc((n + 1) * sizeof(double));
    double *per2 = (double *)malloc((n + 1) * sizeof(double));
    per1[0] = 0;
    per2[0] = 0;
    int ind[n];
    if (k == 0) {
        e = vvod2(n, A, filename, first_row, last_row);
        if (e == -3) {
            free(A);
            free(x);
            free(b);
            MPI_Finalize();
            return -3;
        }
    } else {
        vvod1(n, k, A, first_row, last_row);
    }
    vvod_b(b, A, n, first_row, last_row);
    for (i = 0; i < n; i++) {
        ind[i] = i;
        x[i] = 0.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int h = 0; h < num; h++) {
        if (h == p) {
            for (int i = first_row; i <= last_row && i < m; i++) {
                for (int j = 0; j < m; j++) {
                    if (i != 0 && j == 0) {
                        printf(" ");
                    }
                    printf("%.3e ", A[(i - first_row) * n + j]);
                }
                printf("%.3e", b[(i - first_row)]);
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    double start1 = currentTimeNano1();
    e = function(n, A, b, ind, p, num, per1, per2, &e);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    printf("Time of process ");
    printf("%d", p);
    printf(":");
    printf(" ");
    printf("%f", end - start);
    printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
    if (e == -3) {
        free(A);
        free(x);
        free(b);
        // free(ind);
        free(per1);
        free(per2);
        MPI_Finalize();
        return -4;
    }
    if (k == 0) { //ввод матрицы для подсчета норм
        vvod2(n, A, filename, first_row, last_row);
    } else {
        vvod1(n, k, A, first_row, last_row);
    }
    double a;
    for (int h = 0; h < num; h++) {
        for (int i = 0; i < n; i++) {
            a = x[i];
            MPI_Bcast(&a, 1, MPI_DOUBLE, h, MPI_COMM_WORLD);
            x[i] = a;
        }
        for (int i = first_row; i <= last_row; i++) {
            x[ind[i]] = b[i - first_row];
        }
    }
    vvod_b(b, A, n, first_row, last_row); //ввод матрицы для подсчета норм
    if (p == 0) {
        printf("Solution:\n");
        for (i = 0; i < m; i++) {
            printf("%.3e ", x[i]);
        }
        printf("\n");
    }

    double res = normaM(A, b, x, n, p, num, res1, res2);
    if (p == 0) {
        printf("Residual: "); //вывод нормы невязки
        printf("%.3e", res);
        printf("\n");
    }

    if (p == 0) {
        printf("Error: "); //вывод нормы погрешности
        printf("%.3e", norma(x, n));
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(A);
    free(x);
    free(b);
    free(per1);
    free(per2);
    MPI_Finalize();
    double finish1 = currentTimeNano1();
    if (p == 0) {

        printf("Time: %lf\n", (finish1 - start1) / e93); //вывод времени в секундах
    }

    return 0;
}