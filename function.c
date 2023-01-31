#include <math.h>
#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
double const const15 = 1e-15;
int swap(double *a, double *b)
{
    double c;
    c = *a;
    *a = *b;
    *b = c;
    return 0;
}
int swap1(int *a, int *b)
{
    int c;
    c = *a;
    *a = *b;
    *b = c;
    return 0;
}
double normaM(const double *mat, const double *b, const double *x, int n, int p, int num, double res1, double res2) //норма невязки
{
    int i;
    int first_row = n * p / num;
    int last_row = (n * (p + 1) / num - 1);
    double s = 0;
    double y = 0;
    double t = 0;
    for (i = first_row * n; i <= last_row * n; i++) {
        if (i % n == 0) {
            y = 0;
        }
        y = y + mat[i - n * first_row] * x[i % n];
        if (i % n == n - 1) {
            s = s + (y - b[i / n]) * (y - b[i / n]);
            t = t + b[i / n] * b[i / n];
        }
    }
    MPI_Allreduce(&s, &res1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&t, &res2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    s = sqrt(s);
    t = sqrt(t);
    if (t < const15) {
        return t;
    }
    return s / t;
}
double norma(const double *x, int n) //норма погрешности
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++) {
        if (i % 2 == 0) {
            s = s + (x[i] - 1) * (x[i] - 1);
        }
        if (i % 2 == 1) {
            s = s + (x[i]) * (x[i]);
        }
    }
    s = sqrt(s);
    return s;
}
double normaMat(double *mat, int n, int p, int num, double res1)
{
    //норма матрицы для алгоритма
    double max = 0.0;
    double s = 0.0;
    int i;
    int first_row = n * p / num;
    int last_row = (n * (p + 1) / num - 1);
    for (i = first_row; i <= last_row; i++) {
        for (int j = 0; j < n; j++) {
            s = s + fabs(mat[i - first_row]);
            if (j == n - 1) {
                if (max < s) {
                    max = s;
                }
                s = 0;
            }
        }
    }
    MPI_Allreduce(&max, &res1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return res1;
}
int function(int n, double *mat, double *b, int *ind, int p, int num, double *per1, double *per2, int *err)
{
    int first_row = n * p / num;
    int last_row = (n * (p + 1) / num - 1);
    int temp[3];
    double res1 = 1;
    double eps = fmin(const15, const15 * normaMat(mat, n, p, num, res1));
    // double eps = 1e-12;
    double tp;
    int k;
    int j;
    int sign = 0;
    int i = 0;
    while (i < n) {
        int c_max = i;
        // нахождение максимального элемента, его процесса и строки
        int r_max = fmax(i, first_row);
        for (int k = fmax(i, first_row); k <= last_row; k++) {
            for (int l = i; l < n; l++) {
                if (fabs(mat[(k - first_row) * n + l]) > fabs(mat[(r_max - first_row) * n + c_max])) {
                    r_max = k;
                    c_max = l;
                }
            }
        }
        struct double_int {
            double val;
            int rank;
        } local_max, global_max;
        MPI_Barrier(MPI_COMM_WORLD);
        if (r_max > last_row) {
            local_max.val = 0;
        } else {
            local_max.val = fabs(mat[(r_max - first_row) * n + c_max]);
        }
        local_max.rank = p;
        MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                      MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (p == global_max.rank) {
            temp[0] = r_max;
            temp[1] = c_max;
            if (mat[(r_max - first_row) * n + c_max] > eps) {
                temp[2] = 1;
            } else {
                temp[2] = -1;
            }
        }
        MPI_Bcast(temp, 3, MPI_INT, global_max.rank, MPI_COMM_WORLD);
        r_max = temp[0];
        c_max = temp[1];
        sign = temp[2];
        // cout << global_max.val <<'\n';
        if (global_max.val < eps) // Проверка на вырожденность
        {
            *err = -3;
            return -3;
        }
        // MPI_Barrier(MPI_COMM_WORLD);
        for (int l = first_row; l <= last_row; l++) { // перестановка столбцов
            swap(&mat[(l - first_row) * n + i],
                 &mat[(l - first_row) * n + c_max]);
        }
        swap1(&ind[i], &ind[c_max]);
        int rowI = 0; // поток содержащий i-ую строчку
        while (n * (rowI + 1) / num - 1 < i) {
            rowI++;
        }
        int rowMax = 0; // поток, содержащий строчку с макс элементом
        while (n * (rowMax + 1) / num - 1 < r_max) {
            rowMax++;
        }
        if (p == rowI) {
            for (int l = 0; l < n; l++) {
                per1[l] = mat[(i - first_row) * n + l];
            }
            per1[n] = b[i - first_row];
        }

        MPI_Bcast(per1, n + 1, MPI_DOUBLE, rowI, MPI_COMM_WORLD);
        if (p == rowMax) {

            for (int l = 0; l < n; l++) {
                per2[l] = mat[(r_max - first_row) * n + l];
            }
            per2[n] = b[r_max - first_row];
            for (int l = 0; l < n; l++) {
                mat[(r_max - first_row) * n + l] = per1[l];
            }
            b[r_max - first_row] = per1[n];
        }
        // отправляем строку с главным элементом
        MPI_Bcast(per2, n + 1, MPI_DOUBLE, rowMax, MPI_COMM_WORLD);
        if (p == rowI) { //Нормировка главной строки
            for (int l = 0; l < n; l++) {
                mat[(i - first_row) * n + l] = per2[l] / global_max.val * sign;
            }
            b[i - first_row] = per2[n] / global_max.val * sign;
        }
        for (j = first_row; j <= last_row; j++) { //вычитаем строки
            if (j != i && fabs(mat[(j - first_row) * n + i]) > eps) {
                tp = mat[(j - first_row) * n + i];
                for (k = i; k < n; k++) {
                    mat[(j - first_row) * n + k] -= (per2[k] / global_max.val * sign) * tp;
                }
                b[j - first_row] = b[j - first_row] - per2[n] / global_max.val * sign * tp;
            }
        }
        i++;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}