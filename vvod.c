#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
int const dvad = 20;

double f(int k, int n, int i, int j)
{
    double a = 0;
    switch (k) {
    case 1:
        a = n - fmax(i + 1, j + 1) + 1;
        break;
    case 2:
        a = fmax(i + 1, j + 1);
        break;
    case 3:
        a = abs(i - j);
        break;
    case 4:
        a = (double)1 / (i + j + 1);
        break;
    }
    return a;
}
int vvod1(int n, int k, double *mat, int first_row, int last_row)
{
    for (int i = first_row; i <= last_row; i++) {
        for (int j = 0; j < n; j++) {
            mat[(i - first_row) * n + j] = f(k, n, i, j);
        }
    }
    return 0;
}

int vvod2(int n, double *A, char filename[dvad], int first_row, int last_row)
{
    double u;
    int s1 = 0;
    FILE *fp;
    fp = fopen(filename, "rt");
    if (fp == NULL) {
        return -3;
    }
    fseek(fp, 0, SEEK_END);
    long pos = ftell(fp);
    if (pos > 0) {
        rewind(fp);
    } else {
        fclose(fp);
        return -3;
    }
    int i = 0;
    for (i = 0; i < n * n; i++) {
        s1 = s1 + 1;
        if (!fscanf(fp, "%lf", &u)) {
            fclose(fp);
            return -3;
        }
        if (first_row * n <= i && i <= last_row * n + n - 1) {
            A[(i / n - first_row) * n + i % n] = u;
        }
    }
    if (s1 < n * n) {
        return -3;
    }

    return 0;
}