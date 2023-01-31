#define FUCNTION_H
int function(int n, double *mat, double *b, int *ind, int p, int num,
              double *per1, double *per2, int *err);
double f(int k, int n, int i, int j);
int vvod1(int n, int k, double *mat, int first_row, int last_row);
int vvod2(int n, double *mat, char filename[20], int first_row,
              int last_row);
void vvod_b(double *b, const double *A, int n, int first_row, int last_row);
double normaM(const double *mat, const double *b, const double *x, int n, int p,
              int num, double res1, double res2);
double norma(const double *x, int n);
double normaMat(double *mat, int n, int p, int num, double res1);