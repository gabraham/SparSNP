void crossprod(double *x, double *y, double *z, int m, int n, int p);
void invert2x2(double *y, double *x);
void wcrossprod(double *x, double *y, double *w, double *z,
      int m, int n, int p);
void sqmvprod(double *x, double *y, double *z, int m);
void copyshrink(double *x, double *y, int n, int p, int *active, int m);
void copyshrinkrange(double *x, double *y, int n, int p, int from, int to);
void printmatrix(double *x, int n, int p);
