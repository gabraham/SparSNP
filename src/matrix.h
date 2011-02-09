void crossprod(double *x, double *y, double *z, int m, int n, int p);
void invert2x2(double *y, double *x);
void wcrossprod(double *x, double *y, double *w, double *z,
      int m, int n, int p);
void sqmvprod(double *x, double *y, double *z, int m);
void copyshrink(double *x, double *y, int n, int p, int *active, int m);
void copyshrinkrange(double *x, double *y, int n, int p, int from, int to);
void printmatrix(double *x, int n, int p);
void printmatrix0(double *x, int n, int p);
int cov(double *x, double *S, int n, int p);
void cov2cor(double *S, double *P, int p);
int covmiss(double *x, double *S, int n, int p, int *good);

