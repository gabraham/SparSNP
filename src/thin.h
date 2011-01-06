
#define THIN_WINDOW_SIZE 50
#define THIN_STEP_SIZE 5
#define THIN_COR_MAX 0.8

void cov(double *x, double *S, int n, int p);
void cov2cor(double *S, double *P, int p);
void copyshrink(double *x, double *y, int n, int p, int *active, int m);
int thin(double *x, int n, int p, int *active,
      double cor_max, int windowsize, int stepsize);


