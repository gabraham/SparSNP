
#define THIN_WINDOW_SIZE 50
#define THIN_STEP_SIZE 5
#define THIN_COR_MAX 0.90

void copyshrink(double *x, double *y, int n, int p, int *active, int m);
int thin(double *x, int n, int p, int *active,
      double cor_max, int windowsize, int stepsize);


