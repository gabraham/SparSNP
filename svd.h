#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

void tcrossprod(double *A, double *B, int *m, int *k, int *n, double *C);
int pseudoinverse(double *Aorig, int *m, int *n, double *P);
   

