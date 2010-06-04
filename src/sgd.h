#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "loss.h"
#include "scale.h"
#include "evaluation.h"

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
typedef void (*dloss_pt)(dtype *, double, dtype, int, double *);
typedef void (*predict_gmatrix)(gmatrix *, double *, double *, int *);


double sgd_gmatrix(gmatrix *g,
   dloss_pt, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

void writevectorf(char* file, double* beta, int p);
void writevectorl(char* file, int* beta, int p);

