#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "loss.h"
#include "evaluation.h"

typedef double (*loss_pt)(dtype *, double *, dtype, int);
typedef double (*predict_pt)(sample *, double *, double *, double *, int);
typedef void (*dloss)(dtype *, double *, dtype, int, double *); 

double scd_gmatrix(gmatrix *g,
   dloss, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

void writevectorf(char* file, double* beta, int p);
void writevectorl(char* file, int* beta, int p);

