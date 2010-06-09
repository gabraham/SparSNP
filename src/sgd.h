#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "loss.h"
#include "util.h"
#include "evaluation.h"

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
typedef void (*dloss_pt)(dtype *, double, dtype, int, double *);
typedef void (*predict_gmatrix)(gmatrix *, double *, double *, int *);


typedef double (*optim_gmatrix)(gmatrix *g,
   dloss_pt, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

double sgd_gmatrix(gmatrix *g,
   dloss_pt, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

double scd_gmatrix(gmatrix *g,
   dloss_pt, loss_pt, predict_pt,
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

