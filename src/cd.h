#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "gmatrix.h"


#define MAXLP 20

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
typedef void (*dloss_pt)(dtype *, double, dtype, int, double *);
typedef void (*d2loss_pt)(dtype *, double, int, double *);
typedef double (*d2loss_pt_j)(dtype, double);
typedef void (*predict_gmatrix)(gmatrix *, double *, double *, int *);

typedef double (*optim_gmatrix)(gmatrix *g,
   dloss_pt, d2loss_pt, d2loss_pt_j, loss_pt, predict_pt,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

double cd_gmatrix(gmatrix *g,
   dloss_pt, d2loss_pt, d2loss_pt_j, loss_pt, predict_pt,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc);

double get_lambda1max_gmatrix(gmatrix *g,
      dloss_pt dloss_pt_func,        /* gradient */
      d2loss_pt d2loss_pt_func,        /* 2nd deriv */
      d2loss_pt_j d2loss_pt_j_func,        /* 2nd deriv wrt beta_j */
      loss_pt loss_pt_func,    /* loss for one sample */
      predict_pt predict_pt_func /* prediction for one sample */
);

