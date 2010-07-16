#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "gmatrix.h"


#define MAXLP 20

typedef double (*loss_pt)(double, dtype);
typedef double (*predict_pt)(double);
/*typedef void (*dloss_pt)(dtype *, double, dtype, int, double *);
typedef void (*d2loss_pt)(dtype *, double, int, double *);
typedef double (*d2loss_pt_j)(dtype, double);
typedef void (*predict_gmatrix)(gmatrix *, double *, double *, int *);*/
typedef double (*phi1)(double);
typedef double (*phi2)(double);

int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
      int maxepoch, double *beta, double lambda1, double lambda2,
      double threshold, int verbose, int *trainf, double trunc);

int cd_gmatrix2(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,
      int maxepoch,
      double *beta,
      double lambda1);


double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func
);

