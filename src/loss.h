#include <math.h>
#include <stdlib.h>
#include "common.h"
#include "gmatrix.h"

/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700


void predict_logloss_gmatrix(gmatrix *, double *, double *, int *);
double predict_logloss_pt_gmatrix(sample *, double *, double *, double *, int);
double predict_logloss_pt(double);

void predict_l2loss_gmatrix(gmatrix *, double *, double *, int *);
double predict_l2loss_pt_gmatrix(sample *, double *, double *, double *, int);
double predict_l2loss_pt(double);

double plogis(double);
double dotprod(dtype *, double *, int);

double logloss_pt(double, dtype);
double logloss(double *, dtype *, int);
void logdloss_pt(dtype *, double, dtype, int, double*);

double l2loss_pt(double, dtype);
double l2loss(double *, dtype *, int);
void l2dloss_pt(dtype *, double, dtype, int, double*);
