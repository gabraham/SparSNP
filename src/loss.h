#include <math.h>
#include <stdlib.h>
#include "common.h"
#include "gmatrix.h"

/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700


void predict_logloss(gmatrix *, double *, double *, int *);
double predict_logloss_pt(sample *, double *, double *, double *, int);

void predict_l2loss(gmatrix *, double *, double *, int *);
double predict_l2loss_pt(sample *, double *, double *, double *, int);

double plogis(double);
double dotprod(dtype *, double *, int);

double logloss_pt(dtype *, double *, dtype , int);
double logloss(dtype **, double *, dtype *, int, int);
void logdloss(dtype *, double *, dtype, int, double*);

double l2loss_pt(dtype *, double *, dtype , int);
double l2loss(dtype **, double *, dtype *, int, int);
void l2dloss(dtype *, double *, dtype, int, double*);
