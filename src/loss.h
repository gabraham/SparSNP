#include <math.h>
#include <stdlib.h>
#include "common.h"

/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801

/*inline double exponential(double y);*/

double plogis(double);
double dotprod(dtype *, double *, int);

double logloss_pt(double, dtype);
double logloss(double *, dtype *, int);
void logdloss_pt(dtype *, double, dtype, int, double*);
void logd2loss_pt(dtype *, double, int, double*);
double logd2loss_pt_j(dtype, double);

double linearloss_pt(double, dtype);
double linearloss(double *, dtype *, int);
void lineardloss_pt(dtype *, double, dtype, int, double*);
void lineard2loss_pt(dtype *, double, int, double*);
double lineard2loss_pt_j(dtype, double);

double linearphi1(double lp);
double linearphi2(double lp);
double logphi1(double lp);
double logphi2(double lp);
double linearinv(double lp);
double loginv(double lp);
double sqrhingeinv(double lp);

