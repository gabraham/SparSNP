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

double log_loss_pt(double, dtype);
double log_loss(double *, dtype *, int);
void log_dloss_pt(dtype *, double, dtype, int, double*);
void log_d2loss_pt(dtype *, double, int, double*);
double log_d2loss_pt_j(dtype, double);

double linear_loss_pt(double, dtype);
double linear_loss(double *, dtype *, int);
void linear_dloss_pt(dtype *, double, dtype, int, double*);
void linear_d2loss_pt(dtype *, double, int, double*);
double linear_d2loss_pt_j(dtype, double);

double sqrhinge_loss(double *, dtype *, int);
double sqrhinge_loss_pt(double, dtype);

double linearphi1(double lp);
double linearphi2(double lp);
double logphi1(double lp);
double logphi2(double lp);
double linearinv(double lp);
double loginv(double lp);
double sqrhingeinv(double lp);

typedef double (*loss)(double*, dtype*, int);
typedef double (*loss_pt)(double, dtype);

