#include <math.h>
#include <stdlib.h>
#include "common.h"

/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700

double predict_logloss_pt(double);
double predict_l2loss_pt(double);

double plogis(double);
double dotprod(dtype *, double *, int);

double logloss_pt(double, dtype);
double logloss(double *, dtype *, int);
void logdloss_pt(dtype *, double, dtype, int, double*);
void logd2loss_pt(dtype *, double, int, double*);
double logd2loss_pt_j(dtype, double);

double l2loss_pt(double, dtype);
double l2loss(double *, dtype *, int);
void l2dloss_pt(dtype *, double, dtype, int, double*);
void l2d2loss_pt(dtype *, double, int, double*);
double l2d2loss_pt_j(dtype, double);

inline double l2phi1(double lp);
inline double l2phi2(double lp);
inline double logphi1(double lp);
inline double logphi2(double lp);
inline double l2inv(double lp);
inline double loginv(double lp);

