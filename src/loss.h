/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700

#include "common.h"

double plogis(double);
double dotprod(dtype *, double *, int);
double logloss_pt(dtype *, double *, dtype , int);
double logloss(dtype **, double *, dtype *, int, int);
void logdloss(dtype *, double *, dtype, int, double*);

