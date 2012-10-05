/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <math.h>
#include <stdlib.h>
#include "common.h"

/* Used to truncate exp before infinity, which occurs at ~709 */ 
#define MAXPROD 700

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801

/*inline double exponential(double y);*/

double plogis(double);

double linearphi1(double lp);
double linearphi2(double lp);
double logphi1(double lp);
double logphi2(double lp);
double linearinv(double lp);
double loginv(double lp);
double sqrhingeinv(double lp);

