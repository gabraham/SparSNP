/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <math.h>
#include <stdio.h>
#include "loss.h"

double linearphi1(double lp)
{
   return lp;
}

double linearphi2(double lp)
{
   return 1;
}

double logphi1(double lp)
{
   return 1 / (1 + exp(-lp));
}

double logphi2(double p)
{
   return p * (1 - p);
}

/* linear link function */
double linearinv(double lp)
{
   return lp;
}

/* logistic link function, i.e., logit */
double loginv(double lp)
{
   return log(lp / (1 - lp));
}

double sqrhingeinv(double lp)
{
   return lp;
}

