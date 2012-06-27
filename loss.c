/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <math.h>
#include <stdio.h>
#include "loss.h"

double dotprod(dtype *a, double *b, int m)
{
   int i;
   double s = 0;
   for(i = 0 ; i < m ; i++)
      s += a[i] * b[i];
   return fmin(fmax(s, -MAXPROD), MAXPROD);
}

double log_loss_pt(double d, dtype y)
{
   return -(double)y * d + log(1 + exp(d));
}

double log_loss(double *d, dtype *y, int n)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += log_loss_pt(d[i], y[i]);
   return loss;
}

void log_dloss_pt(dtype *x, double d, dtype y, int p, double* grad)
{
   int j;
   double pr = exp(d);
   for(j = 0 ; j < p ; j++)
      grad[j] = x[j] * (pr / (1 + pr) - (double)y);
}

void log_d2loss_pt(dtype *x, double d, int p, double* d2)
{
   int j;
   double pr = exp(d);
   for(j = 0 ; j < p ; j++)
      d2[j] = log_d2loss_pt_j(x[j], pr / (1 + pr));
}

double log_d2loss_pt_j(dtype x, double P)
{
   return P * (1 - P);
}

double linear_loss_pt(double d, dtype y)
{
   double z = y - d;
   return z * z;
}

double linear_loss(double *d, dtype *y, int n)
{
   int i;
   double loss = 0;
   for(i = n - 1 ; i >= 0 ; --i)
      loss += linear_loss_pt(d[i], y[i]);
   return loss;
}

void linear_dloss_pt(dtype *x, double d, dtype y, int p, double* grad)
{
   int j;
   for(j = 0 ; j < p ; j++)
      grad[j] = x[j] * (d - (double)y);
}

double sqrhinge_loss_pt(double d, dtype y)
{
   double z = 1 - d * y;
   return (z > 0) ? z * z : 0.0;
}

double sqrhinge_loss(double *d, dtype *y, int n)
{
   int i;
   double loss = 0;
   for(i = n - 1 ; i >= 0 ; --i)
      loss += sqrhinge_loss_pt(d[i], y[i]);
   return loss;
}

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
   /*double p = logphi1(lp);*/
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

