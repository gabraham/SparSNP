#include <math.h>

#include "loss.h"

double plogis(double x)
{
   return 1 / (1 + exp(-x));
}

double dotprod(double *a, double *b, int m)
{
   int i;
   double s = 0;
   for(i = 0 ; i < m ; i++)
      s += a[i] * b[i];
   return s;
}

double logloss_pt(double *x, double *beta, int y, int p)
{
   double d = fmin(dotprod(x, beta, p), MAXPROD);
   return -(double)y * d + log(1 + exp(d));
}

double logloss(double **x, double *beta, int *y, int n, int p)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += logloss_pt(x[i], beta, y[i], p);
   return loss;
}

void logdloss(double *x, double *beta, int y, int p, double* grad)
{
   int i;
   double pr = exp(fmin(dotprod(x, beta, p), MAXPROD));
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (pr / (1 + pr) - (double)y);
}

