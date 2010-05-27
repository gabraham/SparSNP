#include <math.h>

#include "loss.h"

double plogis(double x)
{
   return 1 / (1 + exp(-x));
}

double dotprod(dtype *a, double *b, int m)
{
   int i;
   double s = 0;
   for(i = 0 ; i < m ; i++)
      s += a[i] * b[i];
   return fmin(fmax(s, -MAXPROD), MAXPROD);
}

double logloss_pt(dtype *x, double *beta, dtype y, int p)
{
   double d = dotprod(x, beta, p);
   return -(double)y * d + log(1 + exp(d));
}

double logloss(dtype **x, double *beta, dtype *y, int n, int p)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += logloss_pt(x[i], beta, y[i], p);
   return loss;
}

void logdloss(dtype *x, double *beta, dtype y, int p, double* grad)
{
   int i;
   double pr = exp(dotprod(x, beta, p));
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (pr / (1 + pr) - (double)y);
}

double l2loss_pt(dtype *x, double *beta, dtype y, int p)
{
   double d = dotprod(x, beta, p);
   return pow(y - d, 2);
}

double l2loss(dtype **x, double *beta, dtype *y, int n, int p)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += l2loss_pt(x[i], beta, y[i], p);
   return loss;
}

void l2dloss(dtype *x, double *beta, dtype y, int p, double* grad)
{
   int i;
   double pr = dotprod(x, beta, p);
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (pr - (double)y);
}


