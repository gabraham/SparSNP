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

double logloss_pt(double d, dtype y)
{
   return -(double)y * d + log(1 + exp(d));
}

double logloss(double *d, dtype *y, int n)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += logloss_pt(d[i], y[i]);
   return loss;
}

void logdloss_pt(dtype *x, double d, dtype y, int p, double* grad)
{
   int j;
   double pr = exp(d);
   for(j = 0 ; j < p ; j++)
      grad[j] = x[j] * (pr / (1 + pr) - (double)y);
}

void logd2loss_pt(dtype *x, double d, int p, double* d2)
{
   int j;
   double pr = exp(d);
   for(j = 0 ; j < p ; j++)
      d2[j] = logd2loss_pt_j(x[j], pr / (1 + pr));
}

double logd2loss_pt_j(dtype x, double P)
{
   return P * (1 - P);
}

double l2loss_pt(double d, dtype y)
{
   return pow(y - d, 2);
}

double l2loss(double *d, dtype *y, int n)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += l2loss_pt(d[i], y[i]);
   return loss;
}

void l2dloss_pt(dtype *x, double d, dtype y, int p, double* grad)
{
   int j;
   for(j = 0 ; j < p ; j++)
      grad[j] = x[j] * (d - (double)y);
}

double l2phi1(double lp)
{
   return lp;
}

double l2phi2(double lp)
{
   return 1;
}

double logphi1(double lp)
{
   return plogis(lp);
}

double logphi2(double lp)
{
   double p = logphi1(lp);
   return p * (1 - p);
}

double l2inv(double lp)
{
   return lp;
}

/* same as logit */
double loginv(double lp)
{
   return log(lp / (1 - lp));
}

