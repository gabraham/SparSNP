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
   return x * x * P * (1 - P);
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

void l2d2loss_pt(dtype *x, double d, int p, double* d2)
{
   int j;
   for(j = 0 ; j < p ; j++)
      d2[j] = l2d2loss_pt_j(x[j], 0);
}

double l2d2loss_pt_j(dtype x, double pr)
{
   return x * x;
}

double predict_logloss_pt(double d)
{
   return 1 / (1 + exp(-d));
}

double predict_l2loss_pt(double d)
{
   return d;
}

