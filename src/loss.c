#include <math.h>

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
   return -(double)y * dotprod(x, beta, p) + log(1 + exp(dotprod(x, beta, p)));
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
   double pr = exp(dotprod(x, beta, p));
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (pr / (1 + pr) - (double)y);
}

