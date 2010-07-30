#include <math.h>
#include <stdio.h>
#include "loss.h"

#ifndef EXP
#define EXP exp
#endif

/*double plogis(double x)
{
   return 1 / (1 + exp(-x));
}*/

/* Schraudolph 1998, ``A Fast, Compact Approximation of the Exponential Function''
   and
   Cawley 2000, ``On a Fast, Compact Approximation of the Exponential
   Function''
   
   Accuracy depends on y being not too large

 */
inline double exponential(double y)
{
    union
    {
        double d;
#ifdef LITTLE_ENDIAN
        struct { int j, i; } n;
#else
        struct { int i, j; } n;
#endif
    }
    _eco;

    _eco.n.i = (int)(EXP_A*(y)) + (1072693248 - EXP_C);
    _eco.n.j = 0;

    return _eco.d;
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
   double z = y - d;
   return z * z;
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

inline double l2phi1(double lp)
{
   return lp;
}

inline double l2phi2(double lp)
{
   return 1;
}

inline double logphi1(double lp)
{
   return 1 / (1 + EXP(-lp));
}

inline double logphi2(double p)
{
   /*double p = logphi1(lp);*/
   return p * (1 - p);
}

inline double l2inv(double lp)
{
   return lp;
}

/* same as logit */
inline double loginv(double lp)
{
   return log(lp / (1 - lp));
}

