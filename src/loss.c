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


double predict_logloss_pt(sample *s, double *beta, double *mean, double *sd, int p)
{
   int i = 0;
   dtype *x = malloc(sizeof(dtype) * p);
   double yhat = 0;

   x[0] = 1;
   for(i = 0 ; i < p - 1 ; i++)
      x[i+1] = (s->x[i] - mean[i]) / sd[i];
   yhat = 1 / (1 + exp(-dotprod(x, beta, p)));

   free(x);
   return yhat;
}

void predict_logloss(gmatrix *g, double *beta, double *yhat, int *trainf)
{
   int i, k;
   sample sm;

   sample_init(&sm, g->p);

   k = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      g->nextrow(g, &sm);
      if(trainf[i])
      {
	 yhat[k] = predict_logloss_pt(&sm, beta, g->mean, g->sd, g->p + 1);
	 k++;
      }
   } 

   sample_free(&sm);
}

double predict_l2loss_pt(sample *s, double *beta, double *mean, double *sd, int p)
{
   int i = 0;
   dtype *x = malloc(sizeof(dtype) * p);
   double yhat = 0;

   x[0] = 1;
   for(i = 0 ; i < p - 1 ; i++)
      x[i+1] = (s->x[i] - mean[i]) / sd[i];
   yhat = dotprod(x, beta, p);

   free(x);
   return yhat;
}

void predict_l2loss(gmatrix *g, double *beta, double *yhat, int *trainf)
{
   int i, k;
   sample sm;

   sample_init(&sm, g->p);

   k = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      g->nextrow(g, &sm);
      if(trainf[i])
      {
	 yhat[k] = predict_l2loss_pt(&sm, beta, g->mean, g->sd, g->p + 1);
	 k++;
      }
   } 

   sample_free(&sm);
}

