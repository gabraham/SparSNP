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
   int i;
   double pr = exp(d);
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (pr / (1 + pr) - (double)y);
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
   int i;
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * (d - (double)y);
}

double predict_logloss_pt(double d)
{
   return 1 / (1 + exp(-d));
}

double predict_logloss_pt_gmatrix(sample *s, double *beta, int p)
{
   /*int i = 0;
   dtype *x;
   double yhat = 0;

   MALLOCTEST(x, sizeof(dtype) * (p + 1))

   x[0] = 1;
   for(i = 1 ; i < p + 1; i++)
      x[i] = (s->x[i] - mean[i]) / sd[i];
   yhat = predict_logloss_pt(dotprod(s->x, beta, p + 1));

   free(x);
   return yhat; */
   return predict_logloss_pt(dotprod(s->x, beta, p + 1));
}

void predict_logloss_gmatrix(gmatrix *g, double *beta,
      double *yhat, int *trainf)
{
   int i, k;
   sample sm;

   sample_init(&sm, g->inmemory, g->p);

   k = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      g->nextrow(g, &sm);
      if(trainf[i])
      {
	 yhat[k] = predict_logloss_pt_gmatrix(&sm, beta, g->p);
	 k++;
      }
   } 

   sample_free(&sm);
}

double predict_l2loss_pt(double d)
{
   return d;
}

double predict_l2loss_pt_gmatrix(sample *s, double *beta, int p)
{
   /*int i = 0;
   dtype *x;
   double yhat = 0;

   MALLOCTEST(x, sizeof(dtype) * p)

   x[0] = 1;
   for(i = 0 ; i < p - 1 ; i++)
      x[i+1] = (s->x[i] - mean[i]) / sd[i];
   yhat = dotprod(x, beta, p);

   free(x);
   return yhat;*/
   return dotprod(s->x, beta, p + 1);
}

void predict_l2loss_gmatrix(gmatrix *g, double *beta, double *yhat, int *trainf)
{
   int i, k;
   sample sm;

   sample_init(&sm, g->inmemory, g->p);

   k = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      g->nextrow(g, &sm);
      if(trainf[i])
      {
	 yhat[k] = predict_l2loss_pt_gmatrix(&sm, beta, g->p);
	 k++;
      }
   } 

   sample_free(&sm);
}

