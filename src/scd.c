#include "sgd.h"

double soft_threshold(double beta, double gamma)
{
   return sign(beta) * fmax(fabs(beta) - gamma, 0);
}

/* Stochastic coordinate descent */
double scd_gmatrix(gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   loss_pt loss_pt_func,    /* loss for one sample */
   predict_pt predict_pt_func, /* prediction for one sample */
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc)
{
   int i, j, k;
   int epoch = 1;
   sample sm;
   double d, s;
   double loss = 0;
   double *lp;
   double y2 = 0;

   if(!g->inmemory)
   {
      fprintf(stderr, "scd_gmatrix doesn't support disk based gmatrix yet");
      return FAILURE;
   }

   sample_init(&sm, g->n);

   while(epoch <= maxepoch)
   {
      loss = 0;
      CALLOCTEST(lp, g->n, sizeof(double));

      for(j = 0 ; j < g->p + 1; j++)
      {
	 d = 0;

	 for(i = 0 ; i < g->n ; i++)
	 {
	    y2 = 0;
	    for(k = 0 ; k < g->p + 1 ; k++)
	       y2 += g->x[i][k] * beta[k];

	    d += g->x[i][j] * (g->y[i] - y2);
	 }

	 /*beta[j] = soft_threshold(beta[j] + d, lambda1) / (1 + lambda2);*/
	 beta[j] += d;

      }

      for(i = 0 ; i < g->n ; i++)
      {
	 for(j = 0 ; j < g->p + 1 ; j++)
	    lp[i] += g->x[i][j] * beta[j];

	 loss += loss_pt_func(lp[i], g->y[i]) / g->n;
      }

      printf("Epoch %d  training loss: %.5f\n", epoch, loss);

      free(lp);
      epoch++;
   }


   return SUCCESS;
}

