#include "sgd.h"

double soft_threshold(double beta, double step)
{
   return sign(beta) * fmax(fabs(beta) - step, 0);
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
   double q = 0;

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

      for(i = 0 ; i < g->p + 1; i++)
      {
	 g->nextcol(g, &sm);
	 q = 0;
	 for(j = 0 ; j < g->n ; j++)
	 {
	    for(k = 0 ; k <- g->p + 1 ; k++)
	    {
	       if(k != i)
		  q += sm.x[j] * beta[k];
	    }

	    d += sm.x[j] * (g->y[j] - q);
	 }

	 beta[i] = soft_threshold(d, lambda1) / (1 + lambda2);
      }

      for(i = 0 ; i < g->p + 1 ; i++)
      {
	 for(j = 0 ; j < g->n ; j++)
	 {
	    g->nextcol(g, &sm);
	    lp[j] += beta[i] * sm.x[j];
	 }
      }

     /* for(j = 0 ; j < g->n ; j++)
	 loss += loss_pt_func(lp[j], g->y[j]) / g->n;

      printf("Epoch %d  training loss: %.5f\n", epoch, loss); */

      free(lp);
      epoch++;
   }


   return SUCCESS;
}

