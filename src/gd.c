#include "sgd.h"

/* gradient descent */
double gd_gmatrix(gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   loss_pt loss_pt_func,    /* loss for one sample */
   predict_pt predict_pt_func, /* prediction for one sample */
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc)
{
   int epoch = 1, i, j;
   double *grad;
   double stepsize = maxstepsize;
   double dp = 0;
   double *gradsum = NULL;

   MALLOCTEST(grad, sizeof(double) * (g->p + 1))
   CALLOCTEST(gradsum, g->p + 1, sizeof(double))
   
   while(epoch <= maxepoch)
   {
      /* compute gradient */
      for(i = 0 ; i < g->n ; i++)
      { 
	 dp = dotprod(g->x[i], beta, g->p + 1);
	 dloss_pt_func(g->x[i], dp, g->y[i], g->p + 1, grad);
	 for(j = 0 ; j < g->p + 1 ; j++)
	    gradsum[j] += grad[j];
      }

      /* do step */
      for(j = 0 ; j < g->p + 1 ; j++)
      {
	 beta[j] -= stepsize * gradsum[j];
	 gradsum[j] = 0;
      }

      epoch++;
   }

   free(grad);
   free(gradsum);
   return SUCCESS;
}

