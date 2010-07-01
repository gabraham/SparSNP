#include "sgd.h"

/* gradient descent */
double gd_gmatrix(
   gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   d2loss_pt d2loss_pt_func,      /* 2nd deriv */
   d2loss_pt_j d2loss_pt_j_func,      /* 2nd deriv wrt beta_j */
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
   double loss;
   sample sm;

   MALLOCTEST(grad, sizeof(double) * (g->p + 1))
   CALLOCTEST(gradsum, g->p + 1, sizeof(double))
   /*MALLOCTEST(d2, sizeof(double) * (g->p + 1))*/
   /*CALLOCTEST(d2sum, g->p + 1, sizeof(double))*/
   
   if(!sample_init(&sm, g->inmemory, g->p + 1))
      return FAILURE;

   while(epoch <= maxepoch)
   {
      /* compute gradient */
      for(i = 0 ; i < g->n ; i++)
      { 
	 g->nextrow(g, &sm);
	 
	 dp = dotprod(sm.x, beta, g->p + 1);
	 dloss_pt_func(sm.x, dp, sm.y, g->p + 1, grad);
	 /*d2loss_pt_func(sm.x, dp, g->p + 1, d2);*/

	 for(j = 0 ; j < g->p + 1 ; j++)
	 {
	    gradsum[j] += grad[j];
	    /* d2sum[j] += d2[j]; */
	 }
      }

      /* do step */
      /* TODO: don't penalise intercept */
      for(j = 0 ; j < g->p + 1 ; j++)
      {
	 /*s = 0;
	 if(d2sum[j] > 0)
	    s = gradsum[j] / d2sum[j];
	 beta[j] = soft_threshold(beta[j] - s, lambda1) / (1 + lambda2);*/
	 beta[j] = soft_threshold(beta[j] - stepsize * gradsum[j], lambda1) 
	       / (1 + lambda2);
	 gradsum[j] = 0;
	 /*d2sum[j] = 0;*/
      }

      if(verbose)
      {
	 loss = 0;
      	 for(i = 0 ; i < g->n ; i++)
      	 {
	    g->nextrow(g, &sm);
      	    dp = dotprod(sm.x, beta, g->p + 1);
      	    loss += loss_pt_func(dp, sm.y) / g->n;
      	 }
      	 printf("Epoch %d loss=%.5f\n", epoch, loss);
      }

      epoch++;
   }

   free(grad);
   free(gradsum);
   /*free(d2);
   free(d2sum);*/
   return SUCCESS;
}

