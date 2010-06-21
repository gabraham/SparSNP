#include "sgd.h"

/* coordinate descent */
double cd_gmatrix(gmatrix *g,
   dloss_pt dloss_pt_func,        /* gradient */
   d2loss_pt d2loss_pt_func,        /* 2nd deriv */
   d2loss_pt_j d2loss_pt_j_func,        /* 2nd deriv wrt beta_j */
   loss_pt loss_pt_func,    /* loss for one sample */
   predict_pt predict_pt_func, /* prediction for one sample */
   double maxstepsize,
   int maxepoch, double *beta, double lambda1, double lambda2,
   double threshold, int verbose, int *trainf, double trunc)
{
   int i, j, k;
   int epoch = 1;
   double loss = 0;
   double beta_new;
   short *converged = NULL;
   int numconverged = 0;
   double relerr;
   double grad = 0;
   double d2 = 0;
   double *lp = NULL;

   if(!g->inmemory)
   {
      fprintf(stderr, "cd_gmatrix doesn't support disk based gmatrix yet");
      return FAILURE;
   }

   /*sample_init(&sm, g->n);*/

   CALLOCTEST(converged, g->p + 1, sizeof(short));
   /*CALLOCTEST(grad, g->p + 1, sizeof(double));*/
   CALLOCTEST(lp, g->n, sizeof(double));

   while(epoch <= maxepoch && numconverged < g->p + 1)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 if(converged[j])
	    continue;

	 grad = 0;
	 d2 = 0;

	 /* compute gradient */
	 for(i = 0 ; i < g->n ; i++)
	 {
	    grad += g->x[i][j] * (lp[i] - g->y[i]);
	    d2 += pow(g->x[i][j], 2.0);
	 }

	 /* TODO: don't penalise intercept */
	 /*beta_new = soft_threshold(beta[j] + grad / d2, lambda1) / (1 +
	  * lambda2);*/

	 if(d2 != 0)
	    beta_new = beta[j] - grad / d2;
	 else
	    beta_new = beta[j];

	 /* check for convergence */
	 if(epoch > 1)
	 {
	    relerr = fabs(beta[j] - beta_new) / (fabs(beta[j]) + fabs(beta_new));
	    if(relerr < threshold)
	    {
	       converged[j] = TRUE;
	       numconverged++;
	    }
	 }

	 /* update linear predictor */
	 for(i = 0 ; i < g->n ; i++)
	    lp[i] += g->x[i][j] * (beta_new - beta[j]);

	 beta[j] = beta_new;
      }

      if(verbose)
      {
	 loss = 0;
      	 for(i = 0 ; i < g->n ; i++)
      	    loss += loss_pt_func(lp[i], g->y[i]) / g->n;
      	 printf("Epoch %d  training loss: %.5f\n", epoch, loss);
      }

      epoch++;
   }

   free(converged);
   free(lp);

   return SUCCESS;
}

