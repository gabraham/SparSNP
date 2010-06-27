#include "sgd.h"

short convergetest(double a, double b, double threshold)
{
   if(a == 0 && b == 0)
      return TRUE;

   return (fabs(a - b) / (fabs(a) + fabs(b))) < threshold;
}

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
   short *converged2 = NULL;
   int numconverged = 0;
   double relerr;
   double grad = 0;
   double d2 = 0;
   double *lp = NULL;
   double pr;
   double s;
   sample sm;
   double truncl = log((1 - trunc) / trunc);
   double lambda1max = 0;
   short done = FALSE;
   int allconverged = 0;

   sample_init(&sm, g->inmemory, g->n);
   MALLOCTEST(sm.x, sizeof(dtype) * g->n)

   CALLOCTEST(converged, g->p + 1, sizeof(short));
   /*CALLOCTEST(converged2, g->p + 1, sizeof(short));*/
   /*CALLOCTEST(grad, g->p + 1, sizeof(double));*/
   CALLOCTEST(lp, g->n, sizeof(double));

   while(epoch <= maxepoch) /* && numconverged < g->p + 1) */
   /*while(TRUE)*/
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 g->nextcol(g, &sm);
	 
	 if(converged[j])
	    continue;

	 grad = 0;
	 d2 = 0;

	 /* compute gradient */
	 for(i = 0 ; i < g->n ; i++)
	 {
	    if(sm.x[i] == 0)
	       continue;

	    pr = predict_pt_func(lp[i]);
	    grad += sm.x[i] * (pr - g->y[i]);
	    d2 += d2loss_pt_j_func(sm.x[i], pr);
	 }

	 /* don't move if 2nd derivative is zero */
	 s = 0;
	 if(d2 != 0)
	    s = grad / d2;

	 /* don't penalise intercept */
	 if(j == 0)
	    beta_new = beta[j] - s;
	 else
	    beta_new = soft_threshold(beta[j] - s, lambda1) / (1 + lambda2);

	 /* find smallest lambda1 that makes all coefficients zero, by finding
	  * the largest beta[j] - s */
	 if(lambda1max < fabs(beta[j] - s))
	    lambda1max = fabs(beta[j] - s);

	 /* check for convergence */
	 if(epoch > 1 && convergetest(beta[j], beta_new, threshold))
	 {
	    converged[j] = TRUE;
	    numconverged++;
	 }

	 /* update linear predictor */
	 for(i = 0 ; i < g->n ; i++)
	    if(sm.x[i] != 0)
	       lp[i] += sm.x[i] * (beta_new - beta[j]);

	 /* clip very large coefs to prevent divergence */
	 beta[j] = fmin(fmax(beta_new, -truncl), truncl);
      }

      if(verbose)
      {
	 loss = 0;
      	 for(i = 0 ; i < g->n ; i++)
      	    loss += loss_pt_func(lp[i], g->y[i]) / g->n;
      	 printf("Epoch %d  training loss: %.5f  converged: %d\n", epoch, loss,
	 numconverged);
      }

      if(numconverged == g->p + 1)
      {
	 /*if(allconverged == 1)
	    allconverged = 0;*/
	 printf("all converged\n");
	 for(j = 0 ; j < g->p + 1 ; j++)
	    converged[j] = FALSE;
	 numconverged = 0;
      }


      /*printf("lambda1max: %.5f\n", lambda1max);*/

      

      epoch++;
   }

   free(converged);
   /*free(converged2);*/
   free(lp);
   sample_free(&sm);
   free(sm.x);

   return lambda1max;
}

