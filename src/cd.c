#include "sgd.h"

short convergetest(double a, double b, double threshold)
{
   /* absolute convergence */
   if(fabs(a) <= ZERO_THRESH && fabs(b) <= ZERO_THRESH)
      return TRUE;

   /* relative convergence */
   return (fabs(a - b) / (fabs(a) + fabs(b))) < threshold;
}

/* Find smallest lambda1 that makes all coefficients zero (except the intercept)
 */
double get_lambda1max_gmatrix(gmatrix *g,
      dloss_pt dloss_pt_func,        /* gradient */
      d2loss_pt d2loss_pt_func,        /* 2nd deriv */
      d2loss_pt_j d2loss_pt_j_func,        /* 2nd deriv wrt beta_j */
      loss_pt loss_pt_func,    /* loss for one sample */
      predict_pt predict_pt_func) /* prediction for one sample */
{
   int i, j;
   double *lp = NULL;
   double grad, d2, s, pr, zmax = 0;
   sample sm;

   CALLOCTEST(lp, g->n, sizeof(double))
   sample_init(&sm, g->inmemory, g->n);
   MALLOCTEST(sm.x, sizeof(dtype) * g->n)

   for(j = 0 ; j < g->p + 1; j++)
   {
      g->nextcol(g, &sm);

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
	 s = -grad / d2;

      /* find smallest lambda1 that makes all coefficients zero, by finding
       * the largest z, but evaluate the intercept first because
       * it's not penalised.
       */
      if(j == 0)
      {
	 for(i = 0 ; i < g->n ; i++)
	    lp[i] = sm.x[i] * s;
      }
      else if(zmax < fabs(s))
	 zmax = fabs(s);
   } 

   free(lp);
   free(sm.x);

   return zmax;
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
   int i, j;
   int epoch = 1;
   double loss = 0;
   double beta_new;
   short *converged = NULL;
   int numconverged = 0;
   double grad = 0;
   double d2 = 0;
   double *lp = NULL;
   double pr;
   double s;
   sample sm;
   double truncl = log((1 - trunc) / trunc);
   int allconverged = 0;
   int zeros = 0;

   sample_init(&sm, g->inmemory, g->n);
   MALLOCTEST(sm.x, sizeof(dtype) * g->n)

   CALLOCTEST(converged, g->p + 1, sizeof(short));
   CALLOCTEST(lp, g->n, sizeof(double));

   while(epoch <= maxepoch)
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
	    /* skip zeros, they don't change the linear predictor */
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
	    beta_new = soft_threshold(beta[j] - s, lambda1)
		  / (1 + lambda2);

	 /*if(g->p + 1 - numconverged < 100 && beta_new != 0)
	    printf("[%d] %.20f %.20f %.20f %.20f\n", j, beta[j], s, beta[j] - s, beta_new);*/
	    

	 /* check for convergence */
	 if(epoch > 1 && convergetest(beta[j], beta_new, threshold))
	 {
	    converged[j] = TRUE;
	    numconverged++;
	 }

	    

	 /* clip very large coefs to prevent divergence */
	 beta_new = fmin(fmax(beta_new, -truncl), truncl);

	 /* update linear predictor */
	 for(i = 0 ; i < g->n ; i++)
	    if(sm.x[i] != 0)
	       lp[i] += sm.x[i] * (beta_new - beta[j]);

	 beta[j] = beta_new;
      }

      if(verbose)
      {
	 loss = 0;
	 for(i = 0 ; i < g->n ; i++)
	    loss += loss_pt_func(lp[i], g->y[i]) / g->n;
	 
	 zeros = 0;
	 /* don't count intercept */
	 if(!converged[0])
	    printf("not converged: %d (%.20f)\n", 0, beta[0]);
	 for(j = 1 ; j < g->p + 1 ; j++)
	 {
	    if(beta[j] == 0)
	       zeros++;
	    else
	       printf("nonzero: %d (%.20f)\n", j, beta[j]);

	    if(!converged[j])
	       printf("not converged: %d (%.20f)\n", j, beta[j]);
	 }

	 printf("Epoch %d  training loss: %.5f  converged: %d\
  zeros: %d  non-zeros: %d\n",
   epoch, loss, numconverged, zeros, g->p - zeros);

      }

      if(numconverged == g->p + 1)
      {
	 printf("all converged\n");
	 allconverged++;

	 /* converged twice in a row, no need to continue */
	 if(allconverged == 2)
	 {
	    printf("terminating\n");
	    break;
	 }

	 for(j = 0 ; j < g->p + 1 ; j++)
	    converged[j] = FALSE;
	 numconverged = 0;
      }

      epoch++;
   }

   free(converged);
   free(lp);
   sample_free(&sm);
   free(sm.x);

   return SUCCESS;
}

