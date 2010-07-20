#include "common.h"
#include "cd.h"

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
double get_lambda1max_gmatrix(
      gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func)
{
   int i, j;
   double *lp = NULL;
   double grad, d2, s, zmax = 0;
   sample sm;

   CALLOCTEST(lp, g->n, sizeof(double))
   if(!sample_init(&sm, g->n))
      return FAILURE;

   for(j = 0 ; j < g->p + 1; j++)
   {
      gmatrix_disk_nextcol(g, &sm);

      grad = 0;
      d2 = 0;

      /* compute gradient */
      for(i = 0 ; i < g->n ; i++)
      {
	 if(sm.x[i] == 0)
	    continue;

	 grad += sm.x[i] * (phi1_func(lp[i]) - g->y[i]);
	 d2 += sm.x[i] * sm.x[i] * phi2_func(lp[i]);
      }

      /* don't move if 2nd derivative is zero */
      s = 0;
      if(d2 != 0)
	 s = -grad / d2;

      /* find smallest lambda1 that makes all coefficients zero, by finding
       * the largest z, but evaluate the intercept first because
       * it's not penalised. */
      if(j == 0)
	 for(i = 0 ; i < g->n ; i++)
	    lp[i] = sm.x[i] * s;
      else if(zmax < fabs(s))
	 zmax = fabs(s);
   } 

   free(lp);
   free(sm.x);

   return zmax;
}

/*int cd_gmatrix2(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,
      int maxepoch,
      double *beta,
      double lambda1)
{
   int i, j;
   int iter;
   double beta_old = 0, grad, d2, s;
   double *lp = NULL;
   sample sm;

   CALLOCTEST(lp, g->n, sizeof(double))

   if(!sample_init(&sm, g->n))
      return FAILURE;

   for(iter = 0 ; iter < maxepoch ; iter++)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 beta_old = beta[j];
	 gmatrix_disk_nextcol(g, &sm);
	 grad = d2 = 0;
	 for(i = 0 ; i < g->n ; i++)
	 {
	    grad += sm.x[i] * (phi1_func(lp[i]) - g->y[i]);
	    d2 += sm.x[i] * sm.x[i] * phi2_func(lp[i]);
	 }
	 
	 s = 0;
	 if(d2 != 0)
	    s = grad / d2;

	 if(j == 0)
	    beta[j] = soft_threshold(beta[j] - s, lambda1);
	 else
	    beta[j] -= s;

	 for(i = 0 ; i < g->n ; i++)
	    lp[i] += sm.x[i] * (beta[j] - beta_old);
      }
   }

   free(lp);
   sample_free(&sm);

   return SUCCESS;
}*/


/* coordinate descent */
int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
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
   double s;
   sample sm;
   double truncl = log((1 - trunc) / trunc);
   int allconverged = 0;
   int zeros = 0;
   const int CONVERGED = 2;

   if(!sample_init(&sm, g->n))
      return FAILURE;

   CALLOCTEST(converged, g->p + 1, sizeof(short));
   CALLOCTEST(lp, g->n, sizeof(double));

   while(epoch <= maxepoch)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 gmatrix_disk_nextcol(g, &sm);

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

	    grad += sm.x[i] * (phi1_func(lp[i]) - g->y[i]);
	    d2 += sm.x[i] * sm.x[i] * phi2_func(lp[i]);
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

	 /* check for convergence */
	 if(epoch > 1 && convergetest(beta[j], beta_new, threshold))
	 {
	    converged[j] = TRUE;
	    numconverged++;
	 }

	 /* clip very large coefs to limit divergence */
	 beta_new = fmin(fmax(beta_new, -truncl), truncl);

	 /* update linear predictor */
	 for(i = 0 ; i < g->n ; i++)
	    if(sm.x[i] != 0)
	       lp[i] = fmin(fmax(
			lp[i] + sm.x[i] * (beta_new - beta[j]),
		       -MAXLP), MAXLP);

	 beta[j] = beta_new;
      }

      if(verbose > 1)
      {
	 loss = 0;
	 for(i = 0 ; i < g->n ; i++)
	 {
	    if(verbose > 1)
	       printf("\tlp[%d]: %.3f\n", i, lp[i]);
	    loss += loss_pt_func(lp[i], g->y[i]) / g->n;
	 }
      }
	 
      if(epoch > 1)
      {
	 zeros = 0;

	 for(j = 1 ; j < g->p + 1 ; j++)
	 {
	    if(fabs(beta[j]) < ZERO_THRESH)
	    {
	       beta[j] = 0;
	       zeros++;
	    }

	    if(verbose > 1 && g->p - numconverged < 100 && !converged[j])
	    {
	       printf("not converged: %d (%.20f, %.20f, %.20f)\n",
		     j, beta[j], grad, d2);
	    }
	 }
      }

      if(verbose)
      {
	 printf("Epoch %d  converged: %d zeros: %d  non-zeros: %d\n",
	       epoch, numconverged, zeros, g->p - zeros);
      }
      else if(verbose > 1)
      {
	 printf("Epoch %d  training loss: %.5f  converged: %d\
  zeros: %d  non-zeros: %d\n",
   epoch, loss, numconverged, zeros, g->p - zeros);
      }

      if(numconverged == g->p + 1)
      {
	 printf("all converged\n");
	 allconverged++;

	 /* converged twice in a row, no need to continue */
	 if(allconverged == CONVERGED)
	 {
	    printf("terminating with %d non-zero coefs\n\n",
		  g->p - zeros);
	    break;
	 }

	 for(j = 0 ; j < g->p + 1 ; j++)
	    converged[j] = FALSE;
	 numconverged = 0;
      }

      if(verbose)
	 printf("\n");

      epoch++;
   }

   free(converged);
   free(lp);
   sample_free(&sm);

   if(allconverged == CONVERGED)
      return g->p - zeros;
   else
      return FAILURE;
}

