#include "common.h"
#include "cd.h"

#define MAX_SHOW_NOTCONV 100
#define FMAX(a, b) (a < b ? b : a) 
#define FMIN(a, b) (a < b ? a : b)

short convergetest(double a, double b, double threshold)
{
   /* absolute convergence */
   if(fabs(a) <= ZERO_THRESH && fabs(b) <= ZERO_THRESH)
      return TRUE;

   /* relative convergence */
   return (fabs(a - b) / (fabs(a) + fabs(b))) < threshold;
}

/* Find smallest lambda1 that makes all coefficients
 * zero (except the intercept)
 */
double get_lambda1max_gmatrix(
      gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      inv inv_func,
      step step_func
      )
{
   int i, j;
   double *lp = NULL;
   double grad, d2, s, zmax = 0, beta0;
   sample sm;

   CALLOCTEST(lp, g->n, sizeof(double))
   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   g->nextcol(g, &sm);

   /* First compute the intercept. When all other variables
    * are zero, the intercept is just the inv(mean(y)) */
   s = 0;
   for(i = 0 ; i < g->n ; i++)
      s += (double)g->y[i];

   beta0 = inv_func(s / g->n);

   for(i = 0 ; i < g->n ; i++)
      lp[i] = beta0;

   /* find smallest lambda1 that makes all coefficients zero, by
    * finding the largest z, but let the intercept affect lp
    * first because it's not penalised. */
   for(j = 1 ; j < g->p + 1; j++)
   {
      grad = d2 = 0;

      g->nextcol(g, &sm);

      step_func(sm.x, g->y, lp, g->n,
	       phi1_func, phi2_func, &grad, &d2);

      /* don't move if 2nd derivative is zero */
      s = 0;
      if(d2 != 0)
	 s = -grad / d2;

      if(zmax < fabs(s))
	 zmax = fabs(s);
   } 

   free(lp);
   sample_free(&sm);

   return zmax;
}

void step_regular(dtype *x, dtype *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func, double *grad, double *d2)
{
   int i;
   double lphi1;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
   {
      /* skip zeros, they don't change the linear predictor */
      if(x[i] == 0)
	 continue;

      lphi1 = phi1_func(lp[i]);
      (*grad) += x[i] * (lphi1 - y[i]);
      (*d2) += x[i] * x[i] * phi2_func(lphi1);
   }
}

void step_grouped(dtype *x, dtype *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func, double *grad, double *d2)
{
   /*int i;
   double lphi1;*/

   /*for(i = 0 ; i < n ; i++)
   {
      if(x[i] == 0)
	 continue;

      lphi1 = phi1_func(lp[i]);
      (*grad) += x[i] * (lphi1 - y[i]);
      (*d2) += x[i] * x[i] * phi2_func(lphi1);
   }*/
}

/* coordinate descent */
int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
      inv inv_func,
      step step_func,
      int maxepoch, double *beta, double *lp, double lambda1, double lambda2,
      double threshold, int verbose, int *trainf, double trunc)
{
   int i, j, epoch = 1;
   double loss = 0, grad, d2, s, beta_new;
   short *converged = NULL;
   int numconverged = 0;
   /*double *lp = NULL;*/
   sample sm;
   double truncl = log((1 - trunc) / trunc);
   int allconverged = 0, zeros = 0;
   const int CONVERGED = 2;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   CALLOCTEST(converged, g->p + 1, sizeof(short));
   /*CALLOCTEST(lp, g->n, sizeof(double));*/

   while(epoch <= maxepoch)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 g->nextcol(g, &sm);

	 if(converged[j])
	   continue;

	 grad = d2 = 0;

	 step_func(sm.x, g->y, lp, g->n,
	       phi1_func, phi2_func, &grad, &d2);
	 
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
	    {
	       /*lp[i] = fmin(MAXLP, fmax(-MAXLP,
		     lp[i] + sm.x[i] * (beta_new - beta[j])));*/
	       lp[i] += sm.x[i] * (beta_new - beta[j]);
	       if(lp[i] < -MAXLP)
		  lp[i] = -MAXLP;
	       else if(lp[i] > MAXLP)
		  lp[i] = MAXLP;
	    }

	 beta[j] = beta_new;
      }

      if(verbose > 1)
      {
	 loss = 0;
	 for(i = 0 ; i < g->n ; i++)
	 {
	    if(verbose > 3)
	       printf("\tlp[%d]: %.3f\n", i, lp[i]);
	    loss += loss_pt_func(lp[i], g->y[i]) / g->n;
	 }
      }
	 
      if(epoch > 1)
      {
	 zeros = 0;

	 if(verbose > 1 && !converged[0])
	    printf("intercept not converged: %.10f\n", beta[0]); 

	 for(j = 1 ; j < g->p + 1 ; j++)
	 {
	    if(fabs(beta[j]) < ZERO_THRESH)
	    {
	       beta[j] = 0;
	       zeros++;
	    }

	    if(verbose > 1 && !converged[j])
	       printf("beta[%d] not converged: %.10f\n", j, beta[j]);
	 }
      }

      if(verbose == 1)
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

      /* All vars must converge in two consecutive epochs 
       * in order to stop */
      if(numconverged == g->p + 1)
      {
	 printf("all converged\n");
	 allconverged++;
	 
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
      else
	 allconverged = 0;

      if(verbose)
	 printf("\n");

      epoch++;
   }

   free(converged);
 /*  free(lp);*/
   sample_free(&sm);

   if(allconverged == CONVERGED)
      return g->p - zeros + 1;
   else
      return FAILURE;
}

