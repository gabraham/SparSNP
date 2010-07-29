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
   /*double grad, d2;*/
   double s, zmax = 0, beta0;
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
      /*grad = d2 = 0;*/

      g->nextcol(g, &sm);

      /*step_func(&sm, g->y, lp, g->n,
	       phi1_func, phi2_func, &grad, &d2);*/
      s = -step_regular(&sm, g->y, lp, g->n,
	       phi1_func, phi2_func);

      /* don't move if 2nd derivative is zero */
    /*  s = 0;
      if(d2 != 0)
	 s = -grad / d2;*/

      if(s != 0 && zmax < fabs(s))
	 zmax = fabs(s);
   } 

   free(lp);
   sample_free(&sm);

   return zmax;
}

double step_regular(sample *s, double *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func)
{
   int i;
   double lphi1;
   double x, x2;
   double grad = 0, d2 = 0;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
   {
      x = s->x[i];
      /* skip zeros, they don't change the linear predictor */
      if(x == 0)
	 continue;
      x2 = s->x2[i];

      lphi1 = phi1_func(lp[i]);
      grad += x * (lphi1 - y[i]);
      d2 += x2 * phi2_func(lphi1);
   }
   if(d2 == 0)
      return 0;
   return grad / d2;
}

double step_grouped(sample *s, double *y, double *lp, int n,
      phi1 phi1_func, phi2 phi2_func)
{
   /*int i, k;
   double lphi1;

   for(i = 0 ; i < n ; i++)
   {
      lphi1 = phi1_func(lp[i]);
      for(k = 0 ; k < s->nbins ; k++)
      {
	 (*grad) += 	 
      }
   }*/
   return 0;
}

/* coordinate descent */
int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
      inv inv_func,
      step step_func,
      int maxepoch, double *beta, double *lp, double lambda1,
      double lambda2, double threshold, int verbose,
      int *trainf, double trunc)
{
   int i, j, epoch = 1, numconverged = 0, numiter;
   double loss = 0, s, beta_new;
   /*double grad, d2;*/
   short *converged = NULL;
   sample sm;
   double truncl = log2((1 - trunc) / trunc);
   int allconverged = 0, zeros = 0;
   const int CONVERGED = 2;
   double l2recip = 1 / (1 + lambda2);

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   CALLOCTEST(converged, g->p + 1, sizeof(short));

   while(epoch <= maxepoch)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 g->nextcol(g, &sm);
	 numiter = 0;

	 while(!converged[j] && numiter <= maxepoch)
	 {
	    numiter++;
	   /* grad = d2 = 0;*/

	    /*step_func(&sm, g->y, lp, g->n,
		  phi1_func, phi2_func, &grad, &d2);*/
	    s = step_regular(&sm, g->y, lp, g->n,
		  phi1_func, phi2_func);

	    /* don't move if 2nd derivative is zero */
	    /*s = 0;
	    if(d2 != 0)
	       s = grad / d2;*/

	    if(s == 0)
	       continue;

	    /* don't penalise intercept */
	    if(j == 0)
	       beta_new = beta[j] - s;
	    else
	       beta_new = soft_threshold(beta[j] - s, lambda1) * l2recip;

	    /* check for convergence */
	    if(numiter > 1 && convergetest(beta[j], beta_new, threshold))
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
		  lp[i] += sm.x[i] * (beta_new - beta[j]);
		  if(lp[i] < -MAXLP)
		     lp[i] = -MAXLP;
		  else if(lp[i] > MAXLP)
		     lp[i] = MAXLP;
	       }

	    beta[j] = beta_new;
	 }
      }

#ifdef VERBOSE
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
#endif
	 
      if(epoch > 1)
      {
	 zeros = 0;

#ifdef VERBOSE
	 if(verbose > 1 && !converged[0])
	    printf("intercept not converged: %.10f\n", beta[0]); 
#endif

	 for(j = 1 ; j < g->p + 1 ; j++)
	 {
	    if(fabs(beta[j]) < ZERO_THRESH)
	    {
	       beta[j] = 0;
	       zeros++;
	    }

#ifdef VERBOSE
	    if(verbose > 1 && !converged[j])
	       printf("beta[%d] not converged: %.10f\n", j, beta[j]);
#endif
	 }
      }

#ifdef VERBOSE
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
#endif

      /* All vars must converge in two consecutive epochs 
       * in order to stop */
      if(numconverged == g->p + 1)
      {
#ifdef VERBOSE
	 if(verbose)
	    printf("all converged\n");
#endif
	 allconverged++;
	 
	 if(allconverged == CONVERGED)
	 {
#ifdef VERBOSE
	    if(verbose)
	       printf("terminating with %d non-zero coefs\n\n",
		     g->p - zeros);
#endif
	    break;
	 }

	 for(j = 0 ; j < g->p + 1 ; j++)
	    converged[j] = FALSE;
	 numconverged = 0;
      }
      else
	 allconverged = 0;

#ifdef VERBOSE
      if(verbose)
	 printf("\n");
#endif

      epoch++;
   }

   free(converged);
   sample_free(&sm);

   if(allconverged == CONVERGED)
      return g->p - zeros + 1;
   else
      return FAILURE;
}

