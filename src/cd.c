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
   double s, zmax = 0, beta0;
   sample sm;

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
      g->lp[i] = beta0;

   printf("beta0: %.10f\n", beta0);

   /* find smallest lambda1 that makes all coefficients zero, by
    * finding the largest z, but let the intercept affect lp
    * first because it's not penalised. */
   for(j = 1 ; j < g->p + 1; j++)
   {
      g->nextcol(g, &sm);
      if(!g->active[j])
	 continue;

      s = -step_func(&sm, g, phi1_func, phi2_func);
      if(s != 0 && zmax < fabs(s))
	 zmax = fabs(s);
   } 

   sample_free(&sm);

   return zmax;
}

double step_regular(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->n;
   double lphi1;
   double grad = 0, d2 = 0;
   double *restrict x_tmp = s->x,
	  *restrict x2_tmp = s->x2,
	  *restrict lp_tmp = g->lp,
	  *restrict y_tmp = g->y;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
   {
      lphi1 = phi1_func(lp_tmp[i]);
      grad += x_tmp[i] * (lphi1 - y_tmp[i]);
      d2 += x2_tmp[i] * phi2_func(lphi1);
   }
   if(d2 == 0)
      return 0;
   return grad / d2;
}

/* In linear regression, for standardised inputs x, the 2nd derivative is
 * always N since it is the sum of squares \sum_{i=1}^N x_{ij}^2 =
 * \sum_{i=1}^N 1 = N
 */
double step_regular_l2(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->n;
   double grad = 0;
   double *restrict x_tmp = s->x, 
          *restrict lp_tmp = g->lp,
	  *restrict y_tmp = g->y;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
      grad += x_tmp[i] * (lp_tmp[i] - y_tmp[i]);

   return grad / n;
}

double step_regular_logistic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->n;
   double grad = 0, d2 = 0;
   double *restrict y_tmp = g->y,
	  *restrict lp_invlogit_tmp = g->lp_invlogit,
	  *restrict x_tmp = s->x,
	  *restrict x2_tmp = s->x2;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
   {
      grad += x_tmp[i] * (lp_invlogit_tmp[i] - y_tmp[i]);
      d2 += x2_tmp[i] * lp_invlogit_tmp[i] * (1 - lp_invlogit_tmp[i]);
   }
   if(d2 == 0)
      return 0;
   return grad / d2;
}

/*
 * Squared hinge loss, assumes y \in {-1,1}
 */
double step_regular_sqrhinge(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->n;
   double grad = 0;
   double *restrict x_tmp = s->x,
	  *restrict y_tmp = g->y,
	  *restrict ylp_tmp = g->ylp,
	  *restrict ylp_pos_tmp = g->ylp_pos;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
      grad += ylp_pos_tmp[i] * y_tmp[i] * x_tmp[i] * ylp_tmp[i];
   return grad / n;
}

/* coordinate descent */
int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      loss_pt loss_pt_func,    /* loss for one sample */
      inv inv_func,
      step step_func,
      int maxepoch, double *beta, double lambda1,
      double lambda2, double threshold, int verbose,
      int *trainf, double trunc)
{
   const int maxiter = 100;
   int i, j, epoch = 1, numconverged = 0, numiter;
   double s, beta_new;
   short *converged = NULL;
   sample sm;
   double truncl = log2((1 - trunc) / trunc);
   int allconverged = 0, zeros = 0;
   const int CONVERGED = 2;
   double l2recip = 1 / (1 + lambda2);
   double *restrict x;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   CALLOCTEST(converged, g->p + 1, sizeof(short));

   while(epoch <= maxepoch)
   {
      for(j = 0 ; j < g->p + 1; j++)
      {
	 g->nextcol(g, &sm);
	 if(!g->active[j])
	 {
	    if(!converged[j])
	    {
	       converged[j] = TRUE;
	       numconverged++;
	    }
	    continue;
	 }

	 x = sm.x;
	 numiter = 0;

	 while(!converged[j] && numiter <= maxiter)
	 {
	    numiter++;

	    s = step_func(&sm, g, phi1_func, phi2_func);

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
	    beta_new = (beta_new > truncl) ? 
		  truncl : ((beta_new < -truncl) ? -truncl : beta_new);

	    /* update linear predictor */
	    for(i = 0 ; i < g->n ; i++)
	    {
	       g->lp[i] += x[i] * (beta_new - beta[j]);
	       g->lp[i] = (g->lp[i] > MAXLP) ? 
		     MAXLP : ((g->lp[i] < -MAXLP) ? -MAXLP : g->lp[i]); 
	    }

	    /* update functions of linear predictor */
	    if(g->model == MODEL_LOGISTIC)
	    {
	       for(i = 0 ; i < g->n ; i++)
		  g->lp_invlogit[i] = 1 / (1 + exp(-g->lp[i]));
	    }
	    else if(g->model == MODEL_SQRHINGE)
	    {
	       for(i = 0 ; i < g->n ; i++)
	       {
		  g->ylp[i] = g->y[i] * g->lp[i];
		  g->ylp_pos[i] = (g->ylp[i] > 0);
	       }
	    }

	    beta[j] = beta_new;
	 }
	 if(numiter > maxiter)
	    printf("maximum number of internal iterations reached\n");
      }

      /* count number of zero variables */
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
	 }
      }

#ifdef VERBOSE
      if(verbose == 1)
      {
	 printf("Epoch %d  converged: %d zeros: %d  non-zeros: %d\n",
	       epoch, numconverged, zeros, g->p - zeros);
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

	 /* reset convergence for next epoch */
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

