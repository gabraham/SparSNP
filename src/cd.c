#include "common.h"
#include "cd.h"

int convergetest(double a, double b, double threshold)
{
   /* absolute convergence */
   if(fabs(a) <= ZERO_THRESH && fabs(b) <= ZERO_THRESH)
      return TRUE;

   /* relative convergence */
   return (fabs(a - b) / (fabs(a) + fabs(b))) < threshold;
}

double clip(const double x, const double min, const double max)
{
   return (x > max) ? max : ((x < min) ? min : x);
}

double zero(const double x, const double thresh)
{
   return (fabs(x) < thresh) ? 0 : x;
}

/* update linear predictor and clip it if it's too large */
void updatelp(gmatrix *g, const double beta_new, const unsigned int j,
      const double *restrict x)
{
   unsigned int i, n = g->n;
   const double beta = g->beta[j];
   double *restrict lp_invlogit = g->lp_invlogit,
	  *restrict lp = g->lp,
	  *restrict y = g->y,
	  *restrict ylp = g->ylp,
	  *restrict ylp_neg = g->ylp_neg;

   if(x)
   {
      for(i = 0 ; i < n ; i++)
	 g->lp[i] = clip(g->lp[i] + x[i] * (beta_new - beta), -MAXLP, MAXLP);
   }
   else /* update from intercept */
   {
      for(i = 0 ; i < n ; i++)
	 lp[i] = beta_new;
   }

   /* update functions of linear predictor */
   if(g->model == MODEL_LOGISTIC)
   {
      for(i = 0 ; i < n ; i++)
	 lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = 0 ; i < n ; i++)
      {
	 ylp[i] = y[i] * lp[i] - 1;
	 ylp_neg[i] = (ylp[i] < 0);
      }
   }
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
   unsigned int i, j, n = g->n, p1 = g->p + 1;
   double s, zmax = 0, beta_new;
   sample sm;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   g->nextcol(g, &sm);

   /* First compute the intercept. When all other variables
    * are zero, the intercept is just inv(mean(y)) for a suitable inv()
    * function depending on the loss */
   s = 0;
   for(i = 0 ; i < n ; i++)
      s += g->y[i];

   beta_new = inv_func(s / g->n);
   updatelp(g, beta_new, 0, NULL);

   printf("beta[0]: %.10f\n", beta_new);

   /* find smallest lambda1 that makes all coefficients zero, by
    * finding the largest z, but let the intercept affect lp
    * first because it's not penalised. */
   for(j = 1 ; j < p1; j++)
   {
      g->nextcol(g, &sm);
      if(!g->active[j])
	 continue;

      s = fabs(step_func(&sm, g, phi1_func, phi2_func));
      if(zmax < s)
	 zmax = s;
   } 

   sample_free(&sm);

   return zmax;
}

double step_generic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   unsigned int i, n = g->n;
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
double step_regular_linear(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->n;
   double grad = 0;
   double *restrict x_tmp = s->x, 
          *restrict lp_tmp = g->lp,
	  *restrict y_tmp = g->y;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
      grad += x_tmp[i] * (lp_tmp[i] - y_tmp[i]);

   return grad / n;
}

double step_regular_logistic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   unsigned int i, n = g->n;
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
 * Squared hinge loss, assumes y \in {-1,1},
 * and that X is scaled so that the 2nd derivative is always 1
 */
double step_regular_sqrhinge(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   unsigned int i, n = g->n;
   double grad = 0;
   const double *restrict x_tmp = s->x,
	        *restrict y_tmp = g->y,
		*restrict ylp_tmp = g->ylp,
		*restrict ylp_neg_tmp = g->ylp_neg;

   /* compute gradient */
   for(i = 0 ; i < n ; i++)
      grad += ylp_neg_tmp[i] * y_tmp[i] * x_tmp[i] * ylp_tmp[i];
   return grad / n;
}

/* coordinate descent */
int cd_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      step step_func,
      const int maxepochs,
      const int maxiters,
      const double lambda1,
      const double lambda2,
      const double threshold,
      const int verbose,
      const double trunc)
{
   const int CONVERGED = 2;
   unsigned int j, numiter,
       epoch = 1, numconverged = 0,
       allconverged = 0, zeros = 0;
   short *converged = NULL;
   double s, beta_new;
   const double truncl = log2((1 - trunc) / trunc),
	        l2recip = 1 / (1 + lambda2);
   sample sm;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   CALLOCTEST(converged, g->p + 1, sizeof(short));

   while(epoch <= maxepochs)
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

	 numiter = 0;

	 while(!converged[j] && numiter <= maxiters)
	 {
	    numiter++;

	    s = step_func(&sm, g, phi1_func, phi2_func);

	    /* don't penalise intercept */
	    beta_new = g->beta[j] - s;
	    if(j > 0)
	       beta_new = soft_threshold(beta_new, lambda1) * l2recip;

	    /* check for convergence */
	    if(numiter > 1 && convergetest(g->beta[j], beta_new, threshold))
	    {
	       converged[j] = TRUE;
	       numconverged++;
	    }

	    /* clip very large coefs to limit divergence */
	    beta_new = clip(beta_new, -truncl, truncl);
	    beta_new = zero(beta_new, ZERO_THRESH);

	    updatelp(g, beta_new, j, sm.x);
	    g->beta[j] = beta_new;
	 }
	 if(numiter > maxiters)
	    printf("max number of internal iterations (%d) \
reached for variable: %d\n", maxiters, j);
      }

      /* count number of zero variables, excluding the intercept */
      if(epoch > 1)
      {
	 zeros = 0;
	 for(j = 1 ; j < g->p + 1 ; j++)
	 {
	    if(fabs(g->beta[j]) < ZERO_THRESH)
	    {
	       g->beta[j] = 0;
	       zeros++;
	    }
	 }
      }

      if(verbose)
	 printf("Epoch %d  converged: %d zeros: %d  non-zeros: %d\n",
	       epoch, numconverged, zeros, g->p - zeros);

      /* All vars must converge in two consecutive epochs 
       * in order to stop */
      if(numconverged == g->p + 1)
      {
	 if(verbose)
	    printf("all converged\n");
	 allconverged++;
	 
	 if(allconverged == CONVERGED)
	 {
	    if(verbose)
	       printf("terminating with %d non-zero coefs\n\n",
		     g->p - zeros);
	    break;
	 }

	 /* reset convergence for next epoch */
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
   sample_free(&sm);

   if(allconverged == CONVERGED)
      return g->p - zeros + 1;
   return FAILURE;
}

