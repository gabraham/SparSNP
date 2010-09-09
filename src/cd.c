#include "common.h"
#include "cd.h"

#define printfverb(...) if(verbose) printf(__VA_ARGS__)

int convergetest(double a, double b, double threshold)
{
   /* absolute convergence */
   if(fabs(a) <= ZERO_THRESH && fabs(b) <= ZERO_THRESH)
      return TRUE;

   /* relative convergence */
   return (fabs(a - b) / (fabs(a) + fabs(b))) < threshold;
}

inline double clip(const double x, const double min, const double max)
{
   return (x > max) ? max : ((x < min) ? min : x);
}

inline double zero(const double x, const double thresh)
{
   return (fabs(x) < thresh) ? 0 : x;
}

/* update linear predictor and clip it if it's too large */
void updatelp(gmatrix *g, const double beta_new, const int j,
      const double *restrict x)
{
   int i, n = g->ntrain[g->fold];
   const double beta_diff = beta_new - g->beta[j];
   double *restrict lp_invlogit = g->lp_invlogit,
	  *restrict lp = g->lp,
	  *restrict y = g->y,
	  *restrict ylp = g->ylp,
	  *restrict ylp_neg = g->ylp_neg;

   if(x)
   {
      for(i = n - 1 ; i >= 0 ; --i)
	 /*g->lp[i] = clip(g->lp[i] + x[i] * beta_diff, -MAXLP, MAXLP);*/
	 g->lp[i] = lp[i] + x[i] * beta_diff;
   }
   else /* update from intercept */
   {
      for(i = n - 1 ; i >= 0 ; --i)
	 lp[i] = beta_new;
   }

   /* update functions of linear predictor */
   if(g->model == MODEL_LOGISTIC)
   {
      for(i = n - 1 ; i >= 0 ; --i)
	 lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = n - 1 ; i >= 0 ; --i)
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
      step step_func)
{
   int i, j, n = g->ntrain[g->fold], p1 = g->p + 1;
   double s, zmax = 0, beta_new;
   sample sm;

   if(!sample_init(&sm, n, g->inmemory))
      return FAILURE;

   g->nextcol(g, &sm);

   /* First compute the intercept. When all other variables
    * are zero, the intercept is just inv(mean(y)) for a suitable inv()
    * function depending on the loss */
   s = 0;
   for(i = 0 ; i < n ; i++)
      s += g->y[i];

   beta_new = inv_func(s / n);
   updatelp(g, beta_new, 0, NULL);
   printf("intercept: %.10f\n", beta_new);

   /* find smallest lambda1 that makes all coefficients zero, by
    * finding the largest z, but let the intercept affect lp
    * first because it's not penalised. */
   for(j = 1 ; j < p1; j++)
   {
      g->nextcol(g, &sm);
      if(g->ignore[j])
	 continue;

      s = fabs(step_func(&sm, g, phi1_func, phi2_func));
      zmax = (zmax < s) ? s : zmax;
   } 

   sample_free(&sm);

   return zmax;
}

double step_generic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->ntrain[g->fold];
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
   /*int i, n = g->ntrain[g->fold];*/
   int i, n = g->n;
   double grad = 0;
   double *restrict x_tmp = s->x, 
          *restrict lp_tmp = g->lp,
	  *restrict y_tmp = g->y;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
      grad += x_tmp[i] * (lp_tmp[i] - y_tmp[i]);

/*   return grad * g->ntrainrecip[g->fold];*/
   return grad / n;
}

double step_regular_logistic(sample *s, gmatrix *g,
      phi1 phi1_func, phi2 phi2_func)
{
   int i, n = g->ntrain[g->fold];
   double grad = 0, d2 = 0;
   double *restrict y_tmp = g->y,
	  *restrict lp_invlogit_tmp = g->lp_invlogit,
	  *restrict x_tmp = s->x,
	  *restrict x2_tmp = s->x2;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
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
   int i, n = g->ntrain[g->fold];
   double grad = 0;
   const double *restrict x_tmp = s->x,
	        *restrict y_tmp = g->y,
		*restrict ylp_tmp = g->ylp,
		*restrict ylp_neg_tmp = g->ylp_neg;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
      grad += ylp_neg_tmp[i] * y_tmp[i] * x_tmp[i] * ylp_tmp[i];
   return grad * g->ntrainrecip[g->fold];
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
      const double thresh,
      const int verbose,
      const double trunc)
{
   const int n = g->ntrain[g->fold], p = g->p, p1 = g->p + 1;
   int j, iter, allconverged = 0, numactive = 0,
       epoch = 1, numconverged = 0,
       good = FALSE;
   int *converged = NULL,
       *active_old = NULL,
       *active_new = NULL;
   double s, beta_new;
   const double truncl = log2((1 - trunc) / trunc),
	        l2recip = 1 / (1 + lambda2);
   sample sm;

   if(!sample_init(&sm, n, g->inmemory))
      return FAILURE;

   CALLOCTEST(converged, p1, sizeof(int));
   CALLOCTEST(active_new, p1, sizeof(int));
   CALLOCTEST(active_old, p1, sizeof(int));

   /* start off with all variables marked active
    * even though they're all zero, unless they're
    * marked as ignore */
   for(j = p ; j >= 0 ; --j)
      active_new[j] = active_old[j] = !g->ignore[j];

   while(epoch <= maxepochs)
   {
      printf(">>> start epoch %d <<<\n", epoch);
      for(j = 0 ; j < p1; j++)
      {
	 g->nextcol(g, &sm);
	 iter = 0;

	 printf("active[%d]: %d; ", j, active_new[j]);
	 if(active_new[j])
	 {
	    printf("iterating over %d ", j);
	    /* iterate over jth variable */
      	    while(iter < maxiters)
      	    {
	       printf(".");
      	       s = step_func(&sm, g, phi1_func, phi2_func);
      	       beta_new = g->beta[j] - s;
      	       if(j > 0) /* don't penalise intercept */
      	          beta_new = soft_threshold(beta_new, lambda1) * l2recip;
      	       beta_new = clip(beta_new, -truncl, truncl);
      	       beta_new = zero(beta_new, ZERO_THRESH);
      
      	       converged[j] = convergetest(g->beta[j], beta_new, thresh);
      	       updatelp(g, beta_new, j, sm.x);
      	       g->beta[j] = beta_new;
	       if(converged[j])
		  break;
      	       iter++;
      	    }
	 }
	 else
	    printf("skipping %d, ", j);

	 active_new[j] = (g->beta[j] != 0);
	 printf(" ; active[%d]: %d; converged[%d]: %d; beta: %.20f\n",
	       j, active_new[j], j, converged[j], g->beta[j]);


	 if(iter > maxiters)
	    printfverb("max number of internal iterations (%d) \
reached for variable: %d\n", maxiters, j);
      }

      /*printfverb("Epoch %d  converged: %d zeros: %d  \
nonzeros: %d  zero_same: %d\n", epoch, numconverged, zeros,
	    g->p - zeros, ); */

      numconverged = 0;
      for(j = p ; j >= 0; --j)
	 numconverged += converged[j];

      /* state machine for active set convergence */ 
      if(numconverged == p1) 
      {
	 /* prepare for another iteration over all
	  * non-ignored variables, store a copy of the active set
	  * for later */
	 if(allconverged == 0)
	 {
	    printf("all converged\n");
	    for(j = p ; j >= 0 ; --j)
	    {
	       active_old[j] = active_new[j];
	       active_new[j] = !g->ignore[j];
	    }
	    allconverged++;
	 }
	 else /* 2nd iteration done, check
		 whether active set has changed */
	 {
	    numactive = 0;
	    for(j = 0 ; j <= p ; j++)
	    {
	       numactive += active_new[j];
	       printf("%d ", j);
	       if(active_new[j] != active_old[j])
		  break;
	    }
	    printf("\n");

	    /* all equal, terminate */
	    if(j > p)
	    {
	       printf("terminating with %d active vars, active set is the same\n",
		     numactive);
	       good = TRUE;
	       break;
	    }

	    printf("active set changed: %d\n", j);

	    /* active set has changed, copy the new state and
	     * iterate over new active set */
	    for(j = p ; j >= 0 ; --j)
	       active_old[j] = active_new[j];
	    allconverged = 0;
	 }
      }
      else /* reset to first state */
	 allconverged = 0;

      printf(">>> end epoch %d <<<\n\n", epoch);
      epoch++;
   }

   sample_free(&sm);
   free(converged);
   free(active_old);
   free(active_new);

   if(good)
      return numactive;
   return FAILURE;
}

