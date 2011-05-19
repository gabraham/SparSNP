#include "common.h"
#include "cd.h"

#define printfverb(...) if(verbose) printf(__VA_ARGS__)

static double clip(const double x, const double min, const double max);
static double zero(const double x, const double thresh);

inline static double clip(const double x, const double min, const double max)
{
   return (x > max) ? max : ((x < min) ? min : x);
}

inline static double zero(const double x, const double thresh)
{
   return (fabs(x) < thresh) ? 0 : x;
}

/* Update linear predictor and related variables.
 * */
void updatelp(gmatrix *g, const double update,
      const double *restrict x)
{
   int i, n = g->ncurr;
   double tmp;
   double *restrict lp_invlogit = g->lp_invlogit,
	  *restrict lp = g->lp,
	  *restrict y = g->y,
	  *restrict ylp_neg_y_ylp = g->ylp_neg_y_ylp;
   double loss = 0, ylp = 0;

   if(g->model == MODEL_LINEAR)
   {
      if(x)
      {
	 for(i = n - 1 ; i >= 0 ; --i)
	 {
	    lp[i] += x[i] * update;
	    tmp = y[i] - lp[i];
	    loss += tmp * tmp;
	 }
      }
      else
      {
	 for(i = n - 1 ; i >= 0 ; --i)
	 {
	    lp[i] = update;
	    tmp = y[i] - lp[i];
	    loss += tmp * tmp;
	 }
      }
   }
   else if(g->model == MODEL_LOGISTIC)
   {
      if(x)
      {
	 for(i = n - 1 ; i >= 0 ; --i)
	 {
	    lp[i] += x[i] * update;
	    lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	    loss += log(1 + exp(lp[i])) - y[i] * lp[i];
	 }
      }
      else
      {
	 for(i = n - 1 ; i >= 0 ; --i)
	 {
	    lp[i] += update;
	    lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	    loss += log(1 + exp(lp[i])) - y[i] * lp[i];
	 }
      }
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      if(x)
      {
	 for(i = n - 1 ; i >= 0 ; --i)
      	 {
      	    lp[i] += x[i] * update;
      	    ylp = y[i] * lp[i] - 1;
      	    
	    if(ylp < 0)
	    {
      	       ylp_neg_y_ylp[i] = y[i] * ylp;
	       loss += ylp * ylp;
	    }
	    else
	       ylp_neg_y_ylp[i] = 0;
      	 }
      }
      else
      {
	 for(i = n - 1 ; i >= 0 ; --i)
      	 {
      	    lp[i] += update;
      	    ylp = y[i] * lp[i] - 1;

	    if(ylp < 0)
	    {
      	       ylp_neg_y_ylp[i] = y[i] * ylp;
	       loss += ylp * ylp;
	    }
	    else
	       ylp_neg_y_ylp[i] = 0;
      	 }
      }
   }

   loss *= g->ncurr_recip;
   g->loss = loss;
}

/* 
 * Find smallest lambda1 that makes all coefficients
 * zero (except the intercept)
 */
double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func, phi2 phi2_func, inv inv_func, step step_func)
{
   int i, j, n = g->ncurr, n1 = g->ncurr - 1, p1 = g->p + 1;
   double s, zmax = 0, beta_new;
   sample sm;

   if(!sample_init(&sm))
      return FAILURE;

   g->nextcol(g, &sm, 0, NA_ACTION_ZERO);

   /* First compute the intercept. When all other variables
    * are zero, the intercept is just inv(mean(y)) for a suitable inv()
    * function depending on the loss */
   s = 0;
   for(i = n1 ; i >= 0 ; --i)
      s += g->y[i];

   beta_new = inv_func(s / n);
   updatelp(g, beta_new, NULL);
   printf("intercept: %.15f (%d samples)\n", beta_new, n);

   /* find smallest lambda1 that makes all coefficients zero, by
    * finding the largest z, but let the intercept affect lp
    * first because it's not penalised. */
   for(j = 1 ; j < p1; j++)
   {
      if(g->ignore[j])
	 continue;
      g->nextcol(g, &sm, j, NA_ACTION_ZERO);

      s = fabs(step_func(&sm, g));
      zmax = (zmax < s) ? s : zmax;
   } 

   g->beta[0] = beta_new;

   return zmax;
}

/* In linear regression, for standardised inputs x, the 2nd derivative is
 * always N since it is the sum of squares \sum_{i=1}^N x_{ij}^2 =
 * \sum_{i=1}^N 1 = N
 */
double step_regular_linear(sample *s, gmatrix *g)
{
   int i;
   double grad = 0;
   double *restrict x = s->x, 
          *restrict lp = g->lp,
	  *restrict y = g->y;

   /* compute gradient */
   for(i = g->ncurr - 1 ; i >= 0 ; --i)
      grad += x[i] * (lp[i] - y[i]);

   return grad * g->ncurr_recip;
}

double step_regular_logistic(sample *s, gmatrix *g)
{
   int i, n = g->ncurr;
   double grad = 0, d2 = 0;
   double *restrict y = g->y,
	  *restrict lp_invlogit = g->lp_invlogit,
	  *restrict x = s->x;

   /* compute 1st and 2nd derivatives */
   for(i = n - 1 ; i >= 0 ; --i)
   {
      grad += x[i] * (lp_invlogit[i] - y[i]);
      d2 += x[i] * x[i] * lp_invlogit[i] * (1 - lp_invlogit[i]);
   }

   grad *= g->ncurr_recip;
   d2 *= g->ncurr_recip;

   if(d2 == 0)
      return 0;
   return grad / d2;
}

/*
 * Squared hinge loss, assumes y \in {-1,1},
 * and that X is scaled so that the 2nd derivative is always N
 */
double step_regular_sqrhinge(sample *s, gmatrix *g)
{
   int i, n = g->ncurr;
   double grad = 0;
   const double *restrict x = s->x,
		*restrict ylp_neg_y_ylp = g->ylp_neg_y_ylp;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
      grad += ylp_neg_y_ylp[i] * x[i];

   return grad * g->ncurr_recip; /* avoid division */
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
   const int p = g->p, p1 = g->p + 1;
   int j, iter, allconverged = 0, numactive = 0,
       epoch = 1, numconverged = 0,
       good = FALSE;
   int *active_old = NULL;
   double s_old = 0, s = 0, beta_new;
   const double truncl = log((1 - trunc) / trunc),
	        l2recip = 1 / (1 + lambda2);
   double *m = NULL;
   sample sm;
   double old_loss = 0;
   int conv = 0;

   if(!sample_init(&sm))
      return FAILURE;

   /*CALLOCTEST(beta_old, p1, sizeof(double));*/
   CALLOCTEST(active_old, p1, sizeof(int));
   CALLOCTEST(m, p1, sizeof(double));

   /* start off with all variables marked active
    * even though they're all zero, unless they're
    * marked as ignore */
   for(j = p ; j >= 0 ; --j)
   {
      active_old[j] = g->active[j];
      m[j] = 1.0;
   }

   while(epoch <= maxepochs)
   {
      numactive = 0;
      numconverged = 0;
      for(j = 0 ; j < p1; j++)
      {
	 iter = 0;
	 s_old = 0;
	 conv = TRUE;
	 if(g->active[j])
	 {
	    conv = FALSE;
	    g->nextcol(g, &sm, j, NA_ACTION_ZERO);
	    old_loss = g->loss;

	    /* iterate over jth variable */
      	    while(iter < maxiters)
      	    {
      	       s = step_func(&sm, g) * m[j];
      	       beta_new = g->beta[j] - s;
      	       if(j > 0) /* don't penalise intercept */
      	          beta_new = soft_threshold(beta_new, lambda1) * l2recip;
      	       beta_new = clip(beta_new, -truncl, truncl);
      	       beta_new = zero(beta_new, ZERO_THRESH);

	       /* beta_new may have changed by thresholding */
	       s = beta_new - g->beta[j];
	       updatelp(g, s, sm.x);
	       g->beta[j] = beta_new;

	       if(fabs(s) <= ZERO_THRESH
		     || fabs(old_loss - g->loss) / g->loss <= 1e-4
		     || g->loss <= 1e-6
		  )
	       {
		  conv = TRUE;
		  break;
	       }
	       s_old = s;
      	       iter++;
      	    }

	    /*printfverb("%d loss: %.10f\n", j, g->loss);*/
	 }

	 g->active[j] = (g->beta[j] != 0);
	 numactive += g->active[j];
	 numconverged += conv;

	 if(iter >= maxiters)
	    printfverb("max number of internal iterations (%d) \
reached for variable: %d  s: %.15f  old_loss: %.10f  loss: %.10f\n",
	       maxiters, j, s, old_loss, g->loss);
      }

      printfverb("fold: %d  epoch: %d  numactive: %d  numconverged: %d\n", 
	    g->fold, epoch, numactive, numconverged);
      fflush(stdout);

      /* 3-state machine for active set convergence */ 
      if(numconverged == p1) 
      {
	 printfverb("all converged\n");
	 allconverged++;

	 /* prepare for another iteration over *all*
	  * (non-ignored) variables, store a copy of the
	  * current active set for later */
	 if(allconverged == 1)
	 {
	    printfverb("prepare for final epoch\n");
	    for(j = p ; j >= 0 ; --j)
	    {
	       active_old[j] = g->active[j];
	       g->active[j] = !g->ignore[j];
	    }
	 }
	 else /* 2nd iteration over all variables done, check
		 whether active set has changed */
	 {
	    for(j = p ; j >= 0 ; --j)
	       if(g->active[j] != active_old[j])
		  break;

	    /* all equal, terminate */
	    if(j < 0)
	    {
	       printfverb("\n[%ld] terminating at epoch %d \
with %d active vars\n", time(NULL), epoch, numactive);
	       good = TRUE;
	       break;
	    }

	    printfverb("active set changed, %d active vars\n",
		  numactive);

	    /* active set has changed, iterate over
	     * new active variables */
	    for(j = p ; j >= 0 ; --j)
	       active_old[j] = g->active[j];
	    allconverged = 1;
	 }
      }
      else /* reset to first state */
	 allconverged = 0;

      epoch++;
   }
   printfverb("\n");

   /*free(beta_old);*/
   free(active_old);
   free(m);

   return good ? numactive : CDFAILURE;
}

