/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include "common.h"
#include "sparsnp.h"

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

/* 
 * Find smallest lambda1 that makes all coefficients
 * zero (except the intercept)
 */
double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func, phi2 phi2_func, inv inv_func, step step_func)
{
   int i, j, n = g->ncurr, n1 = g->ncurr - 1, p1 = g->p + 1;
   int k, K = g->K;
   double d1, d2, s, zmax = 0, beta_new;
   sample sm;

   if(!sample_init(&sm))
      return FAILURE;

   for(k = 0 ; k < K ; k++)
   {
      g->nextcol(g, &sm, 0, NA_ACTION_RANDOM);
   
      /* First compute the intercept. When all other variables
       * are zero, the intercept is just inv(mean(y)) for a suitable inv()
       * function depending on the loss */
      s = 0;
      for(i = n1 ; i >= 0 ; --i)
         s += g->y[n * k + i];
   
      beta_new = inv_func(s / n);
      updatelp(g, beta_new, sm.x, 0, k);
      printf("[k=%d] intercept: %.15f (%d samples)\n", k, beta_new, n);
   
      /* find smallest lambda1 that makes all coefficients zero, by
       * finding the largest z, but let the intercept affect lp
       * first because it's not penalised. */
      for(j = 1 ; j < p1; j++)
      {
         if(g->ignore[p1 * k + j])
	    continue;
         g->nextcol(g, &sm, j, NA_ACTION_RANDOM);
   
         step_func(&sm, g, k, &d1, &d2);
	 s = fabs(d1 / d2);
         zmax = (zmax < s) ? s : zmax;
      } 
   
      /* intercept */
      g->beta[k * p1] = beta_new;
   }

   return zmax;
}

/* Coordinate descent, multi-task
 *
 * This function will only work for quadratic functions like linear loss and
 * squared hinge loss, since it assumes that the Newton step for each variable
 * is exact. It's easy to modify it to loop over the Newton steps until some
 * convergence is achieved, e.g., for logistic loss.
 * */
int cd_gmatrix(gmatrix *g,
      step step_func,
      const int maxepochs,
      const int maxiters,
      const double *C,
      const double lambda1,
      const double lambda2,
      const double gamma,
      const int verbose,
      const double trunc)
{
   const int p1 = g->p + 1, K = g->K;
   int j, k, allconverged = 0, numactive = 0,
       epoch = 1,
       good = FALSE;
   double beta_new, delta;
   const double truncl = log((1 - trunc) / trunc);
   const double l2recip = 1 / (1 + lambda2);
   sample sm;
   double *beta_old = NULL;
   int *active_old = NULL;
   int pkj, p1K1 = p1 * K - 1;
   double d1, d2, pd1 = 0, pd2 = 0;
   int nE = K * (K - 1) / 2, e;
   double Ckne, Ckne2;
   double beta_pkj;

   if(!sample_init(&sm))
      return FAILURE;

   CALLOCTEST(beta_old, p1 * K, sizeof(double));
   CALLOCTEST(active_old, p1 * K, sizeof(int));

   for(j = p1K1 ; j >= 0 ; --j)
      active_old[j] = g->active[j];

   while(epoch <= maxepochs)
   {
      numactive = 0;

      for(k = 0 ; k < K ; k++)
      {
	 for(j = 0 ; j < p1; j++)
      	 {
	    pkj = p1 * k + j;
	    beta_pkj = g->beta[pkj];
      	    beta_new = beta_old[pkj] = beta_pkj;
      	    
      	    if(g->active[pkj])
      	    {
      	       g->nextcol(g, &sm, j, NA_ACTION_RANDOM);

      	       step_func(&sm, g, k, &d1, &d2);

	       /* Apply inter-task penalty, summing over all the edges */
	       pd1 = 0;
	       pd2 = 0;
	       for(e = 0 ; e < nE ; e++)
	       {
	          Ckne = C[k * nE + e];
		  Ckne2 = Ckne * Ckne;
		  pd1 += Ckne2 * beta_pkj;
		  pd2 += Ckne2;
	       }
	       pd1 *= 2;
	       pd2 *= 2;

      	       beta_new = beta_pkj - (d1 + pd1) / (d2 + pd2);
      	       
	       if(j > 0)
		  beta_new = soft_threshold(beta_new, lambda1) * l2recip;
	       beta_new = clip(beta_new, -truncl, truncl);

      	       delta = beta_new - beta_pkj;
      	       
      	       if(fabs(delta) > ZERO_THRESH)
      	       {
      	          updatelp(g, delta, sm.x, j, k);
      	          g->beta[pkj] = beta_new;
      	       }

      	       g->active[pkj] = (g->beta[pkj] != 0);
      	    }

      	    numactive += g->active[pkj];
      	 }
      }

      printfverb("fold: %d  epoch: %d  numactive: %d\n", 
	    g->fold, epoch, numactive);
      fflush(stdout);

      /* State machine for active set convergence */ 
      allconverged++;

      /* prepare for another iteration over *all*
       * (non-ignored) variables, store a copy of the
       * current active set for later */
      if(allconverged == 1)
      {
	 printfverb("prepare for final epoch\n");
	 for(j = p1K1 ; j >= 0 ; --j)
	 {
	    active_old[j] = g->active[j];
	    g->active[j] = !g->ignore[j];
	 }
      }
      else /* 2nd iteration over all variables done, check
	      whether active set has changed */
      {
	 for(j = p1K1 ; j >= 0 ; --j)
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
	 for(j = p1K1 ; j >= 0 ; --j)
	    active_old[j] = g->active[j];
	 allconverged = 1;
      }
     
      epoch++;
   }
   printfverb("\n");

   FREENULL(beta_old);
   FREENULL(active_old);

   return good ? numactive : CDFAILURE;
}

