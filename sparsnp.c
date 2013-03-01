/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include "common.h"
#include "sparsnp.h"

/* 
 * Find smallest lambda1 that makes all coefficients
 * zero (except the intercept)
 * 
 */
double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func, phi2 phi2_func, inv inv_func, step step_func)
{
   int i, j, n = g->ncurr, n1 = g->ncurr - 1;
   long p1 = g->p + 1;
   int k, K = g->K;
   double d1, d2, s, zmax = 0, beta_new;
   sample sm;

   if(!sample_init(&sm))
      return FAILURE;

   for(k = 0 ; k < K ; k++)
   {
      if(!g->nextcol(g, &sm, 0, NA_ACTION_PROPORTIONAL))
	 return FAILURE;
   
      /* First compute the intercept. When all other variables
       * are zero, the intercept is just inv(mean(y)) for a suitable inv()
       * function depending on the loss */
      s = 0;
      for(i = n1 ; i >= 0 ; --i)
         s += g->y[n * k + i];
   
      beta_new = inv_func(s / n);
      updatelp(g, beta_new, sm.x, 0, k);

      if(g->verbose)
	 printf("[k=%d] intercept: %.15f (%d samples)\n", k, beta_new, n);
   
      /* find smallest lambda1 that makes all coefficients zero, by
       * finding the largest z, but let the intercept affect lp
       * first because it's not penalised. */
      for(j = 1 ; j < p1; j++)
      {
         if(g->ignore[p1 * k + j])
	    continue;
	    
         if(!g->nextcol(g, &sm, j, NA_ACTION_PROPORTIONAL))
	    return FAILURE;
   
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
      const double lambda1,
      const double lambda2,
      const double gamma,
      const double trunc,
      int *numactiveK)
{
   const int p1 = g->p + 1, K = g->K;
   int j, k, allconverged = 0, numactive = 0,
       epoch = 1,
       good = FALSE;
   double beta_new, delta, s;
   const double l2recip = 1 / (1 + lambda2);
   sample sm;
   double *beta_old = NULL;
   int *active_old = NULL;
   long pkj, p1K1 = p1 * K - 1, kK1, l;
   int v1, v2, e, sv;
   int nE = K * (K - 1) / 2;
   double d1, d2, df1, df2;
   double beta_pkj;
   int dofusion = K > 1 && gamma > 0;

   if(!sample_init(&sm))
      return FAILURE;

   CALLOCTEST(beta_old, (long)p1 * K, sizeof(double));
   CALLOCTEST(active_old, (long)p1 * K, sizeof(int));

   for(j = p1K1 ; j >= 0 ; --j)
      active_old[j] = g->active[j];

   while(epoch <= maxepochs)
   {
      numactive = 0;

      for(k = 0 ; k < K ; k++)
      {
	 numactiveK[k] = 0;
	 kK1 = k * (K - 1);

	 for(j = 0 ; j < p1; j++)
      	 {
	    pkj = (long)p1 * k + j;
	    beta_pkj = g->beta[pkj];
      	    beta_new = beta_old[pkj] = beta_pkj;
      	    
      	    if(!g->active[pkj])
	    {
#ifdef DEBUG
	       printf("skipping inactive k=%d j=%d\n", k, j);
#endif
	    }
	    else
      	    {
	       if(!g->nextcol(g, &sm, j, NA_ACTION_PROPORTIONAL))
		  return FAILURE;

      	       step_func(&sm, g, k, &d1, &d2);

	       /* don't penalise intercept */
	       if(j > 0)
	       {
		  //s = beta_pkj - d1 / d2;

		  /* fusion penalty */
		  df1 = 0;
		  df2 = 0;

		  if(dofusion)
		  {
		     for(l = 0 ; l < K - 1 ; l++)
		     {
	       	        e = g->edges[l + kK1];
	       	        v1 = g->pairs[e];
	       	        v2 = g->pairs[e + nE];
	       	        sv = g->beta[j + v1 * p1] * g->C[e + v1 * nE]
	       	           + g->beta[j + v2 * p1] * g->C[e + v2 * nE];
	       	        df1 += sv * g->C[e + k * nE];
		     }

		     df1 *= gamma;

		     /* 2nd derivative of fusion loss */
	             df2 = gamma * g->diagCC[k];
		  }

	          s = beta_pkj - (d1 + df1) / (d2 + df2);

		  beta_new = soft_threshold(s, lambda1) * l2recip;
	       }
	       else
		  s = beta_new = beta_pkj - d1 / d2;

	       /* numerically close enough to zero */
	       if(fabs(beta_new) < ZERO_THRESH)
		  beta_new = 0;

	       /*beta_new = clip(beta_new, -truncl, truncl);*/

      	       delta = beta_new - beta_pkj;
      	       
#ifdef DEBUG
   if(j > 0)
   printf("[k=%d j=%d] d1=%.6f d2=%.6f s=%.6f delta=%.6f\
 beta_old=%.6f beta_new=%.6f\n",
		  k, j, d1, d2, s, delta, beta_pkj, beta_new);
#endif

      	  //     if(fabs(delta) > ZERO_THRESH)
      	       {
      	          updatelp(g, delta, sm.x, j, k);
      	          g->beta[pkj] = beta_new;
      	       }

      	       g->active[pkj] = (g->beta[pkj] != 0);
      	    }

	    numactiveK[k] += g->active[pkj];
	    numactive += g->active[pkj];
      	 }
#ifdef DEBUG
   	 printf("[end of task k=%d] numactive=%d (excl. %d intercept/s)\n",
	    k, numactive - (k + 1), k+1);

	 printf("----------------\n");
#endif
      }

      if(g->verbose)
      {
	 printf("fold: %d  epoch: %d  numactive: %d\n", 
	    g->fold, epoch, numactive);
	 fflush(stdout);
      }

      /* State machine for active set convergence */ 
      allconverged++;

      /* prepare for another iteration over *all*
       * (non-ignored) variables, store a copy of the
       * current active set for later */
      if(allconverged == 1)
      {
	 if(g->verbose)
	    printf("prepare for final epoch\n");
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
	    if(g->verbose)
	       printf("\n[%ld] terminating at epoch %d \
with %d active vars\n", time(NULL), epoch, numactive);
	    good = TRUE;
	    break;
	 }

	 if(g->verbose)
	    printf("active set changed, %d active vars\n",
	       numactive);

	 /* active set has changed, iterate over
	  * new active variables */
	 for(j = p1K1 ; j >= 0 ; --j)
	    active_old[j] = g->active[j];
	 allconverged = 1;
      }
     
      epoch++;
   }
   if(g->verbose)
      printf("\n");

   FREENULL(beta_old);
   FREENULL(active_old);

   return good ? numactive : CDFAILURE;
}

