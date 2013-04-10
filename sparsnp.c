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
#include "util.h"

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
   long jp1k;
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
      //updatelp(g, beta_new, sm.x, 0, k);
      updateloss(g, 0, beta_new, sm.x, 0, k, 0, 0, 0); 

      if(g->verbose)
	 printf("[k=%d] intercept: %.15f (%d samples)\n", k, beta_new, n);
   
      /* find smallest lambda1 that makes all coefficients zero, by
       * finding the largest z, but let the intercept affect lp
       * first because it's not penalised. */
      for(j = 1 ; j < p1; j++)
      {
	 jp1k = j + p1 * k;
         if(g->ignore[jp1k])
	    continue;
	    
         if(!g->nextcol(g, &sm, j, NA_ACTION_PROPORTIONAL))
	    return FAILURE;
   
         step_func(&sm, g, k, &d1, &d2);
	 s = fabs(d1 / d2);
         zmax = (zmax < s) ? s : zmax;
	 g->grad_array[jp1k].value = s;
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
   long pkj, p1K1 = p1 * K - 1, kK1;
   int v1, v2, e, l;
   int nE = K * (K - 1) / 2;
   double d1, d2, df1, df2, sv;
   double beta_pkj;
   double lossold = g->loss;
   int mult = 2, idx;
   int mx;
   double *BCT = NULL;
   int i;
   double tmp;
   int m;

   if(!sample_init(&sm))
      return FAILURE;

   CALLOCTEST(beta_old, (long)p1 * K, sizeof(double));
   CALLOCTEST(active_old, (long)p1 * K, sizeof(int));

   CALLOCTEST(BCT, (long)p1 * nE, sizeof(double));

   //for(j = p1K1 ; j >= 0 ; --j)
     // active_old[j] = g->active[j] = !g->ignore[j];

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
      	    
      	    if(g->active[pkj])
      	    {
	       if(!g->nextcol(g, &sm, j, NA_ACTION_PROPORTIONAL))
		  return FAILURE;
	       
      	       step_func(&sm, g, k, &d1, &d2);

	       /* don't penalise intercept */
	       if(j > 0)
	       {
		  /* fusion penalty */
		  df1 = 0;
		  df2 = 0;

		  if(g->dofusion)
		  {
		     /* for each task k, compute derivative over
		      * the K-1 other tasks, as each task has K-1 pairs
		      */
		     for(l = K - 2 ; l >= 0 ; --l)
		     {
	       	        e = g->edges[l + kK1];
	       	        v1 = g->pairs[e];
	       	        v2 = g->pairs[e + nE];
	       	        sv = g->beta[j + v1 * p1] * g->C[e + v1 * nE]
	       	           + g->beta[j + v2 * p1] * g->C[e + v2 * nE];
	       	        df1 += sv * g->C[e + k * nE];
		     }

		     /* derivatives of fusion loss */
		     df1 *= gamma;
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

      	       delta = beta_new - beta_pkj;

	       if(delta != 0)
		  updateloss(g, beta_pkj, delta, sm.x, j, k, lambda1, gamma, nE); 
      	       
	       /* intercept always deemed active */
	       g->active[pkj] = (j == 0) || (g->beta[pkj] != 0);
      	    }

	    numactiveK[k] += g->active[pkj];
	    numactive += g->active[pkj];
      	 }
      }

      g->loss = 0;
      for(i = 0 ; i < g->ncurr ; i++)
      {
	 for(k = 0 ; k < g->K ; k++)
	 {
	    m = i + k * g->ncurr;
	    tmp = (g->lp[m] - g->y[m]) * (g->lp[m] - g->y[m]);
	    g->loss += tmp;
	 }
      }
      g->loss /= (2.0 * g->ncurr);

      g->l1loss = 0;
      for(j = 1 ; j < p1 ; j++)
	 for(k = 0 ; k < K ; k++)
	    g->l1loss += fabs(g->beta[j + k * p1]);
      g->loss += lambda1 * g->l1loss;
      

      if(g->dofusion)
      {
	 g->floss = 0;
	 for(j = p1 - 1 ; j >= 1 ; --j)
	 {
	    for(e = nE - 1 ; e >= 0 ; --e)
	    {
	       v1 = g->pairs[e];
	       v2 = g->pairs[e + nE];
	       sv = g->beta[j + v1 * p1] * g->C[e + v1 * nE] 
	          + g->beta[j + v2 * p1] * g->C[e + v2 * nE];
	       g->floss += sv * sv;
	    }
	 }

	 g->loss += gamma * 0.5 * g->floss;
      }

      if(fabs(g->loss - lossold) / fabs(lossold) < g->tol)
	 allconverged++;
      else
      {
	 allconverged = 0;
	 //printf("%d loss: %.6f lossold: %.6f\n", epoch, g->loss, lossold);
      }

      if(allconverged == 1)
      {
	 /* reset active-set to contain all (non monomorphic) coordinates, in
	  * order to check whether non-active coordinates become active again 
	  * or vice-versa */
         for(j = p1K1 ; j >= 0 ; --j)
         {
            active_old[j] = g->active[j];
            g->active[j] = !g->ignore[j];
         }
	 if(g->verbose)
	 {
	    timestamp();
	    printf(" resetting activeset at epoch %d, loss: %.6f floss: %.6f\n",
	       epoch, g->loss, g->floss);
	 }
	 mult = 2;
      }
      else if(allconverged == 2)
      {
         for(j = p1K1 ; j >= 0 ; --j)
            if(g->active[j] != active_old[j])
               break;

         if(j < 0)
         {
            if(g->verbose)
	    {
	       timestamp();
               printf(" terminating at epoch %d with %d active vars\n",
		  epoch, numactive);
	    }
            good = TRUE;
            break;
         }

         if(g->verbose)
	 {
	    timestamp();
            printf(" active set changed, %d active vars, mx:", numactive);
	 }

	 /* keep iterating over existing active set, keep track
	  * of the current active set */
         for(j = p1K1 ; j >= 0 ; --j)
            active_old[j] = g->active[j];

	 /* double the size of the active set */
	 for(k = 0 ; k < K ; k++)
	 {
	    mx = fminl(mult * numactiveK[k], p1);
	    printf("%d ", mx);
	    for(j = mx - 1 ; j >= 0 ; --j)
	    {
	       idx = g->grad_array[j + k * p1].index + k * p1;
	       g->active[idx] = !g->ignore[idx];
	    }
	 }

	 printf("\n");

         allconverged = 0;
	 mult *= 2;
      }     

      epoch++;
      lossold = g->loss;
   }

   if(g->verbose)
      printf("\n");

   FREENULL(beta_old);
   FREENULL(active_old);
   FREENULL(BCT);

   return good ? numactive : CDFAILURE;
}

void updateloss(gmatrix *g, double beta_pkj,
   double delta, double *x,
   int j, int k, double lambda1, double gamma, int nE)
{
   double lossoldk = g->lossK[k];
   double l1lossoldk = g->l1lossK[k];
   double beta_new = beta_pkj + delta;
   int p1 = g->p + 1;
   long pkj = (long)p1 * k + j;

   updatelp(g, delta, x, j, k); /* updates g->lossK[k] */
   g->beta[pkj] = beta_new;

   if(j > 0)
      g->l1lossK[k] += fabs(beta_new) - fabs(beta_pkj);

   g->loss += g->lossK[k] - lossoldk + lambda1 * (g->l1lossK[k] - l1lossoldk);
}

