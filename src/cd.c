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
   updatelp(g, beta_new, sm.x, 0);
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

/* coordinate descent */
int cd_gmatrix(gmatrix *g,
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
   int i, j, iter, allconverged = 0, numactive = 0,
       epoch = 1, numconverged = 0,
       good = FALSE;
   double s = 0, beta_new;
   const double truncl = log((1 - trunc) / trunc);
   /*const double l2recip = 1 / (1 + lambda2);*/
   sample sm;
   double old_loss = 0;
   int conv = 0;
   int *zero = NULL;
   double *beta_old = NULL;
   double d = 0;

   if(!sample_init(&sm))
      return FAILURE;

   CALLOCTEST(zero, p1, sizeof(int));
   CALLOCTEST(beta_old, p1, sizeof(double));
   

   while(epoch <= maxepochs)
   {
      numactive = 0;
      numconverged = 0;

      for(j = 0 ; j < p1; j++)
      {
	 iter = 0;
	 
	 if(g->active[j])
	 {
	    g->nextcol(g, &sm, j, NA_ACTION_ZERO);

      	    while(iter < maxiters)
      	    {
	       s = step_func(&sm, g);
	       beta_new = g->beta[j] - s;
	       
      	       if(j > 0)
      	          beta_new = soft_threshold(beta_new, lambda1); /* l2recip;*/
      	       /*beta_new = clip(beta_new, -truncl, truncl); */

	       beta_new = g->beta[j] - s;
      	       if(j > 0) 
      	          beta_new = soft_threshold(beta_new, lambda1);

	       s = beta_new - g->beta[j];
	       
	       if(s == 0)
		  break;

	       updatelp(g, s, sm.x, j);
	       g->beta[j] = beta_new;

	       if(fabs(s) <= 1e-4
		     || fabs(old_loss - g->loss) / g->loss <= 1e-2
		     || g->loss <= 1e-9
		  )
	       {
		  break;
	       }
      	       iter++;
      	    }
	    zero[j] = g->beta[j] == 0;
	    g->active[j] = !zero[j];
	 }

	 conv = fabs(g->beta[j] - beta_old[j]) <= 1e-4;
	 beta_old[j] = g->beta[j];
	 numconverged += conv;
	 numactive += g->active[j];
      }

      printfverb("fold: %d  epoch: %d  numactive: %d  numconverged: %d\n", 
	    g->fold, epoch, numactive, numconverged);
      fflush(stdout);

      if(numconverged == p1) 
      {
	 printfverb("all converged\n");
	 allconverged++;

	 if(allconverged == 1)
	 {
	    printfverb("prepare for final epoch\n");
	    for(j = p ; j >= 0 ; --j)
	       g->active[j] = !g->ignore[j];
	 }
	 else
	 {
	    printfverb("\n[%ld] terminating at epoch %d \
 with %d active vars\n", time(NULL), epoch, numactive);
	    good = TRUE;
	    break;
	 }
      }
      else
      { 
	 allconverged = 0;
	 for(j = p ; j >= 0 ; --j)
	    g->active[j] = !g->ignore[j];
      }

      epoch++;
   }
   printfverb("\n");

   FREENULL(beta_old);
   FREENULL(zero);

   return good ? numactive : CDFAILURE;
}

