#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"
#include "svd.h"
#include "multivariable.h"
#include "matrix.h"
#include "thin.h"

/*
 * Newton's method for logistic regression with l2 penalties
 *
 * \hat{\beta}_{t+1} = \hat{\beta}_t - (X^T W X + \lambda I)^{\dagger} X^T y
 *
 * dagger is pseudoinverse
 *
 * The argument p includes the intercept
 */
int newton(double *x, double *y, double *beta, double *invhessian,
      int n, int p, double lambda2, int verbose)
{
   int i, j, 
       iter = 1, maxiter = 100,
       converged = FALSE, diverged = FALSE,
       ret = NEWTON_SUCCESS;

   double *grad = NULL,
	  *hessian = NULL,
	  *lp = NULL,
	  *lp_invlogit = NULL,
	  *w = NULL,
	  *s = NULL;

   double loss_old = 1e6, loss = 0;

   MALLOCTEST(lp, sizeof(double) * n);
   MALLOCTEST(lp_invlogit, sizeof(double) * n);
   MALLOCTEST(grad, sizeof(double) * p);
   CALLOCTEST(hessian, p * p, sizeof(double));
   CALLOCTEST(w, n, sizeof(double));
   CALLOCTEST(s, p, sizeof(double));

   while(iter <= maxiter) 
   {
      if(verbose)
	 printf("NEWTON iter %d\tdev=%.7f\n", iter, loss);
      
      if(iter > 1)
	 loss_old = loss;
      loss = 0;

      /* setup the linea predictors, and compute the loss */
      for(i = n - 1; i >= 0 ; --i)
      {
	 lp[i] = x[i * p] * beta[0];

	 for(j = 1 ; j < p ; j++)
	 {
	    lp[i] += x[i * p + j] * beta[j];
	    lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	 }

	 loss += log(1 + exp(lp[i])) - y[i] * lp[i];
      }

      loss /= n;

      /* same convergence test as in R's glm.fit */
      if(fabs(loss - loss_old) / (fabs(loss) + 0.1) <= NEWTON_THRESH)
      {
	 converged = TRUE;
	 break;
      }

      /* gradient */
      for(j = 0 ; j < p ; j++)
      {
	 grad[j] = 0.0;
	 for(i = n - 1 ; i >= 0 ; --i)
	    grad[j] += x[i * p + j] * (lp_invlogit[i] - y[i]);
      }

      /* weights */
      for(i = n - 1 ; i >= 0 ; --i)
	 w[i] = lp_invlogit[i] * (1 - lp_invlogit[i]);

      /* hessian */
      wcrossprod(x, x, w, hessian, n, p, p);

      /* Add l2 penalty to diagonal of X^T X, except the intercept */
      /*for(j = 1 ; j < p ; j++)
	 hessian[j * p + j] += lambda2;*/

      pseudoinverse(hessian, &p, &p, invhessian);

      /* Compute the Newton step */
      sqmvprod(invhessian, grad, s, p);

      if(loss <= NEWTON_DEVIANCE_MIN)
      {
      	 diverged = TRUE;
	 break;
      }
      
      diverged = FALSE;
      for(j = 0 ; j < p ; j++)
      {
	 /* Newton step */
	 beta[j] -= s[j];
	 /*printf("beta[%d]: %.5f\n", j, beta[j]);*/
	 if(fabs(beta[j]) >= NEWTON_THRESH_MAX)
	 {
	    diverged = TRUE;
	    break;
	 }
      }

      if(diverged)
	 break;

      /*if(converged)
	 break;*/

      iter++;
   }

   if(iter >= maxiter)
      ret = NEWTON_ERR_NO_CONVERGENCE;
   else if(diverged)
      ret = NEWTON_ERR_DIVERGENCE;

   FREENULL(grad);
   FREENULL(hessian);
   FREENULL(lp);
   FREENULL(lp_invlogit);
   FREENULL(w);
   FREENULL(s);

   return ret;
}

/* 
 * Evaluate the 2 by 2 Hessian of logistic loss at beta
 *
 * in R:
 *
 * Q <- diag(drop(plogis(x %*% beta) * (1 - plogis(x %*% beta)))),
 *    an n by n diagonal matrix
 *
 * H <- crossprod(x, Q) %*% x, a 2 by 2 matrix
 *
 */
int make_hessian(double *hessian, double *x,
      double beta_intercept, double beta, int n)
{
   int k;

   double P, q, z;

   /* row-major ordering */
   for(k = 0 ; k < n ; k++)
   {
      P = 1 / (1 + exp(-beta_intercept - beta * x[k]));
      q = P * (1 - P);
      hessian[0] += q;
      z = x[k] * q;
      hessian[1] += z;
      hessian[2] = hessian[1];
      hessian[3] += x[k] * z;
   }

   return SUCCESS;
}

int multivariable_newton(Opt *opt, gmatrix *g, int nums1,
      int *pselected, int *numselected, int *rets)
{
   int n = g->ncurr, j, k, p1 = g->p + 1;
   double *invhessian = NULL,
	  *beta = NULL;
	  /**se = NULL;*/
   int *activeselected = NULL,
       *activeselected_ind = NULL;

   MALLOCTEST(g->x, sizeof(double) * n * nums1);
   CALLOCTEST(invhessian, nums1 * nums1, sizeof(double));

   /* flags for columns of x that are active, ie a subset of g->active */
   CALLOCTEST(activeselected, nums1, sizeof(int));

   /* which of the columns 0:p1 are active */
   CALLOCTEST(activeselected_ind, nums1, sizeof(int));

   /* read the chosen variables into memory, including the intercept */
   if(!gmatrix_read_matrix(g, g->active, nums1))
      return FAILURE;

   /* Slice the active vector for this window. We start off by
    * making all variables active, then thin() will make some
    * inactive. Ignore intercept for thinning. */
   k = 1;
   for(j = 1 ; j < p1 ; j++)
   {
      if(g->active[j])
      {
	 activeselected[k] = TRUE;
	 activeselected_ind[k] = j;
	 k++;
      }
   }

   /* thin the SNPs based on correlation, but only if there are at
    * least two SNPs (don't forget intercept adds one) */
   if(nums1 > 2 && opt->do_thinning)
   {
      /* activeselected[0] must be FALSE, we don't want to thin the
       * intercept */
      activeselected[0] = FALSE;
      if(!thin(g->x, n, nums1, activeselected, THIN_COR_MAX, nums1, nums1))
	 return FAILURE;

      /* Count the remaining SNPs post thinning, add one for
       * intercept (pselected doesn't include intercept)  */
      activeselected_ind[0] = 0;
      activeselected[0] = TRUE;
      *pselected = 0;
      for(j = 1 ; j < nums1 ; j++)
	 *pselected += activeselected[j];

      printf("After thinning, %d of %d SNPs left\n",
	    *pselected, *numselected);

      MALLOCTEST(g->xthinned, sizeof(double) * n * (*pselected + 1));
      copyshrink(g->x, g->xthinned, n, nums1, activeselected, (*pselected) + 1);
   }
   else /* no thinning */
   {
      g->xthinned = g->x;
      *pselected = *numselected;
   }

   /* coefs only for SNPs that survived thinning */
   CALLOCTEST(beta, (*pselected) + 1, sizeof(double));
   /*CALLOCTEST(se, (*pselected) + 1, sizeof(double));*/

   /* train un-penalised multivariable model on
    * the selected SNPs, with lambda=0 */
   *rets = newton(g->xthinned, g->y, beta, invhessian, n, (*pselected) + 1,
	 opt->lambda2_multivar, TRUE);
   printf("NEWTON returned %d\n", *rets);	    

   if(g->xthinned == g->x)
   {
      FREENULL(g->x);
   }
   else
   {
      FREENULL(g->x);
      FREENULL(g->xthinned);
   }

   /* Copy estimated beta back to the array for all SNPs. Due to
    * thinning, not all columns were used (pselected <= nums1),
    * so check if they were */
   k = 0; /* should run up to pselected */
   activeselected[0] = TRUE;
   for(j = 0 ; j < nums1 ; j++)
   {
      if(activeselected[j])
      {
	 g->beta[activeselected_ind[j]] = beta[k];
	 /*se[j] = sqrt(invhessian[k * nums1 + k]);*/
	 k++;
      }
   }

   FREENULL(activeselected);
   FREENULL(activeselected_ind);
   FREENULL(beta);
   FREENULL(invhessian);
   /*FREENULL(se);*/

   return SUCCESS;
}

int multivariable_lasso(Opt *opt, gmatrix *g, int threshind)
{
   int i, ret;
   char tmp[MAX_STR_LEN];

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    g->ntrain[g->fold], g->ntest[g->fold]);
   
   /* first numnz is always zero by definition */
   CALLOCTEST(g->numnz, opt->nlambda1, sizeof(int));

   /* don't start from zero, getlambda1max already computed that */
   for(i = 1 ; i < opt->nlambda1 ; i++)
   {
      if(opt->verbose)
	 printf("\nFitting with lambda1=%.20f\n", opt->lambda1path[i]);

      /* return value is number of nonzero variables,
       * including the intercept */
      ret = cd_gmatrix(
	    g, opt->step_func,
	    opt->maxepochs, opt->maxiters,
	    opt->lambda1path[i], opt->lambda2,
	    opt->threshold, opt->verbose, opt->trunc);

      if(ret == CDFAILURE)
      {
	 printf("failed to converge after %d epochs\n", opt->maxepochs);
	 break;
      } 

      g->numnz[i] = ret;
      gmatrix_reset(g);

      if(opt->caller == OPTIONS_CALLER_CD)
      {
	 snprintf(tmp, MAX_STR_LEN, "%s.%02d.%02d",
	       opt->beta_files[0], g->fold, i);
      }
      else
      {
	 snprintf(tmp, MAX_STR_LEN, "multivar_%s.%02d.%02d.%02d",
	       opt->beta_files[0], g->fold, threshind, i);
      }

      if(opt->unscale)
      {
	 unscale_beta(g->beta_orig, g->beta, g->mean, g->sd, g->p + 1);
	 if(!writevectorf(tmp, g->beta_orig, g->p + 1))
	    return FAILURE;
      }
      else if(!writevectorf(tmp, g->beta, g->p + 1))
	 return FAILURE;

      if(opt->nzmax != 0 && opt->nzmax <= ret - 1)
      {
	 printf("maximum number of non-zero variables \
reached or exceeded: %d\n", opt->nzmax);
	 i++; /* increment to correct number of models fitted successfully */
	 break;
      }
   }

   /* filename for number of non-zero variables */
   if(opt->caller == OPTIONS_CALLER_CD)
   {
      snprintf(tmp, MAX_STR_LEN, "%s.%02d", opt->numnz_file, g->fold);
   }
   else
   {
      snprintf(tmp, MAX_STR_LEN, "multivar_%s.%02d.%02d.%2d",
	    opt->numnz_file, g->fold, threshind, i);
   }

   /* number of non-zero variables for each successful fit and the all-zero
    * fit, excluding intercept */
   if(!writevectorl(tmp, g->numnz, i))
      return FAILURE;

   FREENULL(g->numnz);

   return SUCCESS;
}

/*
 * Creates a vector of lambda1 penalties
 */
int make_lambda1path(Opt *opt, gmatrix *g, int threshind)
{
   int i;
   double s;
   char tmp[MAX_STR_LEN];

   if(opt->lambda1 >= 0)
   {
      opt->lambda1max = opt->lambda1;
      opt->l1minratio = 1;
      opt->nlambda1 = 1;
   }
   else
   {
      /* create lambda1 path */
      /* get lambda1 max */
      opt->lambda1max = get_lambda1max_gmatrix(g,
	    opt->inv_func, opt->step_func);
      if(opt->verbose)
	 printf("lambda1max: %.20f\n", opt->lambda1max);
      opt->lambda1path[0] = opt->lambda1max;
   }
   
   opt->lambda1min = opt->lambda1max * opt->l1minratio;
   opt->lambda1path[opt->nlambda1 - 1] = opt->lambda1min;
   s = (log2(opt->lambda1max) - log2(opt->lambda1min)) / opt->nlambda1; 
   for(i = 1 ; i < opt->nlambda1 ; i++)
      opt->lambda1path[i] = pow(2, log2(opt->lambda1max) - s * i);

   /* Write the coefs for model with intercept only */
   if(opt->caller == OPTIONS_CALLER_CD)
      snprintf(tmp, MAX_STR_LEN, "%s.%02d.%02d",
	    opt->beta_files[0], g->fold, 0);
   else
      snprintf(tmp, MAX_STR_LEN, "multivar_%s.%02d.%02d.%2d",
	    opt->beta_files[0], g->fold, threshind, i);

   if(opt->unscale)
   {
      unscale_beta(g->beta_orig, g->beta, g->mean, g->sd, g->p + 1);
      if(!writevectorf(tmp, g->beta_orig, g->p + 1))
	 return FAILURE;
   }
   else
      if(!writevectorf(tmp, g->beta, g->p + 1))
	 return FAILURE;

   if(opt->caller == OPTIONS_CALLER_CD)
      snprintf(tmp, MAX_STR_LEN, "%s.%02d", opt->lambda1pathfile, g->fold);
   else
      snprintf(tmp, MAX_STR_LEN, "%s.%02d.%02d", opt->lambda1pathfile,
	    g->fold, threshind);

   return writevectorf(tmp, opt->lambda1path, opt->nlambda1);
}

