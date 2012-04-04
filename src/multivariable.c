#include "common.h"
#include "cd.h"
#include "util.h"
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
   int i, j, ip, 
       iter = 1, maxiter = 100,
       diverged = FALSE,
       ret = NEWTON_SUCCESS;

   double *w = NULL,
	  *grad = NULL,
	  *hessian = NULL,
	  *lp = NULL,
	  *lp_invlogit = NULL,
	  *s = NULL;

   double loss_old = 1e6, loss = 0;

   CALLOCTEST(w, n, sizeof(double));
   MALLOCTEST(grad, sizeof(double) * p);
   CALLOCTEST(hessian, p * p, sizeof(double));
   CALLOCTEST(s, p, sizeof(double));
   MALLOCTEST(lp, sizeof(double) * n);
   MALLOCTEST(lp_invlogit, sizeof(double) * n);

   while(iter <= maxiter) 
   {
      if(verbose)
	 printf("NEWTON iter %d\tdev=%.7f\n", iter, loss);
      
      if(iter > 1)
	 loss_old = loss;
      loss = 0;

      /* setup the linear predictors, and compute the loss */
      for(i = n - 1; i >= 0 ; --i)
      {
	 ip = i * p;
	 lp[i] = x[ip] * beta[0];

	 for(j = 1 ; j < p ; j++)
	 {
	    lp[i] += x[ip + j] * beta[j];
	    lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	 }

	 loss += log(1 + exp(lp[i])) - y[i] * lp[i];
      }

      loss /= n;

      /* same convergence test as in R's glm.fit */
      if(fabs(loss - loss_old) / (fabs(loss) + 0.1) <= NEWTON_THRESH)
      {
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
	 if(fabs(beta[j]) >= NEWTON_THRESH_MAX)
	 {
	    diverged = TRUE;
	    break;
	 }
      }

      if(diverged)
	 break;

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
 * Evaluate the 2 by 2 Hessian at beta
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

      printf("After thinning, %d of %d SNPs left (excluding intercept)\n",
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
    * thinning, not all columns used (pselected <= nums1),
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

/*int multivariable_lasso(Opt *opt, gmatrix *g, int nums1,
      int *pselected, int *numselected, int *rets)
{
   int ret;
   ret = cd_gmatrix(
	    g, opt->phi1_func, opt->phi2_func,
	    opt->step_func,
	    opt->maxepochs, opt->maxiters,
	    opt->lambda1path[i], opt->lambda2,
	    opt->threshold, opt->verbose, opt->trunc);
 
   return SUCCESS;
}*/

