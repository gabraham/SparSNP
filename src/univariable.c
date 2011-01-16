#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"
#include "svd.h"
#include "univariable.h"
#include "matrix.h"
#include "thin.h"

/*
 * Iteratively-Reweighted Least Squares for logistic regression
 * with l2 penalties
 *
 * \hat{\beta}_{t+1} = \hat{\beta}_t - (X^T W X + \lambda I)^{\dagger} X^T y
 *
 * dagger is pseudoinverse
 */
int irls(double *x, double *y, double *beta, double *invhessian,
      int n, int p, double lambda2, int verbose)
{
   int i, j, 
       iter = 1, maxiter = 30,
       converged = FALSE, diverged = FALSE,
       ret = SUCCESS;

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
	 printf("IRLS iter %d\tdev=%.3f\n", iter, loss);
      
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

      /* same convergence test as in R's glm.fit */
      if(fabs(loss - loss_old) / (fabs(loss) + 0.1) <= IRLS_THRESH)
      {
	 converged = TRUE;
	 break;
      }

      /* gradient */
      for(j = 0 ; j < p ; j++)
      {
	 grad[j] = 0.0;
	 for(i = n - 1; i >= 0 ; --i)
	    grad[j] += x[i * p + j] * (lp_invlogit[i] - y[i]);
      }

      /* weights */
      for(i = n - 1; i >= 0 ; --i)
	 w[i] = lp_invlogit[i] * (1 - lp_invlogit[i]);

      /* hessian */
      wcrossprod(x, x, w, hessian, n, p, p);

      /* Add l2 penalty to diagonal of X^T X, except the intercept */
      for(j = 1 ; j < p ; j++)
	 hessian[j * p + j] += lambda2;

      pseudoinverse(hessian, &p, &p, invhessian);

      /* Compute the Newton step */
      sqmvprod(invhessian, grad, s, p);

      if(loss <= IRLS_DEVIANCE_MIN)
      {
      	 diverged = TRUE;
	 break;
      }
      
      diverged = FALSE;
      for(j = 0 ; j < p ; j++)
      {
	 beta[j] -= s[j];
	 if(fabs(beta[j]) >= IRLS_THRESH_MAX)
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
      ret = IRLS_ERR_NO_CONVERGENCE;
   else if(diverged)
      ret = IRLS_ERR_DIVERGENCE;

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

/* Two stages:
 *
 * 1) Select the top k SNPs by univariable p values, using each
 * SNP and the intercept as covariables
 *
 * 2) Fit a multivariable model to the top k SNPs plus intercept
 */
int univar_gmatrix(Opt *opt, gmatrix *g, double *beta, double *zscore)
{
   int i, j, ret,
       p1 = g->p + 1;
   double beta2[2] = {0, 0};
   double *invhessian = NULL,
	  *x = NULL;
   sample sm;

   if(!sample_init(&sm))
      return FAILURE;

   MALLOCTEST(invhessian, sizeof(double) * 4);

   /* get p-values per SNP, skip intercept */
   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d\r", j);
      if(!g->active[j])
      {
	 printf("skipped variable %d\n", j);
	 continue;
      }
      
      g->nextcol(g, &sm, j, NA_ACTION_DELETE);
      CALLOCTEST(x, 2 * sm.n, sizeof(double));
      for(i = sm.n - 1 ; i >= 0 ; --i)
      {
	 x[2 * i] = 1.0;
	 x[2 * i + 1] = sm.x[i];
      }

      beta2[0] = beta2[1] = 0.0;
      ret = irls(x, sm.y, beta2, invhessian, sm.n, 2,
	    opt->lambda2_univar, FALSE);
      FREENULL(x);

      if(ret == FAILURE)
	 return FAILURE;
      else if(ret == IRLS_ERR_NO_CONVERGENCE)
      {
	 printf("IRLS didn't converge for variable %d\n", j);
	 beta[j] = zscore[j] = 0.0;
      }
      else
      {
	 /* z-score for the SNP only, ignore intercept */
	 zscore[j] = beta2[1] / sqrt(invhessian[3]);
	 beta[j] = beta2[1];

	 if(ret == IRLS_ERR_DIVERGENCE)
	    printf("IRLS diverged for variable %d, z=%.3f\n",
		  j, zscore[j]);
      }
   }
   printf("\n");

   FREENULL(invhessian);

   return SUCCESS;
}

int run_train(Opt *opt, gmatrix *g)
{
   int i, j, k,
       n = g->ncurr,
       p1 = g->p + 1,
       nums1 = 0,
       ret = 0,
       *numselected = NULL,
       *pselected = NULL;
   double *x = NULL,
	  *xthinned = NULL,
	  *zscore = NULL,
	  *beta = NULL,
	  *se = NULL,
	  *invhessian = NULL;
   char tmp[MAX_STR_LEN];
   int *activeselected = NULL,
       *activeselected_ind = NULL,
       *rets = NULL;

   CALLOCTEST(zscore, p1, sizeof(double));
   CALLOCTEST(beta, p1, sizeof(double));
   CALLOCTEST(numselected, opt->nzthresh, sizeof(int));
   CALLOCTEST(pselected, opt->nzthresh, sizeof(int));
   CALLOCTEST(rets, opt->nzthresh, sizeof(int));

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    g->ntrain[g->fold], g->ntest[g->fold]);

   /* read existing univar files */
   if(opt->existing_univar)
   {
      snprintf(tmp, MAX_STR_LEN, "univar_beta.00.%02d", g->fold);
      printf("reading %s\n", tmp);
      if(!load_beta(beta, tmp, g->p + 1))
	 return FAILURE;

      snprintf(tmp, MAX_STR_LEN, "univar_zscore.00.%02d", g->fold);
      printf("reading %s\n", tmp);
      if(!load_beta(zscore, tmp, g->p + 1))
	 return FAILURE;
   }
   else /* compute beta and zscore */
   {
      /* Get coefs and z-scores for each SNP using the
       * univariable method */
      if(!univar_gmatrix(opt, g, beta, zscore))
         return FAILURE;

      /* univariable coefs */
      snprintf(tmp, MAX_STR_LEN, "univar_beta.00.%02d", g->fold);
      if(!writevectorf(tmp, beta, g->p + 1))
         return FAILURE;

      /* univariable z scores */
      snprintf(tmp, MAX_STR_LEN, "univar_zscore.00.%02d", g->fold);
      if(!writevectorf(tmp, zscore, g->p + 1))
         return FAILURE;
   }

   printf("univariate selection done\n");

   if(opt->do_multivar)
   {
      for(i = 0 ; i < opt->nzthresh ; i++)
      {
	 gmatrix_zero_model(g); /* reset active variables */

	 numselected[i] = 0; 
	 g->active[0] = TRUE;
	 for(j = 1 ; j < p1 ; j++)
	 {
	    g->active[j] &= (fabs(zscore[j]) >= opt->zthresh[i]);
	    if(g->active[j])
	       numselected[i]++;
	 }
	 nums1 = numselected[i] + 1;
	 printf("total %d SNPs exceeded z-score=%.3f\n", numselected[i],
	       opt->zthresh[i]);

	 FREENULL(se);
	 /* standard error for ALL SNPs 0 to p1 */
	 CALLOCTEST(se, p1, sizeof(double));

	 FREENULL(x);
	 FREENULL(beta);
	 FREENULL(invhessian);

	 if(numselected[i] > 0)
	 {
	    MALLOCTEST(x, sizeof(double) * n * nums1);
	    CALLOCTEST(invhessian, nums1 * nums1, sizeof(double));

	    /* flags for columns of x that are active, ie a subset of g->active */
	    CALLOCTEST(activeselected, nums1, sizeof(int));

	    /* which of the columns 0:p1 are active */
	    CALLOCTEST(activeselected_ind, nums1, sizeof(int));

	    /* read the chosen variables into memory, including the intercept */
	    if(!gmatrix_read_matrix(g, x, g->active, nums1))
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
	    if(nums1 > 2)
	    {
	       /* activeselected[0] must be FALSE, we don't want to thin the
		* intercept */
	       activeselected[0] = FALSE;
	       if(!thin(x, n, nums1, activeselected, THIN_COR_MAX, nums1, nums1))
		  return FAILURE;

	       /* Count the remaining SNPs post thinning, add one for
		* intercept  */
	       activeselected_ind[0] = 0;
	       activeselected[0] = TRUE;
	       pselected[i] = 0;
	       for(j = 1 ; j < nums1 ; j++)
	          pselected[i] += activeselected[j];

	       printf("After thinning, %d of %d SNPs left (excluding intercept)\n",
		     pselected[i] - 1, numselected[i]);

	       MALLOCTEST(xthinned, sizeof(double) * n * (pselected[i] + 1));
	       copyshrink(x, xthinned, n, nums1, activeselected, pselected[i] + 1);
	    }
	    else /* no thinning */
	    {
	       xthinned = x;
	       pselected[i] = numselected[i];
	    }

	    /* coefs only for SNPs that survived thinning */
	    CALLOCTEST(beta, pselected[i] + 1, sizeof(double));

	    /* train un-penalised multivariable model on
	     * the selected SNPs, with lambda=0 */
	    rets[i] = irls(xthinned, g->y, beta, invhessian, n, pselected[i] + 1,
		  opt->lambda2_multivar, TRUE);
	    printf("IRLS returned %d\n", ret);	    

	    if(xthinned == x)
	    {
	       FREENULL(x);
	    }
	    else
	    {
	       FREENULL(x);
	       FREENULL(xthinned);
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
		  se[j] = sqrt(invhessian[k * nums1 + k]);
		  k++;
	       }
	    }
	 }

	 /* estimated coefficients */
	 snprintf(tmp, MAX_STR_LEN, "multivar_%s.%02d.%02d",
	       opt->beta_files[0], i, g->fold);
	 printf("writing %s\n", tmp);
	 if(!writevectorf(tmp, g->beta, g->p + 1))
	    return FAILURE;

	 /* standard errors */
	 snprintf(tmp, MAX_STR_LEN, "multivar_%s_se.%02d.%02d",
	       opt->beta_files[0], i, g->fold);
	 printf("writing %s\n", tmp);
	 if(!writevectorf(tmp, se, g->p + 1))
	    return FAILURE;
      }

      /* number of selected variables, pre thinning */
      snprintf(tmp, MAX_STR_LEN, "multivar_%s.%02d",
	    opt->numnz_file, g->fold);
      printf("writing %s\n", tmp);
      if(!writevectorl(tmp, numselected, opt->nzthresh))
	 return FAILURE;

      /* number of selected variables, post thinning */
      snprintf(tmp, MAX_STR_LEN, "multivar_%s_thinned.%02d",
	    opt->numnz_file, g->fold);
      printf("writing %s\n", tmp);
      if(!writevectorl(tmp, pselected, opt->nzthresh))
	 return FAILURE;

      /* IRLS exit codes */
      snprintf(tmp, MAX_STR_LEN, "multivar_irls.%02d", g->fold);
      printf("writing %s\n", tmp);
      if(!writevectorl(tmp, rets, opt->nzthresh))
	 return FAILURE;
   }

   FREENULL(beta);
   FREENULL(invhessian);
   FREENULL(se);
   FREENULL(zscore);
   FREENULL(activeselected);
   FREENULL(activeselected_ind);
   FREENULL(numselected);
   FREENULL(pselected);
   FREENULL(rets);

   return SUCCESS;
}

int do_train(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, k;

   if(!gmatrix_init(g, opt->filename, opt->n, opt->p,
	    NULL, opt->yformat, opt->model, opt->encoded,
	    opt->binformat, opt->folds_ind_file, opt->mode,
	    opt->loss_pt_func))
      return FAILURE;

   printf("%d CV folds\n", g->nfolds);

   /* cross-validation: training stage */
   if(g->nfolds > 1)
   {
      for(k = 0 ; k < g->nfolds ; k++)
      {
	 printf("CV fold: %d\n", k);
	 /*len = strlen(opt->scalefile) + 1 + 3;
	 snprintf(tmp, len, "%s.%02d", opt->scalefile, k);
	 g->scalefile = tmp;*/
	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 gmatrix_zero_model(g);
	 /*make_lambda1path(opt, g);*/
	 gmatrix_reset(g);

	 if(!(ret &= run_train(opt, g)))
	    break;

	 /*snprintf(tmp, 5, "y.%02d", k);
	 printf("writing y file: %s\n", tmp);
	 if(!(ret &= writevectorf(tmp, g->y, g->ncurr)))
	    break;*/
      }
   }
   else
   {
      /*g->scalefile = opt->scalefile;
      if(!gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;*/
      gmatrix_zero_model(g);
      /*make_lambda1path(opt, g);*/
      gmatrix_reset(g);
      if(!(ret = run_train(opt, g)))
	 return FAILURE;

      /*snprintf(tmp, 5, "y.%02d", 0);
      printf("writing y file: %s\n", tmp);
      if(!writevectorf(tmp, g->y, g->ncurr))
	 return FAILURE;*/
   }
   
   return ret;
}

/*
 * For each beta file, predict outcome using the chosen model
 */
int run_predict_beta(gmatrix *g, predict predict_func,
      char* predict_file)
{
   int i, j, n = g->ncurr, p1 = g->p + 1;
   sample sm;
   double *yhat;
   double *restrict lp = g->lp;
   double *restrict beta = g->beta;

   if(!sample_init(&sm))
      return FAILURE;

   CALLOCTEST(yhat, n, sizeof(double));

   for(j = 0 ; j < p1 ; j++)
   {
      if(beta[j] != 0)
      {
	 g->nextcol(g, &sm, j, NA_ACTION_ZERO);
	 for(i = 0 ; i < n ; i++)
	    lp[i] += sm.x[i] * beta[j];
      }
   }
   
   for(i = 0 ; i < n ; i++)
   {
      yhat[i] = predict_func(lp[i]);
      /*g->loss[ += g->loss_pt(yhat[i], g->y[i]);*/
   }

   printf("writing %s (%d) ... ", predict_file, n);
   if(!writevectorf(predict_file, yhat, n))
      return FAILURE;
   printf("done\n");

   FREENULL(yhat);
   
   return SUCCESS;
}

int run_predict(gmatrix *g, predict predict_func, char **beta_files,
      int n_beta_files)
{
   int i;
   char tmp[MAX_STR_LEN];

   for(i = 0 ; i < n_beta_files ; i++)
   {
      gmatrix_zero_model(g);
      printf("reading %s\n", beta_files[i]);
      if(!load_beta(g->beta_orig, beta_files[i], g->p + 1))
      {
	 printf("skipping %s\n", beta_files[i]);
	 continue;
      }

      /* beta should already be on original scale of data */
      for(int j = 0 ; j < g->p + 1 ; j++)
	 g->beta[j] = g->beta_orig[j];

      snprintf(tmp, MAX_STR_LEN, "%s.pred", beta_files[i]);
      if(!run_predict_beta(g, predict_func, tmp))
	 return FAILURE;
   }

   return SUCCESS;
}

int do_predict(gmatrix *g, Opt *opt, char *tmp)
{
     int ret = SUCCESS, b, k, len;

   if(!gmatrix_init(g, opt->filename, opt->n, opt->p,
	    NULL, opt->yformat, opt->model, opt->encoded,
	    opt->binformat, opt->folds_ind_file, opt->mode,
	    opt->loss_pt_func))
      return FAILURE;

   if(g->nfolds > 1)
   {
      if(!opt->beta_files_fold)
	 MALLOCTEST2(opt->beta_files_fold,
	       sizeof(char*) * opt->n_beta_files);
      for(b = 0 ; b < opt->n_beta_files ; b++)
	 opt->beta_files_fold[b] = NULL;

      /* cross-validation: prediction stage */
      for(k = 0 ; k < g->nfolds ; k++)
      {
	 /*len = strlen(opt->scalefile) + 1 + 3;
	 snprintf(tmp, len, "%s.%02d", opt->scalefile, k);
	 g->scalefile = tmp;
	 printf("reading scale file: %s\n", tmp);*/
	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 /* write y file */
	 snprintf(tmp, 14, "multivar_y.%02d", k);
	 printf("writing y file: %s\n", tmp);
	 if(!(ret &= writevectorf(tmp, g->y, g->ncurr)))
	    break;

	 /*gmatrix_zero_model(&g);*/
	 gmatrix_reset(g);

	 /* set up correct file names */
	 for(b = 0 ; b < opt->n_beta_files ; b++)
	 {
	    len = strlen(opt->beta_files[b]);
	    if(!opt->beta_files_fold[b])
	       MALLOCTEST2(opt->beta_files_fold[b], len + 1 + 3);
	    snprintf(opt->beta_files_fold[b], MAX_STR_LEN, "%s.%02d",
		  opt->beta_files[b], k);
	 }

	 if(!(ret &= run_predict(g, opt->predict_func,
		     opt->beta_files_fold, opt->n_beta_files)))
	    break;
      }
   }
   else
   {
      /*g->scalefile = opt->scalefile;
      if(!gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;*/
      gmatrix_zero_model(g);
      ret = run_predict(g, opt->predict_func,
	    opt->beta_files, opt->n_beta_files);
   }

   return ret;
}

int main(int argc, char* argv[])
{
   int ret = FAILURE;

   Opt opt;
   gmatrix g;
   char tmp[MAX_STR_LEN];

   setbuf(stdout, NULL);

   if(!opt_defaults(&opt, OPTIONS_CALLER_UNIVARIABLE) 
	 || !opt_parse(argc, argv, &opt))
   {
      opt_free(&opt);
      return EXIT_FAILURE;
   }

   if(opt.mode == MODE_TRAIN && !opt.nofit)
      ret = do_train(&g, &opt, tmp);
   else if(opt.mode == MODE_PREDICT)
      ret = do_predict(&g, &opt, tmp);

   gmatrix_free(&g);
   opt_free(&opt);

   return ret == FAILURE ? EXIT_FAILURE : EXIT_SUCCESS;
}

