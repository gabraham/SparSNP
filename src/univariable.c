#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"
#include "univariable.h"
#include "multivariable.h"
#include "matrix.h"
#include "thin.h"

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
      ret = newton(x, sm.y, beta2, invhessian, sm.n, 2,
	    opt->lambda2_univar, FALSE);
      FREENULL(x);

      if(ret == FAILURE)
	 return FAILURE;
      else if(ret == NEWTON_ERR_NO_CONVERGENCE)
      {
	 printf("Newton didn't converge for variable %d\n", j);
	 beta[j] = zscore[j] = 0.0;
      }
      else
      {
	 /* z-score for the SNP only, ignore intercept */
	 zscore[j] = beta2[1] / sqrt(invhessian[3]);
	 beta[j] = beta2[1];

	 if(ret == NEWTON_ERR_DIVERGENCE)
	    printf("NEWTON diverged for variable %d, z=%.3f\n",
		  j, zscore[j]);
      }
   }
   printf("\n");

   FREENULL(invhessian);

   return SUCCESS;
}

int run_train(Opt *opt, gmatrix *g)
{
   int i, j,
       p1 = g->p + 1,
       nums1 = 0,
       *numselected = NULL,
       *pselected = NULL;
   double *x = NULL,
	  *zscore = NULL,
	  *beta = NULL,
	  *se = NULL,
	  *invhessian = NULL;
   char tmp[MAX_STR_LEN];
   int *rets = NULL;

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
	 printf("Threshold %d: %.5f\n", i, opt->zthresh[i]);
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
	    if(opt->multivar == OPTIONS_MULTIVAR_NEWTON)
	       multivariable_newton(opt, g, nums1,
		     pselected + i, numselected + i, rets + i);
	    else
	       multivariable_lasso(opt, g, nums1,
		     pselected + i, numselected + i, rets + i);
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

	 /* no point in testing looser z-scores since they won't converge as
	  * well */
	 if(rets[i] == NEWTON_ERR_NO_CONVERGENCE)
	    break;

	 printf("\n");
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

      /* Newton exit codes */
      snprintf(tmp, MAX_STR_LEN, "multivar_status.%02d", g->fold);
      printf("writing %s\n", tmp);
      if(!writevectorl(tmp, rets, opt->nzthresh))
	 return FAILURE;
   }

   FREENULL(beta);
   FREENULL(invhessian);
   FREENULL(se);
   FREENULL(zscore);
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
	    opt->loss_pt_func, opt->subset_file))
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
	    opt->loss_pt_func, opt->subset_file))
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

