/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <stdlib.h>
#include <libgen.h>
#include "cd.h"
#include "util.h"

#define OPTIONS_CALLER cd

/*
 * Creates a vector of lambda1 penalties
 */
int make_lambda1path(Opt *opt, gmatrix *g)
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
      opt->lambda1max = get_lambda1max_gmatrix(g, opt->phi1_func,
	    opt->phi2_func, opt->inv_func, opt->step_func);
      if(opt->verbose)
	 printf("lambda1max: %.20f\n", opt->lambda1max);
      opt->lambda1path[0] = opt->lambda1max;
   }
   
   opt->lambda1min = opt->lambda1max * opt->l1minratio;
   opt->lambda1path[opt->nlambda1 - 1] = opt->lambda1min;
   s = (log10(opt->lambda1max) - log10(opt->lambda1min)) / (opt->nlambda1 - 1);
   for(i = 1 ; i < opt->nlambda1 ; i++)
      opt->lambda1path[i] = pow(10, log10(opt->lambda1max) - s * i);

   if(opt->verbose)
      printf("lambda1min: %.20f\n", opt->lambda1path[i-1]);

   /* Write the coefs for model with intercept only */
   snprintf(tmp, MAX_STR_LEN, "%s.%02d.%02d",
	 opt->beta_files[0], 0, g->fold);

   unscale_beta(g->beta_orig, g->beta, g->mean, g->sd, g->p + 1);
   if(!write_beta_sparse(tmp, g->beta_orig, g->p + 1))
      return FAILURE;

   snprintf(tmp, MAX_STR_LEN, "%s.%02d", opt->lambda1pathfile, g->fold);
   return writevectorf(tmp, opt->lambda1path, opt->nlambda1);
}

/*
 * Run coordinate descent for each lambda1 penalty
 */
int run_train(Opt *opt, gmatrix *g)
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
	 printf("\n[%d] Fitting with lambda1=%.20f\n", i, opt->lambda1path[i]);

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

      snprintf(tmp, MAX_STR_LEN, "%s.%02d.%02d",
	    opt->beta_files[0], i, g->fold);
      printf("unscaling beta\n");
      unscale_beta(g->beta_orig, g->beta, g->mean, g->sd, g->p + 1);
      if(!write_beta_sparse(tmp, g->beta_orig, g->p + 1))
	 return FAILURE;

      if(opt->nzmax != 0 && opt->nzmax <= ret - 1)
      {
	 printf("maximum number of non-zero variables \
reached or exceeded: %d\n", opt->nzmax);
	 i++; /* increment to correct number of models fitted successfully */
	 break;
      }
   }

   snprintf(tmp, MAX_STR_LEN, "%s.%02d", opt->numnz_file, g->fold);
   /* number of non-zero variables for each successful fit and the all-zero
    * fit, excluding intercept */
   if(!writevectorl(tmp, g->numnz, i))
      return FAILURE;

   FREENULL(g->numnz);

   return SUCCESS;
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
	 g->nextcol(g, &sm, j, NA_ACTION_RANDOM);

	 /* We count up to n, which should be the same as sm.n,
	  * since we're not deleting missing obs */
	 for(i = 0 ; i < n ; i++)
	    lp[i] += sm.x[i] * beta[j];
      }
   }
   
   for(i = 0 ; i < n ; i++)
      yhat[i] = predict_func(lp[i]);

   if(!writevectorf(predict_file, yhat, n))
      return FAILURE;

   FREENULL(yhat);
   
   return SUCCESS;
}

int run_predict(gmatrix *g, predict predict_func,
      char **beta_files, int n_beta_files)
{
   int i;
   char tmp[MAX_STR_LEN];

   for(i = 0 ; i < n_beta_files ; i++)
   {
      gmatrix_zero_model(g);
      if(!load_beta_sparse(g->beta_orig, beta_files[i], g->p + 1))
      {
	 printf("skipping %s\n", beta_files[i]);
	 continue;
      }

      memcpy(g->beta, g->beta_orig, sizeof(double) * (g->p+1));

      snprintf(tmp, MAX_STR_LEN, "%s.pred", basename(beta_files[i]));
      if(!run_predict_beta(g, predict_func, tmp))
	 return FAILURE;
   }

   return SUCCESS;
}

int do_train(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, k, len;

   if(!gmatrix_init(g, opt->filename, opt->n, opt->p,
	    NULL, opt->yformat, opt->model, opt->modeltype, opt->encoded,
	    opt->binformat, opt->folds_ind_file, opt->mode,
	    opt->loss_pt_func, opt->subset_file,
	    opt->famfilename))
      return FAILURE;

   printf("%d CV folds\n", g->nfolds);
   printf("NZmax: %d\n", opt->nzmax);

   /* cross-validation: training stage */
   if(g->nfolds > 1)
   {
      for(k = 0 ; k < g->nfolds ; k++)
      {
	 printf("CV fold: %d\n", k);

	 MALLOCTEST(g->scalefile, sizeof(char) * MAX_STR_LEN);
	 len = strlen(opt->scalefile) + 1 + 3;
	 snprintf(g->scalefile, len, "%s.%02d", opt->scalefile, k);
	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 gmatrix_zero_model(g);
	 make_lambda1path(opt, g);
	 gmatrix_reset(g);

	 if(!(ret &= run_train(opt, g)))
	    break;

	 FREENULL(g->scalefile);
      }
   }
   else
   {
      MALLOCTEST(g->scalefile, sizeof(char) * MAX_STR_LEN);
      strncpy(g->scalefile, opt->scalefile, 
	    FMIN(strlen(opt->scalefile), MAX_STR_LEN));
      if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;
      gmatrix_zero_model(g);
      make_lambda1path(opt, g);
      gmatrix_reset(g);
      ret = run_train(opt, g);
      FREENULL(g->scalefile);
   }
   
   return ret;
}

/* We don't use the scaled inputs during prediction, as the model weights from
 * the training stage have been "unscaled" to original scale */
int do_predict(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, b, k, len;

   if(!gmatrix_init(g, opt->filename, opt->n, opt->p,
	    NULL, opt->yformat, opt->model, opt->modeltype, opt->encoded,
	    opt->binformat, opt->folds_ind_file, opt->mode,
	    opt->loss_pt_func, opt->subset_file,
	    opt->famfilename))
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
	 len = strlen(opt->scalefile) + 1 + 3;
	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 /* write y file */
	 snprintf(tmp, 5, "y.%02d", k);
	 printf("writing y file: %s\n", tmp);
	 if(!(ret &= writevectorf(tmp, g->y, g->ncurr)))
	    break;

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

   if(!opt_defaults(&opt, OPTIONS_CALLER_CD) 
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

