#include <stdlib.h>
#include "cd.h"
#include "loss.h"
#include "util.h"

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
   s = (log(opt->lambda1max) - log(opt->lambda1min)) / opt->nlambda1; 
   for(i = 1 ; i < opt->nlambda1 ; i++)
      opt->lambda1path[i] = exp(log(opt->lambda1max) - s * i);

   /* Write the coefs for model with intercept only */
   snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->beta_files[0], 0);
   if(!writevectorf(tmp, g->beta, opt->p + 1))
      return FAILURE;

   return writevectorf(opt->lambda1pathfile,
	 opt->lambda1path, opt->nlambda1);
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
	    opt->ntrain, opt->n - opt->ntrain);
   
   /* don't start from zero, getlambda1max already computed that */
   for(i = 1 ; i < opt->nlambda1 ; i++)
   {
      if(opt->verbose)
	 printf("\nFitting with lambda1=%.20f\n", opt->lambda1path[i]);

      /* return value is number of nonzero variables,
       * including the intercept */
      ret = cd_gmatrix(
	    g, opt->phi1_func, opt->phi2_func,
	    opt->step_func,
	    opt->maxepochs, opt->maxiters,
	    opt->lambda1path[i], opt->lambda2,
	    opt->threshold, opt->verbose, opt->trunc);

      gmatrix_reset(g);

      if(ret == CDFAILURE)
      {
	 printf("failed to converge after %d epochs\n", opt->maxepochs);
	 break;
      } 

      snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->beta_files[0], i);
      if(!writevectorf(tmp, g->beta, opt->p + 1))
	 return FAILURE;

      if(!opt->warmrestarts)
	 zero_model(g);

      if(opt->nzmax != 0 && opt->nzmax <= ret - 1)
      {
	 printf("maximum number of non-zero variables reached: %d\n", 
	       opt->nzmax);
	 break;
      }
   }

   return SUCCESS;
}

/* zero the lp and adjust the lp-functions */
void zero_model(gmatrix *g)
{
   int i, j, n = g->ntrain[g->fold], p1 = g->p + 1;

   for(j = p1 - 1 ; j >= 0 ; --j)
   {
      g->beta[j] = 0;
      g->active[j] = !g->ignore[j];
   }

   for(i = n - 1 ; i >= 0 ; --i)
      g->lp[i] = 0;
   if(g->model == MODEL_LOGISTIC)
      for(i = n - 1 ; i >= 0 ; --i)
	 g->lp_invlogit[i] = 0.5; /* 0.5 = 1 / (1 + exp(-0)) */
   else if(g->model == MODEL_SQRHINGE)
      for(i = n - 1 ; j >= 0 ; --i)
      {
	 g->ylp[i] = -1;    /* y * 0 - 1 = -1  */
	 g->ylp_neg[i] = 1; /* ylp < 0 => true */
      }
}

/*
 * For each beta file, predict outcome using the chosen model
 */
int run_predict_beta(gmatrix *g, predict predict_func,
      char* predict_file)
{
   int i, j, n = g->ntest[g->fold];
   sample sm;
   double *yhat;
   double *restrict lp = g->lp;
   double *restrict x;
   double *restrict beta = g->beta;

   if(!sample_init(&sm, n, g->inmemory))
      return FAILURE;

   MALLOCTEST(yhat, sizeof(double) * n)

   for(j = 0 ; j < g->p + 1 ; j++)
   {
      g->nextcol(g, &sm, FALSE);
      x = sm.x;
      for(i = 0 ; i < n ; i++)
	 lp[i] += x[i] * beta[j];
   }
   
   for(i = 0 ; i < n ; i++)
      yhat[i] = predict_func(lp[i]);

   printf("writing %s (%d) ... ", predict_file, n);
   if(!writevectorf(predict_file, yhat, n))
      return FAILURE;
   printf("done\n");

   free(yhat);
   sample_free(&sm);
   
   return SUCCESS;
}

int run_predict(gmatrix *g, predict predict_func, char **beta_files,
      int n_beta_files)
{
   int i;
   char tmp[MAX_STR_LEN];

   for(i = 0 ; i < n_beta_files ; i++)
   {
      zero_model(g);
      printf("reading %s\n", beta_files[i]);
      if(!load_beta(g->beta, beta_files[i], g->p + 1))
	 return FAILURE;

      /* scale beta using the scales for this data (beta should already
       * be on original scale, not scaled) */
      scale_beta(g->beta, g->mean, g->sd, g->p + 1);

      snprintf(tmp, MAX_STR_LEN, "%s.pred", beta_files[i]);
      if(!run_predict_beta(g, predict_func, tmp))
	 return FAILURE;
   }

   return SUCCESS;
}

/*int run_pcor(Opt *opt, gmatrix *g)
{
   int i, j, ret;
   double *betahat = NULL;
   double *lp = NULL;
   FILE *out;

   FOPENTEST(out, opt->filename, "wb");


   for(i = 0 ; i < opt->nlambda1 ; i++)
   {
      if(opt->verbose)
	 printf("\nFitting with lambda1=%.20f\n", opt->lambda1path[i]);

      for(j = 0 ; j < opt->p + 1 ; j++)
      {
	 CALLOCTEST(betahat, opt->p + 1, sizeof(double))
	 CALLOCTEST(lp, opt->n, sizeof(double))

	 ret = cd_gmatrix(
      	       g, opt->phi1_func, opt->phi2_func,
	       opt->loss_pt_func, opt->inv_func, opt->step_func,
      	       opt->maxepochs, betahat, opt->lambda1path[i], opt->lambda2,
      	       opt->threshold, opt->verbose, opt->trainf, opt->trunc);

      	 gmatrix_reset(g);

      	 if(ret == FAILURE)
      	 {
      	    printf("failed to converge after %d\n", opt->maxepochs);
      	    break;
      	 }

      	 *snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->betafile, i);
      	 if(!writevectorf(tmp, betahat, opt->p + 1))
      	    return FAILURE;*

      	 * if(!opt->warmrestarts)
      	    for(k = 0 ; k < opt->p + 1 ; k++)
      	       betahat[k] = 0; *

      	 if(opt->nzmax != 0 && opt->nzmax <= ret)
      	 {
      	    printf("maximum number of non-zero variables reached: %d\n", 
      	          opt->nzmax);
      	    break;
      	 }

	 * write output incrementally *
	 FWRITETEST(betahat, sizeof(double), opt->p + 1, out)

	 free(betahat);
	 free(lp);
      }
   }

   fflush(out);
   fclose(out);

   return SUCCESS;
}*/

int main(int argc, char* argv[])
{
   int ret = 0;
   Opt opt;
   gmatrix g;

   setbuf(stdout, NULL);

   if(!opt_defaults(&opt) || !opt_parse(argc, argv, &opt))
   {
      opt_free(&opt);
      return EXIT_FAILURE;
   }

   if(!gmatrix_init(&g, opt.filename, opt.n, opt.p,
	    opt.inmemory, opt.scalefile, opt.yformat, opt.model,
	    opt.encoded, opt.binformat, opt.folds_ind_file, opt.nfolds,
	    opt.mode))
   {
      gmatrix_free(&g);
      opt_free(&opt);
      fflush(stdout);
      return EXIT_FAILURE;
   }
  
   if(opt.mode == MODE_TRAIN && !opt.nofit)
   {
      make_lambda1path(&opt, &g);
      zero_model(&g);
      gmatrix_reset(&g);
      ret = run_train(&opt, &g);
   }
   else if(opt.mode == MODE_PREDICT)
      ret = run_predict(&g, opt.predict_func, opt.beta_files,
	    opt.n_beta_files);

   gmatrix_free(&g);
   opt_free(&opt);

   return ret == FAILURE ? EXIT_FAILURE : EXIT_SUCCESS;
}

