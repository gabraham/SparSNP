#include <stdlib.h>
#include <glob.h>

#include "cd.h"
#include "util.h"
#include "multivariable.h"

#define OPTIONS_CALLER cd

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
      g->nextcol(g, &sm, j);

      /* We count up to n, which should be the same as sm.n,
       * since we're not deleting missing obs */
      for(i = 0 ; i < n ; i++)
	 lp[i] += sm.x[i] * beta[j];
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

      /* Scale beta using the scales for this data (beta
       * should already be on original scale, not scaled to zero-mean and
       * unit-variance) */
      scale_beta(g->beta, g->beta_orig, g->mean, g->sd, g->p + 1);

      snprintf(tmp, MAX_STR_LEN, "%s.pred", beta_files[i]);
      if(!run_predict_beta(g, predict_func, tmp))
	 return FAILURE;
   }

   return SUCCESS;
}

int do_train(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, k, len;

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
	 len = strlen(opt->scalefile) + 1 + 3;
	 snprintf(tmp, len, "%s.%02d", opt->scalefile, k);
	 g->scalefile = tmp;
	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 gmatrix_zero_model(g);
	 make_lambda1path(opt, g, 0);
	 gmatrix_reset(g);

	 /*if(!(ret &= run_train(opt, g)))
	    break;*/
	 if(!(ret &= multivariable_lasso(opt, g, 0)))
	    break;
      }
   }
   else
   {
      g->scalefile = opt->scalefile;
      if(!gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;
      gmatrix_zero_model(g);
      make_lambda1path(opt, g, 0);
      gmatrix_reset(g);
      /*ret = run_train(opt, g);*/
      ret = multivariable_lasso(opt, g, 0);
   }
   
   return ret;
}

int do_predict(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, b, k = 0, len, gret = 0;
   glob_t gl;

   if(!gmatrix_init(g, opt->filename, opt->n, opt->p,
	    NULL, opt->yformat, opt->model, opt->encoded,
	    opt->binformat, opt->folds_ind_file, opt->mode,
	    opt->loss_pt_func, opt->subset_file))
      return FAILURE;

   if(g->nfolds > 1)
   {
      /* cross-validation: prediction stage */
      for(k = 0 ; k < g->nfolds ; k++)
      {
	 len = strlen(opt->scalefile) + 1 + 3;
	 snprintf(tmp, len, "%s.%02d", opt->scalefile, k);
	 g->scalefile = tmp;
	 printf("reading scale file: %s\n", tmp);
	 if(!gmatrix_read_scaling(g, g->scalefile))
	    return FAILURE;

	 if(!(ret &= gmatrix_set_fold(g, k)))
	    break;

	 /* write y file */
	 snprintf(tmp, 5, "y.%02d", k);
	 printf("writing y file: %s\n", tmp);
	 if(!(ret &= writevectorf(tmp, g->y, g->ncurr)))
	    break;

	 snprintf(tmp, MAX_STR_LEN, "%s.%02d.[0-9][0-9]",
	       opt->beta_files[0], k);
	 printf("Searching for glob %s\n", tmp);
      	 if((gret = glob(tmp, 0, NULL, &gl)))
      	 {
      	    fprintf(stderr, "glob error: %d\n", gret);
      	    return FAILURE;
      	 }

	 for(b = 0 ; b < gl.gl_pathc ; b++)
	    printf("found file: %s\n", gl.gl_pathv[b]);

	 gmatrix_reset(g);

	 if(!(ret &= run_predict(g, opt->predict_func,
		     gl.gl_pathv, gl.gl_pathc)))
	    break;

	 printf("\n");
      }
   }
   else
   {
      /* a bit of duplication */
      snprintf(tmp, MAX_STR_LEN, "%s.%02d.[0-9][0-9]",
	    opt->beta_files[0], k);
      printf("Searching for glob %s\n", tmp);
      if((gret = glob(tmp, 0, NULL, &gl)))
      {
	 fprintf(stderr, "glob error: %d\n", gret);
	 return FAILURE;
      }

      for(b = 0 ; b < gl.gl_pathc ; b++)
	 printf("found file: %s\n", gl.gl_pathv[b]);

      g->scalefile = opt->scalefile;
      if(!gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;
      gmatrix_zero_model(g);
      ret = run_predict(g, opt->predict_func,
	    gl.gl_pathv, gl.gl_pathc);
   }

   globfree(&gl);

   return ret;
}

int main(int argc, char* argv[])
{
   int ret = FAILURE;

   Opt opt;
   gmatrix g;
   char tmp[MAX_STR_LEN + 1];

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

