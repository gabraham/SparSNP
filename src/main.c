#include <stdlib.h>
#include "cd.h"
#include "loss.h"
#include "util.h"

/*
 * Split data into training and test set
 */
int cvsplit(Opt *opt)
{
   int i;

   MALLOCTEST2(opt->trainf, sizeof(int) * opt->n)

   for(i = 0 ; i < opt->n ; i++)
   {
      if(opt->cv > 1)
	 opt->trainf[i] = drand48() >= (1.0 / opt->cv);
      else
	 opt->trainf[i] = TRUE;
      opt->ntrain += opt->trainf[i];
   }
   return writevectorl(opt->subsetfile, opt->trainf, opt->n);
}

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
   int i, j, ret;
   char tmp[MAX_STR_LEN];

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    opt->ntrain, opt->n - opt->ntrain);
   
   /* The maximum lambda1 has already been run by make_lambda1path, and the
    * the linear predictor and its related vectors have already been updated */
   for(i = 1 ; i < opt->nlambda1 ; i++)
   {
      if(opt->verbose)
	 printf("\nFitting with lambda1=%.20f\n", opt->lambda1path[i]);

      /* return value is number of nonzero variables,
       * including the intercept */
      ret = cd_gmatrix(
	    g, opt->phi1_func, opt->phi2_func,
	    opt->loss_pt_func, opt->inv_func, opt->step_func,
	    opt->maxepochs, opt->maxiters,
	    opt->lambda1path[i], opt->lambda2,
	    opt->threshold, opt->verbose, opt->trainf, opt->trunc);

      gmatrix_reset(g);

      if(ret == FAILURE)
      {
	 printf("failed to converge after %d epochs\n", opt->maxepochs);
	 break;
      } 

      snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->beta_files[0], i);
      if(!writevectorf(tmp, g->beta, opt->p + 1))
	 return FAILURE;

      if(!opt->warmrestarts)
      {
	 for(j = 0 ; j < opt->p + 1; j++)
	    g->beta[j] = 0;
	 for(j = 0 ; j < opt->n; j++)
	    g->lp[j] = 0;
	 if(opt->model == MODEL_LOGISTIC)
	    for(j = 0 ; j < opt->n; j++)
	       g->lp_invlogit[j] = 0.5; /* 0.5 = 1 / (1 + exp(-0)) */
	 else if(opt->model == MODEL_SQRHINGE)
	    for(j = 0 ; j < opt->n; j++)
	    {
	       g->ylp[j] = -1;    /* y * 0 - 1 = -1  */
	       g->ylp_neg[j] = 1; /* ylp < 0 => true */
	    }
      }

      if(opt->nzmax != 0 && opt->nzmax <= ret - 1)
      {
	 printf("maximum number of non-zero variables reached: %d\n", 
	       opt->nzmax);
	 break;
      }
   }

   return SUCCESS;
}

/*
 * For each beta file, predict outcome using the chosen model
 */
int run_predict_beta(gmatrix *g, predict predict_func,
      char* predict_file)
{
   unsigned int i, j;
   sample sm;
   double *yhat;
   double *restrict lp = g->lp;
   double *restrict x;
   double *restrict beta = g->beta;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   MALLOCTEST(yhat, sizeof(double) * g->n)

   for(j = 0 ; j < g->p + 1 ; j++)
   {
      g->nextcol(g, &sm);
      x = sm.x;
      for(i = 0 ; i < g->n ; i++)
	 lp[i] += x[i] * beta[j];
   }
   
   for(i = 0 ; i < g->n ; i++)
      yhat[i] = predict_func(lp[i]);

   if(!writevectorf(predict_file, yhat, g->n))
      return FAILURE;

   free(yhat);
   sample_free(&sm);
   
   return SUCCESS;
}

/* Assumes ascii, one value per line */
int load_beta(double *beta, char *filename, int p)
{
   int i = 0;
   FILE *in = NULL;
   FOPENTEST(in, filename, "rt");

   while(!feof(in))
   {
      if(fscanf(in, "%lf", beta + i) == EOF)
	 break;
      i++;
   } 
   fclose(in);
   return SUCCESS;
}

int run_predict(gmatrix *g, predict predict_func, char **beta_files,
      int n_beta_files)
{
   int i;
   char tmp[MAX_STR_LEN];

   for(i = 0 ; i < n_beta_files ; i++)
   {
      printf("Reading %s\n", beta_files[i]);
      if(!load_beta(g->beta, beta_files[i], g->p + 1))
	 return FAILURE;

      snprintf(tmp, MAX_STR_LEN, "%s_pred.%d", beta_files[i], i);
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
	    opt.inmemory, opt.scalefile,
	    opt.yformat, opt.model, opt.encoded, opt.binformat))
   {
      gmatrix_free(&g);
      opt_free(&opt);
      fflush(stdout);
      return EXIT_FAILURE;
   }
  
   if(opt.mode == MODE_TRAIN && !opt.nofit)
   {
      make_lambda1path(&opt, &g);
      gmatrix_reset(&g);
      ret = run_train(&opt, &g);
   }
   else if(opt.mode == MODE_PREDICT)
      ret = run_predict(&g, opt.predict_func, opt.beta_files,
	    opt.n_beta_files);

   gmatrix_free(&g);
   opt_free(&opt);

   if(ret == FAILURE)
      return EXIT_FAILURE;
   return EXIT_SUCCESS;
}

