#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"

#define OPTIONS_CALLER univariable

/* 
 * Invert 2x2 matrix x and put it in y, row-major matrix ordering:
 *  0 1
 *  2 3
 */
void invert2x2(double *y, double *x)
{
   double s = x[0] * x[3] - x[1] * x[2];

   y[0] =  x[3] / s;
   y[1] = -x[1] / s;
   y[2] = -x[2] / s;
   y[3] =  x[0] / s;
}

/* Simple coordinate descent for minimising loss, without any penalisation
 *
 * x: an array of length p
 * beta: is an array of length p
 */
int cd_simple(gmatrix *g,
	      sample *sm,
	      double *beta_intercept,
	      double *beta,
	      step step_func, 
	      phi1 phi1_func,
	      phi2 phi2_func,
	      int n,
	      int p)
{
   int epoch = 1, maxepoch = 1000;
   double s1 = 10,
	  s2 = 10;

   sample sm_intercept;
   sm_intercept.x = g->intercept;

   while(epoch <= maxepoch && (fabs(s1) >= 1e-9 || fabs(s2) >= 1e-9))
   {
      /* intercept */
      s1 = step_func(&sm_intercept, g, NULL, NULL);
      *beta_intercept -= s1;
      updatelp(g, -s1, g->intercept);
      /*printf("\t[%d] s1: %.5f\tintercept: %.5f", epoch, s1,
       * *beta_intercept);*/

      /* actual variable */
      s2 = step_func(sm, g, NULL, NULL);
     /* printf("\ts2: %.5f\tbeta: %.5f\n",s2, *beta);*/
      *beta -= s2;
      updatelp(g, -s2, sm->x);

      epoch++;
   }
   /*printf("\n");*/

   return SUCCESS;
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
      double beta_intercept, double beta,
      int n, phi2 phi2_func)
{
   int k;

   double P, q;

   /* row-major ordering */
   for(k = 0 ; k < n ; k++)
   {
      P = 1 / (1 + exp(-beta_intercept - beta * x[k]));
      q = P * (1 - P);
      hessian[0] += q;
      hessian[1] += x[k] * q;
      hessian[2] = hessian[1];
      hessian[3] += x[k] * x[k] * q;
   }

   return SUCCESS;
}

/* Two stages:
 *
 * 1) Select the top k SNPs by univariable p values, using each SNP and the
 * intercept as covariables
 *
 * 2) Fit a multivariable model to the top k SNPs plus intercept
 */
int univar_gmatrix(Opt *opt, gmatrix *g, double *zscore)
{
   int j, k,
       n = g->ncurr,
       p1 = g->p + 1;
   double *hessian = NULL,
	  *invhessian = NULL,
	  *beta = NULL;
   double beta_intercept = 0;

   sample sm;

   if(!sample_init(&sm, n))
      return FAILURE;

   MALLOCTEST(hessian, sizeof(double) * 4);
   MALLOCTEST(invhessian, sizeof(double) * 4);
   MALLOCTEST(beta, sizeof(double) * p1);
   MALLOCTEST(zscore, sizeof(double) * p1);

   /* We don't use this value ever, it's only here for consistency with the
    * multivariable methods such use the intercept */
   beta[0] = 0;

   /* get p-values per SNP */
   for(j = 1 ; j < p1 ; j++)
   {
      beta_intercept = beta[j] = 0.0;
      g->nextcol(g, &sm, j);

      if(!cd_simple(g, &sm, &beta_intercept, beta + j, opt->step_func,
	    opt->phi1_func, opt->phi2_func, n, 2))
	 return FAILURE;

      /* don't let previous estimates affect current ones */
      gmatrix_zero_model(g);

      printf("beta[%d]: %.10f\t", j, beta[j]);

      for(k = 0 ; k < 4 ; k++)
	 hessian[k] = 0.0;

      if(!make_hessian(hessian, sm.x, beta_intercept, beta[j], n, opt->phi2_func))
	 return FAILURE;

      invert2x2(invhessian, hessian);

      /* z-score for the SNP only, ignore intercept */
      zscore[j] = beta[j] / sqrt(invhessian[3]);
      printf("zscore: %.6f\n", zscore[j]);
   }

   FREENULL(hessian);
   FREENULL(invhessian);
   FREENULL(beta);
   FREENULL(zscore);

   return SUCCESS;
}

int run_train(Opt *opt, gmatrix *g)
{
   int ret;
   /*char tmp[MAX_STR_LEN];*/
   double *zscore = NULL;

   MALLOCTEST(zscore, sizeof(double) * g->p);

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    g->ntrain[g->fold], g->ntest[g->fold]);

   ret = univar_gmatrix(opt, g, zscore);

   FREENULL(zscore);

   return SUCCESS;
}

int do_train(gmatrix *g, Opt *opt, char tmp[])
{
   int ret = SUCCESS, k;/*, len;*/

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
      ret = run_train(opt, g);
   }
   
   return ret;
}

int do_predict(gmatrix *g, Opt *opt, char *tmp)
{
   return FAILURE;
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

