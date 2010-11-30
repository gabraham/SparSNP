#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"

#define OPTIONS_CALLER univariable

/* qnorm(abs(0.05 / 185805), lower.tail=FALSE) */
/*#define ZTHRESH 5.012169*/
#define ZTHRESH 1.0

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
   int i, 
	 iter = 1, maxiter = 10000;
   double s1 = 10, /* just a number larger than threshold */
	  s2 = 10;

   sample sm_intercept;
   sm_intercept.x = g->intercept;

   while(iter <= maxiter && (fabs(s1) >= 1e-8 || fabs(s2) >= 1e-8))
   {
      /* intercept */
      s1 = step_func(&sm_intercept, g, NULL, NULL);
      *beta_intercept -= s1;

      /* actual variable */
      s2 = step_func(sm, g, NULL, NULL);
      *beta -= s2;

      for(i = n - 1 ; i >= 0 ; --i)
      {
	 g->lp[i] += -s1 + sm->x[i] * -s2;
	 g->lp_invlogit[i] = 1 / (1 + exp(-g->lp[i]));
      }

      iter++;
   }
   /*printf("\n");*/

   printf("terminated in %d iterations\n", iter);

   return SUCCESS;
}

/*
 * Iteratively-Reweighted Least Squares for logistic regression
 */
int irls(gmatrix *g,
	      sample *sm,
	      double *beta_intercept,
	      double *beta,
	      step step_func, 
	      phi1 phi1_func,
	      phi2 phi2_func,
	      int n,
	      int p)
{
   int i, j, iter = 1, maxiter = 20;

   double grad[2] = {0, 0};
   double hessian[4] = {0, 0, 0, 0},
       invhessian[4] = {0, 0, 0, 0};
   double w, z, s1 = 10, s2 = 10;

   while(iter <= maxiter && (fabs(s1) >= 1e-8 || fabs(s2) >= 1e-8))
   {
      for(j = 3 ; j >= 0 ; --j)
      {
	 grad[j] = 0.0;
	 hessian[j] = 0.0;
      }

      for(i = n - 1; i >= 0 ; --i)
      {
	 g->lp[i] = *beta_intercept + *beta * sm->x[i];
	 g->lp_invlogit[i] = 1 / (1 + exp(-g->lp[i]));

	 /* Gradient */
	 grad[0] +=            (g->lp_invlogit[i] - g->y[i]);
	 grad[1] += sm->x[i] * (g->lp_invlogit[i] - g->y[i]);

	 /* Hessian */
	 w = g->lp_invlogit[i] * (1 - g->lp_invlogit[i]);
	 hessian[0] += w;
	 z = sm->x[i] * w;
	 hessian[1] += z;
	 hessian[2] = hessian[1];
	 hessian[3] += sm->x[i] * z;
      }
      
      invert2x2(invhessian, hessian);
      s1 = invhessian[0] * grad[0] + invhessian[1] * grad[1];
      s2 = invhessian[2] * grad[0] + invhessian[3] * grad[1];
      *beta_intercept -= s1;
      *beta -= s2;

      printf("%d intercept: %.5f beta: %.5f\n", iter, *beta_intercept, *beta);
      iter++;
   }

   if(iter >= maxiter)
      printf("IRLS didn't converge\n");

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
 * 1) Select the top k SNPs by univariable p values, using each SNP and the
 * intercept as covariables
 *
 * 2) Fit a multivariable model to the top k SNPs plus intercept
 */
int univar_gmatrix(Opt *opt, gmatrix *g, double *zscore)
{
   int i, j, k,
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
   CALLOCTEST(beta, p1, sizeof(double));

   /* We don't use this value ever, it's only here for consistency with the
    * multivariable methods such use the intercept */
   beta[0] = 0;

   /* get p-values per SNP */
   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d ", j);
      if(!g->active[j])
	 continue;
      
      beta_intercept = beta[j] = 0.0;
      g->nextcol(g, &sm, j);

      /*if(!cd_simple(g, &sm, &beta_intercept, beta + j, opt->step_func,
	    opt->phi1_func, opt->phi2_func, n, 2))
	 return FAILURE;*/

      if(!irls(g, &sm, &beta_intercept, beta + j, opt->step_func,
	    opt->phi1_func, opt->phi2_func, n, 2))
	 return FAILURE;

      /* don't let previous estimates affect current ones */
      /*gmatrix_zero_model(g);*/
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 g->lp[i] = 0;
	 g->lp_invlogit[i] = 0.5;
      }

      for(k = 0 ; k < 4 ; k++)
	 hessian[k] = 0.0;

      if(!make_hessian(hessian, sm.x, beta_intercept, beta[j], n))
	 return FAILURE;

      invert2x2(invhessian, hessian);

      /* z-score for the SNP only, ignore intercept */
      zscore[j] = beta[j] / sqrt(invhessian[3]);
      printf("z=%.3f\n", zscore[j]);

   }
   printf("\n");

   FREENULL(hessian);
   FREENULL(invhessian);
   FREENULL(beta);

   return SUCCESS;
}

int run_train(Opt *opt, gmatrix *g, double zthresh)
{
   int j, p1 = g->p + 1;
   int ret;
   int numselected = 0;
   double *zscore = NULL;

   CALLOCTEST(zscore, p1, sizeof(double));

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    g->ntrain[g->fold], g->ntest[g->fold]);

   /* select the SNP using univariable method */
   if(!(ret = univar_gmatrix(opt, g, zscore)))
      return FAILURE;

   printf("univariate selection done\n");
   for(j = 0 ; j < p1 ; j++)
   {
      printf("zscore[%d]: %.4f\n", j, zscore[j]);
      g->active[j] &= (fabs(zscore[j]) >= zthresh);
      if(g->active[j])
      {
	 printf("selected var %d\n", j);
	 numselected++;
      }
   }

   if(numselected == 0)
   {
      printf("no SNP exceeded threshold, aborting\n");
   }
   else
   {
      /* train un-penalised multivariable model on
       * the selected SNPs, with lambda=0 */
      ret = cd_gmatrix(
	    g, opt->phi1_func, opt->phi2_func,
	    opt->step_func,
	    opt->maxiters, opt->maxiters,
	    0, opt->lambda2,
	    opt->threshold, opt->verbose, opt->trunc);
   }

   FREENULL(zscore);

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

	 if(!(ret &= run_train(opt, g, ZTHRESH)))
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
      ret = run_train(opt, g, ZTHRESH);
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

