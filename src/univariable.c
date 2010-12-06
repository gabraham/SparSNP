#include "common.h"
#include "cd.h"
#include "util.h"
#include "coder.h"
#include "svd.h"

#define OPTIONS_CALLER univariable

/* qnorm(abs(0.05 / 185805), lower.tail=FALSE) */
/*#define ZTHRESH 5.012169*/
#define ZTHRESH 3.5

#define IRLS_THRESH 1e-9
#define IRLS_THRESH_MAX 10

#define IRLS_ERR_NO_CONVERGENCE 2 /* didn't converge within predefined iterations */
#define IRLS_ERR_DIVERGENCE 3 /* converged but to a very large value */

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

/*
 * Z = X^T Y
 *
 * X: m by n
 * Y: m by p
 * Z: n by p
 *
 * row major ordering
 *
 * This is a bit wasteful for Z = X^T X due to symmetry
 */
void crossprod(double *x, double *y, double *z, int m, int n, int p)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
      for(j = 0 ; j < p ; j++)
	 for(k = 0 ; k < m ; k++)
	    z[i * p + j] += x[k * n + i] * y[k * p + j];
}

/*
 * Weighted cross product
 * Z = X^T W Y
 *
 * W: a diagonal m by m matrix
 * X: m by n
 * Y: m by p
 * Z: n by p
 *
 * Note that w is NOT a matrix, it's an array of length m
 */
void wcrossprod(double *x, double *y, double *w, double *z, int m, int n, int p)
{
   int i, j, k;

   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 k = 0;
	 z[i * p + j] = x[k * n + i] * y[k * p + j] * w[k];
	 for(k = 1 ; k < m ; k++)
	    z[i * p + j] += x[k * n + i] * y[k * p + j] * w[k];
      }
   }
}

/*
 * Square-Matrix-vector product
 *
 * z = X y
 *
 * z: m by 1
 * X: m by m
 * y: m by 1
 *
 */
void sqmvprod(double *x, double *y, double *z, int m)
{
   int i, k;
      
   for(i = 0 ; i < m ; i++)
   {
      k = 0;
      z[i] = x[k * m + i] * y[k];
      for(k = 1 ; k < m ; k++)
	 z[i] += x[k * m + i] * y[k];
   }
}

/*
 * Iteratively-Reweighted Least Squares for logistic regression
 */
int irls(double *x, double *y, double *beta, double *invhessian, int n, int p)
{
   int i, j, 
       iter = 1, maxiter = 50,
       converged = FALSE, diverged = FALSE,
       ret = SUCCESS;

   double *grad = NULL,
	  *hessian = NULL,
	  *lp = NULL,
	  *lp_invlogit = NULL,
	  *w = NULL,
	  *s = NULL;

   MALLOCTEST(lp, sizeof(double) * n);
   MALLOCTEST(lp_invlogit, sizeof(double) * n);
   MALLOCTEST(grad, sizeof(double) * p);
   CALLOCTEST(hessian, p * p, sizeof(double));
   CALLOCTEST(w, n, sizeof(double));
   CALLOCTEST(s, p, sizeof(double));

   while(iter <= maxiter) 
   {
      for(i = n - 1; i >= 0 ; --i)
      {
	 lp[i] = x[i * p] * beta[0];

	 for(j = 1 ; j < p ; j++)
	 {
	    lp[i] += x[i * p + j] * beta[j];
	    lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	 }
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
      pseudoinverse(hessian, &p, &p, invhessian);

      /* Newton step */
      sqmvprod(invhessian, grad, s, p);
      converged = TRUE;
      diverged = FALSE;
      for(j = 0 ; j < p ; j++)
      {
	 beta[j] -= s[j];
	 converged &= (fabs(s[j]) <= IRLS_THRESH);
	 diverged |= (fabs(beta[j]) >= IRLS_THRESH_MAX);
      }

      if(converged || diverged)
	 break;

      iter++;
   }

   if(iter >= maxiter)
   {
      printf("IRLS didn't converge\n");
      ret = IRLS_ERR_NO_CONVERGENCE;
   }
   else if(diverged)
   {
      printf("IRLS diverged\n");
      ret = IRLS_ERR_DIVERGENCE;
   }

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
 * 1) Select the top k SNPs by univariable p values, using each SNP and the
 * intercept as covariables
 *
 * 2) Fit a multivariable model to the top k SNPs plus intercept
 */
int univar_gmatrix(Opt *opt, gmatrix *g, double *zscore)
{
   int i, j, ret,
       n = g->ncurr,
       p1 = g->p + 1;
   double beta[2] = {0, 0};
   double *invhessian = NULL,
	  *x = NULL;
   sample sm;

   if(!sample_init(&sm, n))
      return FAILURE;

   MALLOCTEST(invhessian, sizeof(double) * 4);
   CALLOCTEST(x, 2 * n, sizeof(double));

   /* get p-values per SNP, skip intercept */
   for(j = 1 ; j < p1 ; j++)
   {
      printf("%d ", j);
      if(!g->active[j])
      {
	 printf("skipped variable %d\n", j);
	 continue;
      }
      
      g->nextcol(g, &sm, j);
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 x[2 * i] = 1.0;
	 x[2 * i + 1] = sm.x[i];
      }

      beta[0] = beta[1] = 0.0;
      ret = irls(x, g->y, beta, invhessian, n, 2);

      if(ret == FAILURE)
	 return FAILURE;
      else if(ret == SUCCESS)
      {
	 /* z-score for the SNP only, ignore intercept */
	 zscore[j] = beta[1] / sqrt(invhessian[3]);
      }
      else
	 zscore[j] = 0.0;
      printf("zscore[%d]: %.6f\n", j, zscore[j]);
   }
   printf("\n");

   FREENULL(invhessian);
   FREENULL(x);

   return SUCCESS;
}

int run_train(Opt *opt, gmatrix *g, double zthresh)
{
   int k, j, n = g->ncurr, p1 = g->p + 1;
   int ret;
   int numselected = 0;
   double *zscore = NULL;
   double *x = NULL,
	  *invhessian = NULL,
	  *beta = NULL;
   char tmp[MAX_STR_LEN];

   CALLOCTEST(zscore, p1, sizeof(double));

   if(opt->verbose)
      printf("%d training samples, %d test samples\n",
	    g->ntrain[g->fold], g->ntest[g->fold]);

   /* select the SNP using univariable method */
   if(!(ret = univar_gmatrix(opt, g, zscore)))
      return FAILURE;

   printf("univariate selection done\n");
   numselected = 1;
   g->active[0] = TRUE;
   for(j = 1 ; j < p1 ; j++)
   {
      g->active[j] &= (fabs(zscore[j]) >= zthresh);
      if(g->active[j])
      {
	 printf("selected var %d\n", j);
	 numselected++;
      }
   }

   printf("total %d SNPs exceeded z-score=%.3f\n", numselected, zthresh);
   if(numselected > 0)
   {
      MALLOCTEST(x, sizeof(double) * n * numselected);
      if(!gmatrix_read_matrix(g, x, g->active, numselected))
	 return FAILURE;
      CALLOCTEST(invhessian, numselected * numselected, sizeof(double));
      CALLOCTEST(beta, numselected, sizeof(double));

      /* train un-penalised multivariable model on
       * the selected SNPs, with lambda=0 */
      ret = irls(x, g->y, beta, invhessian, n, numselected);
      printf("IRLS returned %d\n", ret);	    

      /* copy estimated beta to the array for all SNPs */
      k = 0;
      for(j = 0 ; j < p1 ; j++)
      {
	 if(g->active[j])
	 {
	    printf("beta[%d]: %.5f\n", j, beta[k]);
	    g->beta[j] = beta[k];
	    k++;
	 }
      }
   }

   snprintf(tmp, MAX_STR_LEN, "%s.00.%02d", opt->beta_files[0], g->fold);
   if(!writevectorf(tmp, g->beta, g->p + 1))
      return FAILURE;

   FREENULL(zscore);
   FREENULL(x);
   FREENULL(invhessian);
   FREENULL(beta);

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

   if(!sample_init(&sm, n))
      return FAILURE;

   CALLOCTEST(yhat, n, sizeof(double));

   for(j = 0 ; j < p1 ; j++)
   {
      g->nextcol(g, &sm, j);
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
	 snprintf(tmp, 5, "y.%02d", k);
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

