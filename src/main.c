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

void opt_free(Opt *opt)
{
   if(opt->lambda1path)
   {
      free(opt->lambda1path);
      opt->lambda1path = NULL;
   }

   if(opt->trainf)
   {
      free(opt->trainf);
      opt->trainf = NULL;
   }
}

void opt_defaults(Opt *opt)
{
   opt->nlambda1 = 100;
   opt->l1minratio = 1e-3;
   opt->maxepochs = 100;
   opt->lambda1 = -1;
   opt->lambda2 = 0;
   opt->threshold = 1e-4;
   opt->trunc = 1e-15;
   opt->nzmax = 0;
   opt->betafile = "beta.csv";
   opt->n = 0;
   opt->p = 0;
   opt->warmrestarts = FALSE;
   opt->nofit = FALSE;
   opt->filename = NULL;
   opt->lambda1path = NULL;
   opt->verbose = FALSE;
   opt->lambda1max = opt->lambda1min = 1;
   opt->cv = 1;
   opt->seed = time(NULL);
   opt->nzmax = 0;
   opt->trainf = NULL;
   opt->ntrain = opt->n;
   opt->subsetfile = "subset.csv";
   opt->lambda1pathfile = "lambda1path.csv";
   opt->step_func = NULL;
   opt->inmemory = FALSE;
}

int opt_parse(int argc, char* argv[], Opt* opt)
{
   int i;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-f"))
      {
	 i++;
	 opt->filename = argv[i];
      }
      else if(strcmp2(argv[i], "-model"))
      {
	 i++;
	 opt->model = argv[i];
	 if(strcmp2(opt->model, "logistic"))
	 {
	    opt->loss_pt_func = &logloss_pt;
	    opt->phi1_func = &logphi1;
	    opt->phi2_func = &logphi2;
	    opt->inv_func = &loginv;
	    opt->step_func = &step_regular;
	 }
	 else if(strcmp2(opt->model, "grouped"))
	 {
	    opt->loss_pt_func = &logloss_pt;
	    opt->phi1_func = &logphi1;
	    opt->phi2_func = &logphi2;
	    opt->inv_func = &loginv;
	    opt->step_func = &step_grouped;
	 }
	 else if(strcmp2(opt->model, "linear") ||
	       strcmp2(opt->model, "pcor"))
	 {
	    opt->loss_pt_func = &l2loss_pt;
	    opt->phi1_func = &l2phi1;
	    opt->phi2_func = &l2phi2;
	    opt->inv_func = &l2inv;
	    opt->step_func = &step_regular;
	 }
	 else
	 {
	    printf("model not available\n");
	    return FAILURE;
	 }
      }
      else if(strcmp2(argv[i], "-nofit"))
      {
	 opt->nofit = TRUE;
      }
      else if(strcmp2(argv[i], "-n"))
      {
	 i++;
	 opt->n = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-p"))
      {
	 i++;
	 opt->p = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-epochs"))
      {
	 i++;
	 opt->maxepochs = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1"))
      {
	 i++;
	 opt->lambda1 = atof(argv[i]);
	 opt->nlambda1 = 1;
      }
      else if(strcmp2(argv[i], "-l2"))
      {
	 i++;
	 opt->lambda2 = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1min"))
      {
	 i++;
	 opt->l1minratio = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-thresh"))
      {
	 i++;
	 opt->threshold = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-nl1"))
      {
	 i++;
	 opt->nlambda1 = atoi(argv[i]);
      }
      else if(strcmp2(argv[i], "-v"))
      {
	 opt->verbose = TRUE;
      }
      else if(strcmp2(argv[i], "-vv"))
      {
	 opt->verbose = 2;
      }
      else if(strcmp2(argv[i], "-vvv"))
      {
	 opt->verbose = 3;
      }
      else if(strcmp2(argv[i], "-beta"))
      {
	 i++;
	 opt->betafile = argv[i];
      }
      else if(strcmp2(argv[i], "-cv"))
      {
	 i++;
	 opt->cv = atoi(argv[i]);
      }
      else if(strcmp2(argv[i], "-seed"))
      {
	 i++;
	 opt->seed = atol(argv[i]);
      }
      else if(strcmp2(argv[i], "-nzmax"))
      {
	 i++;
	 opt->nzmax = atol(argv[i]);
      }
      else if(strcmp2(argv[i], "-warm"))
      {
	 opt->warmrestarts = TRUE;
      }
      else if(strcmp2(argv[i], "-inmemory"))
      {
	 opt->inmemory = TRUE;
      }
   }


   if(opt->filename == NULL || opt->model == NULL
	 || opt->n == 0 || opt->p == 0)
   {
      printf("usage: cd -model <model> -f <filename> -n <#samples> -p \
<#variables> | -beta <beta filename> -pred <pred filename> -epoch <maxepochs> \
-l1 <lambda1> -l2 <lambda2> -thresh <threshold> \
-pred <prediction file> -cv <cvfolds> -seed <seed> -v -vv\n");
      return FAILURE;
   }

   if(!opt->nzmax)
      opt->nzmax = (int)fmin(opt->n, opt->p); 

   srand(opt->seed);

   CALLOCTEST2(opt->lambda1path, opt->nlambda1, sizeof(double))
   
   cvsplit(opt);

   return SUCCESS; 
}

/*
 * Creates a vector of lambda1 penalties
 */
int make_lambda1path(Opt *opt, gmatrix *g)
{
   int i;
   double s;

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

   return writevectorf(opt->lambda1pathfile, opt->lambda1path, opt->nlambda1);
}

/*
 * Run coordinate descent for each lambda1 penalty
 */
int run(Opt *opt, gmatrix *g)
{
   int i, j, ret;
   double *betahat;
   double *lp;
   char tmp[MAX_STR_LEN];

   if(opt->verbose)
   {
      printf("%d training samples, %d test samples\n",
	    opt->ntrain, opt->n - opt->ntrain);
   }

   CALLOCTEST2(betahat, opt->p + 1, sizeof(double))
   CALLOCTEST2(lp, opt->n, sizeof(double))

   for(i = 0 ; i < opt->nlambda1 ; i++)
   {
      if(opt->verbose)
	 printf("\nFitting with lambda1=%.20f\n", opt->lambda1path[i]);

      /* return value is number of nonzero variables,
       * including the intercept */
      ret = cd_gmatrix(
	    g, opt->phi1_func, opt->phi2_func, opt->loss_pt_func,
	    opt->inv_func, opt->step_func,
	    opt->maxepochs, betahat, lp, opt->lambda1path[i], opt->lambda2,
	    opt->threshold, opt->verbose, opt->trainf, opt->trunc);

      gmatrix_reset(g);

      if(ret == FAILURE)
      {
	 printf("failed to converge after %d\n", opt->maxepochs);
	 break;
      } 

      snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->betafile, i);

      if(!writevectorf(tmp, betahat, opt->p + 1))
	 return FAILURE;

      if(!opt->warmrestarts)
      {
	 for(j = 0 ; j < opt->p + 1; j++)
	    betahat[j] = 0;
	 for(j = 0 ; j < opt->n; j++)
	    lp[j] = 0;
      }

      if(opt->nzmax != 0 && opt->nzmax <= ret - 1)
      {
	 printf("maximum number of non-zero variables reached: %d\n", 
	       opt->nzmax);
	 break;
      }
   }

   free(betahat);
   free(lp);

   return SUCCESS;
}

int run_pcor(Opt *opt, gmatrix *g)
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
      	       g, opt->phi1_func, opt->phi2_func, opt->loss_pt_func,
	       opt->inv_func, opt->step_func,
      	       opt->maxepochs, betahat, lp, opt->lambda1path[i], opt->lambda2,
      	       opt->threshold, opt->verbose, opt->trainf, opt->trunc);

      	 gmatrix_reset(g);

      	 if(ret == FAILURE)
      	 {
      	    printf("failed to converge after %d\n", opt->maxepochs);
      	    break;
      	 }

      	 /*snprintf(tmp, MAX_STR_LEN, "%s.%d", opt->betafile, i);
      	 if(!writevectorf(tmp, betahat, opt->p + 1))
      	    return FAILURE;*/

      	 /* if(!opt->warmrestarts)
      	    for(k = 0 ; k < opt->p + 1 ; k++)
      	       betahat[k] = 0; */

      	 if(opt->nzmax != 0 && opt->nzmax <= ret)
      	 {
      	    printf("maximum number of non-zero variables reached: %d\n", 
      	          opt->nzmax);
      	    break;
      	 }

	 /* write output incrementally */
	 FWRITETEST(betahat, sizeof(double), opt->p + 1, out)

	 free(betahat);
	 free(lp);
      }
   }


   fflush(out);
   fclose(out);


   return SUCCESS;
}

int main(int argc, char* argv[])
{
   Opt opt;
   gmatrix g;

   setbuf(stdout, NULL);

   opt_defaults(&opt);
   if(!opt_parse(argc, argv, &opt))
      return EXIT_FAILURE;

   if(!gmatrix_init(&g, opt.filename, opt.n, opt.p, opt.inmemory))
      return EXIT_FAILURE;
  
   make_lambda1path(&opt, &g);
   gmatrix_reset(&g);

   if(!opt.nofit)
      run(&opt, &g);

   gmatrix_free(&g);
   opt_free(&opt);
   
   return EXIT_SUCCESS;
}

