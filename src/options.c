#include "cd.h"
#include "loss.h"

void opt_free(Opt *opt)
{
   int i;

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

   if(opt->beta_files)
   {
      for(i = 0 ; i < opt->n_beta_files ; i++)
	 free(opt->beta_files[i]);
      free(opt->beta_files);
      opt->beta_files = NULL;
   }
}

int opt_defaults(Opt *opt)
{
   const char *beta_default = "beta.csv";

   opt->mode = MODE_TRAIN;
   opt->model = 0;
   opt->nlambda1 = 100;
   opt->l1minratio = 1e-2;
   opt->maxepochs = 100;
   opt->lambda1 = -1;
   opt->lambda2 = 0;
   opt->threshold = 1e-4;
   opt->trunc = 1e-15;
   opt->nzmax = 0;
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
   opt->scalefile = NULL;
   opt->yformat = YFORMAT01;
   opt->predict_func = NULL;
   opt->predict_file = "predicted.csv";
   opt->encoded = FALSE;

   MALLOCTEST(opt->beta_files, sizeof(char*))
   MALLOCTEST(opt->beta_files[0], sizeof(char) * (strlen(beta_default) + 1))
   strcpy(opt->beta_files[0], beta_default);
   opt->n_beta_files = 1;

   return SUCCESS;
}

int opt_parse(int argc, char* argv[], Opt* opt)
{
   int i, j, k;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-f"))
      {
	 i++;
	 opt->filename = argv[i];
      }
      else if(strcmp2(argv[i], "-train"))
	 opt->mode = MODE_TRAIN;
      else if(strcmp2(argv[i], "-predict"))
	 opt->mode = MODE_PREDICT;
      else if(strcmp2(argv[i], "-model"))
      {
	 i++;
	 if(strcmp2(argv[i], MODEL_NAME_LOGISTIC))
	 {
	    opt->inv_func = &loginv;
	    opt->step_func = &step_regular_logistic;
	    opt->model = MODEL_LOGISTIC;
	    opt->predict_func = &logphi1;
	 }
	 else if(strcmp2(argv[i], MODEL_NAME_SQRHINGE))
	 {
	    opt->inv_func = &sqrhingeinv;
	    opt->step_func = &step_regular_sqrhinge;
	    opt->yformat = YFORMAT11;
	    opt->model = MODEL_SQRHINGE;
	    opt->predict_func = &linearphi1;
	 }
	 else if(strcmp2(argv[i], MODEL_NAME_LINEAR))
	 {
	    opt->inv_func = &linearinv;
	    opt->step_func = &step_regular_linear;
	    opt->model = MODEL_LINEAR;
	    opt->predict_func = &linearphi1;
	 }
	 else if(strcmp2(argv[i], MODEL_NAME_PCOR))
	 {
	    opt->inv_func = &linearinv;
	    opt->step_func = &step_regular_linear;
	    opt->model = MODEL_LINEAR;
	    opt->predict_func = linearphi1;
	 }
	 else
	 {
	    printf("model not available\n");
	    return FAILURE;
	 }
      }
      else if(strcmp2(argv[i], "-nofit"))
	 opt->nofit = TRUE;
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
	 opt->verbose = TRUE;
      else if(strcmp2(argv[i], "-vv"))
	 opt->verbose = 2;
      else if(strcmp2(argv[i], "-vvv"))
	 opt->verbose = 3;
      else if(strcmp2(argv[i], "-scale"))
      {
	 i++;
	 opt->scalefile = argv[i];
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
	 printf("-inmemory not working currently\n");
	 return FAILURE;
      }
      else if(strcmp2(argv[i], "-betafiles"))
      {
	 /* first free the default filename */
	 free(opt->beta_files[0]);
	 free(opt->beta_files);
	 opt->beta_files = NULL;

	 j = ++i;
	 k = 0;

	 /* look ahead to find more tokens */
	 while(j < argc && argv[j][0] != '-')
	 {
	    k++;
	    REALLOCTEST(opt->beta_files, opt->beta_files,
		  sizeof(char*) * k)
	    MALLOCTEST(opt->beta_files[k - 1],
		  sizeof(char) * (strlen(argv[j]) + 1))
	    strcpy(opt->beta_files[k - 1], argv[j]);
	    j++;
	 }
	 opt->n_beta_files = k;

	 if(j < argc && argv[j][0] == '-')
	    i = j - 1;

	 printf(">>>%s\n", opt->beta_files[0]);
      }
      else if(strcmp2(argv[i], "-encoded"))
      {
	 opt->encoded = TRUE;
      }
   }


   if(opt->filename == NULL || opt->model == 0
	 || opt->n == 0 || opt->p == 0)
   {
      printf("usage: cd [-train|-predict] -model <model> \
-f <filename> -n <#samples> -p \
<#variables> | -betafiles <beta filename/s> -pred <pred filename> \
-epoch <maxepochs> -l1 <lambda1> [-encoded] \
-l2 <lambda2> -thresh <threshold> \
-pred <prediction file> -cv <cvfolds> -seed <seed> -v -vv\n");
      return FAILURE;
   }
   else if(opt->n_beta_files > 1 && opt->mode == MODE_TRAIN)
   {
      printf("warning: multiple beta filenames provided in training mode, \
only using the first one\n");
   }

   if(opt->verbose)
      printf("Mode: %s\n",
	    (opt->mode == MODE_TRAIN) ? "train" : "predict");

   if(!opt->nzmax)
      opt->nzmax = (int)fmin(opt->n, opt->p); 

   srand(opt->seed);

   CALLOCTEST2(opt->lambda1path, opt->nlambda1, sizeof(double))
   
   if(opt->cv > 1)
      cvsplit(opt);

   return SUCCESS; 
}
