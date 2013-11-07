/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include "sparsnp.h"

void opt_free(Opt *opt)
{
   int i;

   if(opt->lambda1path)
   {
      free(opt->lambda1path);
      opt->lambda1path = NULL;
   }

   if(opt->beta_files)
   {
      for(i = 0 ; i < opt->n_beta_files ; i++)
	 free(opt->beta_files[i]);
      free(opt->beta_files);
      opt->beta_files = NULL;
   }
   if(opt->beta_files_fold)
   {
      for(i = 0 ; i < opt->n_beta_files ; i++)
	 free(opt->beta_files_fold[i]);
      free(opt->beta_files_fold);
      opt->beta_files_fold = NULL;
   }

   FREENULL(opt->zthresh);
}

int opt_defaults(Opt *opt, short caller)
{
   const char *beta_default = "beta.csv";

   opt->caller = caller;
   opt->mode = MODE_TRAIN;
   opt->model = 0;
   opt->nlambda1 = 100;
   opt->l1minratio = 1e-2;
   opt->l1max = NOT_DEFINED;
   opt->maxepochs = 1e5;
   opt->maxiters = 100;
   opt->lambda1 = NOT_DEFINED;
   opt->lambda2 = 0;
   opt->gamma = 0;
   opt->threshold = 1e-5;
   opt->trunc = 1e-10;
   opt->nzmax = 5000;
   opt->n = 0;
   opt->p = 0;
   opt->warmrestarts = TRUE;
   opt->nofit = FALSE;
   opt->filename = NULL;
   opt->lambda1path = NULL;
   opt->verbose = FALSE;
   opt->lambda1max = opt->lambda1min = 1;
   opt->nfolds = 1;
   opt->folds_ind_file = NULL;
   opt->seed = time(NULL);
   opt->ntrain = opt->n;
   opt->subset_file = NULL;
   opt->lambda1pathfile = "lambda1path.csv";
   opt->step_func = NULL;
   opt->scalefile = "scale.bin";
   opt->yformat = YFORMAT01;
   opt->predict_func = NULL;
   opt->predict_file = "predicted.csv";
   opt->encoded = TRUE;
   opt->beta_files_fold = NULL;
   opt->numnz_file = "nonzero.csv";
   opt->outdir = "";
   opt->scaley = FALSE;

   MALLOCTEST(opt->beta_files, sizeof(char*));
   MALLOCTEST(opt->beta_files[0], sizeof(char) * (strlen(beta_default) + 1));
   strcpy(opt->beta_files[0], beta_default);
   opt->n_beta_files = 1;

   opt->nzthresh = 24;
   MALLOCTEST(opt->zthresh, sizeof(double) * opt->nzthresh);

   opt->zthresh[0] = 30.20559;  /* 1e-200 */
   opt->zthresh[1] = 27.82780;  /* 1e-170 */
   opt->zthresh[2] = 26.12296;  /* 1e-150 */
   opt->zthresh[3] = 23.33408;  /* 1e-120 */
   opt->zthresh[4] = 22.32745;  /* 1e-110 */
   opt->zthresh[5] = 20.16469;  /* 1e-90 */
   opt->zthresh[6] = 18.99164;  /* 1e-80 */
   opt->zthresh[7] = 17.74164;  /* 1e-70 */
   opt->zthresh[8] = 16.39728;  /* 1e-60 */
   opt->zthresh[9] = 14.93334;  /* 1e-50 */
   opt->zthresh[10] = 13.31092; /* 1e-40 */
   opt->zthresh[11] = 11.46403; /* 1e-30 */
   opt->zthresh[12] = 9.262340; /* 1e-20 */
   opt->zthresh[13] = 6.361341; /* 1e-10 */
   opt->zthresh[14] = 5.326724; /* 5e-8  */
   opt->zthresh[15] = 5.199338; /* 1e-7  */
   opt->zthresh[16] = 4.264891; /* 1e-5  */
   opt->zthresh[17] = 3.719016; /* 1e-4  */
   opt->zthresh[18] = 3.570974; /* 1e-4  */
   opt->zthresh[19] = 3.417300; /* 1e-4  */
   opt->zthresh[20] = 3.257323; /* 1e-4  */
   opt->zthresh[21] = 3.090232; /* 1e-3  */
   opt->zthresh[22] = 2.326348; /* 1e-2  */
   opt->zthresh[23] = 1.281552; /* 1e-1  */

   /*opt->lambda2_univar = 1e-3;*/
   opt->lambda2_univar = 0;
   opt->lambda2_multivar = 0;

   opt->do_multivar = TRUE;
   opt->existing_univar = FALSE;
   opt->do_thinning = TRUE;

   opt->multivar = OPTIONS_MULTIVAR_NEWTON;
   opt->famfilename = NULL;
   
   opt->unscale_beta = FALSE;
   opt->cortype = 2;
   opt->corthresh = 0;
   opt->phenoformat = PHENO_FORMAT_FAM;

   /* in bytes */
   opt->maxmem = CACHE_MEM_DEFAULT;

   return SUCCESS;
}

int opt_parse(int argc, char* argv[], Opt* opt)
{
   int i, j, k;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-bin") || strcmp2(argv[i], "-bed"))
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
	    //if(opt->caller == OPTIONS_CALLER_CD)
	    //{
	    //   printf("model `logistic' currently not supported\n");
	    //   return FAILURE;
	    //}
	    //else
	    {
	       opt->inv_func = &loginv;
	       opt->step_func = &step_regular_logistic;
	       opt->model = MODEL_LOGISTIC;
	       opt->predict_func = &logphi1;
	       opt->modeltype = MODELTYPE_CLASSIFICATION;
	    }
	 }
	 else if(strcmp2(argv[i], MODEL_NAME_SQRHINGE))
	 {
	    opt->inv_func = &sqrhingeinv;
	    opt->step_func = &step_regular_sqrhinge;
	    opt->yformat = YFORMAT11;
	    opt->model = MODEL_SQRHINGE;
	    opt->predict_func = &linearphi1;
	    opt->modeltype = MODELTYPE_CLASSIFICATION;
	 }
	 else if(strcmp2(argv[i], MODEL_NAME_LINEAR))
	 {
	    opt->inv_func = &linearinv;
	    opt->step_func = &step_regular_linear;
	    opt->model = MODEL_LINEAR;
	    opt->predict_func = &linearphi1;
	    opt->modeltype = MODELTYPE_REGRESSION;
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
      else if(strcmp2(argv[i], "-maxepochs"))
      {
	 i++;
	 opt->maxepochs = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-maxiters"))
      {
	 i++;
	 opt->maxiters = (int)atof(argv[i]);
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
      else if(strcmp2(argv[i], "-gamma"))
      {
	 i++;
	 opt->gamma = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1min"))
      {
	 i++;
	 opt->l1minratio = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1max"))
      {
	 i++;
	 opt->l1max = atof(argv[i]);
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
      else if(strcmp2(argv[i], "-cold"))
      {
	 opt->warmrestarts = FALSE;
      }
      else if(strcmp2(argv[i], "-outdir"))
      {
	 i++;
	 opt->outdir = argv[i];
      }
      else if(strcmp2(argv[i], "-betafiles"))
      {
	 /* first free the default filename */
	 free(opt->beta_files[0]);
	 free(opt->beta_files);
	 opt->beta_files = NULL;

	 j = ++i;
	 k = 0;

	 printf("outdir: %s\n", opt->outdir);
	 /* look ahead to find more tokens */
	 while(j < argc && argv[j][0] != '-')
	 {
	    k++;
	    REALLOCTEST(opt->beta_files, opt->beta_files,
		  sizeof(char*) * k)
	    MALLOCTEST(opt->beta_files[k - 1],
		  sizeof(char) * (strlen(argv[j]) + strlen(opt->outdir) + 2))
	    /*strcpy(opt->beta_files[k - 1], argv[j]);*/
	    if(strlen(opt->outdir) > 0)
	       sprintf(opt->beta_files[k - 1], "%s/%s", opt->outdir, argv[j]);
	    else
	       sprintf(opt->beta_files[k - 1], "%s", argv[j]);
	    j++;
	 }
	 opt->n_beta_files = k;

	 if(j < argc && argv[j][0] == '-')
	    i = j - 1;
      }
      else if(strcmp2(argv[i], "-notencoded"))
	 opt->encoded = FALSE;
      else if(strcmp2(argv[i], "-foldind"))
      {
	 i++;
	 opt->folds_ind_file = argv[i];
      }
      /*else if(strcmp2(argv[i], "-zthresh"))
      {
	 i++;
	 opt->zthresh = atof(argv[i]);
      }*/
      else if(strcmp2(argv[i], "-nomultivar"))
      {
	 opt->do_multivar = FALSE;
      }
      else if(strcmp2(argv[i], "-existingunivar"))
      {
	 opt->existing_univar = TRUE;
      }
      else if(strcmp2(argv[i], "-nothin"))
      {
	 opt->do_thinning = FALSE;
      }
      else if(strcmp2(argv[i], "-filter"))
      {
	 opt->do_lasso_filter = TRUE;
      }
      else if(strcmp2(argv[i], "-multivarlasso"))
      {
	 opt->multivar = OPTIONS_MULTIVAR_LASSO;
      }
      else if(strcmp2(argv[i], "-subset"))
      {
	 i++;
	 opt->subset_file = argv[i];
      }
      else if(strcmp2(argv[i], "-fam"))
      {
	 i++;
	 opt->famfilename = argv[i];
      }
      else if(strcmp2(argv[i], "-pheno"))
      {
	 i++;
	 opt->phenoformat = PHENO_FORMAT_PHENO;
	 opt->famfilename = argv[i];
      }
      else if(strcmp2(argv[i], "-scaley"))
      {
	 opt->scaley = TRUE;
      }
      else if(strcmp2(argv[i], "-unscale"))
      {
	 opt->unscale_beta = TRUE;
      }
      else if(strcmp2(argv[i], "-cortype"))
      {
	 i++;
	 opt->cortype = atol(argv[i]);
      }
      else if(strcmp2(argv[i], "-corthresh"))
      {
	 i++;
	 opt->corthresh = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-maxmem"))
      {
	 i++;
	 opt->maxmem = atoi(argv[i]);
      }
   }

   if(opt->caller == OPTIONS_CALLER_CD) /* coordinate descent */
   {
      if(opt->filename == NULL || opt->model == 0
            || opt->n == 0 || opt->p == 0 
	    || (opt->mode == MODE_TRAIN && !opt->scalefile))
      {
         printf("usage: sparsnp [-train|-predict] -model <model> \
-bed <filename> -fam <filename> -n <#samples> -p <#variables> -scale <scalefile> \
[-betafiles <beta filename/s>] \
[-maxepochs <maxepochs>] [-maxiters <maxiters>] [-l1 <lambda1>] \
-l2 <lambda2>] [-gamma <gamma>] [-thresh <threshold>] [-foldind <foldsfile>] \
[-pred <prediction file>] [-filter] [-unscale] [-pheno <filename>] \
[-cortype [0/1/2]] [-corthresh <threshold>] [-v] [-vv] [-maxmem]\n");
         return FAILURE;
      }
      else if(opt->n_beta_files > 1 && opt->mode == MODE_TRAIN)
      {
         printf("warning: multiple beta filenames provided in training mode, \
onl   y using the first one\n");
      }
   }
   else /* univariable selection */
   {
      if(opt->filename == NULL || opt->model == 0
            || opt->n == 0 || opt->p == 0)
      {
         printf("usage: univariable [-train|-predict] -model <model> \
-bed <filename> -n <#samples> -p <#variables>  \
[-betafiles <beta filename/s>] \
[-foldind <foldsfile>] [-existingunivar] \
[-pred <prediction file>] [-nomultivar] [-v] [-vv]\n");
         return FAILURE;
      }
      else if(opt->n_beta_files > 1 && opt->mode == MODE_TRAIN)
      {
         printf("warning: multiple beta filenames provided in training mode, \
               only using the first one\n");
      }
   }

   if(opt->mode == MODE_TRAIN && !opt->famfilename)
   {
      printf("Error: you must provide a FAM filename (-fam) or PHENO \
 filename (-pheno) when using plink BED input\n");
      return FAILURE;
   }
	 

   if(!opt->encoded)
   {
      printf("non-encoded (switch -notencoded) data not currently \
supported\n");
      return FAILURE;
   }

   if(opt->verbose)
   {
      printf("Mode: %s\n",
	    (opt->mode == MODE_TRAIN) ? "train" : "predict");
      printf("Model: %s\n",
	    (opt->model == MODEL_LINEAR) ? "linear" :
	       (opt->model == MODEL_SQRHINGE) ? "sqrhinge" : "logistic");
   }

   if(opt->model != MODEL_LINEAR && opt->scaley)
   {
      printf("Refusing to scale phenotypes Y (-scaley) with model: ");
      printf("%s\n",
	    (opt->model == MODEL_LINEAR) ? "linear" :
	    (opt->model == MODEL_SQRHINGE) ? "sqrhinge" : "logistic");
      return FAILURE;
   }

   srand(opt->seed);

   CALLOCTEST2(opt->lambda1path, opt->nlambda1, sizeof(double))
   
   return SUCCESS; 
}

