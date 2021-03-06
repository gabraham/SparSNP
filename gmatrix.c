/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "gmatrix.h"
#include "ind.h"
#include "util.h"
#include "gennetwork.h"

int sample_init(sample *s)
{
   s->n = 0;
   s->x = NULL;
   s->y = NULL;
   s->intercept = FALSE;
   return SUCCESS;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int phenoformat,
      int model, int modeltype,
      short encoded, char *folds_ind_file,
      short mode, char *subset_file, 
      char *famfilename, int scaley, int unscale_beta,
      int cortype, int corthresh, int verbose, long maxmem,
      double gamma)
{
   int i, j, k, p1;
   long seed = 0;
   long p1K;

   g->verbose = verbose;
   g->file = NULL;
   g->map = NULL;
   g->K = 0;
   g->model = model;
   g->modeltype = modeltype;
   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;
   p1 = p + 1;
   g->yidx = 0;
   g->y = NULL;
   g->y_orig = NULL;
   g->sex = NULL;
   g->lossK = NULL;
   g->l1lossK = NULL;
   g->loss = 0;
   g->floss = 0;
   g->err = NULL;
   g->lp = NULL;
   g->err = NULL;
   g->ylp = NULL;
   g->ylp_neg = NULL;
   g->ylp_neg_y = NULL;
   g->ylp_neg_y_ylp = NULL;
   g->lp_invlogit = NULL;
   g->lookup = NULL;
   g->intercept = NULL;
   g->mean = NULL;
   g->sd = NULL;
   g->tmp = NULL;
   g->xtmp = NULL;
   g->x = NULL;
   g->xthinned = NULL;
   g->ignore = NULL;
   g->yformat = yformat;
   g->beta = NULL;
   g->active = NULL;
   g->encoded = encoded;
   g->nencb = (int)ceil((double)n / PACK_DENSITY);
   g->encbuf = NULL;
   g->folds_ind_file = folds_ind_file;
   g->nfolds = 1;
   g->folds = NULL;
   g->fold = 0; 
   g->ntrain = NULL;
   g->ntest = NULL;
   g->ntrainrecip = NULL;
   g->ntestrecip = NULL;
   g->ncurr = 0;
   g->mode = mode;
   g->nseek = sizeof(dtype) * (encoded ? g->nencb : g->n);
   g->beta_orig = NULL;
   g->numnz = NULL;
   g->xcaches = NULL;
   g->ncases = 0;
   g->folds_ind = NULL;
   g->unscale_beta = unscale_beta;
   g->famfilename = famfilename;
   g->C = NULL;
   g->cortype = cortype;
   g->corthresh = corthresh;
   g->phenoformat = phenoformat;
   g->diagCC = NULL;
   g->edges = NULL;
   g->pairs = NULL;
   g->scaley = scaley;
   g->flossK = NULL;
   g->tol = 1e-6;
   g->gamma = gamma;
   g->maxmem = maxmem;

   seed = getpid();
   srand48(seed);

   printf("seed: %ld\n", seed);

   printf("Using PLINK binary format\n");
   printf("FAM file: %s\n", g->famfilename);
   g->offset = 3 - g->nseek; /* TODO: magic number */
   if(g->famfilename)
   {
      printf("pheno format: %s\n",
	 g->phenoformat == PHENO_FORMAT_FAM ? "FAM" : "PHENO");
      if(g->phenoformat == PHENO_FORMAT_FAM)
      {
	 if(!gmatrix_fam_read_y(g))
	    return FAILURE;
	 g->K = 1;
      }
      else if(!gmatrix_pheno_read_y(g))
	 return FAILURE;

      gmatrix_plink_check_pheno(g);
   }

   if(g->mode == MODE_TRAIN && g->modeltype == MODELTYPE_CLASSIFICATION)
   {
      count_cases(g);
      printf("Found %d cases, %d controls\n", g->ncases, g->n - g->ncases);
   }

   MALLOCTEST(g->map, sizeof(mapping));
   mapping_init(g->map);

   srand48(time(NULL));

   p1K = (long)p1 * g->K;
   CALLOCTEST(g->beta_orig, p1K, sizeof(double));
   CALLOCTEST(g->grad_array, p1K, sizeof(pair));

   for(k = 0 ; k < g->K ; k++)
   {
      g->grad_array[k * p1].index = 0;
      g->grad_array[k * p1].value = INFINITY;

      for(j = 1 ; j < p1 ; j++)
      {
	 g->grad_array[j + k * p1].index = j;
	 g->grad_array[j + k * p1].value = 0;
      }
   }

   CALLOCTEST(g->lossK, g->K, sizeof(double));
   CALLOCTEST(g->l1lossK, g->K, sizeof(double));
   
   if(filename)
      FOPENTEST(g->file, filename, "rb");

   /* gmatrix_setup_folds changes g->nfolds */
   if(!gmatrix_setup_folds(g))
      return FAILURE;

   /* We set up an explicit intercept vector so that we don't need to generate
    * it on the fly later  */
   MALLOCTEST(g->intercept, sizeof(double) * g->n);
   for(i = n - 1 ; i >= 0 ; --i)
      g->intercept[i] = 1.0;

   MALLOCTEST(g->tmp, sizeof(dtype) * g->nencb * PACK_DENSITY);

   if(encoded)
      MALLOCTEST(g->encbuf, sizeof(unsigned char) * g->nencb);

   //MALLOCTEST(g->xcaches, g->nfolds * sizeof(cache));
   //for(i = 0 ; i < g->nfolds ; i++)
   //   cache_init(g->xcaches + i, g->n, g->p + 1, maxmem);

   MALLOCTEST(g->xcaches, sizeof(cache));
   g->xcaches->mapping = NULL;
   g->xcaches->revmapping = NULL;
   g->xcaches->counter = NULL;
   g->xcaches->x = NULL;
   g->xcaches->tmp = NULL;

   cache_init(g->xcaches, g->n, g->p + 1, maxmem);

   g->nextcol = gmatrix_disk_nextcol;
   if(!gmatrix_reset(g))
      return FAILURE;

   CALLOCTEST(g->beta, (long)p1 * g->K, sizeof(double));
   CALLOCTEST(g->active, (long)p1 * g->K, sizeof(int));
   CALLOCTEST(g->ignore, (long)p1 * g->K, sizeof(int));

   /* must be called before setting a scalefile */
   if(!gmatrix_init_proportions(g))
      return FAILURE;

   /* don't scale in prediction mode */
   if(g->mode == MODE_TRAIN)
   {
      g->scalefile = scalefile;
      if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
         return FAILURE;
   }
   else
      g->scalefile = NULL;

   for(j = (long)p1 * g->K - 1 ; j >= 0 ; --j)
      g->active[j] = !g->ignore[j];

   gmatrix_set_ncurr(g);
   
   gmatrix_make_y(g);

   if(!gmatrix_init_lp(g))
      return FAILURE;

   return SUCCESS;
}

int gmatrix_setup_folds(gmatrix *g)
{
   int i;
   if(g->folds_ind_file 
	 && !(g->nfolds = ind_getfolds(g->folds_ind_file)))
      return FAILURE;

   MALLOCTEST(g->folds, sizeof(int) * g->n * g->nfolds);
   MALLOCTEST(g->folds_ind, sizeof(int) * g->n * g->nfolds);

   if(g->nfolds > 1 && !ind_read(g->folds_ind_file,
	    g->folds, g->n, g->nfolds))
      return FAILURE;

   MALLOCTEST(g->ntrain, sizeof(int) * g->nfolds);
   MALLOCTEST(g->ntest, sizeof(int) * g->nfolds);
   MALLOCTEST(g->ntrainrecip, sizeof(double) * g->nfolds);
   MALLOCTEST(g->ntestrecip, sizeof(double) * g->nfolds);

   if(g->folds_ind_file && g->nfolds > 1)
   {
      count_fold_samples(g->ntrain, g->ntest,
	    g->ntrainrecip, g->ntestrecip,
	    g->folds, g->nfolds, g->n);
   }
   else /* train and test sets are the same */
   {
      g->ntrain[0] = g->ntest[0] = g->n;
      g->ntrainrecip[0] = g->ntestrecip[0] = 1.0 / g->n;
   }

   for(i = g->n * g->nfolds - 1; i >= 0 ; --i)
      g->folds_ind[i] = (g->mode == MODE_PREDICT) ^ g->folds[i];

   return SUCCESS;
}

void gmatrix_free(gmatrix *g)
{
   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   FREENULL(g->mean);
   FREENULL(g->sd);
   FREENULL(g->y);
   FREENULL(g->y_orig);
   FREENULL(g->sex);
   FREENULL(g->xtmp);
   FREENULL(g->ignore);
   FREENULL(g->tmp);
   FREENULL(g->x);
   FREENULL(g->intercept);
   FREENULL(g->lookup);
   FREENULL(g->lp);
   FREENULL(g->ylp);
   FREENULL(g->err);
   FREENULL(g->ylp_neg);
   FREENULL(g->ylp_neg_y);
   FREENULL(g->ylp_neg_y_ylp);
   FREENULL(g->lp_invlogit);
   FREENULL(g->beta);
   FREENULL(g->beta_orig);
   FREENULL(g->encbuf);
   FREENULL(g->folds);
   FREENULL(g->folds_ind);
   FREENULL(g->ntrain);
   FREENULL(g->ntest);
   FREENULL(g->ntrainrecip);
   FREENULL(g->ntestrecip);
   FREENULL(g->active);
   FREENULL(g->numnz);
   
   if(g->xcaches)
      cache_free(g->xcaches);

   FREENULL(g->xcaches);
   mapping_free(g->map);
   FREENULL(g->map);
   FREENULL(g->C);
   FREENULL(g->proportions);
   FREENULL(g->diagCC);
   FREENULL(g->edges);
   FREENULL(g->pairs);
   FREENULL(g->lossK);
   FREENULL(g->l1lossK);
   FREENULL(g->flossK);
   FREENULL(g->grad_array);
}

/* y is a copy of y_orig as needed for crossval - in training y are the
 * training samples, in prediction y are the testing samples
 */
int gmatrix_split_y(gmatrix *g)
{
   int i, n = g->n, n1 = g->n - 1, k, K1 = g->K - 1;
   int m = g->ncurr - 1;

   if(g->ncurr <= 0) {
      fprintf(stderr, "gmatrix_split_y: ncurr can't be zero");
      return FAILURE;
   }

   FREENULL(g->y);

   MALLOCTEST(g->y, sizeof(double) * g->ncurr * g->K);

   /* represent y as 0/1 or -1/1, allowing for different
    * samples sizes if cross-validation is used */
   if(g->nfolds > 1) /* no cv */
   {
      for(k = K1 ; k >= 0 ; --k)
      {
         m = g->ncurr - 1;
         for(i = n1 ; i >= 0 ; --i)
      	 {
      	    if(g->folds_ind[g->fold * n + i])
            {
      	       g->y[g->ncurr * k + m] = g->y_orig[n * k + i];
               m-- ;
            }
      	 }
      }
   }
   else 
   {
      if(g->n != g->ncurr)
      {
	 fprintf(stderr,
	    "something is wrong in fold splitting, \
g->n != g->ncurr (%d != %d)\n", g->n, g->ncurr);
	 return FAILURE;
      }
      /* simply copy the vector, don't care about folds */
      for(i = n * g->K - 1 ; i >= 0 ; --i)
	 g->y[i] = g->y_orig[i];
   }

   return SUCCESS;
}

/*
 * Read all the data into a preallocated row-major n by m matrix, for n
 * samples and m variables.
 * The variables to be read are determined by the p+1 array ``ind''.
 */
int gmatrix_read_matrix(gmatrix *g, int *ind, int m)
{
   int i, j, k = 0,
       p1 = g->p + 1;
   sample sm;

   if(!sample_init(&sm))
      return FAILURE;

   for(j = 0 ; j < p1 ; j++)
   {
      if(ind[j])
      {
	 /* must not delete obs otherwise matrix structure will break */
	 if(!gmatrix_disk_nextcol(g, &sm, j, NA_ACTION_RANDOM))
	    return FAILURE;
   
         for(i = 0 ; i < sm.n ; i++)
	 {
	    g->x[i * m + k] = sm.x[i];
	    if(g->x[i * m + k] < 0)
	       asm("int3");
	 }
	 k++;
      }
   }

   return SUCCESS;
}

/* Returns one vector of samples from the in-memory matrix. Only used
 * for in-memory lasso.
 *
 * Due to the matrix structure, all columns are of same length,
 * regardless of missing values (missing values are imputed).
 */
int gmatrix_mem_nextcol(gmatrix *g, sample *sm, int j, int na_action)
{
   int i, n = g->ncurr, p1 = g->p + 1;
   sm->y = g->y; 
   for(i = 0 ; i < sm->n ; i++)
      sm->x[i] = g->x[i * n + p1];

   return SUCCESS;
}

long gmatrix_lrand()
{
   return lrand48();
}

double gmatrix_drand()
{
   return drand48();
}

int rand_geno()
{
   return gmatrix_lrand() % 3;
}

/* generate random phenotypes in proportion to their distribution in each SNP
 */
int rand_geno_proportional(gmatrix *g, int j)
{
   double r = gmatrix_drand();
   int p1 = g->p + 1;
   double cumsum[3] = {
      g->proportions[j],                          /* genotype 0 */
      g->proportions[j] + g->proportions[p1 + j], /* genotype 1 */
      1,                                          /* genotype 2 */
   };

   if(r < cumsum[0])
      return 0;
   else if(r < cumsum[1])
      return 1;
   return 2;
}

/* TODO: duplicates code in gmatrix_disk_nextcol */
int gmatrix_disk_nextcol_raw(gmatrix *g, sample *s, int j)
{
   int i, n1 = g->n - 1, n = g->n;
   off_t seek;
 
   /* Get data from disk and unpack, skip y */
   seek = (off_t)j * g->nseek + g->offset;
   FSEEKOTEST(g->file, seek, SEEK_SET);
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   
   decode_plink_mapping(g->map, g->tmp, g->encbuf, g->nencb);

   for(i = n1 ; i >= 0 ; --i)
      s->x[i] = g->tmp[i];
   s->n = n;

   return SUCCESS;
}

/*
 * Reads one column of data from disk (all samples for one variable).
 * Takes into account whether we're doing cross-validation or not
 * and if so which fold we're currently in, and how we choose to 
 * handle missing values.
 *
 * This function does all the heavy lifting in terms of determining
 * which samples belong in the current cross-validation fold etc.
 *
 * na_action: one of NA_ACTION_ZERO, NA_ACTION_DELETE, NA_ACTION_RANDOM/
 */
int gmatrix_disk_nextcol(gmatrix *g, sample *s, int j, int na_action)
{
   int i, l1, n1 = g->n - 1, n = g->n;
   int k = g->ncurr - 1;
   int f = g->fold * n;
   int ngood = 0, ret;
   dtype d;
   off_t seek;
   double *xtmp = NULL;

   if(j == 0)
   {
      /* don't need to worry about cross-validation etc
       * because the intercept is all 1s */
      s->x = g->intercept;

      s->y = g->y;
      s->n = g->ncurr;

      return SUCCESS;
   }

   //ret = cache_get(g->xcaches + g->fold, j, &xtmp);
   ret = cache_get(g->xcaches, j, &xtmp);

   if(ret == SUCCESS)
   {
      s->x = xtmp;
      return SUCCESS;
   }

   /* Get data from disk and unpack, skip y */
   seek = (off_t)j * g->nseek + g->offset;
   FSEEKOTEST(g->file, seek, SEEK_SET);
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   
   decode_plink_mapping(g->map, g->tmp, g->encbuf, g->nencb);

   /* Get the scaled versions of the genotypes */
   if(g->nfolds < 2)
   { /* without crossval */
      if(g->scalefile)
      { /* scale.c doesn't use scalefile so check */
         /* Get the scaled value instead of the original value.
	  * Note that na_action=NA_ACTION_DELETE is NOT supported
	  * for scaled inputs since scaled inputs come from
	  * the lookup table */
	 if(na_action == NA_ACTION_DELETE)
	 {
	    fprintf(stderr, "NA_ACTION_DELETE not supported for scaled \
inputs in gmatrix_disk_nextcol\n");
	    return FAILURE;
	 }

         l1 = j * NUM_X_LEVELS;
         for(i = n1 ; i >= 0 ; --i)
	    xtmp[i] = g->lookup[l1 + g->tmp[i]];
      }
      else
      { /* no scaling */
	 if(na_action == NA_ACTION_ZERO)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       d = g->tmp[i];
	       xtmp[i] = (d == X_LEVEL_NA ? 0 : (double)d);
	    }
	    s->n = n;
	 }
	 else if(na_action == NA_ACTION_RANDOM)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       d = g->tmp[i];
	       xtmp[i] = (d == X_LEVEL_NA ? (double)rand_geno() : (double)d);
	    }
	    s->n = n;
	 }
	 else if(na_action == NA_ACTION_DELETE)
	 {
	    for(i = 0 ; i < n ; ++i)
	    {
	       d = g->tmp[i];
	       if(d != X_LEVEL_NA)
	       {
		  xtmp[ngood] = (double)d;
	       }
	    }
	    s->n = ngood;
	 }
	 else if(na_action == NA_ACTION_PROPORTIONAL)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       d = g->tmp[i];
	       xtmp[i] = (
		  d == X_LEVEL_NA 
		     ? (double)rand_geno_proportional(g, j) : (double)d);
	    }
	    s->n = n;
	 }
	 else return FAILURE;
      }
   }
   else 
   {
      /* with crossval */
      /* Get the scaled value instead of the original value */
      if(g->scalefile)
      {
	 if(na_action == NA_ACTION_DELETE)
	 {
	    fprintf(stderr, "NA_ACTION_DELETE not supported for scaled \
inputs in gmatrix_disk_nextcol\n");
	    return FAILURE;
	 }

	 l1 = j * NUM_X_LEVELS;
	 for(i = n1 ; i >= 0 ; --i)
	 {
	    if(g->folds_ind[f + i])
	    {
	       xtmp[k--] = g->lookup[l1 + g->tmp[i]];
	       ngood++;
	    }
	 }
	 s->n = ngood;
      }
      else 
      { /* no scaling */
	 if(na_action == NA_ACTION_ZERO)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       if(g->folds_ind[f + i])
	       {
	          d = g->tmp[i];
	          xtmp[k--] = (d == X_LEVEL_NA ? 0 : (double)d);
		  ngood++;
	       }
	    }
	    s->n = ngood;
	 } 
	 if(na_action == NA_ACTION_RANDOM)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       if(g->folds_ind[f + i])
	       {
	          d = g->tmp[i];
	          xtmp[k--] = 
		     (d == X_LEVEL_NA ? (double)rand_geno() : (double)d);
		  ngood++;
	       }
	    }
	    s->n = ngood;
	 } 
	 else if(na_action == NA_ACTION_DELETE)
	 {
	    for(i = 0 ; i < n ; ++i)
	    {
	       if(g->folds_ind[f + i])
	       {
	          d = g->tmp[i];
		  if(d != X_LEVEL_NA)
		  {
		     xtmp[ngood] = (double)d;
		  }
	       }
	    }
	    s->n = ngood;
	 }
	 else if(na_action == NA_ACTION_PROPORTIONAL)
	 {
	    for(i = n1 ; i >= 0 ; --i)
	    {
	       if(g->folds_ind[f + i])
	       {
	          d = g->tmp[i];
	          xtmp[k--] = (d == X_LEVEL_NA ?
			(double)rand_geno_proportional(g, j) : (double)d);
		  ngood++;
	       }
	    }
	    s->n = ngood;
	 }
      }
   }

   s->x = xtmp;
   /*s->y = g->ytmp;*/

   return SUCCESS;
}

/* 
 * Reads a plink FAM file and gets the phenotype from the 6th column, using
 * any sort of whitespace separator.
 */
int gmatrix_fam_read_y(gmatrix *g)
{
   int i = 0;
   FILE* famfile = NULL;
   char *famid = NULL,
        *individ = NULL,
	*patid = NULL,
	*matid = NULL,
	*sex = NULL,
	*pheno = NULL;
   int ret;
   int const NUMFIELDS = 6;

   printf("gmatrix_fam_read_y\n");

   /* spaces match tabs as well */
   char f[19] = "%s %s %s %s %s %s";

   CALLOCTEST(g->y_orig, g->n, sizeof(double));

   MALLOCTEST(famid, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(individ, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(patid, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(matid, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(sex, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(pheno, sizeof(char) * FAM_MAX_CHARS);

   FOPENTEST(famfile, g->famfilename, "r");

   while(i < g->n &&
	 EOF != (ret = fscanf(famfile, f,
	    famid, individ, patid, matid, sex, pheno)))
   {
      if(strcmp2(pheno, MISSING_PHENO) || strcmp2(pheno, "NA") 
	 || ret < NUMFIELDS)
      {
	 fprintf(stderr,
	    "Missing phenotypes (%s, NA) currently not supported, remove \
these samples with PLINK; aborting\n",
	    MISSING_PHENO);
	 return FAILURE;
      }
      g->y_orig[i] = atof(pheno);
      //strncpy(g->sex[i], sex, FAM_MAX_CHARS - 1);
      //g->sex[i] = *sex;
      i++;
   }

   fclose(famfile);
   FREENULL(famid);
   FREENULL(individ);
   FREENULL(patid);
   FREENULL(matid);
   FREENULL(sex);
   FREENULL(pheno);

   return SUCCESS;
}

/* Standardise each column of Y to zero mean and unit variance.
 */
int gmatrix_scale_y(gmatrix *g)
{
#ifdef DEBUG
   char buf[100];
#endif
   int i, k, n = g->ncurr, K = g->K;
   double *mean = NULL,
	  *sd = NULL;
   double delta;

   printf("scaling y\n");

#ifdef DEBUG
   sprintf(buf, "y_%s_before_%02d.txt",
      (g->mode == MODE_PREDICT) ? "predict" : "train",
      g->fold);
   writematrixf(g->y, n, K, buf);
#endif

   CALLOCTEST(mean, g->K, sizeof(double));
   CALLOCTEST(sd, g->K, sizeof(double));

   /* first get mean and sd */
   for(k = 0 ; k < K ; k++)
   {
      for(i = 0 ; i < n ; i++)
      {
	 delta = g->y[k * n + i] - mean[k];
	 mean[k] += delta / (i + 1);
	 sd[k] += delta * (g->y[k * n + i] - mean[k]);
      }
      sd[k] = sqrt(sd[k] / (n - 1));
   }

   /* now standardise Y */
   for(k = 0 ; k < K ; k++)
      for(i = 0 ; i < n ; i++)
	 g->y[k * n + i] = (g->y[k * n + i] - mean[k]) / sd[k];

#ifdef DEBUG
   sprintf(buf, "y_%s_after_%02d.txt",
      (g->mode == MODE_PREDICT) ? "predict" : "train",
      g->fold);
   writematrixf(g->y, n, K, buf);
#endif

   FREENULL(mean);
   FREENULL(sd);

   return SUCCESS;
}

/* 
 * Reads a plink pheno file and gets the phenotype from column 3 and onwards,
 * using any sort of whitespace separator.
 *
 * y stores the phenotypes in column-major ordering
 */
int gmatrix_pheno_read_y(gmatrix *g)
{
   int i = 0, k = 0, n = g->n, m;
   FILE* famfile = NULL;
   char *famid = NULL,
        *individ = NULL,
	*pheno = NULL;
   int ret;
   int const NUMFIELDS = 2;
   double *tmp;
   char line[MAX_LINE_CHARS];

   /* spaces match tabs as well */
   char f[7] = "%s %s ";

   printf("gmatrix_pheno_read_y\n");

   /* we don't know how what K is before reading the pheno file */
   CALLOCTEST(tmp, g->n * MAX_NUM_PHENO, sizeof(double));

   MALLOCTEST(famid, sizeof(char) * FAM_MAX_CHARS);
   MALLOCTEST(individ, sizeof(char) * FAM_MAX_CHARS);
   CALLOCTEST(pheno, FAM_MAX_CHARS, sizeof(char));

   FOPENTEST(famfile, g->famfilename, "r");
   printf("g->famfilename %s\n", g->famfilename);

   /* read the family data for each line */
   while(i < g->n && EOF != (ret = fscanf(famfile, f, famid, individ)))
   {
      if(ret < NUMFIELDS)
      {
	 fprintf(stderr,
	    "Error in reading PHENO file '%s', missing phenotypes?",
	    g->famfilename);
	 return FAILURE;
      }

      if(fgets(line, MAX_LINE_CHARS, famfile) == NULL)
      {
      	 if(!feof(famfile))
	 {
	    fprintf(stderr, "error in reading PHENO file '%s'\n",
	       g->famfilename);
	    return FAILURE;
	 }
	 break;
      }
      else
      {
	 k = 0;
	 m = 0;
	 /* now read phenotypes, one at a time */
       	 while(k < MAX_NUM_PHENO && sscanf(&line[m], "%s ", pheno) != EOF)
       	 {
       	    tmp[i + n * k] = atof(pheno);
	    m += strlen(pheno) + 1;
       	    k++;
       	 }
      }

      i++;
   }

   printf("read %d rows and %d phenotype/s from FAM file '%s'\n", i, k,
	 g->famfilename);

   /* now truncate K to what we have observed and copy over */
   g->K = k;
   CALLOCTEST(g->y_orig, g->n * g->K, sizeof(double));
   for(i = n * k - 1 ; i >= 0 ; --i)
      g->y_orig[i] = tmp[i];

   fclose(famfile);
   FREENULL(famid);
   FREENULL(individ);
   FREENULL(pheno);
   FREENULL(tmp);

   return SUCCESS;
}

/*
 * plink phenotypes can be 0/1 and 1/2, so check for 1/2 and convert as
 * necessary
 */
int gmatrix_plink_check_pheno(gmatrix *g)
{
   int i = 0, n = g->n, K = g->K;
   int twofound = FALSE;

   /* don't change inputs for regression */
   if(g->modeltype == MODELTYPE_REGRESSION)
      return SUCCESS;

   for(i = n * K - 1 ; i >= 0 ; --i)
      if((twofound = (g->y_orig[i] == 2)))
	 break;

   if(twofound)
   {
      for(i = n * K - 1 ; i >= 0 ; --i)
      {
	 if(g->yformat == YFORMAT01) 
	    g->y_orig[i] -= 1;
	 else 
	    g->y_orig[i] = 2.0 * g->y_orig[i] - 3.0;
      }
   }
   /* input is 0/1 but we want it to be -1/1 */ 
   else if(g->yformat == YFORMAT11)
   {
      for(i = n * K - 1 ; i >= 0 ; --i)
	 g->y_orig[i] = 2.0 * g->y_orig[i] - 1.0;
   }

   return SUCCESS;
}

int gmatrix_reset(gmatrix *g)
{
   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   FOPENTEST(g->file, g->filename, "rb");

   g->i = g->j = 0;

   return SUCCESS;
}

/* 
 * Populate lookup tables for SNP levels 0, 1, 2 for each SNP.
 *
 * By this point, if we chose to impute data, then all data has been imputed
 * already.
 */
int gmatrix_read_scaling(gmatrix *g, char *file_scale)
{
   int j, k, l1, l2, p1 = g->p + 1;
   FILE *in;

   if(!file_scale)
      return FAILURE;

   printf("reading scale file %s\n", file_scale);

   if(!g->lookup)
      CALLOCTEST(g->lookup, NUM_X_LEVELS * p1, sizeof(double));
   if(!g->mean)
      MALLOCTEST(g->mean, sizeof(double) * p1);
   if(!g->sd)
      MALLOCTEST(g->sd, sizeof(double) * p1);

   FOPENTEST(in, file_scale, "rb");
   FREADTEST(g->mean, sizeof(double), p1, in);
   FREADTEST(g->sd, sizeof(double), p1, in);

   /* intercept */
   for(k = 0 ; k < NUM_X_LEVELS ; k++)
      g->lookup[k] = 1;

   g->ignore[0] = FALSE;

   for(j = 1 ; j < p1 ; j++)
   {
      if(!(g->ignore[j] = (g->sd[j] == 0)))
      {
	 l1 = j * NUM_X_LEVELS;
	 for(k = 0 ; k < NUM_X_LEVELS ; k++)
	 {
	    l2 = l1 + k;
	    /* scale the levels using the mean/sd for each SNP, but
	     * don't scale the NA level */
	    g->lookup[l2] = (k != X_LEVEL_NA) ? 
		  (k - g->mean[j]) / g->sd[j] : 0;
	 }
      }
   }

   fclose(in);

   return SUCCESS;
}

/* We pre-allocate a large beta matrix, then we trim it down if we don't
 * actually use all the K tasks */
int gmatrix_trim_beta(gmatrix *g)
{
   int m = (g->p + 1) * g->K;
   double *tmp = NULL;

   MALLOCTEST(tmp, m * sizeof(double));
   memcpy(tmp, g->beta, m * sizeof(double));

   FREENULL(g->beta);
   g->beta = tmp;

   return SUCCESS;
}

/* number of training and test samples per fold.
 * folds is a vector of binary indicators, for each sample
 * whether it's in or out of the training set in that cv fold.
 */
void count_fold_samples(int *ntrain, int *ntest,
      double *ntrainrecip, double *ntestrecip,
      int *folds, int nfolds, int n)
{
   int i, k;

   for(k = 0 ; k < nfolds ; k++)
   {
      ntrain[k] = 0;
      for(i = 0 ; i < n ; i++)
	 ntrain[k] += folds[n * k + i];
      ntest[k] = n - ntrain[k];
      ntrainrecip[k] = 0;
      if(ntrain[k] > 0)
	 ntrainrecip[k] = 1.0 / ntrain[k];
      ntestrecip[k] = 0;
      if(ntest[k] > 0)
	 ntestrecip[k] = 1.0 / ntest[k];
   }
}

/* Sets the number of samples for actual use in CD */
void gmatrix_set_ncurr(gmatrix *g)
{
   if(g->nfolds > 1)
   {
      g->ncurr = (g->mode == MODE_TRAIN) ? 
	 g->ntrain[g->fold] : g->ntest[g->fold];
   }
   else
   {
      g->ncurr = g->n; 
   }
   g->ncurr_recip = 1.0 / (g->ncurr - 1);
}

int gmatrix_make_y(gmatrix *g)
{
   return g->y_orig 
      && gmatrix_split_y(g) 
      && (!g->scaley || (g->scaley && gmatrix_scale_y(g)))
      && (g->mode == MODE_PREDICT
	    || (g->mode == MODE_TRAIN && gmatrix_make_fusion(g)));
}

int gmatrix_set_fold(gmatrix *g, int fold)
{
   g->fold = fold;
   gmatrix_set_ncurr(g);

   cache_init(g->xcaches, g->n, g->p + 1, g->maxmem);

   if(!gmatrix_init_lp(g))
      return FAILURE;
   if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
      return FAILURE;
   return gmatrix_make_y(g);
 }

/* zero the lp and adjust the lp-functions */
void gmatrix_zero_model(gmatrix *g)
{
   int i, j, k, n = g->ncurr, p1 = g->p + 1, K = g->K;

   for(j = (long)p1 * K - 1 ; j >= 0 ; --j)
   {
      g->beta[j] = 0;
      g->beta_orig[j] = 0;
      g->active[j] = !g->ignore[j];
   }

   for(i = n * K - 1 ; i >= 0 ; --i)
      g->lp[i] = 0;

   if(g->model == MODEL_LINEAR)
   {
      for(i = n * K - 1 ; i >= 0 ; --i)
	 g->lp[i] = 0;
   }
   else if(g->model == MODEL_LOGISTIC)
   {
      for(i = n * K - 1 ; i >= 0 ; --i)
	 g->lp_invlogit[i] = 0.5; /* 0.5 = 1 / (1 + exp(-0)) */
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = n * K - 1 ; i >= 0 ; --i)
      {
	 g->ylp[i] = -1;    /* y * 0 - 1 = -1  */
	 g->ylp_neg[i] = 1; /* ylp < 0 => true */
	 g->ylp_neg_y[i] = g->y[i];
	 g->ylp_neg_y_ylp[i] = g->ylp_neg_y[i] * g->ylp[i];
      }
   }
   
   for(k = 0 ; k < K ; k++)
   {
      g->lossK[k] = 0;
      g->l1lossK[k] = 0;
   }
   
   g->loss = 0;
   g->l1loss = 0;
   g->floss = 0;
}

/* Create the fusion penalty matrix C from each Y matrix */
int gmatrix_make_fusion(gmatrix *g)
{
#ifdef DEBUG
   char buf[100];
#endif
   int nE = g->K * (g->K - 1) / 2, k, e;
   double s, t;
   
   if(g->K == 1 || g->gamma <= 0)
   {
      g->dofusion = FALSE;
      return SUCCESS;
   }

   g->dofusion = TRUE;

   FREENULL(g->C);
   CALLOCTEST(g->C, nE * g->K, sizeof(double));

   FREENULL(g->pairs);
   CALLOCTEST(g->pairs, nE * 2, sizeof(int));

   FREENULL(g->edges);
   CALLOCTEST(g->edges, (g->K - 1) * g->K, sizeof(int));

   FREENULL(g->flossK);
   CALLOCTEST(g->flossK, g->K, sizeof(double));
   
   if(!gennetwork(g->y, g->ncurr, g->K, g->corthresh, g->cortype, g->C,
	    g->pairs, g->edges))
      return FAILURE;

   FREENULL(g->diagCC);
   CALLOCTEST(g->diagCC, g->K, sizeof(double));
   
   for(k = 0 ; k < g->K ; k++)
   {
      s = 0;
      for(e = 0 ; e < nE ; e++)
      {
	 t = g->C[e + k * nE]; 
	 s += t * t; 
      }
      g->diagCC[k] = s;
   }

   // mapping of edges to vertices (two vertices per edge), zero-based index
   // assumes that C is full size, i.e., all K(K-1)/2 edges are in it
   // pairs <- t(apply(C, 1, function(r) which(r != 0))) - 1

   // each kth column represents which edges task k is involved in
   // edges <- matrix(which(C != 0, arr.ind=TRUE)[,1], K - 1) - 1

#ifdef DEBUG
   sprintf(buf, "C_%02d.txt", g->fold);
   if(!writematrixf(g->C, nE, g->K, buf))
      return FAILURE;
#endif

   return SUCCESS;
}

/* Initialises the LP (linear predictor) */
int gmatrix_init_lp(gmatrix *g)
{
   printf("gmatrix_init_lp: ncurr=%d K=%d\n", g->ncurr, g->K);
   FREENULL(g->lp);
   CALLOCTEST(g->lp, g->ncurr * g->K, sizeof(double));

   if(g->model == MODEL_LINEAR)
   {
      FREENULL(g->err);
      CALLOCTEST(g->err, g->ncurr * g->K, sizeof(double));
   }
   else if(g->model == MODEL_LOGISTIC) 
   {
      FREENULL(g->lp_invlogit);
      CALLOCTEST(g->lp_invlogit, g->ncurr * g->K, sizeof(double));
   } 
   else if(g->model == MODEL_SQRHINGE) 
   {
      FREENULL(g->ylp);
      CALLOCTEST(g->ylp, g->ncurr * g->K, sizeof(double));

      FREENULL(g->ylp_neg);
      CALLOCTEST(g->ylp_neg, g->ncurr * g->K, sizeof(double));

      FREENULL(g->ylp_neg_y);
      CALLOCTEST(g->ylp_neg_y, g->ncurr * g->K, sizeof(double));

      FREENULL(g->ylp_neg_y_ylp);
      CALLOCTEST(g->ylp_neg_y_ylp, g->ncurr * g->K, sizeof(double));
   }
   return SUCCESS;
}

/* g->proportions is an (g->p+1) times 3 matrix, we don't use the first row */
int gmatrix_init_proportions(gmatrix *g)
{
   int i, j, p1 = g->p + 1, n = g->n;
   int x, ngood;
   sample sm;
   long sum0, sum1, sum2;
   int gp0 = p1 * 0, /* :-) */
       gp1 = p1 * 1,
       gp2 = p1 * 2;

   if(!sample_init(&sm))
      return FAILURE;

   MALLOCTEST(g->proportions, p1 * 3 * sizeof(double));
   MALLOCTEST(sm.x, sizeof(double) * n);

   for(j = 1 ; j < p1 ; j++)
   {
      if(!gmatrix_disk_nextcol_raw(g, &sm, j))
	 return FAILURE;

      ngood = sum0 = sum1 = sum2 = 0;
      
      for(i = 0 ; i < n ; i++)
      {
	 x = sm.x[i];
	 if(x != X_LEVEL_NA)
	 {
	    sum0 += (x == 0);
	    sum1 += (x == 1);
	    sum2 += (x == 2);
	    ngood++;
	 }
	 g->proportions[gp0 + j] = (double)sum0 / ngood;
	 g->proportions[gp1 + j] = (double)sum1 / ngood;
	 g->proportions[gp2 + j] = (double)sum2 / ngood;
      }
   }

   FREENULL(sm.x);

   return SUCCESS;
}

void step_regular_logistic(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p)
{
   int i, n = g->ncurr, nki = n * k + n - 1;
   double grad = 0, d2 = 0;
   double *y = g->y,
	  *lp_invlogit = g->lp_invlogit,
	  *x = s->x;

   /* compute 1st and 2nd derivatives wrt task k */
   for(i = n - 1 ; i >= 0 ; --i)
   {
      grad += x[i] * (lp_invlogit[nki] - y[nki]);
      d2 += x[i] * x[i] * lp_invlogit[nki] * (1 - lp_invlogit[nki]);
      --nki;
   }

   grad *= g->ncurr_recip;
   d2 *= g->ncurr_recip;

   if(d2 == 0)
   {
      *d1_p = 0;
      *d2_p = 1;
   }
   else
   {
      *d1_p = grad;
      *d2_p = d2;
   }
}

/* In linear regression, for standardised inputs x, the 2nd derivative is
 * always N since it is the sum of squares \sum_{i=1}^N x_{ij}^2 =
 * \sum_{i=1}^N 1 = N
 */
void step_regular_linear(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p)
{
   int i, n = g->ncurr, nki = n * k + n - 1;
   double grad = 0;
   double *err = g->err,
	  *restrict x = s->x;

   /* compute gradient wrt task k*/
   for(i = n - 1 ; i >= 0 ; --i)
   {
      grad += x[i] * err[nki];
      --nki;
   }

   *d1_p = grad * g->ncurr_recip;
   *d2_p = (n - 1) * g->ncurr_recip;
}

/*
 * Squared hinge loss, assumes y \in {-1,1},
 * and that X is scaled so that the 2nd derivative is always <=N
 */
void step_regular_sqrhinge(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p)
{
   int i, n = g->ncurr, nki = n * k + n - 1;
   double grad = 0;
   const double *x = s->x,
		*ylp_neg_y_ylp = g->ylp_neg_y_ylp;

   /* compute gradient wrt task k */
   for(i = n - 1 ; i >= 0 ; --i)
   {
      grad += ylp_neg_y_ylp[nki] * x[i];
      --nki;
   }

   *d1_p = grad * g->ncurr_recip;
   *d2_p = (n - 1) * g->ncurr_recip;
}

/* Update linear predictor and related variables.
 *
 * The updates where x is null are for the get_lambda1max_gmatrix updates
 * i.e. x_i=1 for all i
 */
void updatelp(gmatrix *g, const double update,
      const double *restrict x, int j, int k)
{
   int i, n = g->ncurr, nki;
   double *restrict err = g->err,
	  *restrict lp_invlogit = g->lp_invlogit,
	  *restrict lp = g->lp,
	  *restrict y = g->y,
	  *restrict ylp_neg_y_ylp = g->ylp_neg_y_ylp;
   double loss = 0, ylp = 0, explp = 0;

   nki = n * k + n - 1;
   if(g->model == MODEL_LINEAR)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[nki] += x[i] * update;
	 err[nki] = lp[nki] - y[nki];
	 loss += err[nki] * err[nki];
	 --nki;
      }
   }
   else if(g->model == MODEL_LOGISTIC)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[nki] += x[i] * update;
	 explp = exp(lp[nki]);
	 lp_invlogit[nki] = explp / (1 + explp);
	 loss += log(1 + explp) - y[nki] * lp[nki];
	 --nki;
      }
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[nki] += x[i] * update;
	 ylp = y[nki] * lp[nki] - 1;
	 ylp_neg_y_ylp[nki] = y[nki] * ylp;
	 loss += ylp * ylp;
	 --nki;
      }
   }

   g->lossK[k] = loss * g->ncurr_recip;
}

/*
 * Set x to a pointer to the jth SNP. If it's in the cache, return SUCCESS,
 * otherwise return FAILURE. This way we can read data into the cache by
 * reading it into x.
 */
int cache_get(cache *ca, int j, double **x)
{
   int m = ca->mapping[j];
   int r;
   ca->counter[j]++;

   /* in cache */
   if(m != CACHE_NOT_EXISTS)
   {
      *x = (ca->x + m * ca->n);
      return SUCCESS;
   }
   else /* provide pointer to access the cache for putting data */
   {
      /* check whether we already have data at this spot, and if it is used
       * then reclaim it for another SNP if the new one has been in the cache
       * more often */

      r = ca->revmapping[ca->lastfree]; 

      if(r != CACHE_NOT_EXISTS)
      {
	 /* new data seen more often, replace old with new */ 
	 if(ca->counter[j] > ca->counter[r])
	 {
	    ca->mapping[r] = CACHE_NOT_EXISTS;
	    ca->mapping[j] = ca->lastfree;
	    *x = ca->x + ca->mapping[j] * ca->n;
	    ca->revmapping[ca->lastfree] = j;
	    ca->lastfree++;
	 }
	 else /* provide temporary storage but don't put in cache */
	    *x = ca->tmp;
      }
      else
      {
	 ca->mapping[j] = ca->lastfree;
	 *x = ca->x + ca->mapping[j] * ca->n;
	 ca->revmapping[ca->lastfree] = j;
	 ca->lastfree++;
      }

      /* we've run of out space, go back to beginning of buffer and allow it
       * to be over-written with new data */
      if(ca->lastfree >= ca->nbins)
	 ca->lastfree = 0;

      return FAILURE;
   }
}

/*
 * n: number of samples per SNP
 * p: number of SNPs
 */
int cache_init(cache *ca, int n, int p, long maxmem)
{
   int i;

   cache_free(ca);

   ca->nbins = (int)(maxmem / sizeof(double) / n);
   ca->n = n;
   ca->lastfree = 0;

   printf("[cache_init n:%d p:%d nbins:%d]\n", n, p, ca->nbins);
   fflush(stdout);

   MALLOCTEST(ca->mapping, p * sizeof(int));
   MALLOCTEST(ca->revmapping, ca->nbins * sizeof(int));
   CALLOCTEST(ca->counter, p, sizeof(int));
   MALLOCTEST(ca->x, ca->nbins * n * sizeof(double));
   MALLOCTEST(ca->tmp, n * sizeof(double));

   for(i = ca->nbins * n - 1 ; i >= 0 ; --i)
      ca->x[i] = -0.123456789; /* magic number for debugging */

   for(i = 0 ; i < p ; i++)
      ca->mapping[i] = CACHE_NOT_EXISTS;

   for(i = 0 ; i < ca->nbins ; i++)
      ca->revmapping[i] = CACHE_NOT_EXISTS;

   return SUCCESS;
}

void cache_free(cache *ca)
{
   FREENULL(ca->mapping);
   FREENULL(ca->revmapping);
   FREENULL(ca->tmp);
   FREENULL(ca->x);
   FREENULL(ca->counter);
}

void count_cases(gmatrix *g)
{
   int i = 0;
   for(i = g->n * g->K - 1 ; i >= 0 ; --i)
      g->ncases += g->y_orig[i] == 1;
}

int grad_array_compare(const void *a, const void *b)
{
   double ax = ((const struct pair*)a)->value;
   double bx = ((const struct pair*)b)->value;

   if(ax > bx)
      return -1;
   else if(ax < bx)
      return 1;
   else
      return 0;
}

void gmatrix_sort_grad_array(gmatrix *g)
{
   int k;
   int p1 = g->p + 1;

   for(k = 0 ; k < g->K ; k++)
   {
      qsort(g->grad_array + k * p1, p1, sizeof(pair), grad_array_compare);
   }
}

