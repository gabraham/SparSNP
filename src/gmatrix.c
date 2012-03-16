/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "gmatrix.h"
#include "ind.h"
#include "util.h"

int sample_init(sample *s)
{
   s->n = 0;
   s->x = NULL;
   s->y = NULL;
   s->intercept = FALSE;
   return SUCCESS;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int model, int modeltype,
      short encoded, short binformat, char *folds_ind_file,
      short mode, loss_pt loss_pt_func, char *subset_file, 
      char *famfilename)
{
   int i, j, p1;

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
   g->loss = 0;
   g->lp = NULL;
   g->ylp = NULL;
   g->ylp_neg = NULL;
   g->ylp_neg_y = NULL;
   g->ylp_neg_y_ylp = NULL;
   g->lp_invlogit = NULL;
   g->loss_func = NULL;
   g->loss_pt_func = loss_pt_func;
   g->lookup = NULL;
   g->intercept = NULL;
   g->mean = NULL;
   g->sd = NULL;
   g->tmp = NULL;
   g->xtmp = NULL;
   /*g->ytmp = NULL;*/
   g->x = NULL;
   g->xthinned = NULL;
   g->ignore = NULL;
   g->yformat = yformat;
   g->beta = NULL;
   g->active = NULL;
   g->encoded = encoded;
   g->nencb = (int)ceil((double)n / PACK_DENSITY);
   g->encbuf = NULL;
   /*g->decode = &decode;*/
   g->binformat = binformat;
   g->folds_ind_file = folds_ind_file;
   g->nfolds = 1;
   g->folds = NULL;
   g->fold = 0; 
   g->ntrain = NULL;
   g->ntest = NULL;
   g->ntrainrecip = NULL;
   g->ntestrecip = NULL;
   g->ncurr = 0;
   g->ncurr_j = NULL; /* current number of samples, depending on
			 whether we're in training or testing and
			 which fold we're in, etc. */
   g->ncurr_recip_j = NULL; /* reciprocal of ncurr to
			       allow multiplication instead
			       of division */
   g->mode = mode;
   g->nseek = sizeof(dtype) * (encoded ? g->nencb : g->n);
   g->beta_orig = NULL;
   g->numnz = NULL;
   g->xcaches = NULL;
   g->ncases = 0;
   g->folds_ind = NULL;

   g->famfilename = famfilename;

   MALLOCTEST(g->map, sizeof(mapping));
   mapping_init(g->map);

   srand48(time(NULL));

   CALLOCTEST(g->beta_orig, p1, sizeof(double));
   CALLOCTEST(g->ncurr_j, p1, sizeof(int));
   CALLOCTEST(g->ncurr_recip_j, p1, sizeof(double));
   
   if(filename)
      FOPENTEST(g->file, filename, "rb");

   /* gmatrix_setup_folds changes g->nfolds */
   if(!gmatrix_setup_folds(g))
      return FAILURE;

   MALLOCTEST(g->intercept, sizeof(double) * g->n);
   for(i = n - 1 ; i >= 0 ; --i)
      g->intercept[i] = 1.0;

   MALLOCTEST(g->tmp, sizeof(dtype) * g->nencb * PACK_DENSITY);

   if(encoded)
      MALLOCTEST(g->encbuf, sizeof(unsigned char) * g->nencb);

   MALLOCTEST(g->xcaches, g->nfolds * sizeof(cache));
   for(i = 0 ; i < g->nfolds ; i++)
   {
      /* TODO: uses up g->n cells even though we don't really not all since in
       * crossval there will be fewer training/test samples, but easier than
       * figuring out exactly how many to use */
      cache_init(g->xcaches + i, g->n, g->p + 1);
   }

   g->nextcol = gmatrix_disk_nextcol;
   if(!gmatrix_reset(g))
      return FAILURE;

   if(g->binformat == BINFORMAT_PLINK)
   {
      printf("Using PLINK binary format\n");
      printf("FAM file: %s\n", g->famfilename);
      g->offset = 3 - g->nseek;
      if(g->famfilename)
      {
	 if(!gmatrix_fam_read_y(g))
	    return FAILURE;
	 gmatrix_plink_check_pheno(g);
      }
   }
   else
   {
      g->offset = 0;
      /*g->decode = &decode;*/
      if(!gmatrix_disk_read_y(g))
	 return FAILURE;
   }

   if(g->mode == MODE_TRAIN && g->modeltype == MODELTYPE_CLASSIFICATION)
   {
      count_cases(g);
      printf("Found %d cases, %d controls\n", g->ncases, g->n - g->ncases);
   }

   CALLOCTEST(g->beta, p1, sizeof(double));
   CALLOCTEST(g->active, p1, sizeof(int));
   CALLOCTEST(g->ignore, p1, sizeof(int));

   /* don't scale in prediction mode */
   if(g->mode == MODE_TRAIN)
   {
      g->scalefile = scalefile;
      if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
	 return FAILURE;
   }
   else
      g->scalefile = NULL;

   for(j = g->p ; j >= 0 ; --j)
      g->active[j] = !g->ignore[j];

   gmatrix_set_ncurr(g);
   
   if(g->mode == MODE_TRAIN && g->y_orig && !gmatrix_split_y(g))
      return FAILURE;

   if(!gmatrix_init_lp(g))
      return FAILURE;

   /*MALLOCTEST(g->xtmp, sizeof(double) * g->n);*/
   /*MALLOCTEST(g->ytmp, sizeof(double) * g->n);*/

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
   int i;
   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   FREENULL(g->mean);
   FREENULL(g->sd);
   FREENULL(g->y);
   FREENULL(g->y_orig);
   FREENULL(g->xtmp);
   /*FREENULL(g->ytmp);*/
   FREENULL(g->x);
   FREENULL(g->xthinned);
   FREENULL(g->ignore);
   FREENULL(g->tmp);
   FREENULL(g->intercept);
   FREENULL(g->lookup);
   FREENULL(g->lp);
   FREENULL(g->ylp);
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
   FREENULL(g->ncurr_j);
   FREENULL(g->ncurr_recip_j);
   
   if(g->xcaches)
   {
      for(i = 0 ; i < g->nfolds ; i++)
	 cache_free(g->xcaches + i);
   }
   FREENULL(g->xcaches);

   mapping_free(g->map);
   FREENULL(g->map);
}

/* y_orig is the original vector of labels/responses, it 
 * stays in memory and NEVER changes after reading it from disk */
int gmatrix_disk_read_y(gmatrix *g)
{
   int i, n = g->n, n1 = g->n - 1;

   /* read all of y the first time we see it, then skip it */
   CALLOCTEST(g->y_orig, n, sizeof(double));

   /* The y vector may be byte packed */
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   /*decode_plink(g->tmp, g->encbuf, g->nencb);*/
   decode_plink_mapping(g->map, g->tmp, g->encbuf, g->nencb);
 
   if(g->yformat == YFORMAT01 || g->modeltype == MODELTYPE_REGRESSION) 
   {
      for(i = n1 ; i >= 0 ; --i)
	 g->y_orig[i] = (double)g->tmp[i];
   } 
   else if(g->yformat == YFORMAT11)
   {  /*  -1/1  */
      for(i = n1 ; i >= 0 ; --i)
	 g->y_orig[i] = 2.0 * g->tmp[i] - 1.0;
   }

   return SUCCESS;
}

/* y is a copy of y_orig as needed for crossval - in training y are the
 * training samples, in prediction y are the testing samples
 */
int gmatrix_split_y(gmatrix *g)
{
   int i, n = g->n, n1 = g->n - 1;
   int k = g->ncurr - 1;

   if(g->ncurr <= 0) {
      fprintf(stderr, "gmatrix_split_y: ncurr can't be zero");
      return FAILURE;
   }

   FREENULL(g->y);

   MALLOCTEST(g->y, sizeof(double) * g->ncurr);

   /* represent y as 0/1 or -1/1, allowing for different
    * samples sizes if cross-validation is used */
   if(g->nfolds > 1) /* no cv */
   {
      for(i = n1 ; i >= 0 ; --i)
	 if(g->folds_ind[g->fold * n + i])
	    g->y[k--] = g->y_orig[i];
   }
   else 
   {
      for(i = n1 ; i >= 0 ; --i)
	 g->y[i] = g->y_orig[i];
   }

   return SUCCESS;
}

/*
 * Read all the data into a preallocated row-major n by m matrix, for n
 * samples and m variables.
 * The variables to be read are determined by the p+1 array ``ind''.
 *
 */
int gmatrix_read_matrix(gmatrix *g, int *ind, int m)
{
   int i, j, k = 0,
       p1 = g->p + 1;
   sample sm;

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

int rand_geno()
{
   return lrand48() % 3;
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
      s->x2 = g->intercept;

      s->y = g->y;
      s->n = g->ncurr;

      return SUCCESS;
   }

   ret = cache_get(g->xcaches + g->fold, j, &xtmp);

   if(ret == SUCCESS)
   {
      s->x = xtmp;
      /*s->y = g->ytmp;*/
      return SUCCESS;
   }

   /* Get data from disk and unpack, skip y */
   seek = j * g->nseek + g->offset;
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
		  /*g->ytmp[ngood++] = g->y_orig[i];*/
	       }
	    }
	    s->n = ngood;
	 }
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
	          xtmp[k--] = (d == X_LEVEL_NA ? (double)rand_geno() : (double)d);
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
		     /*g->ytmp[ngood++] = g->y_orig[i];*/
		  }
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

/* Reads a matrix of phenotypes in plain format (no FAM columns)
 *
 * We read the entire dataset into memory
 * */
/*int gmatrix_read_pheno_matrix(gmatrix *g)
{
   
}*/

/* Reads a matrix of phenotypes in plink FAM format */
/*int gmatrix_fam_read_pheno_matrix(gmatrix *g)
{
      
}*/

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

/*
 * plink phenotypes can be 0/1 and 1/2, so check for 1/2 and convert as
 * necessary
 */
int gmatrix_plink_check_pheno(gmatrix *g)
{
   int i = 0;
   int twofound = FALSE;

   /* don't change inputs for regression */
   if(g->modeltype == MODELTYPE_REGRESSION)
      return SUCCESS;

   for(i = 0 ; i < g->n ; i++)
      if((twofound = (g->y_orig[i] == 2)))
	 break;

   if(twofound)
   {
      for(i = 0 ; i < g->n ; i++)
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
      for(i = 0 ; i < g->n ; i++)
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
   /*g->ncurr_recip = 1.0 / g->ncurr;*/
   g->ncurr_recip = 1.0 / (g->ncurr - 1);
}

int gmatrix_set_fold(gmatrix *g, int fold)
{
   g->fold = fold;
   gmatrix_set_ncurr(g);
   if(!gmatrix_init_lp(g))
      return FAILURE;
   if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
      return FAILURE;
   return (g->y_orig) && gmatrix_split_y(g);
}

/* zero the lp and adjust the lp-functions */
void gmatrix_zero_model(gmatrix *g)
{
   int i, j, n = g->ncurr, p1 = g->p + 1;

   for(j = p1 - 1 ; j >= 0 ; --j)
   {
      g->beta[j] = 0;
      g->beta_orig[j] = 0;
      g->active[j] = !g->ignore[j];
   }

   for(i = n - 1 ; i >= 0 ; --i)
      g->lp[i] = 0;

   if(g->model == MODEL_LOGISTIC)
   {
      for(i = n - 1 ; i >= 0 ; --i)
	 g->lp_invlogit[i] = 0.5; /* 0.5 = 1 / (1 + exp(-0)) */
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = n - 1 ; j >= 0 ; --i)
      {
	 g->ylp[i] = -1;    /* y * 0 - 1 = -1  */
	 g->ylp_neg[i] = 1; /* ylp < 0 => true */
	 g->ylp_neg_y[i] = g->y[i];
	 g->ylp_neg_y_ylp[i] = g->ylp_neg_y[i] * g->ylp[i];
      }
   }

   g->loss = 0;
}

int gmatrix_init_lp(gmatrix *g)
{
   FREENULL(g->lp);
   CALLOCTEST(g->lp, g->ncurr, sizeof(double));

   if(g->model == MODEL_LOGISTIC) 
   {
      FREENULL(g->lp_invlogit);
      CALLOCTEST(g->lp_invlogit, g->ncurr, sizeof(double));
   } 
   else if(g->model == MODEL_SQRHINGE) 
   {
      FREENULL(g->ylp);
      CALLOCTEST(g->ylp, g->ncurr, sizeof(double));

      FREENULL(g->ylp_neg);
      CALLOCTEST(g->ylp_neg, g->ncurr, sizeof(double));

      FREENULL(g->ylp_neg_y);
      CALLOCTEST(g->ylp_neg_y, g->ncurr, sizeof(double));

      FREENULL(g->ylp_neg_y_ylp);
      CALLOCTEST(g->ylp_neg_y_ylp, g->ncurr, sizeof(double));
   }
   return SUCCESS;
}

/* In linear regression, for standardised inputs x, the 2nd derivative is
 * always N since it is the sum of squares \sum_{i=1}^N x_{ij}^2 =
 * \sum_{i=1}^N 1 = N
 */
double step_regular_linear(sample *s, gmatrix *g)
{
   int i;
   double grad = 0;
   double *restrict x = s->x, 
          *restrict lp = g->lp,
	  *restrict y = g->y;

   /* compute gradient */
   for(i = g->ncurr - 1 ; i >= 0 ; --i)
      grad += x[i] * (lp[i] - y[i]);

   return grad * g->ncurr_recip;
}

double step_regular_logistic(sample *s, gmatrix *g)
{
   int i, n = g->ncurr;
   double grad = 0, d2 = 0;
   double *restrict y = g->y,
	  *restrict lp_invlogit = g->lp_invlogit,
	  *restrict x = s->x;

   /* compute 1st and 2nd derivatives */
   for(i = n - 1 ; i >= 0 ; --i)
   {
      grad += x[i] * (lp_invlogit[i] - y[i]);
      d2 += x[i] * x[i] * lp_invlogit[i] * (1 - lp_invlogit[i]);
   }

   grad *= g->ncurr_recip;
   d2 *= g->ncurr_recip;

   if(d2 == 0)
      return 0;
   return grad / d2;
}

/*
 * Squared hinge loss, assumes y \in {-1,1},
 * and that X is scaled so that the 2nd derivative is always <=N
 */
double step_regular_sqrhinge(sample *s, gmatrix *g)
{
   int i, n = g->ncurr;
   double grad = 0;
   const double *restrict x = s->x,
		*restrict ylp_neg_y_ylp = g->ylp_neg_y_ylp;

   /* compute gradient */
   for(i = n - 1 ; i >= 0 ; --i)
      grad += ylp_neg_y_ylp[i] * x[i];

   return grad * g->ncurr_recip; /* avoid division */
}

/* Update linear predictor and related variables.
 *
 * The updates where x is null are for the get_lambda1max_gmatrix updates
 * i.e. x_i=1 for all i
 */
void updatelp(gmatrix *g, const double update,
      const double *restrict x, int j)
{
   int i, n = g->ncurr;
   double err; 
   double *restrict lp_invlogit = g->lp_invlogit,
	  *restrict lp = g->lp,
	  *restrict y = g->y,
	  *restrict ylp_neg_y_ylp = g->ylp_neg_y_ylp;
   double loss = 0, ylp = 0;

   if(g->model == MODEL_LINEAR)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[i] += x[i] * update;
	 err = lp[i] - y[i];
	 loss += err * err;
      }
   }
   else if(g->model == MODEL_LOGISTIC)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[i] += x[i] * update;
	 lp_invlogit[i] = 1 / (1 + exp(-lp[i]));
	 loss += log(1 + exp(lp[i])) - y[i] * lp[i];
      }
   }
   else if(g->model == MODEL_SQRHINGE)
   {
      for(i = n - 1 ; i >= 0 ; --i)
      {
	 lp[i] += x[i] * update;
	 ylp = y[i] * lp[i] - 1;

	 if(ylp < 0)
	 {
	    ylp_neg_y_ylp[i] = y[i] * ylp;
	    loss += ylp * ylp;
	 }
	 else
	    ylp_neg_y_ylp[i] = 0;
      }
   }

   g->loss = loss * g->ncurr_recip;
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
int cache_init(cache *ca, int n, int p)
{
   int i;

   ca->nbins = CACHE_MAX_MEM / sizeof(double) / n;
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
      ca->x[i] = -0.123456789;

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
   for(i = 0 ; i < g->n ; i++)
      g->ncases += g->y_orig[i] == 1;
}

