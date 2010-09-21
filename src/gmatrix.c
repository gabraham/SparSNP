#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"
#include "coder.h"
#include "ind.h"

int sample_init(sample *s, int n, short inmemory)
{
   s->inmemory = inmemory;
   s->x = NULL;
   s->cached = FALSE;
   s->intercept = FALSE;
   /*MALLOCTEST(s->x, sizeof(double) * n);*/
   return SUCCESS;
}

void sample_free(sample *s)
{
   if(s->intercept)
      return;

   if(!s->cached && s->x)
   {
      free(s->x);
      s->x = NULL;
   }

   s->x2 = NULL;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      short inmemory,  char *scalefile, short yformat, int model,
      short encoded, short binformat, char *folds_ind_file, int nfolds,
      short mode)
{
   int i, j;

   g->model = model;
   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;
   g->yidx = 0;
   g->y = NULL;
   g->lp = NULL;
   g->ylp = NULL;
   g->ylp_neg = NULL;
   g->lp_invlogit = NULL;
   g->inmemory = inmemory;
   g->scalefile = scalefile;
   g->lookup = NULL;
   g->intercept = NULL;
   g->mean = NULL;
   g->sd = NULL;
   g->tmp = NULL;
   g->ignore = NULL;
   g->yformat = yformat;
   g->beta = NULL;
   g->active = NULL;
   g->encoded = encoded;
   g->nencb = (int)ceil((double)n / PACK_DENSITY);
   g->encbuf = NULL;
   g->decode = &decode;
   g->binformat = binformat;
   g->folds_ind_file = folds_ind_file;
   g->nfolds = nfolds;
   g->folds = NULL;
   g->fold = 0; 
   g->ntrain = NULL;
   g->ntest = NULL;
   g->ntrainrecip = NULL;
   g->ntestrecip = NULL;
   g->ncurr = 0; /* current number of samples, depending on whether we're in
		    training or testing and which fold we're in, etc. */
   g->ncurr_recip = 0.0; /* reciprocal of ncurr to allow multiplication instead
			    of division */
   g->mode = mode;
   g->nseek = sizeof(dtype) * (encoded ? g->nencb : g->n);

   MALLOCTEST(g->ca, sizeof(cache));
   if(!cache_init(g->ca, g->p + 1)) /* +1 not for intercept, which isn't
				       stored but we need to maintain
				       indexing consistency */
      return FAILURE;

   if(filename)
      FOPENTEST(g->file, filename, "rb");

   if(!gmatrix_setup_folds(g))
      return FAILURE;

   MALLOCTEST(g->intercept, sizeof(double) * g->n);
   for(i = 0 ; i < n ; i++)
      g->intercept[i] = 1.0;

   MALLOCTEST(g->tmp, sizeof(dtype) * g->n);

   if(encoded)
      MALLOCTEST(g->encbuf, sizeof(unsigned char) * g->nencb);

   if(g->binformat == BINFORMAT_PLINK)
      g->decode = &decode_plink;

   if(inmemory)
   {
      if(!gmatrix_load(g))
	 return FAILURE;
      g->nextcol = gmatrix_mem_nextcol;
   }
   else
      g->nextcol = gmatrix_disk_nextcol;

   CALLOCTEST(g->beta, g->p + 1, sizeof(double));
   CALLOCTEST(g->active, g->p + 1, sizeof(int));
   CALLOCTEST(g->ignore, g->p + 1, sizeof(int));

   if(g->scalefile && !gmatrix_read_scaling(g, scalefile))
      return FAILURE;

   for(j = g->p ; j >= 0 ; --j)
      g->active[j] = !g->ignore[j];

   gmatrix_set_ncurr(g);

   if(!gmatrix_init_lp(g))
      return FAILURE;

   return SUCCESS;
}

int gmatrix_setup_folds(gmatrix *g)
{
   if(g->folds_ind_file 
	 && !(g->nfolds = ind_getfolds(g->folds_ind_file)))
      return FAILURE;

   MALLOCTEST(g->folds, sizeof(int) * g->n * g->nfolds);
   if(g->nfolds > 1 && !ind_read(g->folds_ind_file, g->folds, g->n, g->nfolds))
      return FAILURE;

   MALLOCTEST(g->ntrain, sizeof(int) * g->nfolds);
   MALLOCTEST(g->ntest, sizeof(int) * g->nfolds);
   MALLOCTEST(g->ntrainrecip, sizeof(double) * g->nfolds);
   MALLOCTEST(g->ntestrecip, sizeof(double) * g->nfolds);

   if(g->folds_ind_file && g->nfolds > 1)
      count_fold_samples(g->ntrain, g->ntest,
	    g->ntrainrecip, g->ntestrecip,
	    g->folds, g->nfolds, g->n);
   else /* train and test sets are the same */
   {
      g->ntrain[0] = g->ntest[0] = g->n;
      g->ntrainrecip[0] = g->ntestrecip[0] = 1.0 / g->n;
   }
   return SUCCESS;
}

void gmatrix_free(gmatrix *g)
{
   int j;

   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   if(g->mean)
      free(g->mean);

   if(g->sd)
      free(g->sd);

   g->mean = g->sd = NULL;

   if(g->y)
      free(g->y);
   g->y = NULL;

   if(g->inmemory && g->x)
   {
      for(j = 0 ; j < g->p + 1 ; j++)
      {
	 if(g->x[j])
	    free(g->x[j]);
      }
      free(g->x);
      g->x = NULL;
   }

   if(g->ignore)
      free(g->ignore);
   g->ignore = NULL;

   if(g->tmp)
      free(g->tmp);
   g->tmp = NULL;

   if(g->intercept)
      free(g->intercept);
   g->intercept = NULL;

   if(g->lookup)
      free(g->lookup);
   g->lookup = NULL;

   if(g->lp)
      free(g->lp);
   g->lp = NULL;

   if(g->ylp)
      free(g->ylp);
   g->ylp = NULL;

   if(g->ylp_neg)
      free(g->ylp_neg);
   g->ylp_neg = NULL;

   if(g->lp_invlogit)
      free(g->lp_invlogit);
   g->lp_invlogit = NULL;

   if(g->beta)
      free(g->beta);
   g->beta = NULL;

   if(g->encbuf)
      free(g->encbuf);
   g->encbuf = NULL;

   if(g->folds)
      free(g->folds);
   g->folds = NULL;

   if(g->ntrain)
      free(g->ntrain);
   g->ntrain = NULL;

   if(g->ntest)
      free(g->ntest);
   g->ntest = NULL;

   if(g->ntrainrecip)
      free(g->ntrainrecip);
   g->ntrainrecip = NULL;

   if(g->ntestrecip)
      free(g->ntestrecip);
   g->ntestrecip = NULL;

   if(g->active)
      free(g->active);
   g->active = NULL;

   if(g->ca)
   {
      cache_free(g->ca);
      free(g->ca);
   }
   g->ca = NULL;
}

int gmatrix_disk_nexty(gmatrix *g)
{
   int i, n = g->n, n1 = g->n - 1;
   int k = g->ncurr - 1;

   if(g->binformat == BINFORMAT_BIN) {
      /* read all of y the first time we see it, then skip it */
      if(!g->y) {
	 CALLOCTEST(g->y, n, sizeof(double));

	 /* The y vector may be byte packed */
	 FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
	 g->decode(g->tmp, g->encbuf, g->nencb);

	 /* represent y as 0/1 or -1/1, allowing for different samples sizes
	  * if cross-validation is used */
	 if(g->nfolds < 2) /* no cv */
	 {
	    if(g->yformat == YFORMAT01)
	       for(i = n1 ; i >= 0 ; --i)
		  g->y[i] = (double)g->tmp[i];
	    else  /*  -1/1  */
	       for(i = n1 ; i >= 0 ; --i)
		  g->y[i] = 2.0 * g->tmp[i] - 1.0;
	 } 
	 else  /* with cv */
	 {
	    if(g->yformat == YFORMAT01)
	    {
	       if(g->mode == MODE_TRAIN) {
		  for(i = n1 ; i >= 0 ; --i)
		     if(g->folds[g->fold * n + i])
			g->y[k--] = (double)g->tmp[i];
	       } else {
		  for(i = n1 ; i >= 0 ; --i)
		     if(!g->folds[g->fold * n + i])
			g->y[k--] = (double)g->tmp[i];
	       }
	    }
	    else /*  y \in -1/1  */
	    {
	       if(g->mode == MODE_TRAIN) {
		  for(i = n1 ; i >= 0 ; --i)
		     if(g->folds[g->fold * n + i])
			g->y[k--] = 2.0 * g->tmp[i] - 1.0;
	       } else {
		  for(i = n1 ; i >= 0 ; --i)
		     if(!g->folds[g->fold * n + i])
			g->y[k--] = 2.0 * g->tmp[i] - 1.0;
	       }
	    }
	 }

      } else { /* don't read y again */
	 FSEEKOTEST(g->file, g->nseek, SEEK_CUR);
      }
   } else { /* skip plink headers */
      FSEEKOTEST(g->file,
	    sizeof(dtype) * PLINK_HEADER_SIZE, SEEK_CUR);
   }

   return SUCCESS;
}

/* big ugly function
 */
int gmatrix_disk_nextcol(gmatrix *g, sample *s, int skip)
{
   int i, l1;
   int n = g->n, n1 = g->n - 1, p1 = g->p + 1;
   int ncurr = g->ncurr;
   int k = g->ncurr - 1;

   if(g->j == p1 && !gmatrix_reset(g))
      return FAILURE;

   /* either read the y vector, or skip rows and/or headers */
   if(g->j == 0)
   {
      s->intercept = TRUE;
      if(!gmatrix_disk_nexty(g))
	 return FAILURE;

      /* No need to copy values, just assign the intercept vector,
       * but first, free any old data we have from previous
       * iterations */
      if(s->x)
	 free(s->x);

      /* don't need to worry about cross-validation etc because the intercept
       * is all 1s */
      s->x = g->intercept;

      g->j++;
      return SUCCESS;
   }

   /* this isn't an intercept, so we need to copy actual values but first make
    * room for them */
   s->intercept = FALSE;
   if(g->j == 1) {
      MALLOCTEST(s->x, sizeof(double) * ncurr);
   }

   if(skip) 
   {
      FSEEKOTEST(g->file, g->nseek, SEEK_CUR);
      g->j++;
      return SUCCESS;
   }

   /* read the data, unpack if necessary */
   /*if(g->encoded) {*/
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   g->decode(g->tmp, g->encbuf, g->nencb);
   /*} else {
     FREADTEST(g->tmp, sizeof(dtype), n, g->file);
     }*/

   /* Get the scaled versions of the genotypes */
   if(g->nfolds < 2) { /* no cv */
      if(g->scalefile) {
         /* Get the scaled value instead of the original value */
         l1 = g->j * NUM_X_LEVELS;
         for(i = n1 ; i >= 0 ; --i)
	    s->x[i] = g->lookup[l1 + g->tmp[i]];
      } else {
	 for(i = n1 ; i >= 0 ; --i)
	    if(g->folds[g->fold * n + i])
	       s->x[k--] = (double)g->tmp[i];
      }
   } else { /* cv */
      /* Get the scaled value instead of the original value */
      if(g->scalefile) {
	 l1 = g->j * NUM_X_LEVELS;
	 for(i = n1 ; i >= 0 ; --i)
	    /* different between train and test */
	    if((g->mode == MODE_PREDICT) ^ g->folds[g->fold * n + i])
	       s->x[k--] = g->lookup[l1 + g->tmp[i]];
      } else {
      	 for(i = n1 ; i >= 0 ; --i)
	    if((g->mode == MODE_PREDICT) ^ g->folds[g->fold * n + i])
	       s->x[k--] = (double)g->tmp[i];
      }
   }

   g->j++;
   return SUCCESS;
}

int gmatrix_mem_nextcol(gmatrix *g, sample *s, int skip)
{
   if(g->j == g->p + 1)
   {
      g->i = g->j = 0;
   }

   s->x = g->x[g->j];
   g->j++;

   return SUCCESS;
}

int gmatrix_load(gmatrix *g)
{
   int i, j;
   sample s;

   if(!sample_init(&s, g->n, FALSE))
      return FAILURE;

   MALLOCTEST(g->x, sizeof(double*) * (g->p + 1))

      for(j = 0 ; j < g->p + 1 ; j++)
      {
	 MALLOCTEST(g->x[j], sizeof(double) * g->n)
	    if(!gmatrix_disk_nextcol(g, &s, FALSE))
	       return FAILURE;

	 for(i = 0 ; i < g->n ; i++)
	    g->x[j][i] = (double)s.x[i];
      }

   sample_free(&s);

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

/* Populate lookup tables for SNP levels 0, 1, 2
 */
int gmatrix_read_scaling(gmatrix *g, char *file_scale)
{
   int j, k, l1, l2, p1 = g->p + 1;
   FILE *in;

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
	    g->lookup[l2] = (k - g->mean[j]) / g->sd[j];
	 }
      }
   }

   fclose(in);

   return SUCCESS;
}

/* number of training and test samples per fold */
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

int cache_init(cache *ca, int nkeys)
{
   int i;
   ca->nkeys = nkeys;
   ca->size = HASH_SIZE;
   ca->active = 0;

   MALLOCTEST(ca->buckets, sizeof(bucket) * ca->size);
   CALLOCTEST(ca->weights, nkeys, sizeof(ca->weights));
   CALLOCTEST(ca->keys, nkeys, sizeof(ca->keys));

   for(i = 0 ; i < ca->size ; i++)
   {
      ca->buckets[i].active = FALSE;
      ca->buckets[i].key = 0.0;
   }

   return SUCCESS;
}

void cache_free(cache *ca)
{
   if(ca->buckets)
      free(ca->buckets);
   ca->buckets = NULL;
   ca->active = 0;

   if(ca->weights)
      free(ca->weights);
   ca->weights = NULL;

   if(ca->keys)
      free(ca->keys);
   ca->keys = NULL;
}

int cache_put(cache *ca, int key, double *value)
{
   bucket *bk = NULL;
   int hval = hash(key);

   bk = &ca->buckets[hval];

   if(!bk->active || ca->weights[key] > ca->weights[bk->key])
   {
      ca->keys[bk->key] = FALSE;
      bk->key = key;
      bk->value = value;
      ca->keys[key] = TRUE;
      if(!bk->active)
      {
	 bk->active = TRUE;
	 ca->active++;
      }
      return SUCCESS;
   }

   return FAILURE;
}

double* cache_get(cache *ca, int key)
{
   bucket *bk;
   int hval = hash(key);

   /* ca->weights[key]++; */

   bk = &ca->buckets[hval];
   if(bk && bk->active)
      return bk->value;

   return NULL;
}

int hash(int key)
{
   return key >> 28;
}

/* Sets the number of samples for actual use in CD */
void gmatrix_set_ncurr(gmatrix *g)
{
   if(g->mode == MODE_TRAIN)
      g->ncurr = g->ntrain[g->fold];
   else
      g->ncurr = g->ntest[g->fold];

   g->ncurr_recip = 1.0 / g->ncurr;
}

int gmatrix_set_fold(gmatrix *g, int fold)
{
   g->fold = fold;
   gmatrix_set_ncurr(g);
   if(!gmatrix_init_lp(g))
      return FAILURE;
   if(g->scalefile)
      return gmatrix_read_scaling(g, g->scalefile);
   return SUCCESS;
}

/* zero the lp and adjust the lp-functions */
void gmatrix_zero_model(gmatrix *g)
{
   int i, j, n = g->ncurr, p1 = g->p + 1;

   for(j = p1 - 1 ; j >= 0 ; --j)
   {
      g->beta[j] = 0;
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
      }
   }
}

int gmatrix_init_lp(gmatrix *g)
{
   if(g->lp)
      free(g->lp);
   CALLOCTEST(g->lp, g->ncurr, sizeof(double));
   if(g->model == MODEL_LOGISTIC) {
      if(g->lp_invlogit)
	 free(g->lp_invlogit);
      CALLOCTEST(g->lp_invlogit, g->ncurr, sizeof(double));
   } else if(g->model == MODEL_SQRHINGE) {
      if(g->ylp)
	 free(g->ylp);
      CALLOCTEST(g->ylp, g->ncurr, sizeof(double));
      if(g->ylp_neg)
	 free(g->ylp_neg);
      CALLOCTEST(g->ylp_neg, g->ncurr, sizeof(double));
   }
   return SUCCESS;
}
