#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"
#include "coder.h"
#include "ind.h"

static inline int hash(int key);

int sample_init(sample *s, int n)
{
   s->x = NULL;
   s->cached = FALSE;
   s->intercept = FALSE;
   return SUCCESS;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int model,
      short encoded, short binformat, char *folds_ind_file,
      short mode, loss_pt loss_pt_func)
{
   int i, j, p1;

   g->model = model;
   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;
   p1 = p + 1;
   g->yidx = 0;
   g->y = NULL;
   g->y_orig = NULL;
   g->lp = NULL;
   g->ylp = NULL;
   g->ylp_neg = NULL;
   g->lp_invlogit = NULL;
   g->loss_func = NULL;
   g->loss_pt_func = loss_pt_func;
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

   CALLOCTEST(g->beta_orig, p1, sizeof(double));

   CALLOCTEST(g->ncurr_j, p1, sizeof(int));
   CALLOCTEST(g->ncurr_recip_j, p1, sizeof(double));
   
   MALLOCTEST(g->ca, sizeof(cache));
   if(!cache_init(g->ca, p1)) /* +1 not for intercept, which isn't
				       stored in cache, but we need to maintain
				       indexing consistency */
      return FAILURE;

   if(filename)
      FOPENTEST(g->file, filename, "rb");

   if(!gmatrix_setup_folds(g))
      return FAILURE;

   MALLOCTEST(g->intercept, sizeof(double) * g->n);
   for(i = n - 1 ; i >= 0 ; --i)
      g->intercept[i] = 1.0;

   MALLOCTEST(g->tmp, sizeof(dtype) * g->n);

   if(encoded)
      MALLOCTEST(g->encbuf, sizeof(unsigned char) * g->nencb);

   if(g->binformat == BINFORMAT_PLINK)
      g->decode = &decode_plink;

   g->nextcol = gmatrix_disk_nextcol;

   if(!gmatrix_reset(g))
      return FAILURE;

   if(!gmatrix_disk_read_y(g))
      return FAILURE;

   CALLOCTEST(g->beta, p1, sizeof(double));
   CALLOCTEST(g->active, p1, sizeof(int));
   CALLOCTEST(g->ignore, p1, sizeof(int));

   if(g->scalefile && !gmatrix_read_scaling(g, scalefile))
      return FAILURE;

   for(j = g->p ; j >= 0 ; --j)
      g->active[j] = !g->ignore[j];

   gmatrix_set_ncurr(g);
   
   if(!gmatrix_split_y(g))
      return FAILURE;

   if(!gmatrix_init_lp(g))
      return FAILURE;

   MALLOCTEST(g->xtmp, sizeof(double) * g->n);

   return SUCCESS;
}

int gmatrix_setup_folds(gmatrix *g)
{
   if(g->folds_ind_file 
	 && !(g->nfolds = ind_getfolds(g->folds_ind_file)))
      return FAILURE;

   MALLOCTEST(g->folds, sizeof(int) * g->n * g->nfolds);
   if(g->nfolds > 1 && !ind_read(g->folds_ind_file,
	    g->folds, g->n, g->nfolds))
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
   FREENULL(g->ignore);
   FREENULL(g->tmp);
   FREENULL(g->intercept);
   FREENULL(g->lookup);
   FREENULL(g->lp);
   FREENULL(g->ylp);
   FREENULL(g->ylp_neg);
   FREENULL(g->lp_invlogit);
   FREENULL(g->beta);
   FREENULL(g->beta_orig);
   FREENULL(g->encbuf);
   FREENULL(g->folds);
   FREENULL(g->ntrain);
   FREENULL(g->ntest);
   FREENULL(g->ntrainrecip);
   FREENULL(g->ntestrecip);
   FREENULL(g->active);
   FREENULL(g->numnz);
   FREENULL(g->ncurr_j);
   FREENULL(g->ncurr_recip_j);

   if(g->ca)
   {
      cache_free(g->ca);
      free(g->ca);
   }
   g->ca = NULL;

}

/* y_orig stays in memory and never changes */
int gmatrix_disk_read_y(gmatrix *g)
{
   int i, n = g->n, n1 = g->n - 1;

   /* read all of y the first time we see it, then skip it */
   CALLOCTEST(g->y_orig, n, sizeof(double));

   /* The y vector may be byte packed */
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   g->decode(g->tmp, g->encbuf, g->nencb);
 
   if(g->yformat == YFORMAT01) {
      for(i = n1 ; i >= 0 ; --i)
	 g->y_orig[i] = (double)g->tmp[i];
   } else {  /*  -1/1  */
      for(i = n1 ; i >= 0 ; --i)
	 g->y_orig[i] = 2.0 * g->tmp[i] - 1.0;
   }

   return SUCCESS;
}

/* y is a copy of y_orig as needed for crossval */
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
	 if((g->mode == MODE_PREDICT) ^ g->folds[g->fold * n + i])
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
 * Read all the data into a preallocated row-major n by k matrix, which
 * variables are read is determined by the p+1 array ind
 */
int gmatrix_read_matrix(gmatrix *g, double *x, int *ind)
{
   int i, j,
       p1 = g->p + 1, 
       n = g->n;
   double prev = -1;
   sample sm;

   for(j = 0 ; j < p1 ; j++)
   {
      g->active[j] = TRUE;
      if(ind[j])
      {
	 if(!g->nextcol(g, &sm, j))
	    return FAILURE;
   
         for(i = 0 ; i < n ; i++)
	 {
	    x[i * p1 + j] = sm.x[i];
	    g->active[j] &= (sm.x[i] != prev);
	    prev = sm.x[i];
	 }
      }
   }

   return SUCCESS;
}

/*
 * reads one column of data from disk (all samples for one variable)
 */
int gmatrix_disk_nextcol(gmatrix *g, sample *s, int j)
{
   int i, l1, n1 = g->n - 1, n = g->n;
   int k = g->ncurr - 1;
   int f = g->fold * n;

   if(j == 0)
   {
      /* don't need to worry about cross-validation etc
       * because the intercept is all 1s */
      s->x = g->intercept;
      return SUCCESS;
   }

   /* column is in cache */
   /*if((s->x = cache_get(g->ca, j)))
      return SUCCESS;*/

   /* Get data from disk and unpack.
    * Skip y. Try to put data in cache. */
   FSEEKOTEST(g->file, j * g->nseek, SEEK_SET);
   FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
   g->decode(g->tmp, g->encbuf, g->nencb);

   /* Get the scaled versions of the genotypes */
   if(g->nfolds < 2) { /* without crossval */
      if(g->scalefile) { /* scale.c doesn't use scalefile so check */
         /* Get the scaled value instead of the original value */
         l1 = j * NUM_X_LEVELS;
         for(i = n1 ; i >= 0 ; --i)
	    g->xtmp[i] = g->lookup[l1 + g->tmp[i]];
      } else { /* no scaling */
	 for(i = n1 ; i >= 0 ; --i)
	    g->xtmp[i] = (double)g->tmp[i];
      }
   } else {
      /* with crossval */
      /* Get the scaled value instead of the original value */
      if(g->scalefile) {
	 l1 = j * NUM_X_LEVELS;
	 for(i = n1 ; i >= 0 ; --i)
	    /* different between train and test */
	    if((g->mode == MODE_PREDICT) ^ g->folds[f + i])
	       g->xtmp[k--] = g->lookup[l1 + g->tmp[i]];
      } else { /* no scaling */
	 for(i = n1 ; i >= 0 ; --i)
	    if((g->mode == MODE_PREDICT) ^ g->folds[f + i])
	       g->xtmp[k--] = (double)g->tmp[i];
      }
   }

   /*cache_put(g->ca, j, g->xtmp, g->ncurr);*/
   s->x = g->xtmp;

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

/* Populate lookup tables for SNP levels 0, 1, 2, 3
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
	    /* Scale unless missing obs. If missing, make the observed a zero
	     * so that it contributes nothing to the gradient later. */
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

/* Hashtable as cache. Each variable has an associated weight,
 * i.e. the number of requests for it, and putting a variable
 * only works if it has higher weight than the variable in
 * the same bucket (if they both hash to the same bucket).*/
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
      ca->buckets[i].n = 0;
      ca->buckets[i].value = NULL;
   }

   return SUCCESS;
}

void cache_free(cache *ca)
{
   int i;

   for(i = 0 ; i < ca->size ; i++)
      FREENULL(ca->buckets[i].value);

   FREENULL(ca->buckets);
   FREENULL(ca->weights);
   FREENULL(ca->keys);
   ca->active = 0;
}

int cache_put(cache *ca, int key, double *value, int n)
{
   int i;
   bucket *bk = NULL;
   int hval = hash(key);

   bk = &ca->buckets[hval];

   /* accept new value only if bucket is empty and if new weight is higher
    * than old weight */
   if(!bk->active || ca->weights[key] > ca->weights[bk->key])
   {
      ca->keys[bk->key] = FALSE;
      bk->key = key;
      bk->n = n;
      FREENULL(bk->value);
      MALLOCTEST(bk->value, sizeof(double) * n);
      for(i = n - 1 ; i >= 0 ; --i)
	 bk->value[i] = value[i];
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
   int hval = 0;
   
   if(!ca->keys[key])
      return NULL;

   hval = hash(key);
   ca->weights[key]++;

   bk = &ca->buckets[hval];
   if(bk && bk->active && bk->key == key)
      return bk->value;

   return NULL;
}

/* assumes 32 bit ints */
static inline int hash(int key)
{
   int s = 0;
   s += key & 1;
   s += key & 2;
   s += key & 4;
   s += key & 8;
   s += key & 16;
   s += key & 32;
   return s;
}

/* Sets the number of samples for actual use in CD */
void gmatrix_set_ncurr(gmatrix *g)
{
   g->ncurr = (g->mode == MODE_TRAIN) ? 
      g->ntrain[g->fold] : g->ntest[g->fold];

   g->ncurr_recip = 1.0 / g->ncurr;
}

int gmatrix_set_fold(gmatrix *g, int fold)
{
   g->fold = fold;
   gmatrix_set_ncurr(g);
   cache_free(g->ca);
   if(!cache_init(g->ca, g->p + 1))
      return FAILURE;
   if(!gmatrix_init_lp(g))
      return FAILURE;
   if(g->scalefile && !gmatrix_read_scaling(g, g->scalefile))
      return FAILURE;
   return gmatrix_split_y(g);
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
   FREENULL(g->lp);
   CALLOCTEST(g->lp, g->ncurr, sizeof(double));
   if(g->model == MODEL_LOGISTIC) {
      FREENULL(g->lp_invlogit);
      CALLOCTEST(g->lp_invlogit, g->ncurr, sizeof(double));
   } else if(g->model == MODEL_SQRHINGE) {
      FREENULL(g->ylp);
      CALLOCTEST(g->ylp, g->ncurr, sizeof(double));
      FREENULL(g->ylp_neg);
      CALLOCTEST(g->ylp_neg, g->ncurr, sizeof(double));
   }
   return SUCCESS;
}

