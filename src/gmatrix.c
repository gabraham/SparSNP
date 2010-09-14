#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"
#include "coder.h"
#include "ind.h"

int sample_init(sample *s, int n, short inmemory)
{
   s->inmemory = inmemory;
   s->x = NULL;
   /*s->x2 = NULL;*/
   s->intercept = FALSE;
   return SUCCESS;
}

void sample_free(sample *s)
{
   if(s->intercept)
      return;
   
   if(s->x)
      free(s->x);
   s->x = NULL;

   /*if(s->x2)
      free(s->x2);*/
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
   g->lookup2 = NULL;
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
   g->mode = mode;

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
   CALLOCTEST(g->lp, g->n, sizeof(double));
   CALLOCTEST(g->active, g->p + 1, sizeof(int));
   CALLOCTEST(g->ignore, g->p + 1, sizeof(int));

   if(scalefile && !gmatrix_read_scaling(g, scalefile))
      return FAILURE;

   for(j = g->p ; j >= 0 ; --j)
      g->active[j] = !g->ignore[j];

   if(g->model == MODEL_LOGISTIC) {
      CALLOCTEST(g->lp_invlogit, g->n, sizeof(double));
   } else if(g->model == MODEL_SQRHINGE) {
      CALLOCTEST(g->ylp, g->n, sizeof(double));
      CALLOCTEST(g->ylp_neg, g->n, sizeof(double));
   }

   return SUCCESS;
}

int gmatrix_setup_folds(gmatrix *g)
{
   MALLOCTEST(g->folds, sizeof(int) * g->n * g->nfolds);
   if(g->folds_ind_file && g->nfolds > 1 && 
	 !read_ind(g->folds_ind_file, g->folds, g->n, g->nfolds))
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

   if(g->lookup2)
      free(g->lookup2);
   g->lookup2 = NULL;

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

/* big ugly function
 */
int gmatrix_disk_nextcol(gmatrix *g, sample *s)
{
   int i, l1;/*, l2;*/
   int n = g->n, n1 = g->n - 1;
   
   if(g->j == g->p + 1 && !gmatrix_reset(g))
      return FAILURE;

   /* either read the y vector, or skip rows and/or headers */
   if(g->j == 0)
   {
      s->intercept = TRUE;

      if(g->binformat == BINFORMAT_BIN) {
	 /* read y the first time we see it */
      	 if(!g->y) {
      	    MALLOCTEST(g->y, sizeof(double) * n);

      	    /* The y vector may be byte packed */
      	    if(g->encoded) {
      	       FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
      	       g->decode(g->tmp, g->encbuf, g->nencb);
      	    } else {
      	       FREADTEST(g->tmp, sizeof(dtype), n, g->file);
      	    }

	    /* represent y as 0/1 or -1/1 */
      	    if(g->yformat == YFORMAT01) {
      	       for(i = n1 ; i >= 0 ; --i)
      	          g->y[i] = (double)g->tmp[i];
	    } else { /*  -1/1  */
      	       for(i = n1 ; i >= 0 ; --i)
      	          g->y[i] = 2.0 * g->tmp[i] - 1.0;
	    }

      	 } else if(g->encoded) { /* don't read y again */
	    FSEEKOTEST(g->file, sizeof(dtype) * g->nencb, SEEK_CUR);
      	 } else {
	    FSEEKOTEST(g->file, sizeof(dtype) * n, SEEK_CUR);
      	 }
      } else { /* skip plink headers */
	 FSEEKOTEST(g->file,
	       sizeof(dtype) * PLINK_HEADER_SIZE, SEEK_CUR);
      }
      
      /* No need to copy values, just assign the intercept vector,
       * but first, free any old data we have from previous
       * iterations*/
      if(s->x)
	 free(s->x);
      /*if(s->x2)
	 free(s->x2);*/
      s->x = g->intercept;
      /*s->x2 = g->intercept;*/
      
      g->j++;
      return SUCCESS;
   }

   /* this isn't an intercept, so we need to copy actual values */
   s->intercept = FALSE;
   if(g->j == 1) {
      MALLOCTEST(s->x, sizeof(double) * n);
   }

   /* Get the scaled versions of the genotypes */
   if(g->scalefile)
   {
      s->intercept = FALSE;

      /* check cache */
      if(g->ca->keys[g->j]) {
	 if(!(s->x = cache_get(g->ca, g->j)))
	 {
	    fprintf(stderr, "j=%d in keys but not in cache\n", g->j);
	    return FAILURE;
	 }
	 printf("cache hit: %d\n", g->j);
      } else {
	 printf("cache miss: %d\n", g->j);
	 if(g->encoded) {
      	    FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
      	    g->decode(g->tmp, g->encbuf, g->nencb);
      	 } else {
      	    FREADTEST(g->tmp, sizeof(dtype), n, g->file);
      	 }

      	 /* Get the scaled value instead of the original value */
      	 l1 = g->j * NUM_X_LEVELS;
      	 for(i = n1 ; i >= 0 ; --i)
      	    s->x[i] = g->lookup[l1 + g->tmp[i]];

	 cache_put(g->ca, g->j, s->x);
      }
      
      g->j++;
      return SUCCESS;
   }

   if(g->encoded) {
      FREADTEST(g->encbuf, sizeof(dtype), g->nencb, g->file);
      g->decode(g->tmp, g->encbuf, g->nencb);
   } else {
      FREADTEST(g->tmp, sizeof(dtype), n, g->file);
   }

   for(i = n1 ; i >= 0 ; --i)
      s->x[i] = (double)g->tmp[i];

   g->j++;
   return SUCCESS;
}

int gmatrix_mem_nextcol(gmatrix *g, sample *s)
{
   if(g->j == g->p + 1)
   {
      g->i = g->j = 0;
   }
      
   s->x = g->x[g->j];
   g->j++;

   return SUCCESS;
}

/* Expects the binary data column to be x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextcol_no_y(gmatrix *g, sample *s)
{
   /*if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   if(g->j == g->yidx) {
      FREADTEST(g->y, sizeof(dtype), g->n, g->file)
   } else {
      FREADTEST(s->x, sizeof(dtype), g->n, g->file)
   }

   g->j++;*/

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
      if(!gmatrix_disk_nextcol(g, &s))
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
   
   FOPENTEST(g->file, g->filename, "rb")

   g->i = g->j = 0;

   return SUCCESS;
}

/* Populate lookup tables for SNP levels 0, 1, 2
 */
int gmatrix_read_scaling(gmatrix *g, char *file_scale)
{
   int j, k, l1, l2, p1 = g->p + 1;
   FILE *in;

   CALLOCTEST(g->lookup, NUM_X_LEVELS * p1, sizeof(double));
   CALLOCTEST(g->lookup2, NUM_X_LEVELS * p1, sizeof(double));
   MALLOCTEST(g->mean, sizeof(double) * p1);
   MALLOCTEST(g->sd, sizeof(double) * p1);

   FOPENTEST(in, file_scale, "rb")

   FREADTEST(g->mean, sizeof(double), p1, in);
   FREADTEST(g->sd, sizeof(double), p1, in);

   /* intercept */
   for(k = 0 ; k < NUM_X_LEVELS ; k++)
   {
      g->lookup[k] = 1;
      g->lookup2[k] = 1;
   }
   
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
	    g->lookup2[l2] = g->lookup[l2] * g->lookup[l2];
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
   int i, j;

   for(i = 0 ; i < nfolds ; i++)
   {
      ntrain[i] = 0;
      for(j = 0 ; j < n ; j++)
	 ntrain[i] += (folds[n * i + j] > 0);
      ntest[i] = n - ntrain[i];
      ntrainrecip[i] = 0;
      if(ntrain[i] > 0)
	 ntrainrecip[i] = 1.0 / ntrain[i];
      ntestrecip[i] = 0;
      if(ntest[i] > 0)
	 ntestrecip[i] = 1.0 / ntest[i];
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
  
   printf("key: %d\n", ca->weights[key]);
   printf("bk->key: %d\n", ca->weights[bk->key]);
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
   }

   return SUCCESS;
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

