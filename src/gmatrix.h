#include <stdio.h>
#include <math.h>

#include "common.h"

#define BUFSIZE 10

/* categorical x inputs: 0, 1, 2 */
#define NUM_X_LEVELS 3
#define X_LEVELS {0, 1, 2}

#define YFORMAT01 1
#define YFORMAT11 2

#define HASH_SIZE 16

typedef struct bucket {
   int key;
   double *value;
   unsigned int active;
   struct bucket *next;
} bucket;

typedef struct cache {
   unsigned int size;
   unsigned int active;
   bucket *buckets;
   int *keys;
   int nkeys;
   unsigned int *weights;
} cache;

typedef struct sample {
   double *x;
   double *x2;
   short inmemory;
   short intercept;
   short cached;
   int nbins;
   int *counts;
   double *values;
} sample;

typedef struct gmatrix {
   char* filename;
   FILE* file;
   double *y;
   double **x;
   int n;
   int *ntrain;
   int *ntest;
   double *ntrainrecip;
   double *ntestrecip;
   int p;
   int i;
   int j;
   double *mean, *sd;
   double *lookup, *lookup2;
   double *lp, *ylp, *ylp_neg, *lp_invlogit;
   double *beta;
   int *active;
   double *intercept;
   int *ignore;
   dtype *buffer;
   int bufsize;
   int bufidx;
   int yidx; /* only used for pcor */
   short inmemory;
   int (*nextcol)(struct gmatrix*, sample*, int skip);
   char *scalefile;
   dtype *tmp;
   unsigned char *encbuf;
   short yformat;
   int model;
   short encoded;
   int nencb;
   short binformat;
   void (*decode)(unsigned char *out,
	 const unsigned char *in,
	 const int n);
   char *folds_ind_file;
   int nfolds;
   int *folds;
   int fold;
   short mode;
   cache *ca;
} gmatrix;

int sample_init(sample *, int, short);
void sample_free(sample *);
int gmatrix_init(gmatrix *, char *, int, int, short, char*, short,
      int, short, short, char *folds_ind_file, int nfolds, short);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *, int skip);
int gmatrix_disk_read_y(gmatrix *g);
int gmatrix_disk_skipcol(gmatrix *g);
int gmatrix_mem_nextcol(gmatrix *, sample *, int skip);
int gmatrix_load(gmatrix *g);
int gmatrix_disk_skipcol(gmatrix *g);
int gmatrix_read_scaling(gmatrix *g, char *file_scale);
void count_fold_samples(int *ntrain, int *ntest, double *ntrainrecip,
      double *ntestrecip, int *folds, int nfolds, int n);
int gmatrix_setup_folds(gmatrix *g);

int cache_init(cache *ht, int nkeys);
void cache_free(cache *ht);
int cache_put(cache *ht, int key, double *value);
double* cache_get(cache *ht, int key);
int hash(int key);

