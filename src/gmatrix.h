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
   short intercept;
   short cached;
   int nbins;
   int *counts;
   double *values;
} sample;

typedef struct gmatrix {
   char* filename;
   FILE* file;
   double *y_orig;
   double *y;
   double *xtmp;
   int n;
   int *ntrain;
   int *ntest;
   double *ntrainrecip;
   double *ntestrecip;
   int ncurr;
   double ncurr_recip;
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
   int nseek;
   double *beta_orig;
} gmatrix;

int sample_init(sample *, int);
void sample_free(sample *);
int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int model,
      short encoded, short binformat, char *folds_ind_file,
      short mode);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *, int skip);
int gmatrix_disk_read_y(gmatrix *g);
int gmatrix_disk_skipcol(gmatrix *g);
int gmatrix_disk_skipcol(gmatrix *g);
int gmatrix_read_scaling(gmatrix *g, char *file_scale);
void count_fold_samples(int *ntrain, int *ntest, double *ntrainrecip,
      double *ntestrecip, int *folds, int nfolds, int n);
int gmatrix_setup_folds(gmatrix *g);
void gmatrix_set_ncurr(gmatrix *g);
int gmatrix_set_fold(gmatrix *g, int fold);
void gmatrix_zero_model(gmatrix *g);
int gmatrix_init_lp(gmatrix *g);
int gmatrix_split_y(gmatrix *g);
int gmatrix_disk_read_y(gmatrix *g);

int cache_init(cache *ht, int nkeys);
void cache_free(cache *ht);
int cache_put(cache *ht, int key, double *value);
double* cache_get(cache *ht, int key);
int hash(int key);

