/*
 * Copyright (c) 2011, National ICT Australia (NICTA)
 * All rights reserved.
 */

#include <stdio.h>
#include <math.h>

#include "common.h"
#include "loss.h"

/* categorical x inputs: 0, 1, 2, 3 */
#define NUM_X_LEVELS 4
#define X_LEVELS {0, 1, 2, 3}
#define X_LEVEL_NA 3
#define PLINK_X_LEVEL_NA 2

#define YFORMAT01 1
#define YFORMAT11 2

#define HASH_SIZE 64

#define NA_ACTION_DELETE 1
#define NA_ACTION_ZERO 2
#define NA_ACTION_RANDOM 3

/* Size of cache itself, excluding the counters and mappings
 * Remember: there are g->folds caches, not just one, so total memory required
 * is CACHE_MAX_MEM * g->nfolds.
 * */
/*#define CACHE_MAX_MEM 134217728 */ /* 2^27=128MB */
#define CACHE_MAX_MEM 268435456 /* 2^27=128MB */

#define CACHE_NOT_EXISTS -1

typedef struct cache {
   int nbins;
   int n;
   int *mapping;
   int *revmapping;
   int *counter;
   int lastfree;
   double *x;
   double *tmp;
} cache;

typedef struct sample {
   int n;
   int ngood;
   int *good;
   double *x;
   double *y;
   double *x2;
   short intercept;
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
   double *ytmp;
   int *good;
   int n;
   int *ntrain;
   int *ntest;
   double *ntrainrecip;
   double *ntestrecip;
   int ncurr;
   double ncurr_recip;
   int *ncurr_j;
   double *ncurr_recip_j;
   int p;
   int i;
   int j;
   double *mean, *sd;
   double *lookup, *lookup2;
   double loss;
   double *lp, *ylp, *ylp_neg, *ylp_neg_y, *ylp_neg_y_ylp, *lp_invlogit;
   double *newtonstep;
   double *beta;
   int *active;
   double *intercept;
   int *ignore;
   dtype *buffer;
   int bufsize;
   int bufidx;
   int yidx; /* only used for pcor */
   int (*nextcol)(struct gmatrix*, sample*, int skip, int na_action);
   char *scalefile;
   dtype *tmp;
   unsigned char *encbuf;
   short yformat;
   int model;
   int modeltype;
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
   int mode;
   int nseek;
   double *beta_orig;
   int *numnz;
   loss loss_func;
   loss_pt loss_pt_func;
   int *ngood;
   double *x;
   double *xthinned;
   char *subset_file;
   int *subsets;
   int nsubsets;
   int offset;
   char *famfilename;
   cache *xcaches;
   int ncases;
} gmatrix;

int sample_init(sample *);
int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int model, int modeltype,
      short encoded, short binformat, char *folds_ind_file,
      short mode, loss_pt, char *subsample_file,
      char *famfilename);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *g, sample *sm, int skip, int na_action);
int gmatrix_mem_nextcol(gmatrix *g, sample *sm, int j, int na_action);
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
int gmatrix_fam_read_y(gmatrix *g);
int gmatrix_read_matrix(gmatrix *g, int *ind, int m);

int gmatrix_load_subsets(gmatrix *g);
int gmatrix_plink_check_pheno(gmatrix *g);

double step_regular_linear(sample *s, gmatrix *g);
double step_regular_sqrhinge(sample *s, gmatrix *g);
double step_regular_logistic(sample *s, gmatrix *g);
int init_newton(gmatrix *g);
void updatelp(gmatrix *g, const double update,
      const double *restrict x, int j);

int cache_get(cache *ca, int j, double **x);
int cache_init(cache *ca, int n, int p);
void cache_free(cache *ca);
void count_cases(gmatrix* g);

int rand_geno();

