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
#include <math.h>
#include <unistd.h>

#include "common.h"
#include "link.h"
#include "coder.h"

/* categorical x inputs: 0, 1, 2, 3 */
#define NUM_X_LEVELS 4
#define X_LEVELS {0, 1, 2, 3}
#define X_LEVEL_NA 3
#define PLINK_X_LEVEL_NA 2

#define MISSING_PHENO "-9"
#define FAM_MAX_CHARS 100

#define YFORMAT01 1
#define YFORMAT11 2

#define HASH_SIZE 64

#define NA_ACTION_NONE 0
#define NA_ACTION_DELETE 1
#define NA_ACTION_ZERO 2
#define NA_ACTION_RANDOM 3
#define NA_ACTION_PROPORTIONAL 4

/* Size of cache itself, excluding the counters and mappings
 * Remember: there are g->folds caches, not just one, so total memory required
 * is CACHE_MAX_MEM * g->nfolds.
 * */
#define CACHE_MAX_MEM 134217728  /* 2^27=128MB */
/*#define CACHE_MAX_MEM   268435456*/ /* 2^27=128MB */
/*#define CACHE_MAX_MEM   357913941*/ /* 2^27=128MB */

#define CACHE_NOT_EXISTS -1

/* PLINK FAM (6 columns) or PHENO (2 FID/IID + the phenotype columns) */
#define PHENO_FORMAT_FAM 1
#define PHENO_FORMAT_PHENO 2

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
   int loss;
   int p;
   int K;
   int i;
   int j;
   double *mean, *sd;
   double *lookup, *lookup2;
   double *lp, *err, *ylp, *ylp_neg, *ylp_neg_y;
   double *ylp_neg_y_ylp, *lp_invlogit;
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
   dtype *tmp2;
   unsigned char *encbuf;
   short yformat;
   int model;
   int modeltype;
   short encoded;
   int nencb;
   /*void (*decode)(unsigned char *out,
	 const unsigned char *in,
	 const int n);*/
   char *folds_ind_file;
   int nfolds;
   int *folds;
   int fold;
   int mode;
   int nseek;
   double *beta_orig;
   int *numnz;
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
   mapping *map;
   int *folds_ind; /* xor of folds with mode == MODE_PREDICT */
   int unscale_beta;
   double *C;
   int cortype;
   int corthresh;
   int phenoformat;
   int verbose;
   double *proportions;
} gmatrix;

int sample_init(sample *);
int gmatrix_init(gmatrix *g, char *filename, int n, int p,
      char *scalefile, short yformat, int phenoformat,
      int model, int modeltype,
      short encoded, char *folds_ind_file,
      short mode, char *subsample_file,
      char *famfilename, int scaley, int unscale_beta,
      int cortype, int corthresh, int verbose);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *g, sample *sm, int skip, int na_action);
int gmatrix_mem_nextcol(gmatrix *g, sample *sm, int j, int na_action);
int gmatrix_disk_read_y(gmatrix *g);
int gmatrix_pheno_read_y(gmatrix *g);
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
int gmatrix_scale_y(gmatrix *g);
int gmatrix_trim_beta(gmatrix *g);
int gmatrix_make_fusion(gmatrix *g);
int gmatrix_load_subsets(gmatrix *g);
int gmatrix_plink_check_pheno(gmatrix *g);
int gmatrix_disk_nextcol_raw(gmatrix *g, sample *s, int j);
int gmatrix_init_proportions(gmatrix *g);
int rand_geno_proportional(gmatrix *g, int j);

void step_regular_linear(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p);
void step_regular_sqrhinge(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p);
void step_regular_logistic(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p);
int init_newton(gmatrix *g);
void updatelp(gmatrix *g, const double update,
      const double *restrict x, int j, int k);

int cache_get(cache *ca, int j, double **x);
int cache_init(cache *ca, int n, int p);
void cache_free(cache *ca);
void count_cases(gmatrix* g);

int rand_geno();
long gmatrix_lrand();
double gmatrix_drand();

