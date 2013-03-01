/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include "gmatrix.h"

#define MAXLP 7

#define MAX_SHOW_NOTCONV 100

#define OPTIONS_CALLER_CD 1
#define OPTIONS_CALLER_UNIVARIABLE 2

#define OPTIONS_MULTIVAR_NEWTON 1
#define OPTIONS_MULTIVAR_LASSO 2

#ifndef PATH_MAX
#define PATH_MAX 8192
#endif

typedef double (*predict_pt)(double);
typedef double (*phi1)(double);
typedef double (*phi2)(double);
typedef double (*inv)(double);
typedef void (*step)(sample *s, gmatrix *g, int k,
   double *restrict d1_p, double *restrict d2_p);

typedef double (*predict)(double x);

typedef struct Opt {
   short caller;
   short mode;
   short model;
   short modeltype;
   int maxepochs;
   int maxiters;
   double lambda1;
   double lambda2;
   double gamma;
   double threshold;
   double l1minratio;
   double l1max;
   int nlambda1;
   double trunc;
   phi1 phi1_func;
   phi2 phi2_func;
   inv inv_func;
   short nofit;
   int n, p;
   short warmrestarts;
   char *filename;
   double *lambda1path;
   double lambda1max;
   double lambda1min;
   int verbose;
   int nfolds;
   long seed;
   int nzmax;
   int ntrain;
   int *trainf;
   char *lambda1pathfile;
   step step_func;
   short inmemory;
   char *scalefile;
   short yformat;
   predict predict_func;
   char **beta_files;
   char **beta_files_fold;
   int n_beta_files;
   char *predict_file;
   short encoded;
   char *folds_ind_file;
   char *numnz_file;
   double *zthresh;
   int nzthresh;
   double lambda2_univar;
   double lambda2_multivar;
   int do_multivar;
   int existing_univar;
   int do_thinning;
   int do_lasso_filter;
   int multivar;
   char *subset_file;
   char *famfilename;
   char *outdir;
   int scaley;
   int unscale_beta;
   int cortype;
   int corthresh;
   int phenoformat;
   long maxmem;
} Opt;

int cd_gmatrix(gmatrix *g,
      step step_func,
      const int maxepochs, const int maxiters, 
      const double lambda1, const double lambda2,
      const double gamma,
      const double trunc,
      int *numactiveK);

double get_lambda1max_gmatrix(gmatrix *g,
      phi1 phi1_func,
      phi2 phi2_func,
      inv inv_func,
      step step_func
);

int cvsplit(Opt *opt);
void opt_free(Opt *opt);
int opt_defaults(Opt *opt, short caller);
int opt_parse(int argc, char* argv[], Opt* opt);
int make_lambda1path(Opt *opt, gmatrix *g);
int run(Opt *opt, gmatrix *g);
void zero_model(gmatrix *g);

//double step_regular_linear(sample *s, gmatrix *g);
//double step_regular_logistic(sample *s, gmatrix *g);
//double step_regular_sqrhinge(sample *s, gmatrix *g);


