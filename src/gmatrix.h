#include <stdio.h>
#include <math.h>
#include "common.h"

#define BUFSIZE 10

/* categorical x inputs: 0, 1, 2 */
#define NUM_X_LEVELS 3
#define X_LEVELS {0, 1, 2}

#define YFORMAT01 1
#define YFORMAT11 2

typedef struct sample {
   double *x;
   double *x2;
   short inmemory;
   short intercept;

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
   int p;
   int i;
   int j;
   double *mean, *sd;
   double *lookup, *lookup2;
   double *lp, *ylp, *ylp_neg, *lp_invlogit;
   double *beta;
   double *intercept;
   int *active;
   dtype *buffer;
   int bufsize;
   int bufidx;
   int yidx; /* only used for pcor */
   short inmemory;
   int (*nextcol)(struct gmatrix*, sample*);
   char *scalefile;
   dtype *tmp;
   short yformat;
   int model;
   short encoded;
} gmatrix;

int sample_init(sample *, int, short);
void sample_free(sample *);
int gmatrix_init(gmatrix *, char *, int, int, short, char*, short,
      int, short);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *);
int gmatrix_mem_nextcol(gmatrix *, sample *);
int gmatrix_disk_nextcol2(gmatrix *, sample *);
int gmatrix_load(gmatrix *g);
int gmatrix_disk_skipcol(gmatrix *g);
int gmatrix_read_scaling(gmatrix *g, char *file_scale);

