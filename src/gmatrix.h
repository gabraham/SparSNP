#include <stdio.h>
#include <math.h>
#include "common.h"

typedef struct sample {
   dtype *x;
   dtype y;
   int p;
   short inmemory;
} sample;

/*
 * When the data is in memory, it is *always* stored in row major ordering, i.e.
 * the x[0] is a pointer of length p+1 referring to the first sample across
 * all variables. This is regardless of whether we later read the data in rows
 * or columns.
 */
typedef struct gmatrix {
   char* filename;
   FILE* file;
   dtype *y;
   int n;
   int p;
   int i;
   int j;
   double *mean;
   double *sd;
   short inmemory;
   short pcor;
   short rowmajor;
   dtype **x;
   int skip;

   int (*nextrow)(struct gmatrix*, sample*);
   dtype (*next_y)(struct gmatrix*);
   int (*nextcol)(struct gmatrix*, sample*);
} gmatrix;


int gmatrix_init(gmatrix *, short, short, short,
      char *, dtype **, dtype *, int, int);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);

int gmatrix_disk_nextrow(gmatrix *, sample *);
int gmatrix_disk_nextcol(gmatrix *, sample *);
dtype gmatrix_disk_next_y(gmatrix *);

int gmatrix_mem_nextrow(gmatrix *, sample *);
int gmatrix_mem_nextcol(gmatrix *, sample *);
dtype gmatrix_mem_next_y(gmatrix *);

int gmatrix_load_rowmajor(gmatrix *);
int gmatrix_load_colmajor(gmatrix *);
int gmatrix_load_pcor(gmatrix *);
int gmatrix_scale(gmatrix *);

int sample_init(sample *, short, int);
void sample_free(sample *);

int gmatrix_mem_nextcol(gmatrix *, sample *);

