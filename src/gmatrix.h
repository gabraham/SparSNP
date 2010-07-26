#include <stdio.h>
#include <math.h>
#include "common.h"

#define BUFSIZE 10

typedef struct sample {
   dtype *x;
   short inmemory;
} sample;

typedef struct gmatrix {
   char* filename;
   FILE* file;
   dtype *y;
   dtype **x;
   int n;
   int p;
   int i;
   int j;
   double *mean;
   double *sd;
   dtype *buffer;
   int bufsize;
   int bufidx;
   int yidx; /* only used for pcor */
   short inmemory;
   int (*nextcol)(struct gmatrix*, sample*);
} gmatrix;

int sample_init(sample *, int, short);
void sample_free(sample *);
int gmatrix_init(gmatrix *, char *, int, int, short);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *);
int gmatrix_mem_nextcol(gmatrix *, sample *);
int gmatrix_disk_nextcol2(gmatrix *, sample *);
int gmatrix_load(gmatrix *g);

