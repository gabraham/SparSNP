#include <stdio.h>
#include <math.h>
#include "common.h"

typedef struct sample {
   dtype *x;
   dtype y;
   int p;
   short inmemory;
} sample;

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
   dtype **x;
   int skip;

   int (*nextrow)(struct gmatrix*, sample*);
   dtype (*next_y)(struct gmatrix*);
   int (*nextcol)(struct gmatrix*, sample*);
} gmatrix;


int gmatrix_init(gmatrix *, short, short, char *, dtype **, dtype *, int, int);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);

int gmatrix_disk_nextrow(gmatrix *, sample *);
dtype gmatrix_disk_next_y(gmatrix *);

int gmatrix_mem_nextrow(gmatrix *, sample *);
dtype gmatrix_mem_next_y(gmatrix *);

int gmatrix_load(gmatrix *);
int gmatrix_load_pcor(gmatrix *);
int gmatrix_scale(gmatrix *);

int sample_init(sample *, int);
void sample_free(sample *);

int gmatrix_mem_nextcol(gmatrix *, sample *);

