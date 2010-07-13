#include <stdio.h>
#include <math.h>
#include "common.h"

typedef struct sample {
   dtype *x;
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
} gmatrix;

int sample_init(sample *, int);
void sample_free(sample *);
int gmatrix_init(gmatrix *, char *, int, int);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *);

