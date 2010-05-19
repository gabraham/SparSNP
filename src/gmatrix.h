#include <stdio.h>
#include "common.h"

typedef struct {
   char* filename;
   FILE* file;
   dtype *y;
   int n;
   int p;
   int i;
   double *mean;
   double *sd;
} gmatrix;

typedef struct {
   dtype *x;
   dtype *x1; /* keep the old pointer so that we can reset x to it after
	        * incrementing x, to avoid copying data */
   dtype y;
   int p;
} sample;


void gmatrix_init(gmatrix *, char *, int, int);
void gmatrix_nextrow(gmatrix *, sample *);
void gmatrix_free(gmatrix *);
dtype gmatrix_next_y(gmatrix *);
void gmatrix_reset(gmatrix *);
void sample_init(sample *, int);
void sample_free(sample *);

