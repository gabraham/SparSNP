#include <stdio.h>

typedef struct {
   char* filename;
   FILE* file;
   double *y;
   int n;
   int p;
   int i;
   double *mean;
   double *sd;
} gmatrix;

typedef struct {
   double *x;
   double *x1; /* keep the old pointer so that we can reset x to it after
	        * incrementing x, to avoid copying data */
   double y;
   int p;
} sample;


void gmatrix_init(gmatrix *, char *, int, int);
void gmatrix_nextrow(gmatrix *, sample *);
void gmatrix_free(gmatrix *);
void gmatrix_reset(gmatrix *);
void sample_init(sample *, int);
void sample_free(sample *);

