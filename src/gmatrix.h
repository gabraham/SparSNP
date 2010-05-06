#include <stdio.h>

typedef struct {
   char* filename;
   FILE* file;
   int *y;
   int n;
   int p;
   int i;
} gmatrix;

typedef struct {
   double *x;
   int y;
   int p;
} sample;


void gmatrix_init(gmatrix *, char *, int, int, int*);
void gmatrix_nextrow(gmatrix *, sample *, int);
void gmatrix_free(gmatrix *);
void sample_init(sample *, int);
void sample_free(sample *);

