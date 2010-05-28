#include <stdio.h>
#include "common.h"

/*typedef struct gmatrix_struct* gmatrix_ptr;*/

typedef struct sample {
   dtype *x;
   dtype *x1; /* keep the old pointer so that we can reset x to it after
	       * incrementing x, to avoid copying data */
   dtype y;
   int p;
} sample;


typedef struct gmatrix {
   char* filename;
   FILE* file;
   dtype *y;
   int n;
   int p;
   int i;
   double *mean;
   double *sd;
   short inmemory;
   short pcor;
   dtype **x;
   int skip;

   /*void (*init)(gmatrix*, short, short, char*, int, int);
   void (*reset)(gmatrix*);*/
   void (*nextrow)(struct gmatrix*, sample*);
   dtype (*next_y)(struct gmatrix*);
} gmatrix;


void gmatrix_init(gmatrix *, short, short, char *, dtype **, dtype *, int, int);
void gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);

void gmatrix_disk_nextrow(gmatrix *, sample *);
dtype gmatrix_disk_next_y(gmatrix *);

void gmatrix_mem_nextrow(gmatrix *, sample *);
dtype gmatrix_mem_next_y(gmatrix *);

void sample_init(sample *, int);
void sample_free(sample *);

