#include <stdio.h>
#include <math.h>
#include "common.h"

#define BUFSIZE 10

typedef struct tabulation {
   int p;
   int nbins;
   int **counts;
} tabulation;

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
   tabulation *tab;
} gmatrix;

int sample_init(sample *, int, short);
void sample_free(sample *);
int gmatrix_init(gmatrix *, char *, int, int, short, short);
int gmatrix_reset(gmatrix *);
void gmatrix_free(gmatrix *);
int gmatrix_disk_nextcol(gmatrix *, sample *);
int gmatrix_mem_nextcol(gmatrix *, sample *);
int gmatrix_disk_nextcol2(gmatrix *, sample *);
int gmatrix_load(gmatrix *g);

int tabulation_init(tabulation *t, int p, int nbins);
void tabulation_free(tabulation *t);
int tabulate(gmatrix *g, char *fileout, short doscale,
      int xlevels, int ylevels, short verbose);
int tabulation_write(tabulation *t, char* fileout);
int tabulation_read(tabulation *t, char* filein);

int gmatrix_tabulation_load(gmatrix *g);
int gmatrix_tabulation_nextcol(gmatrix *g, sample *s);

