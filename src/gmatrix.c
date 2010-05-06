#include <stdio.h>

#include "gmatrix.h"

void gmatrix_init(gmatrix *g, char *filename, int n, int p, int *y)
{
   g->file = fopen(filename, "rb");
   g->i = 1;
   g->n = n;
   g->p = p;
   g->y = y;
}

void gmatrix_nextrow(gmatrix *g, double *x)
{
   if(g->i < g->n)
   {
      fread(x, sizeof(double), g->p, g->file);
      g->i ++;
   }
   else
   {
      fclose(g->file);
      g->file = fopen(g->filename, "rb");
      g->i = 1;
   }
}

