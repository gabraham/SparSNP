#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

void sample_init(sample *s, int p)
{
   s->p = p;
   s->x = malloc(sizeof(double) * p);
}

void sample_free(sample *s)
{
   free(s->x);
}

void gmatrix_init(gmatrix *g, char *filename, int n, int p, int *y)
{
   g->file = fopen(filename, "rb");
   g->filename = filename;
   g->i = 0;
   g->n = n;
   g->p = p;
   g->y = y;
}

void gmatrix_free(gmatrix *g)
{
   /*free(g->y);*/
   fclose(g->file);
   /*free(g->filename); */
}

void gmatrix_nextrow(gmatrix *g, sample *s, int loop)
{
   if(g->i < g->n - 1)
   {
      fread(s->x, sizeof(double), g->p, g->file);
      s->y = g->y[g->i];
      g->i ++;
   }
   else if(loop)
   {
      gmatrix_reset(g);
      fread(s->x, sizeof(double), g->p, g->file);
      s->y = g->y[g->i];
   }
}

void gmatrix_reset(gmatrix *g)
{
   fclose(g->file);
   g->file = fopen(g->filename, "rb");
   g->i = 0;
}

