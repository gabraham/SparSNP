#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

void sample_init(sample *s, int p)
{
   s->p = p;
   s->x1 = malloc(sizeof(double) * (p + 1));
   s->x = s->x1;
}

void sample_free(sample *s)
{
   free(s->x1);
}

void gmatrix_init(gmatrix *g, char *filename, int n, int p)
{
   int i;

   g->file = fopen(filename, "rb");
   g->filename = filename;
   g->i = 0;
   g->n = n;
   g->p = p;
   g->mean = calloc(p, sizeof(double));
   g->sd = malloc(sizeof(double) * p);

   for(i = 0 ; i < g->p ; i++)
      g->sd[i] = 1;
}

void gmatrix_free(gmatrix *g)
{
   /*free(g->y);*/
   fclose(g->file);
   /*free(g->filename); */
   free(g->mean);
   free(g->sd);
}

/* Expects the binary data row to be y, x_1, x_2, x_3, ..., x_p */
void gmatrix_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   /* reset to old pointer so we don't increment beyond
    * the allocated vector later
    */
   s->x = s->x1;
   fread(s->x, sizeof(double), g->p + 1, g->file);
   s->y = s->x[0];
   s->x++;
   g->i++;
}

dtype gmatrix_next_y(gmatrix *g)
{
   dtype y;
   if(g->i == g->n)
      gmatrix_reset(g);

   fread(&y, sizeof(dtype), 1, g->file);
   g->i++;

   /* ignore the x vector*/
   fseek(g->file, sizeof(dtype) * (g->p), SEEK_CUR);
   return y;
}

void gmatrix_reset(gmatrix *g)
{
   fclose(g->file);
   g->file = fopen(g->filename, "rb");
   g->i = 0;
}

