#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

void sample_init(sample *s, int p)
{
   s->p = p;
   s->x1 = malloc(sizeof(dtype) * (p + 1));
   s->x = s->x1;
}

void sample_free(sample *s)
{
   free(s->x1);
}

void gmatrix_init(gmatrix *g, short inmemory, short pcor,
      char *filename, dtype **x, dtype *y, int n, int p)
{
   int i;

   if(filename && !inmemory)
      g->file = fopen(filename, "rb");

   g->filename = filename;
   g->i = 0;
   g->n = n;
   g->p = p;
   g->inmemory = inmemory;
   g->pcor = pcor;
   g->x = x;
   g->y = y;
   g->skip = -1;


   /* TODO: this code isn't needed for discrete inputs, but sgd_gmatrix will
    * need to be fixed too  */
   g->mean = calloc(p, sizeof(double));
   g->sd = malloc(sizeof(double) * p);

   for(i = 0 ; i < g->p ; i++)
      g->sd[i] = 1;

   g->nextrow = gmatrix_disk_nextrow;
   g->next_y = gmatrix_disk_next_y;

   if(inmemory)
   {
      g->nextrow = gmatrix_mem_nextrow;
      g->next_y = gmatrix_mem_next_y;
      
      if(filename != NULL && x == NULL)
	 gmatrix_load(g, filename, n, p);
      /*if(pcor)
      {
	 g->nextrow = gmatrix_mem_pcor_nextrow;
	 g->next_y = gmatrix_mem_pcor_next_y;
      }*/
   }
}

void gmatrix_free(gmatrix *g)
{
   if(!g->inmemory && g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }
   free(g->mean);
   free(g->sd);
   g->mean = g->sd = NULL;
}

/* Expects the binary data row to be y, x_1, x_2, x_3, ..., x_p */
void gmatrix_disk_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   /* reset to old pointer so we don't increment beyond
    * the allocated vector later */
   s->x = s->x1;

   fread(s->x, sizeof(dtype),  g->p + 1, g->file);
   s->y = s->x[0];
   s->x++;
   g->i++;
}

void gmatrix_load(gmatrix *g, char *filename, int n, int p)
{
   int i;
   FILE* fin = fopen(filename, "rb");
   g->x = malloc(sizeof(dtype*) * n);
   g->y = malloc(sizeof(dtype) * n);

   for(i = 0 ; i < n ; i++)
   {
      g->x[i] = malloc(sizeof(dtype) * (p + 1));
      fread(g->x[i], sizeof(dtype),  g->p + 1, fin);
      g->y[i] = g->x[i][0];
      g->x[i]++;
   }

   fclose(fin);
}

void gmatrix_mem_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   s->x = g->x[g->i];
   s->y = g->y[g->i];
   g->i++;
}

/*void gmatrix_mem_pcor_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   if(!g->skip)
   {
      fprintf(stderr,
	 "skip undefined in gmatrix_mem_pcor_nextrow(), aborting");
      return FAILURE;
   }

   s->x = g->x[g->i];
   s->y = g->x[g->i];
   g->i++;
}*/

dtype gmatrix_disk_next_y(gmatrix *g)
{
   dtype y;
   if(g->i == g->n)
      gmatrix_reset(g);

   fread(&y, sizeof(dtype), 1, g->file);
   g->i++;

   /* ignore the x vector*/
   fseek(g->file, sizeof(dtype) * g->p, SEEK_CUR);

   return y;
}

dtype gmatrix_mem_next_y(gmatrix *g)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   return g->y[g->i++];
}


void gmatrix_reset(gmatrix *g)
{
   if(!g->inmemory)
   {
      if(g->file)
      {
         fclose(g->file);
         g->file = NULL;
      }
      g->file = fopen(g->filename, "rb");
   }
   g->i = 0;
}

