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
   g->sd = malloc(sizeof(double) * (p + 1));

   for(i = 1 ; i < g->p + 1 ; i++)
      g->sd[i] = 1;

   g->nextrow = gmatrix_disk_nextrow;
   g->next_y = gmatrix_disk_next_y;

   if(inmemory)
   {
      g->nextrow = gmatrix_mem_nextrow;
      g->next_y = gmatrix_mem_next_y;
      
      if(filename != NULL && x == NULL)
	 gmatrix_load(g);
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
   int i;
   dtype *tmp = malloc(sizeof(intype) * (g->p + 1));
   if(g->i == g->n)
      gmatrix_reset(g);

   /* reset to old pointer so we don't increment beyond
    * the allocated vector later */
   /*s->x = s->x1;*/

   fread(tmp, sizeof(intype), g->p + 1, g->file);
   s->y = tmp[0];
   s->x[0] = 1.0; /* intercept */
   for(i = 1 ; i < g->p + 1; i++)
      s->x[i] = ((dtype)tmp[i] - g->mean[i]) / g->sd[i];
   /*s->x++;*/
   g->i++;
   free(tmp);
}

void gmatrix_load(gmatrix *g)
{
   int i, j;
   FILE* fin = fopen(g->filename, "rb");
   intype *tmp = malloc(sizeof(intype) * (g->p + 1));
   g->x = malloc(sizeof(dtype*) * g->n);
   g->y = malloc(sizeof(dtype) * g->n);

   for(i = 0 ; i < g->n ; i++)
   {
      g->x[i] = malloc(sizeof(dtype) * (g->p + 1));
      /*fread(g->x[i], sizeof(intype), g->p + 1, fin); */
      fread(tmp, sizeof(intype), g->p + 1, fin);

     /* g->y[i] = g->x[i][0]; */

      /* intercept */
      /* g->x[i][0] = 1.0; */

      g->y[i] = (dtype)tmp[0];

      g->x[i][0] = 1.0;
      for(j = 1 ; j < g->p + 1 ; j++)
	 g->x[i][j] = (dtype)tmp[j];
   }

   fclose(fin);
   free(tmp);
}

/* Only applicable for memory-based matrices
 */
void gmatrix_scale(gmatrix *g)
{
   int i, j;
   if(g->inmemory && g->mean && g->sd)
   {
      for(i = 0 ; i < g->n ; i++)
	 for(j = 1 ; j < g->p + 1; j++)
	    g->x[i][j] = (g->x[i][j] - g->mean[j]) / g->sd[j];
   }
}

void gmatrix_mem_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      gmatrix_reset(g);

   s->x = g->x[g->i];
   s->y = g->y[g->i];
   g->i++;
}

/*void gmatrix_mem_nextcol(gmatrix *g, var *v)
{
}*/

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

   fread(&y, sizeof(intype), 1, g->file);
   g->i++;

   /* ignore the x vector*/
   fseek(g->file, sizeof(intype) * (g->p + 1), SEEK_CUR);

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

