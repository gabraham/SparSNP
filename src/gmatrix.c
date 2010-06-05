#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

int sample_init(sample *s, int p)
{
   s->p = p;
   MALLOCTEST(s->x1, sizeof(dtype) * (p + 1))
   s->x = s->x1;
   return SUCCESS;
}

void sample_free(sample *s)
{
   free(s->x1);
}

int gmatrix_init(gmatrix *g, short inmemory, short pcor,
      char *filename, dtype **x, dtype *y, int n, int p)
{
   int i;

   if(filename && !inmemory)
      FOPENTEST(g->file, filename, "rb")

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
   CALLOCTEST(g->mean, p + 1, sizeof(double))
   MALLOCTEST(g->sd, sizeof(double) * (p + 1))

   for(i = 0 ; i < g->p + 1 ; i++)
      g->sd[i] = 1;

   g->nextrow = gmatrix_disk_nextrow;
   g->next_y = gmatrix_disk_next_y;

   if(inmemory)
   {
      g->nextrow = gmatrix_mem_nextrow;
      g->next_y = gmatrix_mem_next_y;
      
      if(filename != NULL && x == NULL)
	 if(!gmatrix_load(g))
	    return FAILURE;
      /*if(pcor)
      {
	 g->nextrow = gmatrix_mem_pcor_nextrow;
	 g->next_y = gmatrix_mem_pcor_next_y;
      }*/
   }

   return SUCCESS;
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
int gmatrix_disk_nextrow(gmatrix *g, sample *s)
{
   int i;
   dtype *tmp;
   
   MALLOCTEST(tmp, sizeof(intype) * (g->p + 1))

   if(g->i == g->n)
      if(!gmatrix_reset(g))
	 return FAILURE;

   /* reset to old pointer so we don't increment beyond
    * the allocated vector later */
   /*s->x = s->x1;*/

   FREADTEST(tmp, sizeof(intype), g->p + 1, g->file)

   s->y = tmp[0];
   s->x[0] = 1.0; /* intercept */
   for(i = 1 ; i < g->p + 1; i++)
      s->x[i] = ((dtype)tmp[i] - g->mean[i]) / g->sd[i];
   /*s->x++;*/
   g->i++;
   free(tmp);
   return SUCCESS;
}

int gmatrix_load(gmatrix *g)
{
   int i, j;
   FILE* fin;
   intype *tmp;
   
   MALLOCTEST(tmp, sizeof(intype) * (g->p + 1))
   MALLOCTEST(g->x, sizeof(dtype*) * g->n)
   MALLOCTEST(g->y, sizeof(dtype) * g->n)

   FOPENTEST(fin, g->filename, "rb")

   for(i = 0 ; i < g->n ; i++)
   {
      MALLOCTEST(g->x[i], sizeof(dtype) * (g->p + 1))

      FREADTEST(tmp, sizeof(intype), g->p + 1, fin)

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

   return SUCCESS;
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

int gmatrix_mem_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      if(!gmatrix_reset(g))
	 return FAILURE;

   s->x = g->x[g->i];
   s->y = g->y[g->i];
   g->i++;

   return SUCCESS;
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

int gmatrix_reset(gmatrix *g)
{
   if(!g->inmemory)
   {
      if(g->file)
      {
         fclose(g->file);
         g->file = NULL;
      }
      FOPENTEST(g->file, g->filename, "rb")
   }
   g->i = 0;
   return SUCCESS;
}

