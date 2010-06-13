#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

int sample_init(sample *s, int p)
{
   s->p = p;
   /*MALLOCTEST(s->x, sizeof(dtype) * (p + 1))*/
   return SUCCESS;
}

void sample_free(sample *s)
{
   /* if data is loaded into memory, s->x is just a pointer to that data
    * and we don't want to free it */
   if(s->x && !s->inmemory)
      free(s->x);
}

int gmatrix_init(gmatrix *g, short inmemory, short pcor,
      char *filename, dtype **x, dtype *y, int n, int p)
{
   int i;

   if(filename && !inmemory)
      FOPENTEST(g->file, filename, "rb")

   g->filename = filename;
   g->i = 0;
   g->j = 0;
   g->n = n;
   g->p = p;
   g->inmemory = inmemory;
   g->pcor = pcor;
   g->x = x;
   g->y = y;
   /*g->skip = -1;*/


   /* TODO: this code isn't needed for discrete inputs, but sgd_gmatrix will
    * need to be fixed too  */
   CALLOCTEST(g->mean, p + 1, sizeof(double))
   MALLOCTEST(g->sd, sizeof(double) * (p + 1))

   for(i = 0 ; i < g->p + 1 ; i++)
      g->sd[i] = 1;

   g->nextrow = gmatrix_disk_nextrow;
   g->next_y = gmatrix_disk_next_y;
   g->nextcol = NULL;

   if(inmemory)
   {
      g->nextrow = gmatrix_mem_nextrow;
      g->next_y = gmatrix_mem_next_y;
      g->nextcol = gmatrix_mem_nextcol;
      
      if(filename != NULL && x == NULL)
      {
	 if(g->pcor)
	    return gmatrix_load_pcor(g);
	 return gmatrix_load(g);
      }
   }

   return SUCCESS;
}

void gmatrix_free(gmatrix *g)
{
   int i;

   if(!g->inmemory && g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   if(g->mean)
      free(g->mean);

   if(g->sd)
      free(g->sd);

   g->mean = g->sd = NULL;

   /* if in memory, the sample struct will contain a
    * pointer to x and be freed later */
   if(g->x)
   {
      for(i = 0 ; i < g->n ; i++)
      {
	 if(g->x[i])
	    free(g->x[i]);
      }
      free(g->x);
      g->x = NULL;
   }

   if(g->y)
   {
      free(g->y);
      g->y = NULL;
   }
}

/* Expects the binary data row to be y, x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextrow(gmatrix *g, sample *s)
{
   int i;
   intype *tmp;
   
   s->inmemory = g->inmemory;

   MALLOCTEST(s->x, sizeof(dtype) * (g->p + 1))
   MALLOCTEST(tmp, sizeof(intype) * (g->p + 1))

   if(g->i == g->n)
      if(!gmatrix_reset(g))
	 return FAILURE;

   FREADTEST(tmp, sizeof(intype), g->p + 1, g->file)

   s->y = (dtype)tmp[0];
   s->x[0] = 1.0; /* intercept */
   for(i = 1 ; i < g->p + 1; i++)
      /*s->x[i] = ((dtype)tmp[i] - g->mean[i]) / g->sd[i];*/
      s->x[i] = (dtype)tmp[i];
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

      g->y[i] = (dtype)tmp[0];

      g->x[i][0] = 1.0; /* intercept */
      for(j = 1 ; j < g->p + 1 ; j++)
	 g->x[i][j] = (dtype)tmp[j];
   }

   fclose(fin);
   free(tmp);

   return SUCCESS;
}

/* No y variable
 * There is a useless intercept (zero due to scaling), but necessary to avoid
 * off-by-one problems in code that expects an intercept
 * Only p columns
 * The first x column is stored in g->y
 * */
int gmatrix_load_pcor(gmatrix *g)
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

      g->y[i] = 1.0;

      for(j = 0 ; j < g->p + 1 ; j++)
	 g->x[i][j] = (dtype)tmp[j];
   }

   fclose(fin);
   free(tmp);

   return SUCCESS;
}

/* Only applicable for memory-based matrices
 */
int gmatrix_scale(gmatrix *g)
{
   int i, j;
   double delta;
   sample sm;
   double *mean, *sd;

   MALLOCTEST(mean, sizeof(double) * (g->p + 1))
   MALLOCTEST(sd, sizeof(double) * (g->p + 1))

   sample_init(&sm, g->p + 1);

   /*mean[0] = 0;
   sd[0] = 1;*/

   /* sd is really the sum of squares, not the SD, but we
    * use the same variable to save memory */

   for(i = 0 ; i < g->n ; i++)
   {
      g->nextrow(g, &sm);

      for(j = 0 ; j < g->p + 1; j++)
      {
	 if(i == 0)
	    mean[j] = sd[j] = 0;

	 delta = (double)sm.x[j] - mean[j];
	 mean[j] += delta / (i + 1);
	 sd[j] += delta * ((double)sm.x[j] - mean[j]);
      }
   }

   for(j = 0 ; j < g->p + 1 ; j++)
      sd[j] = sqrt(sd[j] / (g->n - 1));

   sample_free(&sm);

   free(g->mean);
   free(g->sd);

   g->mean = mean;
   g->sd = sd;

   if(g->inmemory)
   {
      for(i = 0 ; i < g->n ; i++)
      {
	 for(j = 0 ; j < g->p + 1; j++)
	 {
	    g->x[i][j] = g->x[i][j] - g->mean[j];
	    if(g->sd[j] > SDTHRESH)
	       g->x[i][j] /= g->sd[j];
	 }
      }
   }

   return SUCCESS;
}

int gmatrix_mem_nextrow(gmatrix *g, sample *s)
{
   if(g->i == g->n)
      if(!gmatrix_reset(g))
	 return FAILURE;

   s->inmemory = g->inmemory;

   s->x = g->x[g->i];
   s->y = g->y[g->i];
   g->i++;

   return SUCCESS;
}

int gmatrix_mem_nextcol(gmatrix *g, sample *s)
{
   int i;

   if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   for(i = 0 ; i < g->n ; i++)
   {
      s->x[i] = g->x[i][g->j];
      /*s->y[j] = g->y[j];*/
   }

   g->j++;

   return SUCCESS;
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
   intype tmp;
   if(g->i == g->n)
      gmatrix_reset(g);

   fread(&tmp, sizeof(intype), 1, g->file);

   g->i++;

   /* ignore the x vector*/
   fseek(g->file, sizeof(intype) * g->p, SEEK_CUR);

   return (dtype)tmp;
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
   g->j = 0;
   return SUCCESS;
}

