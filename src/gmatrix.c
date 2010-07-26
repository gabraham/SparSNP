#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

int sample_init(sample *s, int n, short inmemory)
{
   s->inmemory = inmemory;
   if(!inmemory)
      MALLOCTEST(s->x, sizeof(dtype) * n)
   return SUCCESS;
}

void sample_free(sample *s)
{
   if(s->inmemory)
      return;
   
   if(s->x)
      free(s->x);
   s->x = NULL;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p, short inmemory)
{
   int i;

   if(filename)
      FOPENTEST(g->file, filename, "rb")

   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;
   g->yidx = 0;
   g->inmemory = inmemory;

   MALLOCTEST(g->y, sizeof(dtype) * g->n)

   if(inmemory)
   {
      if(!gmatrix_load(g))
	 return FAILURE;
      g->nextcol = gmatrix_mem_nextcol;
   }
   else
      g->nextcol = gmatrix_disk_nextcol;

   CALLOCTEST(g->mean, p + 1, sizeof(double))
   MALLOCTEST(g->sd, sizeof(double) * (p + 1))

   for(i = 0 ; i < g->p + 1 ; i++)
      g->sd[i] = 1;


   return SUCCESS;
}

int gmatrix_init2(gmatrix *g, char *filename, int n, int p)
{
   int i;

   if(filename)
      FOPENTEST(g->file, filename, "rb")

   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;

   CALLOCTEST(g->mean, p + 1, sizeof(double))
   MALLOCTEST(g->sd, sizeof(double) * (p + 1))

   for(i = 0 ; i < g->p + 1 ; i++)
      g->sd[i] = 1;

   MALLOCTEST(g->y, sizeof(dtype) * g->n)

   g->bufidx = 0;
   g->bufsize = BUFSIZE;
   MALLOCTEST(g->buffer, sizeof(dtype*) * g->n * g->bufsize)
   FREADTEST(g->y, sizeof(dtype), g->n, g->file)
   FREADTEST(g->buffer, sizeof(dtype), g->n * (int)fmin(g->bufsize, g->p), g->file)

   return SUCCESS;
}

void gmatrix_free(gmatrix *g)
{
   int j;

   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }

   if(g->mean)
      free(g->mean);

   if(g->sd)
      free(g->sd);

   g->mean = g->sd = NULL;

   if(g->y)
      free(g->y);
   g->y = NULL;

   if(g->inmemory && g->x)
   {
      for(j = 0 ; j < g->p + 1 ; j++)
      {
	 if(g->x[j])
	    free(g->x[j]);
      }
      free(g->x);
      g->x = NULL;
   }

   /*if(g->buffer)
      free(g->buffer);
   g->buffer = NULL;*/
}

/* Expects the binary data column to be y, x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextcol(gmatrix *g, sample *s)
{
   int i;
   
   if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   if(g->j == 0)
   {
      FREADTEST(g->y, sizeof(dtype), g->n, g->file)
      
      /* intercept */
      for(i = 0 ; i < g->n ; i++)
	 s->x[i] = 1.0;
   }
   else
      FREADTEST(s->x, sizeof(dtype), g->n, g->file)

   g->j++;

   return SUCCESS;
}

/* Expects the binary data column to be x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextcol_no_y(gmatrix *g, sample *s)
{
   if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   if(g->j == g->yidx) {
      FREADTEST(g->y, sizeof(dtype), g->n, g->file)
   } else {
      FREADTEST(s->x, sizeof(dtype), g->n, g->file)
   }

   g->j++;

   return SUCCESS;
}

int gmatrix_mem_nextcol(gmatrix *g, sample *s)
{
   if(g->j == g->p + 1)
   {
      g->i = g->j = 0;
   }
      
   s->x = g->x[g->j];

   g->j++;

   return SUCCESS;
 
}

int gmatrix_load(gmatrix *g)
{
   int i, j;
   sample s;

   if(!sample_init(&s, g->n, FALSE))
      return FAILURE;

   MALLOCTEST(g->x, sizeof(dtype*) * (g->p + 1))

   for(j = 0 ; j < g->p + 1 ; j++)
   {
      MALLOCTEST(g->x[j], sizeof(dtype) * g->n)
      if(!gmatrix_disk_nextcol(g, &s))
	 return FAILURE;

      for(i = 0 ; i < g->n ; i++)
	 g->x[j][i] = s.x[i];
   }

   sample_free(&s);

   return SUCCESS;
}

/* Expects the binary data column to be y, x_1, x_2, x_3, ..., x_p */
/*int gmatrix_disk_nextcol2(gmatrix *g, sample *s)
{
   int i;
   
   if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   if(g->j == 0)
   {
      FREADTEST(g->y, sizeof(dtype), g->n, g->file)
      
      * intercept *
      for(i = 0 ; i < g->n ; i++)
	 s->x[i] = 1.0;
   }
   else
   {
      if(g->bufidx == g->bufsize)
      {
	 FREADTEST(g->x, sizeof(dtype), g->n * g->bufsize, g->file)
	 g->bufidx = 0;
      }

      s->x = g->x[g->bufidx * g->n];
      
      g->bufidx++;
   }

   g->j++;

   return SUCCESS;
}*/

int gmatrix_reset(gmatrix *g)
{
   if(g->file)
   {
      fclose(g->file);
      g->file = NULL;
   }
   
   FOPENTEST(g->file, g->filename, "rb")

   g->i = g->j = 0;

   return SUCCESS;
}

