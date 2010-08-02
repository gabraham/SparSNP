#include <stdio.h>
#include <stdlib.h>

#include "gmatrix.h"

#define YFORMAT01 0
#define YFORMAT11 1

int sample_init(sample *s, int n, short inmemory)
{
   s->inmemory = inmemory;
   s->x = NULL;
   s->x2 = NULL;
   s->intercept = FALSE;
   return SUCCESS;
}

void sample_free(sample *s)
{
   if(s->intercept)
      return;
   
   if(s->x)
      free(s->x);
   s->x = NULL;

   if(s->x2)
      free(s->x2);
   s->x2 = NULL;
}

int gmatrix_init(gmatrix *g, char *filename, int n, int p, short inmemory,
      short tabulate, char *scalefile, short yformat)
{
   int i;

   if(filename)
      FOPENTEST(g->file, filename, "rb")

   g->filename = filename;
   g->i = g-> j = 0;
   g->n = n;
   g->p = p;
   g->yidx = 0;
   g->y = NULL;
   g->inmemory = inmemory;
   g->tab = NULL;
   g->scalefile = scalefile;
   g->lookup = NULL;
   g->lookup2 = NULL;
   g->intercept = NULL;
   g->mean = NULL;
   g->sd = NULL;
   g->tmp = NULL;
   g->active = NULL;
   g->yformat = yformat;

   MALLOCTEST(g->intercept, sizeof(double) * g->n)
   for(i = 0 ; i < n ; i++)
      g->intercept[i] = 1.0;

   MALLOCTEST(g->tmp, sizeof(dtype) * g->n)

   if(tabulate)
   {
      /*if(!gmatrix_tabulation_load(g))
	 return FAILURE;
      g->nextcol = gmatrix_tabulation_nextcol;*/
   }
   else if(inmemory)
   {
      if(!gmatrix_load(g))
	 return FAILURE;
      g->nextcol = gmatrix_mem_nextcol;
   }
   else
      g->nextcol = gmatrix_disk_nextcol;

   if(scalefile && !gmatrix_read_scaling(g, scalefile))
      return FAILURE;

   return SUCCESS;
}

/*int gmatrix_init2(gmatrix *g, char *filename, int n, int p)
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
}*/

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

   if(g->active)
      free(g->active);
   g->active = NULL;

   if(g->tmp)
      free(g->tmp);
   g->tmp = NULL;

   if(g->intercept)
      free(g->intercept);
   g->intercept = NULL;

   if(g->lookup)
      free(g->lookup);
   g->lookup = NULL;

   if(g->lookup2)
      free(g->lookup2);
   g->lookup2 = NULL;

   /*if(g->tab)
      tabulation_free(g->tab);
   free(g->tab);
   g->tab = NULL;

   *if(g->buffer)
      free(g->buffer);
   g->buffer = NULL;*/
}

int gmatrix_disk_skipcol(gmatrix *g)
{
   if(g->j == 0)
   {
   }
   FSEEKTEST(g->file, sizeof(dtype) * g->n, SEEK_CUR)
   return SUCCESS;
}

/* Expects the binary data column to be y, x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextcol(gmatrix *g, sample *s)
{
   int i, l1, l2;
   
   if(g->j == g->p + 1 && !gmatrix_reset(g))
      return FAILURE;

   /* intercept treated separately */
   if(g->j == 0)
   {
      s->intercept = TRUE;

      /* read y the first time we see it */
      if(!g->y)
      {
	 MALLOCTEST(g->y, sizeof(double) * g->n)
	 FREADTEST(g->tmp, sizeof(dtype), g->n, g->file)
	 if(g->yformat == YFORMAT01)
	    for(i = 0 ; i < g->n ; i++)
	       g->y[i] = (double)g->tmp[i];
	 else
	    for(i = 0 ; i < g->n ; i++)
	       g->y[i] = 2.0 * g->tmp[i] - 1.0;
      }
      else /* don't read y again */
	 FSEEKTEST(g->file, sizeof(dtype) * g->n, SEEK_CUR)
      
      /* No need to copy values, just assign the intercept vector,
       * but first, free any old data we have from previous
       * iterations*/
      if(s->x)
	 free(s->x);
      if(s->x2)
	 free(s->x2);
      s->x = g->intercept;
      s->x2 = g->intercept;
   }
   else
   {
      /* not intercept, need to copy values */
      s->intercept = FALSE;
      if(g->j == 1)
      {
	 MALLOCTEST(s->x, sizeof(double) * g->n)
	 MALLOCTEST(s->x2, sizeof(double) * g->n)
      }
      
      if(g->scalefile)
      {
         s->intercept = FALSE;
         FREADTEST(g->tmp, sizeof(dtype), g->n, g->file)

         /* Get the scaled value instead of the original value */
	 l1 = g->j * NUM_X_LEVELS;
         for(i = 0 ; i < g->n ; i++)
	 {
	    l2 = l1 + g->tmp[i];
            s->x[i] = g->lookup[l2];
            s->x2[i] = g->lookup2[l2];
	 }
      }
      else
      {
         FREADTEST(g->tmp, sizeof(dtype), g->n, g->file)
         for(i = 0 ; i < g->n ; i++)
	 {
            s->x[i] = (double)g->tmp[i];
            s->x2[i] = s->x[i] * s->x[i];
	 }
      }
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

/* Expects the binary data column to be x_1, x_2, x_3, ..., x_p */
int gmatrix_disk_nextcol_no_y(gmatrix *g, sample *s)
{
   /*if(g->j == g->p + 1)
      if(!gmatrix_reset(g))
	 return FAILURE;

   if(g->j == g->yidx) {
      FREADTEST(g->y, sizeof(dtype), g->n, g->file)
   } else {
      FREADTEST(s->x, sizeof(dtype), g->n, g->file)
   }

   g->j++;*/

   return SUCCESS;
}

/*int gmatrix_tabulation_nextcol(gmatrix *g, sample *s)
{
   if(g->j == g->p + 1)
   {
      g->i = g->j = 0;
   }
      
   s->counts = g->tab->counts[g->j];
   s->values = g->tab->values;
   s->nbins = g->tab->nbins;
   g->j++;

   return SUCCESS;
}*/

int gmatrix_load(gmatrix *g)
{
   int i, j;
   sample s;

   if(!sample_init(&s, g->n, FALSE))
      return FAILURE;

   MALLOCTEST(g->x, sizeof(double*) * (g->p + 1))

   for(j = 0 ; j < g->p + 1 ; j++)
   {
      MALLOCTEST(g->x[j], sizeof(double) * g->n)
      if(!gmatrix_disk_nextcol(g, &s))
	 return FAILURE;

      for(i = 0 ; i < g->n ; i++)
	 g->x[j][i] = (double)s.x[i];
   }

   sample_free(&s);

   return SUCCESS;
}

/*int gmatrix_tabulation_load(gmatrix *g)
{
   MALLOCTEST(g->tab, sizeof(tabulation))

   * TODO: magic number, replace with xlevels * ylevels *
   if(!tabulation_init(g->tab, g->p + 1, 6))
      return FAILURE;

   if(!tabulation_read(g->tab, g->filename))
      return FAILURE;

   return SUCCESS;
}*/

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

/* Populate lookup tables for SNP levels 0, 1, 2
 */
int gmatrix_read_scaling(gmatrix *g, char *file_scale)
{
   int j, k, l1, l2;
   FILE *in;

   CALLOCTEST(g->lookup, NUM_X_LEVELS * (g->p + 1), sizeof(double))
   CALLOCTEST(g->lookup2, NUM_X_LEVELS * (g->p + 1), sizeof(double))
   MALLOCTEST(g->mean, sizeof(double) * (g->p + 1))
   MALLOCTEST(g->sd, sizeof(double) * (g->p + 1))
   MALLOCTEST(g->active, sizeof(int) * (g->p + 1))

   FOPENTEST(in, file_scale, "rb")

   FREADTEST(g->mean, sizeof(double), g->p + 1, in);
   FREADTEST(g->sd, sizeof(double), g->p + 1, in);

   /* intercept */
   for(k = 0 ; k < NUM_X_LEVELS ; k++)
   {
      g->lookup[k] = 1;
      g->lookup2[k] = 1;
   }
   
   g->active[0] = TRUE;

   for(j = 1 ; j < g->p + 1; j++)
   {
      if((g->active[j] = (g->sd[j] != 0)))
      {
	 l1 = j * NUM_X_LEVELS;
	 for(k = 0 ; k < NUM_X_LEVELS ; k++)
	 {
	    l2 = l1 + k;
	    g->lookup[l2] = (k - g->mean[j]) / g->sd[j];
	    g->lookup2[l2] = g->lookup[l2] * g->lookup[l2];
	 }
      }
   }
   
   fclose(in);

   return SUCCESS;
}

