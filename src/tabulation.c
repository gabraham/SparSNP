#include "cd.h"
#include "util.h"

int tabulation_init(tabulation *t , int p, int nbins)
{
   int j;
   t->nbins = nbins;
   t->p = p;

   MALLOCTEST(t->values, sizeof(double) * nbins)

   MALLOCTEST(t->counts, sizeof(int*) * p)
   for(j = 0 ; j < p ; j++)
      CALLOCTEST(t->counts[j], nbins, sizeof(int))

   return SUCCESS;
}

void tabulation_free(tabulation *t)
{
   int j;
   if(t->counts)
      for(j = 0 ; j < t->p ; j++)
	 free(t->counts[j]);
   free(t->counts);
   t->counts = NULL;

   if(t->values)
      free(t->values);
   t->values = NULL;
}

int tabulation_write(tabulation *t, char* fileout)
{
   int j;
   FILE* out;

   FOPENTEST(out, fileout, "wb")

   for(j = 0 ; j < t->p ; j++)
      FWRITETEST(t->counts[j], sizeof(int), t->nbins, out);

   fclose(out);

   return SUCCESS;
}

int tabulation_read(tabulation *t, char* filein)
{
   int j;
   FILE* in;

   if(!t->counts)
      return FAILURE;

   FOPENTEST(in, filein, "rb")

   for(j = 0 ; j < t->p ; j++)
   {
      if(!t->counts[j])
	 return FAILURE;
      FREADTEST(t->counts[j], sizeof(int), t->nbins, in);
   }

   fclose(in);

   return SUCCESS;
}

/* For each variable j, produce a crosstabulation of its values
 * with those of y
 */
int tabulate(gmatrix *g, char *fileout, short doscale,
      int xlevels, int ylevels, short verbose)
{
   int i, j, k;
   tabulation t;
   sample sm;

   if(!tabulation_init(&t, g->p + 1, xlevels * ylevels))
      return FAILURE;

   if(!sample_init(&sm, g->n, g->inmemory))
      return FAILURE;

   for(j = 0 ; j < g->p + 1 ; j++)
   {
      g->nextcol(g, &sm);

      for(i = 0 ; i < g->n ; i++)
      {
	 for(k = 0 ; k < xlevels ; k++)
      	 {
      	    if(sm.x[i] == k)
	       t.counts[j][xlevels * (int)g->y[i] + k]++;
      	 }
      }

      if(verbose)
      {
	 for(k = 0 ; k < xlevels * ylevels ; k++)
      	    printf("%d ", t.counts[j][k]);
      	 printf("\n");
      }
   }

   if(!tabulation_write(&t, fileout))
      return FAILURE;

   tabulation_free(&t);
   sample_free(&sm);

   return SUCCESS;
}

