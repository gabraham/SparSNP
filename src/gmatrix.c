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

void gmatrix_test()
{
   int n = 1e2,
       p = 5;
   int i, j;
   double **x = malloc(n * sizeof(double*));
   double s;
   int *y;
   
   srand48(12345);

   y = malloc(n * sizeof(int));
   

   for(i = 0 ; i < n ; i++)
   {
      x[i] = malloc(p * sizeof(double));
      s = 1;
      for(j = 0 ; j < p ; j++)
      {
	 x[i][j] = drand48() - 0.5;
	 s += x[i][j];
      }
      y[i] = drand48() <= plogis(s) ? 1 : 0;
      
   }

}

