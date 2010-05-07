#include <stdlib.h>
#include <stdio.h>

#include "gmatrix.h"

double accuracy(double *yhat, int *y, int n)
{
   int s = 0;
   int i;
   for(i = 0 ; i < n ; i++)
      s += (yhat[i] >= 0.5) == y[i];
   return (double)s / n;
}

double auc(double *yhat, int *y, int n)
{
   int i, j, k;
   double s = 0;
   int m1 = 0, m2;
   double *y1, *y2;

   for(i = 0 ; i < n ; i++)
      m1 += y[i];

   m2 = n - m1;

   y1 = malloc(m1 * sizeof(double));
   y2 = malloc(m2 * sizeof(double));

   i = j = k = 0;
   while(i < n)
   {
      if(y[i] == 1)
      {
	 y1[j] = yhat[i];
	 j++;
      }
      else
      {
	 y2[k] = yhat[i];
	 k++;
      }

      i++;
   }

   for(i = 0 ; i < m1 ; i++)
      for(j = 0 ; j < m2 ; j++)
	 s += (y1[i] > y2[j]) + 0.5 * (y1[i] == y2[j]);

   free(y1);
   free(y2);

   return s / (m1 * m2);
}

double gmatrix_auc(double *yhat, gmatrix *g)
{
   int i, j, k;
   double s = 0;
   int m1 = 0, m2;
   double *y1, *y2;
   double y;

   for(i = 0 ; i < g->n ; i++)
      m1 += gmatrix_next_y(g);

   m2 = g->n - m1;

   y1 = malloc(m1 * sizeof(double));
   y2 = malloc(m2 * sizeof(double));

   gmatrix_reset(g);

   i = j = k = 0;
   while(i < g->n)
   {
      y = gmatrix_next_y(g); 
      if(y == 1)
      {
	 y1[j] = yhat[i];
	 j++;
      }
      else
      {
	 y2[k] = yhat[i];
	 k++;
      }

      i++;
   }

   for(i = 0 ; i < m1 ; i++)
      for(j = 0 ; j < m2 ; j++)
	 s += (y1[i] > y2[j]) + 0.5 * (y1[i] == y2[j]);

   free(y1);
   free(y2);

   return s / (m1 * m2);
}


