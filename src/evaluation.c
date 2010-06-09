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

   MALLOCTEST(y1, m1 * sizeof(double))
   MALLOCTEST(y2, m2 * sizeof(double))

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

double gmatrix_auc(double *yhat, gmatrix *g, int *trainf, int ntrain)
{
   int i, j, k, l;
   double s = 0;
   int m1 = 0, m2;
   double *y1, *y2;
   dtype z = 0;

   /*for(i = 0 ; i < g->n ; i++)
      if(g->next_y(g) && trainf[i])
	 m1++;*/
   for(i = 0 ; i < g->n ; i++)
   {
      z = g->next_y(g);
      if(z && trainf[i])
	 m1++;
   }

   m2 = ntrain - m1;

   if(!m1 || !m2)
   {
      printf("cannot evaluate AUC, too few observations: positives=%d\
 negatives=%d\n", m1, m2);
      return -1;
   }

   printf("%d %d\n", m1, m2);

   MALLOCTEST(y1, m1 * sizeof(double))
   MALLOCTEST(y2, m2 * sizeof(double))

   printf("Positives: %d  Negatives: %d\n", m1, m2);

   gmatrix_reset(g);

   j = k = l = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      z = g->next_y(g);
      if(trainf[i])
      {
	 if(z == ONE)
      	 {
      	    y1[j] = yhat[l];
      	    j++;
      	 }
      	 else
      	 {
      	    y2[k] = yhat[l];
      	    k++;
      	 }
	 l++;
      }
   }

   for(i = 0 ; i < m1 ; i++)
      for(j = 0 ; j < m2 ; j++)
	 s += (y1[i] > y2[j]) + 0.5 * (y1[i] == y2[j]);

   free(y1);
   free(y2);

   return s / (m1 * m2);
}

double gmatrix_accuracy(double *yhat, gmatrix *g, double threshold,
      int *trainf, int ntrain)
{
   int i, k = 0;
   long s = 0;
   int y;

   for(i = 0 ; i < g->n ; i++)
   {
      y = (int)(g->next_y(g));

      if(trainf[i])
      {
	 s += (yhat[k] >= threshold) == y;
	 k++;
      }
   }
   return (double)s / ntrain;
}


