#include <math.h>

#include "gmatrix.h"

#include "scale.h"

void scale(gmatrix *g, double *mean, double *sd)
{
   int i, j;
   int n = g->n;
   int p = g->p;
   double delta;
   sample sm;

   sample_init(&sm, p);
 
   /* sd is really the sum of squares, not the SD, but we
    * use the same variable to save memory */

   for(i = 0 ; i < n ; i++)
   {
      g->nextrow(g, &sm);

      for(j = 0 ; j < p ; j++)
      {
          if(i == 0)
             mean[j] = sd[j] = 0;

          delta = sm.x[j] - mean[j];
          mean[j] += delta / (i + 1);
          sd[j] += delta * (sm.x[j] - mean[j]);
       }
   }

   for(j = 0 ; j < p ; j++)
      sd[j] = sqrt(sd[j] / (n - 1));

   sample_free(&sm);
}

