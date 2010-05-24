#include <stdlib.h>
#include <stdio.h>
#include "evaluation.h"

void auc_test(long seed)
{
   int i, n = 1000;
   double *yhat = malloc(sizeof(double) * n);
   int *y = malloc(sizeof(int) * n);

   srand48(seed);

   for(i = 0 ; i < n ; i++)
   {
      y[i] = drand48() > 0.5 ? 1 : 0;
      yhat[i] = drand48();
   }
   printf("AUC: %.5f\n", auc(yhat, y, n));

   for(i = 0 ; i < n ; i++)
      yhat[i] = y[i];

   printf("AUC: %.5f\n", auc(yhat, y, n));

   for(i = 0 ; i < n ; i++)
      yhat[i] = y[i] ? 0 : 1;

   printf("AUC: %.5f\n", auc(yhat, y, n));
}

int main()
{
   auc_test(123);

   return EXIT_SUCCESS;
}

