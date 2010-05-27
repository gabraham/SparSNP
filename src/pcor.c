#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "gmatrix.h"
#include "sgd.h"
#include "loss.h"
#include "evaluation.h"

int main()
{
   int i, j;
   int n = 1000;
   int p = 100;
   long seed = time(NULL);
   double **x = malloc(sizeof(double*) * n);

   srand48(seed);

   for(i = 0 ; i < n ; i++)
   {
      x[i] = malloc(sizeof(double) * p);
      for(j = 0 ; j < p ; j++)
      {
	 x[i][j] = drand48();
      }
   }


   
   return EXIT_SUCCESS;
}

