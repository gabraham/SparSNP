#include <stdlib.h>

#include "gmatrix.h"


int main()
{
   int i, j, n = 100, p = 10;
   gmatrix g;
   sample sm;

   gmatrix_init(&g, "../R/out.bin", n, p);
   sample_init(&sm, p);

   for(i = 0 ; i < n ; i++)
   {
      gmatrix_nextrow(&g, &sm);
      printf("[%d]\ty=%d\t", i, (int)sm.y);
      for(j = 0 ; j < p ; j++)
	 printf("%.5f ", sm.x[j]);
      printf("\n");
   }

   sample_free(&sm);
   gmatrix_free(&g);

   return EXIT_SUCCESS;
}

