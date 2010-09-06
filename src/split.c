#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "util.h"
#include "ind.h"

/*
 * Split data into training and test set
 */
void split(int *folds, int n, int nfolds)
{
   int i;

   for(i = 0 ; i < n ; i++)
      folds[i] = (int)(drand48() * nfolds);
}

void make_ind(int *ind, int *folds, int n, int nfolds)
{
   int i, j;

   for(i = 0 ; i < nfolds ; i++)
      for(j = 0 ; j < n ; j++)
	 ind[i * n + j] = (folds[j] == i);
}

int main(int argc, char *argv[])
{
   long seed = time(NULL);
   int ret, i, nfolds = 10, n = 0;
   char *file_folds = NULL,
        *file_ind = NULL;
   int *folds = NULL,
       *ind = NULL;

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-n"))
      {
	 i++;
	 n = (int)atof(argv[i]); 
      }
      else if(strcmp2(argv[i], "-nfolds"))
      {
	 i++;
	 nfolds = (int)atof(argv[i]); 
      }
      else if(strcmp2(argv[i], "-folds"))
      {
	 i++;
	 file_folds = argv[i];
      }
      else if(strcmp2(argv[i], "-ind"))
      {
	 i++;
	 file_ind = argv[i];
      }
      else if(strcmp2(argv[i], "-seed"))
      {
	 i++;
	 seed = atol(argv[i]);
      }
   }

   if(n == 0 || file_folds == NULL)
   {
      printf("usage: split -folds <filename> -ind <filename> \
-nfolds <#nfolds> -n <#n> [-seed <seed>]\n");
      return EXIT_FAILURE;
   }

   srand48(seed);
   
   MALLOCTEST(folds, sizeof(int) * n);
   MALLOCTEST(ind, sizeof(int) * n * nfolds);

   split(folds, n, nfolds);

   if(!writevectorl(file_folds, folds, n))
   {
      free(folds);
      return EXIT_FAILURE;
   }

   make_ind(ind, folds, n, nfolds);

   ret = write_ind(file_ind, ind, n, nfolds);

   free(folds);
   free(ind);

   return ret && EXIT_SUCCESS;
}

