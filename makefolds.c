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

/*
 * indicator vector of length nfolds * n, each group of n int represents
 * whether the samples in the kth fold are training (1) or testing (0)
 */
void make_ind(int *ind, int *folds, int n, int nfolds)
{
   int i, k;

   for(k = 0 ; k < nfolds ; k++)
      for(i = 0 ; i < n ; i++)
	 ind[k * n + i] = (folds[i] != k);
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

   /* write folds file (which fold each sample belongs to) */
   if(!writevectorl(file_folds, folds, n))
   {
      free(folds);
      return EXIT_FAILURE;
   }

   /* write indicator file (in each fold, which samples are training and which
    * are testing) */
   make_ind(ind, folds, n, nfolds);
   ret = ind_write(file_ind, ind, n, nfolds);

   free(folds);
   free(ind);

   return ret && EXIT_SUCCESS;
}

