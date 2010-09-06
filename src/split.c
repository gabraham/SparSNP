#include <stdio.h>
#include <stdlib.h>

/*
 * Split data into training and test set
 */
unsigned int split(int *trainf, int n, int nfolds)
{
   unsigned int i, ntrain = 0;
   double prop = 1.0 / nfolds;

   for(i = 0 ; i < n ; i++)
   {
      trainf[i] = drand48() >= prop;
      ntrain += trainf[i];
   }

   return ntrain;
}

int main(int argc, char *argv[])
{
   
   unsigned int ntrain, n = 0;
   int *trainf = NULL;
   
   MALLOCTEST(trainf, sizeof(int) * n);

   ntrain = split(trainf, n, nfolds);

   return writevectorl(opt->subsetfile, opt->trainf, opt->n);

   free(trainf);
}

