#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "ind.h"

/* Writes a binary indicator file. The first int is the number of cv folds.
 * The following ints are indicator variables for each of the k folds, of
 * length n, with 1 when the ith sample is training and 0 otherwise. 
 */
int ind_write(char *file, int *folds, int n, int nfolds)
{
   FILE *out = NULL;
   FOPENTEST(out, file, "wb");
   FWRITETEST(&nfolds, sizeof(int), 1, out);
   FWRITETEST(folds, sizeof(int), n * nfolds, out);
   fclose(out);
   return SUCCESS;
}

int ind_getfolds(char *file)
{
   int nfolds;
   FILE *in = NULL;
   FOPENTEST(in, file, "rb")
   FREADTEST(&nfolds, sizeof(int), 1, in);
   fclose(in);
   return nfolds;
}

int ind_read(char *file, int *folds, const int n, const int nfolds)
{
   int tmp;
   FILE *in = NULL;
   FOPENTEST(in, file, "rb")
   FREADTEST(&tmp, sizeof(int), 1, in); /* ignore it */
   FREADTEST(folds, sizeof(int), n * nfolds, in);
   fclose(in);
   return SUCCESS;
}

