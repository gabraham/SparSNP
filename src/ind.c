#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "ind.h"

int write_ind(char *file, int *folds, int n, int nfolds)
{
   FILE *out = NULL;
   FOPENTEST(out, file, "wb");
   FWRITETEST(&nfolds, sizeof(int), 1, out);
   FWRITETEST(folds, sizeof(int), n * nfolds, out);
   fclose(out);
   return SUCCESS;
}

int read_ind(char *file, int *folds, int n, int nfolds)
{
   FILE *in = NULL;
   FOPENTEST(in, file, "rb")
   FREADTEST(&nfolds, sizeof(int), 1, in);
   FREADTEST(folds, sizeof(int), n * nfolds, in);
   fclose(in);
   return SUCCESS;
}

