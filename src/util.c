#include <stdio.h>
#include "util.h"

void writevectorf(char* file, double* beta, int p)
{
   int i;
   FILE* out = fopen(file, "w");
   for(i = 0 ; i < p ; i++)
      fprintf(out, "%.20f\n", beta[i]);
   fflush(out);
   fclose(out);
}

void writevectorl(char* file, int* beta, int p)
{
   int i;
   FILE* out = fopen(file, "w");
   for(i = 0 ; i < p ; i++)
      fprintf(out, "%d\n", beta[i]);
   fflush(out);
   fclose(out);
}


