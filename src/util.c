#include <stdio.h>
#include "util.h"
#include "common.h"

int writevectorf(char* file, double* beta, int p)
{
   int i;
   FILE* out;
   
   FOPENTEST(out, file, "w")
   for(i = 0 ; i < p ; i++)
      fprintf(out, "%.20f\n", beta[i]);
   fflush(out);
   fclose(out);

   return SUCCESS;
}

int writevectorl(char* file, int* beta, int p)
{
   int i;
   FILE* out;
   
   FOPENTEST(out, file, "w")
   for(i = 0 ; i < p ; i++)
      fprintf(out, "%d\n", beta[i]);
   fflush(out);
   fclose(out);

   return SUCCESS;
}

int writematrixf(double **x, int n, int p, char* file)
{
   int i, j;
   FILE* out;
   
   FOPENTEST(out, file, "wt")
   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 fprintf(out, "%.20f", x[i][j]);
	 if(j < p - 1)
	    fprintf(out, ",");
      }
      fprintf(out, "\n");
   }

   fflush(out);
   fclose(out);

   return SUCCESS;
}

