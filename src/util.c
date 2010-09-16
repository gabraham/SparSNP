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

/* Assumes ascii, one value per line */
int load_beta(double *beta, char *filename, int p)
{
   int i = 0;
   FILE *in = NULL;
   FOPENTEST(in, filename, "rt");

   while(!feof(in))
   {
      if(fscanf(in, "%lf", beta + i) == EOF)
	 break;
      i++;
   } 
   fclose(in);
   return SUCCESS;
}

void scale_beta(double *beta, double *mean, double *sd, int p)
{
   int j;

   beta[0] = 1.0;

   for(j = 1 ; j < p ; j++)
   {
      beta[j] = beta[j] - mean[j];
      if(sd[j] != 0)
	 beta[j] = (beta[j] - mean[j]) / sd[j];
   }
}

