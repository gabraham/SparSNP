#include <stdio.h>
#include "util.h"
#include "common.h"

int writebinvectorf(char* file, double* x, int n)
{
   FILE* out = NULL;
   FOPENTEST(out, file, "wb");
   FWRITETEST(x, sizeof(double), n, out);
   fclose(out);
   return SUCCESS;
}

int writebinvectorl(char* file, int* x, int n)
{
   FILE* out = NULL;
   FOPENTEST(out, file, "wb");
   FWRITETEST(x, sizeof(int), n, out);
   fclose(out);
   return SUCCESS;
}

int writevectorf(char* file, double* beta, int p)
{
   int i;
   FILE* out = NULL;
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
   FILE* out = NULL;
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

int readvectorl(char *filename, int *x, int n)
{
   int i = 0;
   FILE *in = NULL;
   FOPENTEST(in, filename, "rt");

   while(!feof(in))
   {
      if(fscanf(in, "%d", x + i) == EOF)
	 break;
      i++;
   } 
   fclose(in);
   return SUCCESS;
}

/* takes beta on original scale and puts it on zero-mean unit-variance scale 
 * of new data */
void scale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p)
{
   int j;
   double t, s = 0;
   for(j = p - 1 ; j >= 0 ; --j)
   {
      t = beta1[j] * mean[j];
      if(sd[j] != 0)
	 t /= sd[j];

      beta2[j] = beta1[j] * sd[j];
      s += t;
   }
   beta2[0] += s;
}

/* assumes beta0 is intercept, i.e. beta runs from 0 to p (inclusive).
 * 
 * beta_j = beta_j^* / sd_j
 *
 * beta_0 = beta_0^* - \sum_{j=1}^p beta_j^* mean_j / sd_j
 *
 * Note that zero beta^* remains zero in beta
 * */
void unscale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p)
{
   int j;
   double t, s = 0;
   for(j = p - 1 ; j >= 0 ; --j)
   {
      t = beta1[j] * mean[j];
      beta2[j] = beta1[j];
      if(sd[j] != 0)
      {
	 t /= sd[j];
	 beta2[j] /= sd[j];
      }
      s += t;
   }
   beta2[0] -= s;
}

