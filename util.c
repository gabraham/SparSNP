/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
 * All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "common.h"
#include "matrix.h"

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
      fprintf(out, "%.10f\n", beta[i]);
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

/* Writes an n by p matrix to file, in row-major ordering
 * assumes that the matrix is in column-major ordering
 */
int writematrixf(double *x, int n, int p, char* file)
{
   int i, j;
   FILE* out;

   FOPENTEST(out, file, "wt")
   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 fprintf(out, "%.10f", x[j * n + i]);
	 if(j < p - 1)
	    fprintf(out, ",");
	 else
	    fprintf(out, "\n");
      }
   }

   fflush(out);
   fclose(out);

   return SUCCESS;
}

/* Writes an n by p matrix to file, in row-major ordering
 * assumes that the matrix is in column-major ordering
 */
int writematrixl(int *x, int n, int p, char* file)
{
   int i, j;
   FILE* out;

   FOPENTEST(out, file, "wt")
   for(i = 0 ; i < n ; i++)
   {
      for(j = 0 ; j < p ; j++)
      {
	 fprintf(out, "%d", x[j * n + i]);
	 if(j < p - 1)
	    fprintf(out, ",");
	 else
	    fprintf(out, "\n");
      }
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
 * beta_j = beta_j^* / sd_j for j = 1,...,p
 *
 * beta_0 = beta_0^* - \sum_{j=1}^p beta_j^* mean_j / sd_j
 *
 * Note that zero in beta^* remains zero in beta
 * */
void unscale_beta(double *beta2, double *beta1,
      double *mean, double *sd, int p, int K)
{
   int j, k;
   double t, s;

   for(k = 0 ; k < K ; k++)
   {
      s = 0;
      for(j = 1 ; j < p ; j++)
      {
         t = beta1[p * k + j] * mean[j];
         beta2[p * k + j] = beta1[p * k + j];
         if(sd[j] != 0) /* TODO: better checking for zero */
         {
            t /= sd[j];
            beta2[p * k + j] /= sd[j];
         }
         s += t;
      }
      beta2[p * k] = beta1[p * k] - s;
   }
}

int write_beta_sparse(char* file, double* beta, int p, int K)
{
   int j = 0, k;
   FILE* out = NULL;
   FOPENTEST(out, file, "w");
   
   /* always write the intercept */
   fprintf(out, "%d,", j);
   for(k = 0 ; k < K ; k++)
   {
      fprintf(out, "%.10f", beta[p * k + j]);
      if(k == K - 1)
	 fprintf(out, "\n");
      else
	 fprintf(out, ",");
   }

   /* now the other variables, don't print out a row with weights for 
    * all K tasks that are zero */
   for(j = 1 ; j < p ; j++)
   {
      /* check whether row is all zero */
      for(k = 0 ; k < K ; k++)
	 if(beta[p * k + j] != 0)
	    break;

      /* at least one non-zero in row */
      if(k < K)
      {
	 fprintf(out, "%d,", j);
	 for(k = 0 ; k < K ; k++)
	 {
	    if(beta[p * k + j] != 0)
	       fprintf(out, "%.10f", beta[p * k + j]);
	    else
	       fprintf(out, "0"); /* less space, TODO: use a sparser format */
	    if(k == K - 1)
	       fprintf(out, "\n");
	    else
	       fprintf(out, ",");
	 }
      }
   }

   fflush(out);
   fclose(out);
   return SUCCESS;
}

/* expects format:
 * <variable index>:<variable value>
 * %d:%lf
 * one row per non-zero variable
 *  
 * index 0 is intercept
 *
 * Assumes that we don't know what the number of tasks K is beforehand, we
 * just use an upper bound
 *
 * also assumes that each row of beta is stored as a *dense* vector when any
 * of the tasks have non-zero weights
 */
int load_beta_sparse(double *beta, char *filename, int p)
{
   int j = 0, k, K = 0, m;
   FILE *in = NULL;
   char line[MAX_LINE_CHARS];
   char *tmp = NULL;

   CALLOCTEST(tmp, MAX_STR_LEN, sizeof(char));

   FOPENTEST(in, filename, "rt");

   while(!feof(in))
   {
      /* read the variable id */
      if(fscanf(in, "%d,", &j) < 0)
	 break;
	 
      if(fgets(line, MAX_LINE_CHARS, in) == NULL)
      {
      	 if(!feof(in))
	 {
	    fprintf(stderr, "error in reading FAM file '%s'\n", filename);
	    return FAILURE;
	 }
	 break;
      }
      else
      {
         m = 0;
	 k = 0;
         while(k < MAX_NUM_PHENO && sscanf(&line[m], "%[^,\n]", tmp) != EOF)
         {
	    beta[k * p + j] = atof(tmp);
	    m += strlen(tmp) + 1;
	    k++;
         }
	 K = k;
      }
   } 

   fclose(in);

   FREENULL(tmp);

   /* we don't check whether the phenotype file is well formed */
   return K;
}

