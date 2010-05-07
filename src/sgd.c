#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "gmatrix.h"
#include "loss.h"
#include "evaluation.h"

void writebeta(char*, double*, int);
double sgd_gmatrix(gmatrix *, double,
      int, double *, double, double,
      double, int);
void predict_logloss(gmatrix *, double *, double *);
void scale(gmatrix *, double *, double *);

double sgd_gmatrix(gmatrix *g, double maxstepsize,
      int maxepoch, double *beta, double lambda1, double lambda2,
      double threshold, int verbose)
{
   int epoch = 1, i, j;
   double *grad = malloc((g->p + 1) * sizeof(double));
   double *x = malloc((g->p + 1) * sizeof(double));
   double prevloss = 0, loss;
   double stepsize = maxstepsize;
   sample sm;

   sample_init(&sm, g->p);

   while(epoch <= maxepoch)
   {
      loss = 0;
      for(i = 0 ; i < g->n ; i++)
      { 
	 gmatrix_nextrow(g, &sm);

	 /* Intercept */
	 x[0] = 1;

	 /* Scale parameters except the intercept */
	 for(j = 0 ; j < g->p ; j++)
	    x[j+1] = (sm.x[j] - g->mean[j]) / g->sd[j];

	 logdloss(x, beta, sm.y, g->p + 1, grad);
	 loss += logloss_pt(x, beta, sm.y, g->p + 1);

	 /* Update weights */
	 for(j = 0 ; j < g->p + 1; j++)
	 {
	    beta[j] -= stepsize * (grad[j] 
	       + lambda1 * sign(grad[j]) 
	       + lambda2 * grad[j] * grad[j]);
	 }
      }
      stepsize = stepsize / (1 + maxstepsize);

      if(verbose)
      {
	 printf("Epoch %d Loss: %.5f stepsize: %.15f\n",
	       epoch, loss, stepsize);
      }

      if(fabs(prevloss - loss) <= threshold)
      {
	 if(verbose)
	    printf("Termination condition met\n");
	 break;
      }

      prevloss = loss;

      epoch++;
   }

   sample_free(&sm);
   free(grad);
   free(x);
   return loss;
}

void predict_logloss(gmatrix *g, double *beta, double *yhat)
{
   int i, j;
   sample sm;
   double *x = malloc(sizeof(double) * (g->p + 1));

   sample_init(&sm, g->p);

   for(i = 0 ; i < g->n ; i++)
   {
      gmatrix_nextrow(g, &sm);
      x[0] = 1;
      for(j = 0 ; j < g->p ; j++)
	 x[j+1] = (sm.x[j] - g->mean[j]) / g->sd[j];
      yhat[i] = 1 / (1 + exp(-dotprod(x, beta, g->p + 1)));
   } 

   sample_free(&sm);
   free(x);
}

void writeout(char* file, double **x, int *y, double n, double p)
{
   int i, j;
   FILE* out = fopen(file, "w+");

   for(i = 0 ; i < n ; i++)
   {
      fprintf(out, "%d,", y[i]);
      for(j = 0 ; j < p ; j++)
      {
	 if(j < p - 1)
	    fprintf(out, "%.20f,", x[i][j]);
	 else
	    fprintf(out, "%.20f\n", x[i][j]);
      }
   }

   fflush(out);
   fclose(out);
}

void scale(gmatrix *g, double *mean, double *sd)
{
   int i, j;
   int n = g->n;
   int p = g->p;
   double delta;
   sample sm;

   sample_init(&sm, p);
 
   /* sd is really the sum of squares, not the SD, but we
    * use the same variable to save memory */

   for(i = 0 ; i < n ; i++)
   {
      gmatrix_nextrow(g, &sm);

     for(j = 0 ; j < p ; j++)
     {
         if(i == 0)
            mean[j] = sd[j] = 0;

         delta = sm.x[j] - mean[j];
         mean[j] += delta / (i + 1);
         sd[j] += delta * (sm.x[j] - mean[j]);
      }
   }

   for(j = 0 ; j < p ; j++)
      sd[j] = sqrt(sd[j] / (n - 1));

   sample_free(&sm);
}

/* void scale_test()
{
   int i, j;
   int n = 1000, p = 5;
   double **s = malloc(sizeof(double*) * n);
   double **x = malloc(sizeof(double*) * n);
   double mean, sumsq;

   for(i = 0 ; i < n ; i++)
   {
      x[i] = malloc(sizeof(double) * p);
      s[i] = malloc(sizeof(double) * p);
      for(j = 0 ; j < p ; j++)
	 x[i][j] = drand48();
   }

   scale(x, n, p, s);
   for(i = 0 ; i < n ; i++)
      mean += x[i][0] / n;
   printf("Mean x[,1]: %.5f\n", mean);

   mean = 0;
   sumsq = 0;
   for(i = 0 ; i < n ; i++)
      mean += s[i][0] / n;

   for(i = 0 ; i < n ; i++)
      sumsq += pow(s[i][0] - mean, 2);
   printf("Mean s[,1]: %.5f SD s[,1]: %.5f\n", mean, sqrt(sumsq / (n - 1)));
      
} */

/*double test()
{
   int n = 1e2,
       p = 5; * not including intercept *
   int i, j;
   double **xtmp = malloc(n * sizeof(double*));
   double **x = malloc(n * sizeof(double*));
   double *betahat = malloc((p + 1) * sizeof(double));
   double *tmp = malloc(p * sizeof(double));
   int *y;
   double s, loss, err;
   double const Rloss = 68.76173782419018;
   
   srand48(12345);

   for(i = 0 ; i < n ; i++)
   {
      xtmp[i] = malloc(p * sizeof(double));
      x[i] = malloc(p * sizeof(double));

      for(j = 0 ; j < p ; j++)
	 xtmp[i][j] = drand48() - 0.5;
   }

   scale(xtmp, n, p, x);

   * Add intercept term to the scaled data *
   for(i = 0 ; i < n ; i++)
   {
      memcpy(tmp, x[i], sizeof(double) * p);
      free(x[i]);
      x[i] = malloc((p + 1) * sizeof(double));
      x[i][0] = 1;
      for(j = 1 ; j < p + 1; j++)
	 x[i][j] = tmp[j - 1];
   }
   free(tmp);

   y = malloc(n * sizeof(int));
   for(i = 0 ; i < n ; i++)
   {
      s = 1; * intercept  *
      for(j = 0 ; j < p + 1 ; j++)
	 s += x[i][j];
      y[i] = drand48() <= plogis(s) ? 1 : 0;
   }

   loss = sgd_mem(x, y, n, p + 1, 1e-3, 1, betahat, 0, 0, FALSE);
   err = pow(loss - Rloss, 2);
   return err;
}*/

void writebeta(char* file, double* beta, int p)
{
   int i;
   FILE* out = fopen(file, "w");
   for(i = 0 ; i < p ; i++)
      fprintf(out, "%.20f\n", beta[i]);
   fflush(out);
   fclose(out);
}

int main(int argc, char* argv[])
{
   int i;
   double *betahat;
   double *yhat;
   gmatrix g;
   char *filename = NULL;
   char *model = NULL;
   char *betafile = "beta.csv";
   char *predfile = "pred.csv";
   int n = 0, p = 0;
   int verbose = FALSE;

   /* Parameters */
   int maxepochs = 20;
   double stepsize = 1e-4;
   double lambda1 = 0;
   double lambda2 = 0;
   double threshold = 1e-6;
   /* double alpha = 0; */

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp(argv[i], "-f") == 0)
      {
	 i++;
	 filename = argv[i];
      }
      else if(strcmp(argv[i], "-m") == 0)
      {
	 i++;
	 model = argv[i];
      }
      else if(strcmp(argv[i], "-n") == 0)
      {
	 i++;
	 n = atoi(argv[i]);
      }
      else if(strcmp(argv[i], "-p") == 0)
      {
	 i++;
	 p = atoi(argv[i]);
      }
      else if(strcmp(argv[i], "-e") == 0)
      {
	 i++;
	 maxepochs = atoi(argv[i]);
      }
      else if(strcmp(argv[i], "-s") == 0)
      {
	 i++;
	 stepsize = atof(argv[i]);
      }
      else if(strcmp(argv[i], "-l1") == 0)
      {
	 i++;
	 lambda1 = atof(argv[i]);
      }
      else if(strcmp(argv[i], "-l2") == 0)
      {
	 i++;
	 lambda2 = atof(argv[i]);
      }
      else if(strcmp(argv[i], "-t") == 0)
      {
	 i++;
	 threshold = atof(argv[i]);
      }
      else if(strcmp(argv[i], "-v") == 0)
      {
	 verbose = TRUE;
      }
      else if(strcmp(argv[i], "-vv") == 0)
      {
	 verbose = 2;
      }
      else if(strcmp(argv[i], "-b") == 0)
      {
	 i++;
	 betafile = argv[i];
      }
      else if(strcmp(argv[i], "-p") == 0)
      {
	 i++;
	 predfile = argv[i];
      }
   }

   if(filename == NULL || model == NULL || n == 0 || p == 0)
   {
      printf("usage: sgd -m <model> -f <filename> -n <#samples> -p \
<#variables> | -b <beta filename> -p <pred filename> -e <maxepochs> \
-s <stepsize> -l1 <lambda1> -l2 <lambda2> -t <threshold> -v -vv\n");
      return EXIT_FAILURE;
   }

   betahat = calloc(p + 1, sizeof(double));
   gmatrix_init(&g, filename, n, p);

   if(verbose)
      printf("Scaling ... ");
   scale(&g, g.mean, g.sd);
   if(verbose)
      printf("done\n");

   gmatrix_reset(&g);

   if(verbose)
      printf("Starting SGD ...\n");
   sgd_gmatrix(&g, stepsize, maxepochs, betahat,
	 lambda1, lambda2, threshold, verbose);

   gmatrix_reset(&g);
   yhat = malloc(n * sizeof(double));
   predict_logloss(&g, betahat, yhat);

   writebeta(betafile, betahat, p + 1);
   writebeta(predfile, yhat, n);

   printf("AUC: %.5f\n", gmatrix_auc(yhat, &g));

   gmatrix_free(&g);
   free(betahat);
   free(yhat);
   
   return EXIT_SUCCESS;
}

