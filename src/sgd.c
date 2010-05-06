#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "gmatrix.h"
#include "loss.h"


double accuracy(double *p, int *y, int n)
{
   int s = 0;
   int i;
   for(i = 0 ; i < n ; i++)
      s += (p[i] >= 0.5) == y[i];
   return (double)s / n;
}

double auc(double *p, int *y, int n)
{
   int i, j, k;
   double s = 0;
   int m1 = 0, m2;
   double *y1, *y2;

   for(i = 0 ; i < n ; i++)
      m1 += y[i];

   m2 = n - m1;

   y1 = malloc(m1 * sizeof(double));
   y2 = malloc(m2 * sizeof(double));

   i = j = k = 0;
   while(i < n)
   {
      if(y[i] == 1)
      {
	 y1[j] = y[i];
	 j++;
      }
      else
      {
	 y2[k] = y[i];
	 k++;
      }

      i++;
   }

   for(i = 0 ; i < m1 ; i++)
      for(j = 0 ; j < m2 ; j++)
	 s += y1[i] > y2[j] + 0.5 * (y1[i] == y2[k]);

   return s / (m1 * m2);
}

double sgd_gmatrix(gmatrix *g, int n, int p, double maxstepsize,
      int maxepoch, double *beta, double lambda1, double lambda2, int verbose)
{
   int epoch,
       i, j;
   double *grad = malloc(p * sizeof(double));
   double loss;
   double stepsize = maxstepsize;
   double l;
   sample sm;

   sample_init(&sm, p);

   for(epoch = 1 ; epoch <= maxepoch ; epoch++)
   {
      loss = 0;
      for(i = 0 ; i < n ; i++)
      { 
	 gmatrix_nextrow(g, &sm, TRUE);
	 logdloss(sm.x, beta, sm.y, p, grad);
	 l = logloss_pt(sm.x, beta, sm.y, p);
	 loss += l;
	 for(j = 0 ; j < p ; j++)
	 {
	    beta[j] -= stepsize * (grad[j] 
	       + lambda1 * sign(grad[j]) 
	       + lambda2 * grad[j] * grad[j]);
	 }

      }
      /* stepsize = maxstepsize / (1 + epoch); */
      /*stepsize = stepsize / (1 + maxstepsize);*/
      if(verbose)
	 printf("Epoch %d Loss: %.5f stepsize: %.15f\n", epoch, loss, stepsize);
   }

   sample_free(&sm);

   free(grad);
   return loss;
 
}

double sgd_mem(double **x, int *y, int n, int p, double maxstepsize,
      int maxepoch, double *beta, double lambda1, double lambda2, int verbose)
{
   int epoch,
       i, j;
   double *grad = malloc(p * sizeof(double));
   double loss;
   double stepsize = maxstepsize;
   double l;

   for(epoch = 1 ; epoch <= maxepoch ; epoch++)
   {
      loss = 0;
      for(i = 0 ; i < n ; i++)
      { 
	 logdloss(x[i], beta, y[i], p, grad);
	 l = logloss_pt(x[i], beta, y[i], p);
	 loss += l;
	 for(j = 0 ; j < p ; j++)
	 {
	    beta[j] -= stepsize * (grad[j] 
	       + lambda1 * sign(grad[j]) 
	       + lambda2 * grad[j] * grad[j]);
	 }

      }
      /* stepsize = maxstepsize / (1 + epoch); */
      /*stepsize = stepsize / (1 + maxstepsize);*/
      if(verbose)
	 printf("Epoch %d Loss: %.5f stepsize: %.15f\n", epoch, loss, stepsize);
   }

   free(grad);
   return loss;
}

void predict_logloss(double **x, double *beta, int n, int p, double *yhat)
{
    int i;
    for(i = 0 ; i < n ; i++)
       yhat[i] = 1 / (1 + exp(-dotprod(x[i], beta, p)));
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

void scale(double **x, int n, int p, double **s)
{
   int i, j;
   double *mean = calloc(p, sizeof(double));
   double *sumsq = calloc(p, sizeof(double));
   double delta;
   
   for(j = 0 ; j < p ; j++)
   {
      for(i = 0 ; i < n ; i++)
      {
	 delta = x[i][j] - mean[j];
	 mean[j] = mean[j] + delta / (i + 1);
	 sumsq[j] = sumsq[j] + delta * (x[i][j] - mean[j]);
      }

      for(i = 0 ; i < n ; i++)
	 s[i][j] = (x[i][j] - mean[j]) / sqrt(sumsq[j] / (n - 1));
   }

   free(mean);
   free(sumsq);
}

void scale_test()
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
      
}

double test()
{
   int n = 1e2,
       p = 5; /* not including intercept */
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

   /* Add intercept term to the scaled data */
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
      s = 1; /* intercept  */
      for(j = 0 ; j < p + 1 ; j++)
	 s += x[i][j];
      y[i] = drand48() <= plogis(s) ? 1 : 0;
   }

   loss = sgd_mem(x, y, n, p + 1, 1e-3, 1, betahat, 0, 0, FALSE);
   err = pow(loss - Rloss, 2);
   /*printf("Square-error: %.20f\n", err);*/
   return err;
}

int main()
{
   int n = 1e4,
       p = 1e4; /* not including intercept */
   int i, j;
   /*double **xtmp = malloc(n * sizeof(double*));*/
   double *xtmp = malloc(p * sizeof(double));
   /*double *x;*/
   double *betahat = calloc(p + 1, sizeof(double));
   /*double *tmp = malloc(p * sizeof(double));*/
   int *y;
   double *yhat;
   /* double acc, a; */
   double s;
   int maxepochs = 10;
   FILE* out;
   double const ONE = 1;
   gmatrix g;

   /*assert(test() <= 1e-9);*/

   srand48(12345);

   out = fopen("x.bin", "w+");
   

   /*for(i = 0 ; i < n ; i++)
   {
      xtmp[i] = malloc(p * sizeof(double));
      x[i] = malloc(p * sizeof(double));

      for(j = 0 ; j < p ; j++)
	 xtmp[i] = drand48() - 0.5;
   }*/

   /*printf("Scaling ... ");
   scale(xtmp, n, p, x);
   printf("done\n");*/

   y = malloc(n * sizeof(int));
   /* Add intercept term to the scaled data */
   for(i = 0 ; i < n ; i++)
   {
      
      for(j = 0 ; j < p ; j++)
	 xtmp[j] = drand48() - 0.5;

      /*memcpy(tmp, x[i], sizeof(double) * p);*/
      /*free(x[i]);*/
      /*x[i] = malloc((p + 1) * sizeof(double));*/
      /*x[i][0] = 1;
      for(j = 1 ; j < p + 1; j++)
	 x[i][j] = tmp[j - 1];*/

      fwrite(&ONE, sizeof(double), 1, out);
      fwrite(xtmp, sizeof(double), p, out);

      /* all betas = 1*/
      s = 1;
      for(j = 0 ; j < p ; j++)
	 s += xtmp[j];
      y[i] = drand48() <= plogis(s) ? 1 : 0;
   }
   /*free(tmp);*/
   fflush(out);
   fclose(out);

   /*yhat = malloc(n * sizeof(double));*/
   /*printf("Writing out data ... ");
   writeout("out.csv", x, y, n, p + 1);
   printf("done\n");*/

   gmatrix_init(&g, "x.bin", n, p + 1, y);

   printf("Starting SGD ...\n");
   sgd_gmatrix(&g, n, p + 1, 1e-3, maxepochs, betahat, 0, 0, TRUE);

   printf("gmatrix_free ... ");
   gmatrix_free(&g);
   printf("done\n");

   /*predict_logloss(x, betahat, n, p + 1, yhat);
   acc = accuracy(yhat, y, n);
   a = auc(yhat, y, n);
   printf("Accuracy: %.3f AUC: %.3f\n", acc, a);

   printf("betahat: ");
   for(i = 0 ; i < fmin(p + 1, 10) ; i++)
      printf("%.7f ", betahat[i]);
   printf("\n");*/
   
   /*for(i = 0 ; i < n ; i++)
   {
      free(xtmp[i]);
   }*/
   /*free(x);*/
   free(xtmp);
   free(y);
   /*free(yhat);*/
   free(betahat);
   
   return EXIT_SUCCESS;
}

