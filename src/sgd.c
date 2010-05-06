#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>

/*typedef struct gmatrix */

#define sign(x) ((x > 0) - (x < 0))

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
   {
      for(j = 0 ; j < m2 ; j++)
      {
	 s += y1[i] > y2[j] + 0.5 * (y1[i] == y2[k]);
      }
   }

   return s / (m1 * m2);
}

double plogis(double x)
{
   return 1 / (1 + exp(-x));
}

double dotprod(double *a, double *b, int m)
{
   int i;
   double s = 0;
   for(i = 0 ; i < m ; i++)
      s += a[i] * b[i];
   return s;
}

double logloss(double **x, double *beta, int *y, int n, int p)
{
   int i;
   double loss = 0;
   for(i = 0 ; i < n ; i++)
      loss += -(double)y[i] * dotprod(x[i], beta, p) 
	    + log(1 + exp(dotprod(x[i], beta, p)));
   return loss;
}

void logdloss(double *x, double *beta, int y, int p, double* grad)
{
   int i;
   double pr = exp(dotprod(x, beta, p));
   for(i = 0 ; i < p ; i++)
      grad[i] = x[i] * pr / (1 + pr) - (double)y;
}

void sgd(double **x, int *y, int n, int p, double maxstepsize,
      int maxepoch, double *beta, double lambda1, double lambda2)
{
   int epoch,
       i, j;
   double *grad = malloc((p + 1)* sizeof(double));
   double loss = 0.0;
   double stepsize = maxstepsize;

   for(epoch = 1 ; epoch <= maxepoch ; epoch++)
   {
      for(i = 0 ; i < n ; i++)
      { 
	 logdloss(x[i], beta, y[i], p, grad);
	 for(j = 0 ; j < p ; j++)
	 {
	    beta[j] -= stepsize * (grad[j] 
	       + lambda1 * sign(grad[j]) 
	       + lambda2 * grad[j] * grad[j]);
	 }
      }
      /*stepsize = maxstepsize / (1 + epoch); */
      stepsize = stepsize / (1 + maxstepsize); 
      loss = logloss(x, beta, y, n, p);
      printf("Epoch %d Loss: %.4f stepsize: %.15f\n", epoch, loss, stepsize);
   }

   free(grad);
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
	    fprintf(out, "%.9f,", x[i][j]);
	 else
	    fprintf(out, "%.9f\n", x[i][j]);
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
      printf("%.5f ", mean[j]);

      for(i = 0 ; i < n ; i++)
	 s[i][j] = (x[i][j] - mean[j]) / sqrt(sumsq[j] / (n - 1));
   }
   printf("\n");

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

int main()
{
   int n = 1e1,
       p = 3; /* not including intercept */
   int i, j;
   double **xtmp = malloc(n * sizeof(double*));
   double **x = malloc(n * sizeof(double*));
   double *betahat = malloc((p + 1) * sizeof(double));
   double *tmp = malloc(p * sizeof(double));
   int *y;
   double *yhat;
   double acc;
   double s;
   int epoch;
   double a;

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
   yhat = malloc(n * sizeof(double));
   for(i = 0 ; i < n ; i++)
   {
      s = 1; /* intercept  */
      for(j = 0 ; j < p + 1 ; j++)
	 s += x[i][j];
      y[i] = drand48() <= plogis(s) ? 1 : 0;
   }

   writeout("out.csv", x, y, n, p + 1);

   /*for(epoch = 1 ; epoch <= 5 ; epoch++) */
   epoch = 10;
   {
      for(j = 0 ; j < p + 1; j++)
         betahat[j] = 0.0;
      sgd(x, y, n, p + 1, 1e-4, epoch, betahat, 0, 0);
      predict_logloss(x, betahat, n, p + 1, yhat);
      /* for(i = 0 ; i < n ; i++)
	 printf("%d %.5f %d\n", i, yhat[i], y[i]);  */
      acc = accuracy(yhat, y, n);
      a = auc(yhat, y, n);
      printf("Accuracy: %.3f AUC: %.3f\n\n", acc, a);
   }

   for(i = 0 ; i < p + 1 ; i++)
      printf("%.4f ", betahat[i]);
   printf("\n");
   
   for(i = 0 ; i < n ; i++)
   {
      free(x[i]);
      free(xtmp[i]);
   }
   free(x);
   free(xtmp);
   free(y);
   free(yhat);
   free(betahat);
   
   return EXIT_SUCCESS;
}

