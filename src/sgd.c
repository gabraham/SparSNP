#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "gmatrix.h"
#include "sgd.h"
#include "loss.h"
#include "evaluation.h"


double softthreshold(double beta, double step)
{
   double d = beta - step;
   int s = sign(beta);

   if(!s)
      return d;

   if(s == sign(d))
      return d;

   return 0.0;
}

/* Stochastic gradient descent */
double sgd_gmatrix(gmatrix *g, double maxstepsize,
      int maxepoch, double *beta, double lambda1, double lambda2,
      double threshold, int verbose, int *trainf, double trunc)
{
   int epoch = 1, i, j;
   double *grad = malloc((g->p + 1) * sizeof(double));
   dtype *x = malloc((g->p + 1) * sizeof(dtype));
   double prevloss = 0, loss = 0, bestloss = 1e9;
   double stepsize = maxstepsize;
   sample sm;
   double yhat;
   double trainacc = 0;
   double testacc = 0;
   double testloss = 0;
   double diff = 1e9;
   int ntests = 0;
   double ptloss = 0;

   sample_init(&sm, g->p);

   while(epoch <= maxepoch)
   {
      loss = 0;
      testloss = 0;
      testacc = 0;
      trainacc = 0;
      ntests = 0;

      for(i = 0 ; i < g->n ; i++)
      { 
	 gmatrix_nextrow(g, &sm);

	 /* Intercept */
	 x[0] = 1;

	 /* Scale parameters except the intercept */
	 for(j = 0 ; j < g->p ; j++)
	    x[j+1] = (sm.x[j] - g->mean[j]) / g->sd[j];

	 ptloss = logloss_pt(x, beta, sm.y, g->p + 1); 
	 yhat = predict_logloss_pt(&sm, beta, g->mean, g->sd, g->p + 1);

	 /* train */
	 if(trainf[i])
	 {
	    logdloss(x, beta, sm.y, g->p + 1, grad);
	    loss += ptloss;
	    trainacc += (double)((yhat >= 0.5) == (int)sm.y);

	    /* Update weights */
	    for(j = 0 ; j < g->p + 1; j++)
	    {
	       beta[j] -= stepsize * (grad[j] 
	          + lambda1 * sign(beta[j]) 
	          + lambda2 * beta[j] * beta[j]);
	    }
	 }
	 /* test */
	 else
	 {
	    testloss += ptloss;
	    testacc += (double)((yhat >= 0.5) == (int)sm.y);
	    ntests++;
	 }
      }

      trainacc = trainacc / (g->n - ntests);

      if(ntests > 0)
      {
	 testacc = testacc / ntests;
	 testloss = testloss / ntests;
      }
      /*printf("total loss: %.10f over %d samples\n", loss, g->n - ntests);*/
      loss = loss / (g->n - ntests);
      if(bestloss > loss)
	 bestloss = loss;

      /* truncate small weights when lasso is active */
      if(lambda1 > 0)
         for(j = 0 ; j < g->p + 1; j++)
	    if(fabs(beta[j]) <= trunc)
	       beta[j] = 0;

      diff = prevloss - loss;

      if(verbose)
      {
	 printf("Epoch %d  training loss: %.5f diff: %.5f stepsize: %.15f\
 training accuracy: %.8f", epoch, loss, diff, stepsize, trainacc);
	 if(ntests > 0)
	    printf(" test accuracy: %.8f test loss: %.8f", testacc, testloss);
	 printf("\n");
      }
 
      if(epoch > 1 && diff < -threshold)
	 stepsize = fmax(stepsize / 2, 1e-20);
      else if(fabs(diff) <= threshold)
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

double predict_logloss_pt(sample *s, double *beta, double *mean, double *sd, int p)
{
   int i = 0;
   dtype *x = malloc(sizeof(dtype) * p);
   double yhat = 0;

   x[0] = 1;
   for(i = 0 ; i < p - 1 ; i++)
      x[i+1] = (s->x[i] - mean[i]) / sd[i];
   yhat = 1 / (1 + exp(-dotprod(x, beta, p)));

   free(x);
   return yhat;
}

void predict_logloss(gmatrix *g, double *beta, double *yhat, int *trainf)
{
   int i, k;
   sample sm;

   sample_init(&sm, g->p);

   k = 0;
   for(i = 0 ; i < g->n ; i++)
   {
      gmatrix_nextrow(g, &sm);
      if(trainf[i])
      {
	 yhat[k] = predict_logloss_pt(&sm, beta, g->mean, g->sd, g->p + 1);
	 k++;
      }
   } 

   sample_free(&sm);
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

int main(int argc, char* argv[])
{
   int i;
   double *betahat;
   double *yhat_train = NULL, *yhat_test = NULL;
   gmatrix g;
   char *filename = NULL;
   char *model = NULL;
   char *betafile = "beta.csv";
   char *predfile = "pred.csv";
   char *subsetfile = "subset.csv";
   int n = 0, p = 0;
   int verbose = FALSE;
   int *trainf, *testf;
   int ntrain = 0, ntest = 0;
   int cv = 1;
   long seed = time(NULL);
   sample sm;

   /* Parameters */
   int maxepochs = 20;
   double stepsize = 1e-4;
   double lambda1 = 0;
   double lambda2 = 0;
   double threshold = 1e-3;
   double trunc = 1e-9;
   /* double alpha = 0; */

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-f"))
      {
	 i++;
	 filename = argv[i];
      }
      else if(strcmp2(argv[i], "-m"))
      {
	 i++;
	 model = argv[i];
      }
      else if(strcmp2(argv[i], "-n"))
      {
	 i++;
	 n = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-p"))
      {
	 i++;
	 p = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-e"))
      {
	 i++;
	 maxepochs = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-s"))
      {
	 i++;
	 stepsize = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l1"))
      {
	 i++;
	 lambda1 = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-l2"))
      {
	 i++;
	 lambda2 = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-t"))
      {
	 i++;
	 threshold = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-v"))
      {
	 verbose = TRUE;
      }
      else if(strcmp2(argv[i], "-vv"))
      {
	 verbose = 2;
      }
      else if(strcmp2(argv[i], "-b"))
      {
	 i++;
	 betafile = argv[i];
      }
      else if(strcmp2(argv[i], "-pr"))
      {
	 i++;
	 predfile = argv[i];
      }
      else if(strcmp2(argv[i], "-cv"))
      {
	 i++;
	 cv = atoi(argv[i]);
      }
      else if(strcmp2(argv[i], "-seed"))
      {
	 i++;
	 seed = atol(argv[i]);
      }
   }

   if(filename == NULL || model == NULL || n == 0 || p == 0)
   {
      printf("usage: sgd -m <model> -f <filename> -n <#samples> -p \
<#variables> | -b <beta filename> -pr <pred filename> -e <maxepochs> \
-s <stepsize> -l1 <lambda1> -l2 <lambda2> -t <threshold> \
-pr <prediction file> -cv <cvfolds> -v -vv\n");
      return EXIT_FAILURE;
   }

   srand48(seed);
   betahat = calloc(p + 1, sizeof(double));
   gmatrix_init(&g, filename, n, p);
 
   /*if(verbose)
      printf("Scaling ... ");
   scale(&g, g.mean, g.sd);
   if(verbose)
      printf("done\n");*/

   gmatrix_reset(&g);

   trainf = malloc(sizeof(int) * g.n);
   testf = malloc(sizeof(int) * g.n);

   for(i = 0 ; i < g.n ; i++)
   {
      if(cv > 1)
	 trainf[i] = drand48() >= (1.0 / cv);
      else
	 trainf[i] = TRUE;
      ntrain += trainf[i];
      testf[i] = 1 - trainf[i];
   }
   ntest = g.n - ntrain;

   if(verbose)
   {
      printf("Parameters: dtype=%s maxepochs=%d stepsize=%.9f \
lambda1=%.9f lambda2=%.9f \n",
	    type, maxepochs, stepsize, lambda1, lambda2);
      printf("%d training samples, %d test samples\n", ntrain, g.n - ntrain);
   }

   if(verbose)
      printf("Starting SGD ...\n");
   sgd_gmatrix(&g, stepsize, maxepochs, betahat,
	 lambda1, lambda2, threshold, verbose, trainf, trunc);

   gmatrix_reset(&g);
   yhat_train = malloc(ntrain * sizeof(double));
   predict_logloss(&g, betahat, yhat_train, trainf);

   if(ntest > 0)
   {
      gmatrix_reset(&g);
      yhat_test = malloc(ntest * sizeof(double));
      predict_logloss(&g, betahat, yhat_test, testf);
   }

   writevectorf(betafile, betahat, p + 1);
   writevectorf(predfile, yhat_train, ntrain);
   writevectorl(subsetfile, trainf, g.n);

   printf("###############################\n");

   /*gmatrix_reset(&g);
   for(i = 0 ; i < g.n ; i++)
      printf("Y1=%d\n", gmatrix_next_y(&g));

   gmatrix_reset(&g);
   sample_init(&sm, p);
   for(i = 0 ; i < g.n ; i++)
   {
      gmatrix_nextrow(&g, &sm);
      printf("Y2=%d\n", sm.y);
   }*/

   printf("Training AUC (fixed beta): %.5f\n",
	 gmatrix_auc(yhat_train, &g, trainf, ntrain));

   gmatrix_reset(&g);
   printf("Training Accuracy (fixed beta): %.8f\n",
	 gmatrix_accuracy(yhat_train, &g, 0.5, trainf, ntrain));

   printf("\n");

   printf("###############################\n");

   if(ntest > 0)
   {
      gmatrix_reset(&g);
      printf("Test AUC (fixed beta): %.5f\n",
	    gmatrix_auc(yhat_test, &g, testf, ntest));
   
      gmatrix_reset(&g);
      printf("Test Accuracy (fixed beta): %.8f\n",
	    gmatrix_accuracy(yhat_test, &g, 0.5, testf, ntest));
   }

   printf("\n");

   gmatrix_free(&g);
   free(betahat);
   free(yhat_train);
   free(yhat_test);
   
   return EXIT_SUCCESS;
}

