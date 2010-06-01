#include "sgd.h"

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
   loss_pt loss_pt_func = NULL;
   predict_pt predict_pt_func = NULL;
   dloss dloss_func = NULL;
   short inmemory = FALSE;

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
	 if(strcmp2(model, "logistic"))
	 {
	    loss_pt_func = &logloss_pt;
	    predict_pt_func = &predict_logloss_pt;
	    dloss_func = &logdloss;
	 }
	 else if(strcmp2(model, "linear"))
	 {
	    loss_pt_func = &l2loss_pt;
	    predict_pt_func = &predict_l2loss_pt;
	    dloss_func = &l2dloss;
	 }
	 else
	 {
	    printf("model not available\n");
	    return EXIT_FAILURE;
	 }
	 /*else if(strcmp2(model, "hinge"))
	 {
	 }*/
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
      else if(strcmp2(argv[i], "-inmemory"))
      {
	 inmemory = TRUE;
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
   gmatrix_init(&g, inmemory, FALSE, filename, NULL, NULL, n, p);
 
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
      printf("Parameters: model=%s, dtype=%s inmemory=%d, maxepochs=%d stepsize=%.9f \
lambda1=%.9f lambda2=%.9f \n",
      model, type, inmemory, maxepochs, stepsize, lambda1, lambda2);
      printf("%d training samples, %d test samples\n", ntrain, g.n - ntrain);
   }

   if(verbose)
      printf("Starting SGD ...\n");
   sgd_gmatrix(&g, dloss_func, loss_pt_func, predict_pt_func,
	 stepsize, maxepochs, betahat, lambda1, lambda2, threshold,
	 verbose, trainf, trunc);

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

