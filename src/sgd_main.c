#include "sgd.h"

int main(int argc, char* argv[])
{
   int i;
   double *betahat, *betahat_unsc;
   double *yhat_train = NULL, *yhat_test = NULL;
   gmatrix g;
   char *filename = NULL;
   char *model = NULL;
   char *betafile = "beta.csv";
   char *betaunscfile = "beta_unsc.csv";
   char *predtrainfile = "pred_train.csv";
   char *predtestfile = "pred_test.csv";
   char *subsetfile = "subset.csv";
   int n = 0, p = 0;
   int verbose = FALSE;
   int *trainf = NULL, *testf = NULL;
   int ntrain = 0, ntest = 0;
   int cv = 1;
   long seed = time(NULL);
   loss_pt loss_pt_func = NULL;
   predict_pt predict_pt_func = NULL;
   dloss_pt dloss_pt_func = NULL;
   predict_gmatrix predict_gmatrix_func = NULL;
   short inmemory = FALSE;
   short scaleflag = FALSE;
   optim_gmatrix optim_gmatrix_func = sgd_gmatrix;

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
      if(strcmp2(argv[i], "-optim"))
      {
	 i++;
	 if(strcmp2(argv[i], "sgd"))
	    optim_gmatrix_func = sgd_gmatrix;
	 else if(strcmp2(argv[i], "cd"))
	    optim_gmatrix_func = cd_gmatrix;
	 else if(strcmp2(argv[i], "gd"))
	    optim_gmatrix_func = gd_gmatrix;
      }
      else if(strcmp2(argv[i], "-f"))
      {
	 i++;
	 filename = argv[i];
      }
      else if(strcmp2(argv[i], "-model"))
      {
	 i++;
	 model = argv[i];
	 if(strcmp2(model, "logistic"))
	 {
	    loss_pt_func = &logloss_pt;
	    predict_pt_func = &predict_logloss_pt;
	    dloss_pt_func = &logdloss_pt;
	    predict_gmatrix_func = &predict_logloss_gmatrix;
	 }
	 else if(strcmp2(model, "linear"))
	 {
	    loss_pt_func = &l2loss_pt;
	    predict_pt_func = &predict_l2loss_pt;
	    dloss_pt_func = &l2dloss_pt;
	    predict_gmatrix_func = &predict_l2loss_gmatrix;
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
      else if(strcmp2(argv[i], "-epochs"))
      {
	 i++;
	 maxepochs = (int)atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-step"))
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
      else if(strcmp2(argv[i], "-thresh"))
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
      else if(strcmp2(argv[i], "-beta"))
      {
	 i++;
	 betafile = argv[i];
      }
      else if(strcmp2(argv[i], "-betaunsc"))
      {
	 i++;
	 betaunscfile = argv[i];
      }
      else if(strcmp2(argv[i], "-predtrain"))
      {
	 i++;
	 predtrainfile = argv[i];
      }
      else if(strcmp2(argv[i], "-predtest"))
      {
	 i++;
	 predtestfile = argv[i];
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
      else if(strcmp2(argv[i], "-scale"))
      {
	 scaleflag = TRUE;
      }
   }

   if(filename == NULL || model == NULL || n == 0 || p == 0)
   {
      printf("usage: sgd -model <model> -f <filename> -n <#samples> -p \
<#variables> | -beta <beta filename> -pred <pred filename> -epoch <maxepochs> \
-step <stepsize> -l1 <lambda1> -l2 <lambda2> -thresh <threshold> \
-pred <prediction file> -cv <cvfolds> -inmemory -scale -seed <seed> -v -vv\n");
      return EXIT_FAILURE;
   }

   srand48(seed);
   CALLOCTEST2(betahat, p + 1, sizeof(double))
   CALLOCTEST2(betahat_unsc, p + 1, sizeof(double))

   if(!gmatrix_init(&g, inmemory, FALSE, filename, NULL, NULL, n, p))
      return EXIT_FAILURE;
 
   if(scaleflag)
   {
      if(verbose)
	 printf("Scaling ... ");
   
      gmatrix_scale(&g);

      writevectorf("mean.csv", g.mean, p + 1);
      writevectorf("sd.csv", g.sd, p + 1);

      if(verbose)
	 printf("done\n");
   }

   
   gmatrix_reset(&g);

   MALLOCTEST2(trainf, sizeof(int) * g.n)
   MALLOCTEST2(testf, sizeof(int) * g.n)

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

   optim_gmatrix_func(&g, dloss_pt_func, loss_pt_func, predict_pt_func,
	 stepsize, maxepochs, betahat, lambda1, lambda2, threshold,
	 verbose, trainf, trunc);

   
   gmatrix_reset(&g);
   MALLOCTEST2(yhat_train, ntrain * sizeof(double))
   predict_gmatrix_func(&g, betahat, yhat_train, trainf);

   if(ntest > 0)
   {
      gmatrix_reset(&g);
      MALLOCTEST2(yhat_test, ntest * sizeof(double))
      predict_gmatrix_func(&g, betahat, yhat_test, testf);
   }

   /* unscale, return beta to original scale */
   for(i = 1 ; i < p + 1 ; i++)
   {
      betahat_unsc[i] = g.sd[i] * betahat[i] + g.mean[i];
   }

   writevectorf(betafile, betahat, p + 1);
   writevectorf(betaunscfile, betahat_unsc, p + 1);
   writevectorf(predtrainfile, yhat_train, ntrain);
   writevectorl(subsetfile, trainf, g.n);


   printf("###############################\n");

   gmatrix_reset(&g);
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
   
      writevectorf(predtestfile, yhat_test, ntest);
   }

   printf("\n");

   gmatrix_free(&g);
   free(betahat);
   free(betahat_unsc);
   free(yhat_train);
   free(yhat_test);
   free(trainf);
   if(testf)
      free(testf);
   
   return EXIT_SUCCESS;
}

