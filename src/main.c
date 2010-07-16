#include "cd.h"
#include "loss.h"
#include "util.h"

int main(int argc, char* argv[])
{
   int i, j;
   double *betahat, *betahat_unsc;
/*   double *yhat_train = NULL, *yhat_test = NULL;*/
   gmatrix g;
   char *filename = NULL;
   char *model = NULL;
   char tmp[100];
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
   /*predict_pt predict_pt_func = NULL;*/
   phi1 phi1_func = NULL;
   phi2 phi2_func = NULL;
   double lambda1max = 1, lambda1min = 1;
   double *lambda1path = NULL;
   int nlambda1 = 100;
   double s;
   double l1minratio = 1e-3;
   short nofit = FALSE;

   /* Parameters */
   int maxepochs = 1000;
   double lambda1 = -1;
   double lambda2 = 0;
   double threshold = 1e-4;
   double trunc = 1e-9;
   int nzmax = 0;
   /* double alpha = 0; */

   for(i = 1 ; i < argc ; i++)
   {
      if(strcmp2(argv[i], "-f"))
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
	    phi1_func = &logphi1;
	    phi2_func = &logphi2;
	 }
	 else if(strcmp2(model, "linear"))
	 {
	    loss_pt_func = &l2loss_pt;
	    phi1_func = &l2phi1;
	    phi2_func = &l2phi2;
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
      else if(strcmp2(argv[i], "-nofit"))
      {
	 nofit = TRUE;
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
      else if(strcmp2(argv[i], "-l1min"))
      {
	 i++;
	 l1minratio = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-thresh"))
      {
	 i++;
	 threshold = atof(argv[i]);
      }
      else if(strcmp2(argv[i], "-nl1"))
      {
	 i++;
	 nlambda1 = atoi(argv[i]);
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
      else if(strcmp2(argv[i], "-nzmax"))
      {
	 i++;
	 nzmax = atol(argv[i]);
      }
   }

   if(filename == NULL || model == NULL || n == 0 || p == 0)
   {
      printf("usage: cd -model <model> -f <filename> -n <#samples> -p \
<#variables> | -beta <beta filename> -pred <pred filename> -epoch <maxepochs> \
-l1 <lambda1> -l2 <lambda2> -thresh <threshold> \
-pred <prediction file> -cv <cvfolds> -seed <seed> -v -vv\n");
      return EXIT_FAILURE;
   }

   nzmax = (int)fmin(n, p);
  
   srand48(seed);
   CALLOCTEST2(betahat, p + 1, sizeof(double))
   CALLOCTEST2(betahat_unsc, p + 1, sizeof(double))

   CALLOCTEST2(lambda1path, nlambda1, sizeof(double))

   if(!gmatrix_init(&g, filename, n, p))
      return EXIT_FAILURE;
 
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
   writevectorl(subsetfile, trainf, g.n);

   if(verbose)
   {
      printf("Parameters: model=%s, maxepochs=%d lambda1=%.9f lambda2=%.9f \n",
      model, maxepochs, lambda1, lambda2);
      printf("%d training samples, %d test samples\n", ntrain, g.n - ntrain);
   }



   if(lambda1 >= 0)
   {
      lambda1max = lambda1;
      l1minratio = 1;
      nlambda1 = 1;
   }
   else
   {
      /* create lambda1 path */
      /* get lambda1 max */
      lambda1max = get_lambda1max_gmatrix(&g, phi1_func, phi2_func);
      if(verbose)
	 printf("lambda1max: %.5f\n", lambda1max);
      lambda1path[0] = lambda1max;
   }
   
   lambda1min = lambda1max * l1minratio;
   lambda1path[nlambda1 - 1] = lambda1min;
   s = (log(lambda1max) - log(lambda1min)) / nlambda1; 
   for(i = 1 ; i < nlambda1 ; i++)
      lambda1path[i] = exp(log(lambda1max) - s * i);

   writevectorf("lambda1path.csv", lambda1path, nlambda1);

   if(!nofit)
   {
      for(i = 0 ; i < nlambda1 ; i++)
      {
         if(verbose)
            printf("\nFitting with lambda1=%.20f\n", lambda1path[i]);
         if(nzmax != 0 && nzmax < cd_gmatrix(
		  &g, phi1_func, phi2_func, loss_pt_func, maxepochs,
		  betahat, lambda1path[i], lambda2, threshold, verbose,
		  trainf, trunc))
	    break;
         /*if(nzmax != 0 && nzmax < cd_gmatrix2(
		  &g, phi1_func, phi2_func, loss_pt_func,
		  maxepochs, betahat, lambda1path[i]))
	    break;*/
         snprintf(tmp, 100, "%s.%d", betafile, i);
         writevectorf(tmp, betahat, p + 1);

	 for(j = 0 ; j < p + 1; j++)
	    betahat[j] = 0;
      }
   }

   /*gmatrix_reset(&g);
   MALLOCTEST2(yhat_train, ntrain * sizeof(double))*/
   /*predict_gmatrix_func(&g, betahat, yhat_train, trainf);

   if(ntest > 0)
   {
      gmatrix_reset(&g);
      MALLOCTEST2(yhat_test, ntest * sizeof(double))
      predict_gmatrix_func(&g, betahat, yhat_test, testf);
   }*/

   /* unscale, return beta to original scale */
   /*for(i = 1 ; i < p + 1 ; i++)
      betahat_unsc[i] = g.sd[i] * betahat[i] + g.mean[i];*/

   /*writevectorf(betafile, betahat, p + 1);
   writevectorf(betaunscfile, betahat_unsc, p + 1);*/
   /*writevectorf(predtrainfile, yhat_train, ntrain);*/


   printf("###############################\n");

   /*gmatrix_reset(&g);
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

   */

   gmatrix_free(&g);
   free(betahat);
   free(betahat_unsc);
   /*free(yhat_train);*/
/*   free(yhat_test);*/
   free(trainf);
   free(lambda1path);
   if(testf)
      free(testf);
   
   return EXIT_SUCCESS;
}

